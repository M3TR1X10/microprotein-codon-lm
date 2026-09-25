"""Freeze or run the bounded secondary matrix, one isolated CPU worker at a time."""
import argparse
import copy
import ctypes
from datetime import datetime, timezone
import math
import os
from pathlib import Path
import subprocess
import sys
import time

from microprotein_lm.io import now, read_json, sha256
from microprotein_lm.secondary_io import write_json


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix('.partial.json')
    write_json(temporary, value)
    temporary.replace(path)


def cohort_identity(root, dataset):
    directory = root/f'data/processed/secondary/{dataset}'
    manifest = read_json(directory/'manifest.json')
    actual = sha256(directory/'cohort.jsonl')
    if manifest['cohort_sha256'] != actual:
        raise ValueError(f'Cohort and manifest disagree: {dataset}')
    return {'cohort_sha256': actual, 'manifest_sha256': sha256(directory/'manifest.json')}


def verify_cohort(root, job):
    identity = cohort_identity(root, job['dataset'])
    if any(identity[key] != job[key] for key in identity):
        raise ValueError(f'Frozen cohort or manifest changed: {job["dataset"]}')


def verify_inputs(root, plan):
    from microprotein_lm.secondary_train import engine_source_hashes
    if engine_source_hashes() != plan['engine_source_hashes']:
        raise ValueError('Engine source changed after plan freeze')
    for name, checksum in plan['input_hashes'].items():
        if sha256(root/name) != checksum:
            raise ValueError(f'Plan input changed: {name}')


def verify_summary(summary, job, plan):
    expected = {'config': job['config'], 'mode': job['mode'], 'seed': job['seed'],
                'dataset': job['dataset'], 'cohort_sha256': job['cohort_sha256'],
                'source_hashes': plan['engine_source_hashes']}
    if any(summary.get(k) != v for k, v in expected.items()):
        raise ValueError(f'Existing summary does not match frozen job: {job["arm"]}')
    if any(summary.get('data_identity', {}).get(k) != job[k]
           for k in ('cohort_sha256', 'manifest_sha256')):
        raise ValueError(f'Existing summary metadata differs: {job["arm"]}')
    updates = summary.get('updates')
    if (type(summary.get('complete')) is not bool or type(updates) is not int
            or not 0 < updates <= job['config']['updates']
            or summary['complete'] != (updates == job['config']['updates'])):
        raise ValueError(f'Inconsistent completion state: {job["arm"]}')


def recovered_elapsed(prior, plan_hash):
    if not prior:
        return 0.0, None
    if prior.get('plan_sha256') != plan_hash:
        raise ValueError('Execution state belongs to a different frozen plan')
    elapsed = prior.get('elapsed_seconds', 0)
    if not isinstance(elapsed, (int, float)) or not math.isfinite(elapsed) or elapsed < 0:
        raise ValueError('Invalid accumulated execution time')
    active = prior.get('in_progress')
    if not active:
        return float(elapsed), prior.get('recovery_accounting')
    # After interruption, charge downtime too rather than erase possibly spent
    # training time. This conservative rule never silently replenishes a budget.
    began = datetime.fromisoformat(active['started_utc'])
    if began.tzinfo is None:
        raise ValueError('In-progress execution timestamp must include its timezone')
    charged = max(0.0, (datetime.now(timezone.utc)-began).total_seconds())
    return float(elapsed)+charged, {'interrupted_job': active['job_index'],
        'extra_wall_seconds_charged': charged,
        'policy': 'Conservatively charge elapsed wall time since interrupted worker launch, including downtime'}


def make_plan(root):
    from microprotein_lm.secondary_train import engine_source_hashes, _configuration
    base = read_json(root / 'configs/secondary/training.json')
    definitions = [
        ('human_atp8_base', 'human_atp8', 'base', {}),
        ('human_atp8_codon', 'human_atp8', 'codon', {}),
        ('human_complex_codon', 'human_complex', 'codon', {}),
        ('mammal_complex_codon', 'mammal_complex', 'codon', {}),
        ('vertebrate_complex_codon', 'vertebrate_complex', 'codon', {}),
        ('vertebrate_small', 'vertebrate_complex', 'codon',
         {'model': {'d_model': 48, 'n_heads': 4, 'n_layers': 2, 'dropout': 0.1}}),
        ('vertebrate_family_macro', 'vertebrate_complex', 'codon', {'loss': 'family_macro'}),
        ('vertebrate_no_position', 'vertebrate_complex', 'codon', {'position_encoding': False}),
    ]
    jobs = []
    for seed in base['seeds']:
        for arm, dataset, mode, overrides in definitions:
            config = _configuration({**copy.deepcopy(base), **overrides, 'dataset': dataset})
            jobs.append({'arm': arm, 'dataset': dataset, 'mode': mode, 'seed': seed,
                         'config': config, 'role': 'matched',
                         **cohort_identity(root, dataset)})
    config = _configuration({**copy.deepcopy(base), 'dataset': 'full_vertebrate',
              'model': {'d_model': 384, 'n_heads': 8, 'n_layers': 6, 'dropout': 0.1}})
    jobs.append({'arm': 'full_vertebrate_capacity', 'dataset': 'full_vertebrate', 'mode': 'codon',
                 'seed': 17, 'config': config, 'role': 'maximum_pool_engineering',
                 **cohort_identity(root, 'full_vertebrate')})
    plan = {'created_utc': now(), 'scope': 'training_only', 'aggregate_budget_seconds': 14400,
            'shutdown_reserve_seconds': 120, 'engine_source_hashes': engine_source_hashes(),
            'input_hashes': {p: sha256(root/p) for p in ['configs/secondary/biology.json',
                'configs/secondary/training.json', 'docs/secondary-protocol.md',
                'reports/secondary/cohorts.json', 'reports/secondary/hardware-benchmark.json',
                'reports/secondary/acquisition.json', 'reports/secondary/genome-acquisition.json',
                'src/microprotein_lm/secondary_acquire.py', 'src/microprotein_lm/secondary_data.py',
                'src/microprotein_lm/secondary_genomes.py',
                'scripts/build_secondary_data.py', 'scripts/run_secondary.py']},
            'jobs': jobs,
            'budget_basis': '2000 updates, measured short CPU profiles. Matched matrix ~3 hours including 25% overhead; larger full-pool arm adds approximately 22 minutes. All training serialized, four threads, one interop thread; deadline preserves partial runs without claiming completion.',
            'selection_warning': 'Largest profiled practical matched body, not proof of an absolute hardware maximum. Full-pool engineering arm changes both capacity and dataset and is excluded from matched causal comparisons.'}
    target = root/'reports/secondary/plan.json'
    if target.exists():
        old = read_json(target)
        for key in ['jobs', 'input_hashes', 'engine_source_hashes', 'aggregate_budget_seconds',
                    'shutdown_reserve_seconds', 'scope']:
            if old[key] != plan[key]:
                raise ValueError('Frozen plan exists with different inputs; use a new experiment version')
        return old
    atomic_json(target, plan)
    return plan


def lower_priority():
    if os.name == 'nt':
        from ctypes import wintypes
        kernel = ctypes.WinDLL('kernel32', use_last_error=True)
        kernel.GetCurrentProcess.restype = wintypes.HANDLE
        kernel.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
        if not kernel.SetPriorityClass(kernel.GetCurrentProcess(), 0x4000):
            raise OSError('Cannot set below-normal worker priority')


def worker(args, root, plan):
    import torch
    from microprotein_lm.secondary_train import train_secondary
    verify_inputs(root, plan)
    if args.worker < 0 or args.worker >= len(plan['jobs']):
        raise ValueError('Worker index is outside the frozen plan')
    if (args.remaining_seconds is None or not math.isfinite(args.remaining_seconds)
            or args.remaining_seconds <= 0):
        raise ValueError('Worker requires a positive finite remaining budget')
    job = plan['jobs'][args.worker]
    verify_cohort(root, job)
    lower_priority()
    torch.set_num_interop_threads(1)
    output = root/'runs/secondary'/f"{job['arm']}-seed{job['seed']}"
    resume = (output/'checkpoint.pt').exists()
    train_secondary(job['config'], job['mode'], job['seed'],
                    root/'data/processed/secondary'/job['dataset'], output,
                    resume=resume, max_seconds=args.remaining_seconds)


def execute(root, plan):
    verify_inputs(root, plan)
    progress_path = root/'reports/secondary/execution.json'
    results_path = root/'reports/secondary/results.json'
    plan_hash = sha256(root/'reports/secondary/plan.json')
    prior = read_json(progress_path) if progress_path.exists() else {}
    elapsed_before, recovery = recovered_elapsed(prior, plan_hash)
    started = time.monotonic()
    summaries = []

    def persist(status, **extra):
        elapsed = elapsed_before + time.monotonic() - started
        atomic_json(results_path, summaries)
        atomic_json(progress_path, {'updated_utc': now(), 'elapsed_seconds': elapsed,
            'status': status, 'completed_runs': sum(s['complete'] for s in summaries),
            'planned_runs': len(plan['jobs']), 'complete': status == 'complete',
            'plan_sha256': plan_hash, 'recovery_accounting': recovery, 'in_progress': None, **extra})
        return elapsed

    for index, job in enumerate(plan['jobs']):
        verify_inputs(root, plan)
        verify_cohort(root, job)
        output = root/'runs/secondary'/f"{job['arm']}-seed{job['seed']}"
        summary_file = output/'summary.json'
        if summary_file.exists():
            existing = read_json(summary_file)
            verify_summary(existing, job, plan)
            if existing['complete']:
                summaries.append({'arm': job['arm'], 'role': job['role'], **existing})
                continue
        elapsed = elapsed_before + time.monotonic() - started
        remaining = plan['aggregate_budget_seconds'] - elapsed - plan['shutdown_reserve_seconds']
        if remaining <= 0:
            break
        print(f"START {index+1}/{len(plan['jobs'])} {job['arm']} seed={job['seed']}; remaining budget {remaining/60:.1f} min", flush=True)
        persist('running', last_job=index,
                in_progress={'job_index': index, 'started_utc': now()})
        hard_remaining = plan['aggregate_budget_seconds']-(elapsed_before+time.monotonic()-started)
        if hard_remaining <= 0:
            break
        try:
            run = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--root', str(root),
                                  '--worker', str(index), '--remaining-seconds', str(remaining)],
                                 cwd=root, timeout=hard_remaining)
        except subprocess.TimeoutExpired as error:
            # subprocess.run kills and waits for its child on timeout. Atomic
            # checkpoint writes leave the preceding checkpoint intact.
            persist('timed_out', last_job=index, worker_timeout_seconds=hard_remaining,
                    failure='Worker exceeded the remaining aggregate wall-time budget; latest checkpoint retained')
            raise RuntimeError(f'Aggregate deadline reached in {job["arm"]}; preserved outputs are incomplete') from error
        if run.returncode:
            persist('worker_failed', last_job=index, worker_returncode=run.returncode)
            raise RuntimeError(f"Worker failed: {job['arm']} seed {job['seed']}; inspect preserved outputs")
        try:
            summary = read_json(summary_file)
            verify_summary(summary, job, plan)
        except (OSError, ValueError, KeyError) as error:
            persist('invalid_output', last_job=index, failure=str(error))
            raise
        summaries.append({'arm': job['arm'], 'role': job['role'], **summary})
        persist('running' if summary['complete'] else 'budget_stopped', last_job=index,
                worker_returncode=run.returncode)
        if not summary['complete']:
            break
    complete = len(summaries) == len(plan['jobs']) and all(s['complete'] for s in summaries)
    elapsed = persist('complete' if complete else 'budget_stopped')
    print(f"FINISHED {sum(s['complete'] for s in summaries)}/{len(plan['jobs'])} complete runs; {elapsed/60:.1f} minutes", flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', default='.')
    parser.add_argument('--freeze-only', action='store_true')
    parser.add_argument('--worker', type=int)
    parser.add_argument('--remaining-seconds', type=float)
    args = parser.parse_args()
    root = Path(args.root).resolve()
    if args.worker is not None:
        worker(args, root, read_json(root/'reports/secondary/plan.json'))
    else:
        plan = make_plan(root)
        if not args.freeze_only:
            execute(root, plan)
        else:
            print(f"Frozen {len(plan['jobs'])} runs before fitting")
