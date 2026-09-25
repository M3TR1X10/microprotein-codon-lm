"""Sequential, CPU-only checkpoint integrity audit after secondary training.

No fitting, inference, loss calculation or new biological diagnostics occur.
The first recorded file hashes establish future artifact identity; they are not
independent evidence that training was scientifically successful or reproduced.
"""
import argparse
from datetime import datetime, timezone
import gc
import hashlib
import importlib.util
from pathlib import Path
import pickle
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))

import torch

from microprotein_lm.io import now, read_json, sha256
from microprotein_lm.secondary_io import write_json
from microprotein_lm.secondary_train import SecondaryCorpus, build_secondary_model


def load_runner(root):
    spec = importlib.util.spec_from_file_location('secondary_runner_checkpoint_audit', root/'scripts/run_secondary.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def tensor_inventory(value, prefix='checkpoint'):
    """Check every serialized floating/complex tensor, including optimizer state."""
    count, elements = 0, 0
    if isinstance(value, torch.Tensor):
        count, elements = 1, value.numel()
        if (value.is_floating_point() or value.is_complex()) and not bool(torch.isfinite(value).all()):
            raise ValueError(f'Nonfinite tensor: {prefix}')
    elif isinstance(value, dict):
        for key, child in value.items():
            n, size = tensor_inventory(child, f'{prefix}.{key}')
            count, elements = count+n, elements+size
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            n, size = tensor_inventory(child, f'{prefix}[{index}]')
            count, elements = count+n, elements+size
    return count, elements


def audit_file(path, job, summary, corpus, plan):
    before = path.stat()  # Missing required completed checkpoints fail explicitly.
    checksum = sha256(path)
    checkpoint = torch.load(path, map_location='cpu', weights_only=True)
    if not isinstance(checkpoint, dict):
        raise ValueError(f'Checkpoint payload must be a mapping: {path.name}')
    expected_identity = {'config': job['config'], 'mode': job['mode'], 'seed': job['seed'],
        'data_identity': corpus.identity, 'source_hashes': plan['engine_source_hashes']}
    if checkpoint.get('identity') != expected_identity:
        raise ValueError(f'Checkpoint identity disagrees with frozen plan/cohort: {path.name}')
    if checkpoint.get('update') != summary['updates'] or checkpoint['update'] != job['config']['updates']:
        raise ValueError(f'Checkpoint update disagrees with completed budget: {path.name}')
    history = checkpoint.get('history')
    if not isinstance(history, list) or not history or history[0] != summary['initial'] or history[-1] != summary['final']:
        raise ValueError(f'Checkpoint history endpoints disagree with summary: {path.name}')
    if history[0].get('update') != 0 or history[-1].get('update') != checkpoint['update']:
        raise ValueError(f'Checkpoint history has inconsistent update endpoints: {path.name}')
    for key in ('bases_seen', 'sequences_seen', 'per_family_exposures', 'sampling_trace_sha256'):
        if checkpoint.get(key) != summary[key]:
            raise ValueError(f'Checkpoint {key} disagrees with summary: {path.name}')
    if summary['data_identity'] != corpus.identity:
        raise ValueError('Summary biological metadata does not match the frozen cohort')
    tensors, elements = tensor_inventory(checkpoint)
    model = build_secondary_model(job['config'], corpus)
    expected_state = model.state_dict()
    state = checkpoint['model']
    # strict=True enforces names/shapes, but can silently cast a supplied dtype.
    if not isinstance(state, dict) or set(state) != set(expected_state) or any(
            not isinstance(state[key], torch.Tensor) or state[key].dtype != expected_state[key].dtype
            for key in expected_state):
        raise ValueError(f'Checkpoint model keys/dtypes disagree with frozen architecture: {path.name}')
    if not torch.equal(state['embedding.weight'], state['head.weight']):
        raise ValueError(f'Tied embedding/output weights disagree: {path.name}')
    if not job['config'].get('position_encoding', True) and bool(state['position'].count_nonzero()):
        raise ValueError(f'No-position checkpoint contains explicit positional values: {path.name}')
    model.load_state_dict(state, strict=True)
    parameters = sum(p.numel() for p in model.parameters())
    if parameters != summary['parameters']:
        raise ValueError(f'Checkpoint architecture parameter count disagrees with summary: {path.name}')
    state_hash = hashlib.sha256()
    for key in sorted(state):
        tensor = state[key].detach().contiguous()
        state_hash.update(f'{key}|{tensor.dtype}|{tuple(tensor.shape)}\n'.encode())
        state_hash.update(tensor.numpy().tobytes())
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f'Checkpoint changed during audit: {path.name}')
    return {'sha256': checksum, 'bytes': after.st_size, 'parameters': parameters,
            'model_state_sha256': state_hash.hexdigest(),
            'tensor_entries_checked': tensors, 'tensor_elements_checked': elements,
            'update': checkpoint['update'], 'strict_model_load': True, 'finite_tensors': True,
            'identity_and_summary_match': True}


def verify_established_files(files, previous):
    """Never replace previously established checkpoint byte/state identities silently."""
    old_files = previous['files']
    if set(old_files) != {'final.pt', 'checkpoint.pt'}:
        raise ValueError('Established inventory has an unexpected checkpoint-file set')
    for name in old_files:
        for field in ('file', 'sha256', 'bytes', 'model_state_sha256'):
            if files[name].get(field) != old_files[name].get(field):
                raise ValueError(f'Established checkpoint identity drift: {name} {field}')


def audit(root):
    root = Path(root).resolve()
    torch.set_num_threads(1)
    # CLI runs in a new process; imports above execute no tensor operations.
    torch.set_num_interop_threads(1)
    destination = root/'reports/secondary/checkpoints.json'
    prior_exists = destination.exists()
    report = {'created_utc': now(), 'scope': 'CHECKPOINT INTEGRITY ONLY; no fitting, inference or new training diagnostics',
        'audit_script_sha256': sha256(Path(__file__)), 'device': 'cpu', 'cpu_threads': 1,
        'status': 'failed', 'verified_runs': [], 'failed_runs': [], 'excluded_runs': [],
        'limitations': [
            'Checks loadability, finite serialized tensors, frozen architecture and metadata consistency; does not reproduce optimization.',
            'Hashes establish identity at audit time and cannot independently prove earlier unrecorded file contents.',
            'Re-audits enforce established file and model-state hashes for the same frozen plan; failure preserves that inventory.',
            'No inference, model selection, held-out evaluation, biological-function or structural claim is made.',
            'Only completed runs present in the frozen-plan results are audited; excluded runs are listed explicitly.',
            'An outer timeout may leave an older rolling checkpoint for an incomplete run; such runs are not relabeled complete.']}
    try:
        plan_path, results_path = root/'reports/secondary/plan.json', root/'reports/secondary/results.json'
        plan, results = read_json(plan_path), read_json(results_path)
        report.update(plan_sha256=sha256(plan_path), results_sha256=sha256(results_path),
                      engine_source_hashes=plan['engine_source_hashes'], planned_runs=len(plan['jobs']))
        prior_verified = {}
        if prior_exists:
            previous = read_json(destination)
            report['previous_inventory_sha256'] = sha256(destination)
            if previous.get('plan_sha256') != report['plan_sha256']:
                raise ValueError('Established checkpoint inventory belongs to a different frozen plan')
            for verified in previous.get('verified_runs', []):
                key = verified['arm'], verified['seed']
                if key in prior_verified:
                    raise ValueError('Duplicate verified run in established checkpoint inventory')
                prior_verified[key] = verified
            report['previous_verified_run_count'] = len(prior_verified)
        if plan.get('scope') != 'training_only':
            raise ValueError('Expected a frozen training-only plan')
        if sha256(root/'scripts/run_secondary.py') != plan['input_hashes']['scripts/run_secondary.py']:
            raise ValueError('Frozen runner source changed; refusing to import it for verification')
        runner = load_runner(root)
        runner.verify_inputs(root, plan)
        jobs = {(job['arm'], job['seed']): job for job in plan['jobs']}
        if len(jobs) != len(plan['jobs']):
            raise ValueError('Duplicate jobs in frozen plan')
        if set(prior_verified)-set(jobs):
            raise ValueError('Established inventory contains a run outside the frozen plan')
        indexed = {}
        for summary in results:
            key = summary['arm'], summary['seed']
            if key not in jobs or key in indexed:
                raise ValueError(f'Unplanned or duplicate result: {key}')
            indexed[key] = summary
        for key, job in jobs.items():
            label = {'arm': key[0], 'seed': key[1], 'dataset': job['dataset']}
            if key not in indexed:
                if key in prior_verified:
                    report['failed_runs'].append({**label, 'error_type': 'ValueError',
                        'error': 'Previously verified completed run is no longer present in reported results'})
                else:
                    report['excluded_runs'].append({**label, 'reason': 'no_reported_result'})
                continue
            summary = indexed[key]
            try:
                runner.verify_cohort(root, job)
                runner.verify_summary(summary, job, plan)
                if summary.get('role') != job['role']:
                    raise ValueError('Result role disagrees with frozen plan')
                if not summary['complete']:
                    if key in prior_verified:
                        raise ValueError('Previously verified completed run is now reported incomplete')
                    report['excluded_runs'].append({**label, 'reason': 'incomplete', 'updates': summary['updates'],
                                                    'stop_reason': summary.get('stop_reason')})
                    continue
                directory = root/'runs/secondary'/f'{key[0]}-seed{key[1]}'
                local_summary = read_json(directory/'summary.json')
                runner.verify_summary(local_summary, job, plan)
                if any(summary.get(k) != v for k, v in local_summary.items()):
                    raise ValueError('Local summary differs from reported result')
                corpus = SecondaryCorpus(root/'data/processed/secondary'/job['dataset'], job['mode'])
                files = {}
                # Load, check and release each payload/model before the next file.
                for name in ('final.pt', 'checkpoint.pt'):
                    path = directory/name
                    files[name] = {'file': path.relative_to(root).as_posix(),
                                   **audit_file(path, job, summary, corpus, plan)}
                    gc.collect()
                if files['final.pt']['model_state_sha256'] != files['checkpoint.pt']['model_state_sha256']:
                    raise ValueError('Final and rolling checkpoint model states disagree')
                if key in prior_verified:
                    verify_established_files(files, prior_verified[key])
                report['verified_runs'].append({**label, 'summary_sha256': sha256(directory/'summary.json'),
                    'cohort_sha256': job['cohort_sha256'], 'manifest_sha256': job['manifest_sha256'], 'files': files})
            except (OSError, ValueError, RuntimeError, KeyError, TypeError, EOFError, pickle.UnpicklingError) as error:
                report['failed_runs'].append({**label, 'error_type': type(error).__name__, 'error': str(error)})
                gc.collect()
        if sha256(plan_path) != report['plan_sha256'] or sha256(results_path) != report['results_sha256']:
            raise ValueError('Frozen plan or reported results changed during audit')
        report['status'] = 'passed' if not report['failed_runs'] else 'failed'
    except (OSError, ValueError, RuntimeError, KeyError, TypeError, EOFError, pickle.UnpicklingError) as error:
        report['global_error'] = {'error_type': type(error).__name__, 'error': str(error)}
    report['verified_run_count'] = len(report['verified_runs'])
    report['failed_run_count'] = len(report['failed_runs'])
    report['excluded_run_count'] = len(report['excluded_runs'])
    report['all_planned_runs_verified'] = (report['status'] == 'passed' and
        report['verified_run_count'] == report.get('planned_runs', -1))
    if prior_exists and report['status'] != 'passed':
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        output = destination.with_name(f'checkpoints-attempt-{stamp}.json')
        report['established_inventory_preserved'] = True
    else:
        output = destination
        report['established_inventory_preserved'] = False
    report['output_file'] = output.relative_to(root).as_posix()
    temporary = output.with_suffix('.partial.json')
    write_json(temporary, report)
    temporary.replace(output)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default=str(ROOT))
    args = parser.parse_args()
    result = audit(args.root)
    print(f"Checkpoint audit {result['status']}: {result['verified_run_count']} verified, "
          f"{result['failed_run_count']} failed, {result['excluded_run_count']} excluded", flush=True)
    print(f"Audit record: {result['output_file']}", flush=True)
    raise SystemExit(0 if result['status'] == 'passed' else 1)
