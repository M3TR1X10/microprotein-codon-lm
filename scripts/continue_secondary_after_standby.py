"""Archive and apply the explicitly documented one-time standby budget amendment.

No model, checkpoint, frozen plan, engine or scientific setting is changed.
"""
import argparse
from datetime import datetime, timezone
import importlib.util
import math
from pathlib import Path
import shutil

from microprotein_lm.io import now, read_json, sha256
from microprotein_lm.secondary_io import write_json

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = 'reports/secondary/interruption'


def utc(value):
    parsed = datetime.fromisoformat(value.replace('Z', '+00:00'))
    if parsed.tzinfo is None:
        raise ValueError('Power evidence timestamps must have explicit time zones')
    return parsed.astimezone(timezone.utc)


def archive_pairs():
    pairs = [(f'reports/secondary/{name}.json', f'{ARCHIVE}/original-{name}.json')
             for name in ('execution', 'results')]
    pairs.extend((f'runs/secondary/human_atp8_base-seed29/{name}',
                  f'runs/secondary-interruption/20260925-human-atp8-base-seed29/{name}')
                 for name in ('checkpoint.pt', 'summary.json', 'metrics.json', 'run.json', 'baselines.json'))
    return pairs


def credit_from_events(evidence, worker_start, execution_stop):
    if evidence['provider'] != 'Microsoft-Windows-Kernel-Power' or evidence['channel'] != 'System':
        raise ValueError('Unexpected power evidence source')
    events = evidence['events']
    if len(events) != 4 or [e['record_id'] for e in events] != [52915, 52920, 52922, 52931]:
        raise ValueError('This amendment accepts only the reviewed interruption records')
    total_us, previous_end = 0, worker_start
    for start, end in zip(events[::2], events[1::2]):
        began, ended = utc(start['utc']), utc(end['utc'])
        if (start['event_id'] != 506 or end['event_id'] != 507 or start['scenario'] != end['scenario']
                or not end['sleep_entered'] or not previous_end <= began < ended <= execution_stop):
            raise ValueError('Standby windows must be paired, nonoverlapping and wholly inside the worker lifetime')
        duration = (ended-began).total_seconds()
        hardware = end['hardware_drips_us']
        if type(hardware) is not int or not 0 < hardware <= end['software_drips_us'] <= end['duration_us']:
            raise ValueError('Invalid hardware deep-idle residency')
        if abs(duration-end['duration_us']/1e6) > .1 or hardware/1e6 > duration+.01:
            raise ValueError('Event residency and measured interval disagree')
        total_us += hardware
        previous_end = ended
    return total_us/1e6


def checked_copy(root, source, destination):
    source_path, destination_path = root/source, root/destination
    if destination_path.exists():
        raise ValueError(f'Archive already exists: {destination}')
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    checksum = sha256(source_path)
    shutil.copy2(source_path, destination_path)
    if sha256(destination_path) != checksum or sha256(source_path) != checksum:
        raise ValueError(f'Archive copy changed: {source}')
    return {'source': source, 'archive': destination, 'sha256': checksum, 'bytes': source_path.stat().st_size}


def verify_frozen(root, plan):
    runner_path = root/'scripts/run_secondary.py'
    if sha256(runner_path) != plan['input_hashes']['scripts/run_secondary.py']:
        raise ValueError('Frozen runner changed')
    spec = importlib.util.spec_from_file_location('standby_original_runner', runner_path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    runner.verify_inputs(root, plan)
    return runner


def prepare(root):
    destination = root/ARCHIVE/'amendment.json'
    if destination.exists():
        raise ValueError('Amendment is already prepared; never overwrite its original evidence')
    plan = read_json(root/'reports/secondary/plan.json')
    execution = read_json(root/'reports/secondary/execution.json')
    results = read_json(root/'reports/secondary/results.json')
    runner = verify_frozen(root, plan)
    if (execution['status'] != 'budget_stopped' or execution['in_progress'] is not None
            or execution['plan_sha256'] != sha256(root/'reports/secondary/plan.json')):
        raise ValueError('Expected the stopped original execution')
    partial = [s for s in results if not s['complete']]
    if len(partial) != 1 or (partial[0]['arm'], partial[0]['seed'], partial[0]['updates']) != ('human_atp8_base', 29, 268):
        raise ValueError('Unexpected interrupted run; this amendment is incident-specific')
    summary = partial[0]
    job = next(j for j in plan['jobs'] if (j['arm'], j['seed']) == (summary['arm'], summary['seed']))
    runner.verify_summary(summary, job, plan)
    runner.verify_cohort(root, job)
    run_directory = 'runs/secondary/human_atp8_base-seed29'
    run = read_json(root/run_directory/'run.json')
    evidence_path = root/ARCHIVE/'power-evidence.json'
    credit = credit_from_events(read_json(evidence_path), utc(run['created_utc']), utc(execution['updated_utc']))
    original = execution['elapsed_seconds']
    charged = original-credit
    if not 0 < credit < summary['elapsed_seconds'] or not 0 < charged < plan['aggregate_budget_seconds']:
        raise ValueError('Correction does not produce a valid conservative remaining budget')
    archives = [checked_copy(root, source, destination) for source, destination in archive_pairs()]
    if any(sha256(root/item['source']) != item['sha256'] for item in archives):
        raise ValueError('Original execution changed while preparing its archive')
    amendment = {'created_utc': now(), 'kind': 'standby_budget_amendment',
        'scope': 'One-time time accounting only; original model/data/seed/update/RNG plan unchanged',
        'plan_sha256': sha256(root/'reports/secondary/plan.json'),
        'script_sha256': sha256(Path(__file__)), 'evidence_sha256': sha256(evidence_path),
        'document_sha256': sha256(root/'docs/secondary-runtime-amendment.md'),
        'original_elapsed_wall_seconds': original, 'excluded_hardware_deep_idle_seconds': credit,
        'charged_elapsed_before_resume_seconds': charged,
        'remaining_budget_seconds': plan['aggregate_budget_seconds']-charged,
        'aggregate_budget_seconds': plan['aggregate_budget_seconds'],
        'affected_arm': summary['arm'], 'affected_seed': summary['seed'], 'resume_update': summary['updates'],
        'exclude_from_hardware_speed_comparisons': True, 'archives': archives,
        'future_interruptions': 'Original conservative policy remains; no automatic further credit'}
    write_json(destination, amendment)
    return amendment


def apply(root):
    directory = root/ARCHIVE
    receipt_path = directory/'application.json'
    if receipt_path.exists():
        raise ValueError('Standby credit has already been applied')
    amendment = read_json(directory/'amendment.json')
    for path, expected in [('reports/secondary/plan.json', amendment['plan_sha256']),
                           ('scripts/continue_secondary_after_standby.py', amendment['script_sha256']),
                           (f'{ARCHIVE}/power-evidence.json', amendment['evidence_sha256']),
                           ('docs/secondary-runtime-amendment.md', amendment['document_sha256'])]:
        if sha256(root/path) != expected:
            raise ValueError(f'Amendment input changed: {path}')
    plan = read_json(root/'reports/secondary/plan.json')
    runner = verify_frozen(root, plan)
    pairs = [(item['source'], item['archive']) for item in amendment['archives']]
    if len(pairs) != len(archive_pairs()) or set(pairs) != set(archive_pairs()):
        raise ValueError('Amendment requires the exact complete reviewed archive paths')
    for item in amendment['archives']:
        if (type(item['bytes']) is not int or item['bytes'] < 0
                or any((root/item[k]).stat().st_size != item['bytes'] or sha256(root/item[k]) != item['sha256']
                       for k in ('source', 'archive'))):
            raise ValueError(f'Original or archive changed before one-time application: {item["source"]}')
    original = read_json(directory/'original-execution.json')
    originals = read_json(directory/'original-results.json')
    partials = [s for s in originals if not s['complete']]
    if (original['status'] != 'budget_stopped' or original['in_progress'] is not None
            or original['plan_sha256'] != amendment['plan_sha256'] or len(partials) != 1
            or (partials[0]['arm'], partials[0]['seed'], partials[0]['updates']) != ('human_atp8_base', 29, 268)):
        raise ValueError('Archived interruption state does not match this amendment')
    partial = partials[0]
    job = next(j for j in plan['jobs'] if (j['arm'], j['seed']) == (partial['arm'], partial['seed']))
    runner.verify_summary(partial, job, plan)
    runner.verify_cohort(root, job)
    archived_run = read_json(root/'runs/secondary-interruption/20260925-human-atp8-base-seed29/run.json')
    credit = credit_from_events(read_json(directory/'power-evidence.json'), utc(archived_run['created_utc']), utc(original['updated_utc']))
    recomputed = {'original_elapsed_wall_seconds': original['elapsed_seconds'],
        'excluded_hardware_deep_idle_seconds': credit,
        'charged_elapsed_before_resume_seconds': original['elapsed_seconds']-credit,
        'remaining_budget_seconds': plan['aggregate_budget_seconds']-original['elapsed_seconds']+credit,
        'aggregate_budget_seconds': plan['aggregate_budget_seconds']}
    if any(type(amendment[k]) not in (int, float) or not math.isfinite(amendment[k])
           or abs(amendment[k]-value) > 1e-8 for k, value in recomputed.items()):
        raise ValueError('Amendment accounting differs from independently recomputed original evidence')
    if (amendment['kind'] != 'standby_budget_amendment'
            or (amendment['affected_arm'], amendment['affected_seed'], amendment['resume_update']) !=
               (partial['arm'], partial['seed'], partial['updates'])
            or not 0 < credit < partial['elapsed_seconds']
            or not 0 < recomputed['charged_elapsed_before_resume_seconds'] < plan['aggregate_budget_seconds']):
        raise ValueError('Amendment target or conservative budget bounds changed')
    accounting = {k: amendment[k] for k in ('kind', 'original_elapsed_wall_seconds',
        'excluded_hardware_deep_idle_seconds', 'charged_elapsed_before_resume_seconds', 'affected_arm', 'affected_seed')}
    accounting.update(amendment_path=f'{ARCHIVE}/amendment.json', amendment_sha256=sha256(directory/'amendment.json'))
    updated = {**original, 'updated_utc': now(), 'elapsed_seconds': amendment['charged_elapsed_before_resume_seconds'],
               'recovery_accounting': accounting}
    temporary = root/'reports/secondary/execution.partial.json'
    write_json(temporary, updated)
    temporary.replace(root/'reports/secondary/execution.json')
    write_json(receipt_path, {'applied_utc': now(), 'amendment_sha256': accounting['amendment_sha256'],
        'original_execution_sha256': sha256(directory/'original-execution.json'),
        'adjusted_execution_sha256': sha256(root/'reports/secondary/execution.json')})
    return accounting


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'apply'))
    parser.add_argument('--root', default=str(ROOT))
    args = parser.parse_args()
    result = prepare(Path(args.root).resolve()) if args.action == 'prepare' else apply(Path(args.root).resolve())
    print(result)
