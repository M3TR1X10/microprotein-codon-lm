"""Mocked controller checks; fixture files never enter research training."""
from datetime import datetime, timedelta, timezone
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

from microprotein_lm.io import read_json, sha256, write_json
import microprotein_lm.secondary_train as engine


@pytest.fixture
def runner():
    path = Path(__file__).resolve().parents[1]/'scripts/run_secondary.py'
    spec = importlib.util.spec_from_file_location('secondary_runner_software_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def frozen(tmp_path, monkeypatch, runner):
    root = tmp_path/'SOFTWARE_FIXTURES_ONLY'
    config = {'stage': 'training_only', 'model': {'d_model': 8, 'n_heads': 2, 'n_layers': 1},
              'max_nt_context': 12, 'batch_size': 2, 'gradient_accumulation': 1, 'updates': 4,
              'diagnostic_interval': 2, 'checkpoint_interval': 2, 'cpu_threads': 1,
              'gradient_clip': 1., 'warmup_updates': 1, 'learning_rate': .001,
              'min_learning_rate': .0001, 'seeds': [17, 29, 43]}
    write_json(root/'configs/secondary/training.json', config)
    inputs = ['configs/secondary/biology.json', 'docs/secondary-protocol.md',
              'reports/secondary/cohorts.json', 'reports/secondary/hardware-benchmark.json',
              'reports/secondary/acquisition.json', 'reports/secondary/genome-acquisition.json',
              'src/microprotein_lm/secondary_acquire.py', 'src/microprotein_lm/secondary_data.py',
              'src/microprotein_lm/secondary_genomes.py', 'scripts/build_secondary_data.py',
              'scripts/run_secondary.py']
    for name in inputs:
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('SOFTWARE TEST FIXTURE ONLY\n', encoding='utf-8')
    for dataset in ['human_atp8', 'human_complex', 'mammal_complex', 'vertebrate_complex', 'full_vertebrate']:
        directory = root/f'data/processed/secondary/{dataset}'
        directory.mkdir(parents=True)
        (directory/'cohort.jsonl').write_text('{"software_fixture_only":true}\n')
        write_json(directory/'manifest.json', {'cohort_sha256': sha256(directory/'cohort.jsonl')})
    monkeypatch.setattr(engine, 'engine_source_hashes', lambda: {'software_fixture_engine': 'frozen'})
    plan = runner.make_plan(root)
    # A single mock job makes budget/error tests independent of training costs.
    plan['jobs'] = plan['jobs'][:1]
    plan['aggregate_budget_seconds'] = 100
    plan['shutdown_reserve_seconds'] = 10
    write_json(root/'reports/secondary/plan.json', plan)
    clock = [1000.]
    monkeypatch.setattr(runner, 'time', SimpleNamespace(monotonic=lambda: clock[0]))
    return root, plan, clock


def output_for(root, job):
    return root/'runs/secondary'/f"{job['arm']}-seed{job['seed']}"


def valid_summary(job, plan):
    return {'config': job['config'], 'dataset': job['dataset'], 'mode': job['mode'],
            'seed': job['seed'], 'cohort_sha256': job['cohort_sha256'],
            'source_hashes': plan['engine_source_hashes'],
            'data_identity': {k: job[k] for k in ('cohort_sha256', 'manifest_sha256')},
            'updates': job['config']['updates'], 'complete': True}


def test_frozen_manifest_and_runner_hashes_exist_and_are_checked(frozen, runner):
    root, plan, _ = frozen
    assert plan['jobs'][0]['manifest_sha256']
    assert 'scripts/run_secondary.py' in plan['input_hashes']
    runner.verify_inputs(root, plan)
    (root/'scripts/run_secondary.py').write_text('CHANGED SOFTWARE FIXTURE')
    with pytest.raises(ValueError, match='Plan input changed'):
        runner.verify_inputs(root, plan)


@pytest.mark.parametrize('entry', ['controller', 'worker'])
def test_both_entry_points_reject_cohort_and_manifest_refrozen_after_plan(frozen, runner, monkeypatch, entry):
    root, plan, _ = frozen
    job = plan['jobs'][0]
    directory = root/f'data/processed/secondary/{job["dataset"]}'
    (directory/'cohort.jsonl').write_text('{"software_fixture_only":"changed"}\n')
    # A self-consistent new manifest must still fail the earlier frozen identity.
    write_json(directory/'manifest.json', {'cohort_sha256': sha256(directory/'cohort.jsonl')})
    monkeypatch.setattr(runner.subprocess, 'run', lambda *a, **k: pytest.fail('Unexpected worker'))
    monkeypatch.setattr(runner, 'lower_priority', lambda: pytest.fail('Worker reached training setup'))
    with pytest.raises(ValueError, match='Frozen cohort or manifest changed'):
        if entry == 'controller':
            runner.execute(root, plan)
        else:
            runner.worker(SimpleNamespace(worker=0, remaining_seconds=30), root, plan)


def test_manifest_only_change_is_rejected(frozen, runner):
    root, plan, _ = frozen
    job = plan['jobs'][0]
    manifest_path = root/f'data/processed/secondary/{job["dataset"]}/manifest.json'
    manifest = read_json(manifest_path)
    manifest['changed_metadata'] = True
    write_json(manifest_path, manifest)
    with pytest.raises(ValueError, match='Frozen cohort or manifest changed'):
        runner.verify_cohort(root, job)


def test_success_persists_in_progress_before_launch_and_enforces_timeout(frozen, runner, monkeypatch):
    root, plan, clock = frozen
    job = plan['jobs'][0]
    def fake_run(command, cwd, timeout):
        active = read_json(root/'reports/secondary/execution.json')
        assert active['in_progress']['job_index'] == 0
        assert active['in_progress']['started_utc']
        assert active['status'] == 'running'
        assert active['plan_sha256'] == sha256(root/'reports/secondary/plan.json')
        assert float(command[-1]) == 90 and timeout == 100
        clock[0] += 5
        write_json(output_for(root, job)/'summary.json', valid_summary(job, plan))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(runner.subprocess, 'run', fake_run)
    runner.execute(root, plan)
    state = read_json(root/'reports/secondary/execution.json')
    assert state['complete'] and state['elapsed_seconds'] == 5
    assert b'\r\n' not in (root/'reports/secondary/execution.json').read_bytes()
    assert state['in_progress'] is None
    assert len(read_json(root/'reports/secondary/results.json')) == 1


def test_outer_timeout_records_failure_and_preserves_checkpoint(frozen, runner, monkeypatch):
    root, plan, clock = frozen
    checkpoint = output_for(root, plan['jobs'][0])/'checkpoint.pt'
    checkpoint.parent.mkdir(parents=True)
    checkpoint.write_bytes(b'SOFTWARE FIXTURE CHECKPOINT')
    def timeout_run(command, cwd, timeout):
        clock[0] += timeout
        raise runner.subprocess.TimeoutExpired(command, timeout)
    monkeypatch.setattr(runner.subprocess, 'run', timeout_run)
    with pytest.raises(RuntimeError, match='Aggregate deadline'):
        runner.execute(root, plan)
    state = read_json(root/'reports/secondary/execution.json')
    assert state['status'] == 'timed_out' and not state['complete']
    assert state['elapsed_seconds'] == 100 and state['in_progress'] is None
    assert checkpoint.read_bytes() == b'SOFTWARE FIXTURE CHECKPOINT'


@pytest.mark.parametrize('field,value', [('seed', 99), ('cohort_sha256', 'different'),
                                       ('source_hashes', {'changed': 'source'}), ('updates', 1)])
def test_completed_summary_must_match_frozen_job(frozen, runner, monkeypatch, field, value):
    root, plan, _ = frozen
    job = plan['jobs'][0]
    summary = valid_summary(job, plan)
    summary[field] = value
    write_json(output_for(root, job)/'summary.json', summary)
    monkeypatch.setattr(runner.subprocess, 'run', lambda *a, **k: pytest.fail('Unexpected worker'))
    with pytest.raises(ValueError):
        runner.execute(root, plan)


def test_existing_complete_run_is_verified_and_skipped(frozen, runner, monkeypatch):
    root, plan, _ = frozen
    job = plan['jobs'][0]
    write_json(output_for(root, job)/'summary.json', valid_summary(job, plan))
    monkeypatch.setattr(runner.subprocess, 'run', lambda *a, **k: pytest.fail('Unexpected worker'))
    runner.execute(root, plan)
    assert read_json(root/'reports/secondary/execution.json')['complete']


def test_execution_state_cannot_be_reused_for_different_plan(frozen, runner, monkeypatch):
    root, plan, _ = frozen
    write_json(root/'reports/secondary/execution.json', {'plan_sha256': 'another-plan', 'elapsed_seconds': 4})
    monkeypatch.setattr(runner.subprocess, 'run', lambda *a, **k: pytest.fail('Unexpected worker'))
    with pytest.raises(ValueError, match='different frozen plan'):
        runner.execute(root, plan)


def test_interrupted_worker_time_is_conservatively_charged_on_recovery(frozen, runner, monkeypatch):
    root, plan, _ = frozen
    write_json(root/'reports/secondary/execution.json', {
        'plan_sha256': sha256(root/'reports/secondary/plan.json'), 'elapsed_seconds': 10,
        'in_progress': {'job_index': 0, 'started_utc': (datetime.now(timezone.utc)-timedelta(seconds=120)).isoformat()}})
    monkeypatch.setattr(runner.subprocess, 'run', lambda *a, **k: pytest.fail('Budget must not be replenished'))
    runner.execute(root, plan)
    state = read_json(root/'reports/secondary/execution.json')
    assert state['elapsed_seconds'] >= 130 and not state['complete']
    assert state['recovery_accounting']['extra_wall_seconds_charged'] >= 120
    assert state['status'] == 'budget_stopped'


def test_worker_failure_records_time_without_erasing_prior_checkpoint(frozen, runner, monkeypatch):
    root, plan, clock = frozen
    def failed_run(*args, **kwargs):
        clock[0] += 7
        return SimpleNamespace(returncode=3)
    monkeypatch.setattr(runner.subprocess, 'run', failed_run)
    with pytest.raises(RuntimeError, match='Worker failed'):
        runner.execute(root, plan)
    state = read_json(root/'reports/secondary/execution.json')
    assert state['status'] == 'worker_failed' and state['worker_returncode'] == 3
    assert state['elapsed_seconds'] == 7 and state['in_progress'] is None
