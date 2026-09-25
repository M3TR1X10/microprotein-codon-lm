"""Offline synthetic standby-accounting fixtures; never read research checkpoints."""
from copy import deepcopy
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

from microprotein_lm.io import read_json, sha256
from microprotein_lm.secondary_io import write_json


def fixture_evidence():
    return {'provider': 'Microsoft-Windows-Kernel-Power', 'channel': 'System',
        'scenario_field': 'ScenarioInstanceIdV2', 'events': [
            {'record_id': 52915, 'event_id': 506, 'utc': '2026-09-25T03:11:53.9776625Z', 'scenario': 22, 'boot_id': 98},
            {'record_id': 52920, 'event_id': 507, 'utc': '2026-09-25T07:09:54.5582076Z', 'scenario': 22, 'boot_id': 98,
             'sleep_entered': True, 'duration_us': 14280582985, 'software_drips_us': 14005354549, 'hardware_drips_us': 13390759541},
            {'record_id': 52922, 'event_id': 506, 'utc': '2026-09-25T07:09:54.5593153Z', 'scenario': 24, 'boot_id': 98},
            {'record_id': 52931, 'event_id': 507, 'utc': '2026-09-25T07:52:17.8827133Z', 'scenario': 24, 'boot_id': 98,
             'sleep_entered': True, 'duration_us': 2543323574, 'software_drips_us': 2436965574, 'hardware_drips_us': 2329413931}]}


@pytest.fixture
def prepared_inputs(tmp_path, monkeypatch):
    script_path = Path(__file__).resolve().parents[1]/'scripts/continue_secondary_after_standby.py'
    spec = importlib.util.spec_from_file_location('standby_software_fixture', script_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    root = tmp_path/'SOFTWARE_STANDBY_FIXTURES_ONLY'
    (root/'scripts').mkdir(parents=True)
    (root/'scripts/continue_secondary_after_standby.py').write_bytes(script_path.read_bytes())
    monkeypatch.setattr(module, '__file__', str(root/'scripts/continue_secondary_after_standby.py'))
    # Frozen runner identity/cohort helpers have their own tests. Isolate this
    # amendment's archive and accounting logic from any actual training artifacts.
    runner = SimpleNamespace(verify_summary=lambda *args: None, verify_cohort=lambda *args: None)
    monkeypatch.setattr(module, 'verify_frozen', lambda *args: runner)
    jobs = [{'arm': f'SOFTWARE_COMPLETE_{i}', 'seed': 17, 'config': {'updates': 2000}} for i in range(8)]
    jobs += [{'arm': 'human_atp8_base', 'seed': 29, 'config': {'updates': 2000}}]
    plan = {'scope': 'training_only', 'aggregate_budget_seconds': 14400, 'shutdown_reserve_seconds': 120,
            'jobs': jobs, 'engine_source_hashes': {}, 'input_hashes': {}}
    write_json(root/'reports/secondary/plan.json', plan)
    execution = {'status': 'budget_stopped', 'in_progress': None, 'elapsed_seconds': 20466.25,
        'plan_sha256': sha256(root/'reports/secondary/plan.json'), 'complete': False,
        'completed_runs': 8, 'planned_runs': 25, 'recovery_accounting': None,
        'updated_utc': '2026-09-25T07:52:27.636731+00:00'}
    write_json(root/'reports/secondary/execution.json', execution)
    summary = {'arm': 'human_atp8_base', 'seed': 29, 'updates': 268, 'complete': False,
               'elapsed_seconds': 17016.391, 'stop_reason': 'max_seconds'}
    results = [{**job, 'complete': True, 'updates': 2000} for job in jobs[:8]] + [summary]
    write_json(root/'reports/secondary/results.json', results)
    run = root/'runs/secondary/human_atp8_base-seed29'
    run.mkdir(parents=True)
    (run/'checkpoint.pt').write_bytes(b'SOFTWARE FIXTURE CHECKPOINT; NOT RESEARCH DATA')
    write_json(run/'summary.json', summary)
    write_json(run/'metrics.json', [{'update': 0, 'fixture_only': True}, {'update': 268, 'fixture_only': True}])
    write_json(run/'baselines.json', {'software_fixture_only': True})
    write_json(run/'run.json', {'created_utc': '2026-09-25T03:08:47.253861+00:00', 'fixture_only': True})
    write_json(root/module.ARCHIVE/'power-evidence.json', fixture_evidence())
    (root/'docs').mkdir()
    (root/'docs/secondary-runtime-amendment.md').write_text('SOFTWARE FIXTURE DOCUMENT ONLY\n')
    return root, module


def test_evidence_credits_only_reviewed_nonoverlapping_hardware_idle(prepared_inputs):
    _, module = prepared_inputs
    credit = module.credit_from_events(fixture_evidence(), module.utc('2026-09-25T03:08:47Z'),
                                       module.utc('2026-09-25T07:52:27Z'))
    assert credit == pytest.approx(15720.173472)
    assert 20466.25-credit == pytest.approx(4746.076528)
    assert 14400-(20466.25-credit) == pytest.approx(9653.923472)


@pytest.mark.parametrize('mutation', ['record', 'overlap', 'outside_worker', 'scenario', 'hardware', 'duration'])
def test_invalid_event_pairing_or_residency_is_rejected(prepared_inputs, mutation):
    _, module = prepared_inputs
    evidence = fixture_evidence()
    start, stop = module.utc('2026-09-25T03:08:47Z'), module.utc('2026-09-25T07:52:27Z')
    if mutation == 'record':
        evidence['events'][0]['record_id'] = 999
    elif mutation == 'overlap':
        evidence['events'][2]['utc'] = '2026-09-25T07:09:50Z'
    elif mutation == 'outside_worker':
        start = module.utc('2026-09-25T04:00:00Z')
    elif mutation == 'scenario':
        evidence['events'][1]['scenario'] = 999
    elif mutation == 'hardware':
        evidence['events'][1]['hardware_drips_us'] = evidence['events'][1]['duration_us']+1
    else:
        evidence['events'][1]['duration_us'] += 1000000
    with pytest.raises(ValueError):
        module.credit_from_events(evidence, start, stop)


def test_prepare_archives_exact_inputs_then_apply_is_one_time(prepared_inputs):
    root, module = prepared_inputs
    amendment = module.prepare(root)
    assert len(amendment['archives']) == 7
    for item in amendment['archives']:
        assert sha256(root/item['source']) == sha256(root/item['archive']) == item['sha256']
        assert (root/item['archive']).stat().st_size == item['bytes']
    old_checkpoint = (root/'runs/secondary/human_atp8_base-seed29/checkpoint.pt').read_bytes()
    with pytest.raises(ValueError):
        module.prepare(root)
    accounting = module.apply(root)
    assert accounting['charged_elapsed_before_resume_seconds'] == pytest.approx(4746.076528)
    updated = read_json(root/'reports/secondary/execution.json')
    assert updated['elapsed_seconds'] == pytest.approx(4746.076528)
    assert updated['recovery_accounting']['amendment_sha256'] == sha256(root/module.ARCHIVE/'amendment.json')
    assert read_json(root/module.ARCHIVE/'original-execution.json')['elapsed_seconds'] == 20466.25
    assert (root/'runs/secondary/human_atp8_base-seed29/checkpoint.pt').read_bytes() == old_checkpoint
    before = sha256(root/'reports/secondary/execution.json')
    with pytest.raises(ValueError):
        module.apply(root)
    assert sha256(root/'reports/secondary/execution.json') == before


@pytest.mark.parametrize('target', ['source', 'archive', 'evidence', 'script', 'document'])
def test_apply_rejects_changed_archived_or_identity_inputs(prepared_inputs, target):
    root, module = prepared_inputs
    amendment = module.prepare(root)
    paths = {'source': amendment['archives'][0]['source'], 'archive': amendment['archives'][0]['archive'],
             'evidence': f'{module.ARCHIVE}/power-evidence.json',
             'script': 'scripts/continue_secondary_after_standby.py', 'document': 'docs/secondary-runtime-amendment.md'}
    (root/paths[target]).write_bytes(b'CHANGED SOFTWARE FIXTURE')
    with pytest.raises((ValueError, KeyError)):
        module.apply(root)
    assert not (root/module.ARCHIVE/'application.json').exists()


@pytest.mark.parametrize('field', ['original_elapsed_wall_seconds', 'excluded_hardware_deep_idle_seconds',
    'charged_elapsed_before_resume_seconds', 'remaining_budget_seconds', 'aggregate_budget_seconds', 'resume_update'])
def test_apply_recomputes_amendment_math_and_incident_identity(prepared_inputs, field):
    root, module = prepared_inputs
    module.prepare(root)
    path = root/module.ARCHIVE/'amendment.json'
    amendment = read_json(path)
    amendment[field] += 1
    write_json(path, amendment)
    before = sha256(root/'reports/secondary/execution.json')
    with pytest.raises(ValueError):
        module.apply(root)
    assert sha256(root/'reports/secondary/execution.json') == before


@pytest.mark.parametrize('mutation', ['omitted', 'duplicate', 'traversal', 'absolute', 'wrong_size'])
def test_apply_requires_exact_safe_archive_set(prepared_inputs, mutation):
    root, module = prepared_inputs
    module.prepare(root)
    path = root/module.ARCHIVE/'amendment.json'
    amendment = read_json(path)
    if mutation == 'omitted':
        amendment['archives'].pop(0)
    elif mutation == 'duplicate':
        amendment['archives'].append(deepcopy(amendment['archives'][0]))
    elif mutation == 'wrong_size':
        amendment['archives'][0]['bytes'] += 1
    else:
        original = root/amendment['archives'][0]['source']
        outside = root.parent/'outside-software-fixture.json'
        outside.write_bytes(original.read_bytes())
        amendment['archives'][0]['source'] = str(outside) if mutation == 'absolute' else '../outside-software-fixture.json'
    write_json(path, amendment)
    before = sha256(root/'reports/secondary/execution.json')
    with pytest.raises(ValueError):
        module.apply(root)
    assert sha256(root/'reports/secondary/execution.json') == before


def test_crash_after_adjustment_before_receipt_cannot_double_credit(prepared_inputs, monkeypatch):
    root, module = prepared_inputs
    module.prepare(root)
    original_write = module.write_json
    def receipt_failure(path, value):
        if Path(path).name == 'application.json':
            raise OSError('SOFTWARE FIXTURE simulated receipt-write interruption')
        original_write(path, value)
    monkeypatch.setattr(module, 'write_json', receipt_failure)
    with pytest.raises(OSError):
        module.apply(root)
    adjusted = sha256(root/'reports/secondary/execution.json')
    monkeypatch.setattr(module, 'write_json', original_write)
    with pytest.raises(ValueError):
        module.apply(root)
    assert sha256(root/'reports/secondary/execution.json') == adjusted
