"""Tiny serialized SOFTWARE FIXTURES; no model fitting or research diagnostics."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest
import torch

from microprotein_lm.io import read_json, sha256
from microprotein_lm.secondary_io import write_json
from microprotein_lm.secondary_train import (
    SecondaryCorpus, _configuration, build_secondary_model, engine_source_hashes,
)


@pytest.fixture
def auditor(monkeypatch):
    path = Path(__file__).resolve().parents[1]/'scripts/audit_secondary_checkpoints.py'
    spec = importlib.util.spec_from_file_location('secondary_checkpoint_software_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Pytest shares a process with prior tests. The actual CLI uses a fresh process.
    monkeypatch.setattr(module.torch, 'set_num_interop_threads', lambda threads: None)
    return module


@pytest.fixture
def checkpoint_fixture(tmp_path):
    root = tmp_path/'SOFTWARE_CHECKPOINT_FIXTURE_ONLY'
    source_root = Path(__file__).resolve().parents[1]
    (root/'scripts').mkdir(parents=True)
    (root/'scripts/run_secondary.py').write_bytes((source_root/'scripts/run_secondary.py').read_bytes())
    directory = root/'data/processed/secondary/SOFTWARE_FIXTURE'
    directory.mkdir(parents=True)
    records = [{'rna': rna, 'sequence_sha256': hashlib.sha256(rna.encode()).hexdigest(),
                'family': 'SOFTWARE_FIXTURE', 'tax_id': 9606, 'translation_table': 1}
               for rna in ('AUGAAAUAA', 'AUGAACUAA')]
    (directory/'cohort.jsonl').write_bytes((''.join(json.dumps(r)+'\n' for r in records)).encode())
    write_json(directory/'manifest.json', {'cohort_sha256': sha256(directory/'cohort.jsonl'), 'stage': 'training_only'})
    config = _configuration({'stage': 'training_only', 'dataset': 'SOFTWARE_FIXTURE',
        'model': {'d_model': 8, 'n_heads': 2, 'n_layers': 1, 'dropout': 0}, 'max_nt_context': 12,
        'batch_size': 2, 'gradient_accumulation': 1, 'updates': 4, 'diagnostic_interval': 2,
        'checkpoint_interval': 2, 'cpu_threads': 1, 'gradient_clip': 1., 'warmup_updates': 1,
        'learning_rate': .001, 'min_learning_rate': .0001, 'position_encoding': False})
    corpus = SecondaryCorpus(directory, 'codon')
    model = build_secondary_model(config, corpus)
    job = {'arm': 'SOFTWARE_ARM', 'seed': 17, 'dataset': config['dataset'], 'config': config,
           'mode': 'codon', 'role': 'matched', 'cohort_sha256': corpus.identity['cohort_sha256'],
           'manifest_sha256': corpus.identity['manifest_sha256']}
    plan = {'scope': 'training_only', 'jobs': [job], 'engine_source_hashes': engine_source_hashes(),
            'input_hashes': {'scripts/run_secondary.py': sha256(root/'scripts/run_secondary.py')}}
    # This is a constructed serialization contract, deliberately not trained weights.
    history = [{'update': 0, 'fixture_only': True}, {'update': 4, 'fixture_only': True}]
    summary = {'config': config, 'dataset': job['dataset'], 'mode': 'codon', 'seed': 17,
        'complete': True, 'updates': 4, 'cohort_sha256': job['cohort_sha256'],
        'data_identity': corpus.identity, 'source_hashes': plan['engine_source_hashes'],
        'parameters': sum(p.numel() for p in model.parameters()), 'initial': history[0], 'final': history[-1],
        'sequences_seen': 8, 'bases_seen': 48, 'sampling_trace_sha256': 'SOFTWARE_FIXTURE_TRACE',
        'per_family_exposures': {'SOFTWARE_FIXTURE': {'sequences': 8, 'target_bases': 48}}}
    payload = {'identity': {'config': config, 'mode': 'codon', 'seed': 17,
                 'data_identity': corpus.identity, 'source_hashes': plan['engine_source_hashes']},
        'update': 4, 'model': model.state_dict(), 'history': history, 'optimizer': {'state': {}},
        **{k: summary[k] for k in ('sequences_seen', 'bases_seen', 'sampling_trace_sha256', 'per_family_exposures')}}
    output = root/'runs/secondary/SOFTWARE_ARM-seed17'
    output.mkdir(parents=True)
    for name in ('final.pt', 'checkpoint.pt'):
        torch.save(payload, output/name)
    write_json(output/'summary.json', summary)
    write_json(root/'reports/secondary/plan.json', plan)
    write_json(root/'reports/secondary/results.json', [{'arm': job['arm'], 'role': job['role'], **summary}])
    return root, output, payload


def test_checkpoint_audit_hashes_loads_cpu_weights_only_and_preserves_files(checkpoint_fixture, auditor, monkeypatch):
    root, output, _ = checkpoint_fixture
    before = {name: sha256(output/name) for name in ('final.pt', 'checkpoint.pt')}
    original_load, calls = torch.load, []
    def checked_load(*args, **kwargs):
        calls.append(kwargs)
        return original_load(*args, **kwargs)
    monkeypatch.setattr(auditor.torch, 'load', checked_load)
    result = auditor.audit(root)
    assert result['status'] == 'passed' and result['verified_run_count'] == 1
    assert result['all_planned_runs_verified']
    assert len(calls) == 2 and all(c['weights_only'] is True and c['map_location'] == 'cpu' for c in calls)
    for name, metadata in result['verified_runs'][0]['files'].items():
        assert metadata['sha256'] == before[name] == sha256(output/name)
        assert metadata['bytes'] == (output/name).stat().st_size and metadata['strict_model_load']
    assert b'\r\n' not in (root/'reports/secondary/checkpoints.json').read_bytes()


@pytest.mark.parametrize('damage', ['nonfinite', 'identity', 'exposure', 'shape', 'position', 'unreadable', 'missing'])
def test_checkpoint_failures_are_recorded_without_silently_excluding_complete_run(checkpoint_fixture, auditor, damage):
    root, output, payload = checkpoint_fixture
    changed = deepcopy(payload)
    if damage == 'nonfinite':
        changed['optimizer']['state']['software_tensor'] = torch.tensor([float('nan')])
    elif damage == 'identity':
        changed['identity']['seed'] = 999
    elif damage == 'exposure':
        changed['bases_seen'] += 3
    elif damage == 'shape':
        changed['model']['norm.weight'] = torch.zeros(99)
    elif damage == 'position':
        changed['model']['position'].fill_(1)
    if damage == 'unreadable':
        (output/'final.pt').write_bytes(b'NOT A PYTORCH CHECKPOINT; SOFTWARE FIXTURE')
    elif damage == 'missing':
        (output/'final.pt').unlink()
    else:
        torch.save(changed, output/'final.pt')
    result = auditor.audit(root)
    assert result['status'] == 'failed' and result['failed_run_count'] == 1
    assert result['verified_run_count'] == 0 and result['excluded_run_count'] == 0
    assert result['failed_runs'][0]['error']


def test_missing_and_incomplete_expected_runs_are_explicit_exclusions(checkpoint_fixture, auditor):
    root, _, _ = checkpoint_fixture
    plan_path = root/'reports/secondary/plan.json'
    plan = read_json(plan_path)
    first = plan['jobs'][0]
    plan['jobs'] += [{**deepcopy(first), 'arm': 'UNREPORTED_SOFTWARE_ARM'},
                     {**deepcopy(first), 'arm': 'PARTIAL_SOFTWARE_ARM'}]
    write_json(plan_path, plan)
    results_path = root/'reports/secondary/results.json'
    results = read_json(results_path)
    results.append({**deepcopy(results[0]), 'arm': 'PARTIAL_SOFTWARE_ARM', 'complete': False,
                    'updates': 2, 'stop_reason': 'max_seconds'})
    write_json(results_path, results)
    result = auditor.audit(root)
    assert result['status'] == 'passed' and result['verified_run_count'] == 1
    assert not result['all_planned_runs_verified']
    assert {r['reason'] for r in result['excluded_runs']} == {'incomplete', 'no_reported_result'}


def test_global_plan_input_failure_is_recorded(checkpoint_fixture, auditor):
    root, _, _ = checkpoint_fixture
    (root/'scripts/run_secondary.py').write_text('CHANGED SOFTWARE FIXTURE')
    result = auditor.audit(root)
    assert result['status'] == 'failed' and result['global_error']


def test_finite_but_different_rolling_model_is_rejected(checkpoint_fixture, auditor):
    root, output, payload = checkpoint_fixture
    changed = deepcopy(payload)
    changed['model']['norm.weight'].add_(.25)
    torch.save(changed, output/'checkpoint.pt')
    result = auditor.audit(root)
    assert result['status'] == 'failed'
    assert 'model states disagree' in result['failed_runs'][0]['error']


def test_reaudit_rejects_finite_consistent_modification_and_preserves_established_inventory(checkpoint_fixture, auditor):
    root, output, payload = checkpoint_fixture
    assert auditor.audit(root)['status'] == 'passed'
    inventory_path = root/'reports/secondary/checkpoints.json'
    established = inventory_path.read_bytes()
    changed = deepcopy(payload)
    changed['model']['norm.weight'].add_(.125)
    # Both files agree and are finite/loadable; the earlier hashes still catch drift.
    for name in ('final.pt', 'checkpoint.pt'):
        torch.save(changed, output/name)
    attempted = auditor.audit(root)
    assert attempted['status'] == 'failed' and attempted['established_inventory_preserved']
    assert 'Established checkpoint identity drift' in attempted['failed_runs'][0]['error']
    assert inventory_path.read_bytes() == established
    assert attempted['output_file'] != 'reports/secondary/checkpoints.json'
    assert read_json(root/attempted['output_file'])['status'] == 'failed'


@pytest.mark.parametrize('field', ['sha256', 'bytes', 'model_state_sha256'])
def test_each_established_fingerprint_field_is_enforced(auditor, field):
    files = {name: {'file': name, 'sha256': 'OLD', 'bytes': 100, 'model_state_sha256': 'STATE'}
             for name in ('final.pt', 'checkpoint.pt')}
    changed = deepcopy(files)
    changed['final.pt'][field] = 200 if field == 'bytes' else 'CHANGED'
    with pytest.raises(ValueError, match=field):
        auditor.verify_established_files(changed, {'files': files})


def test_reaudit_allows_newly_completed_previously_excluded_job(checkpoint_fixture, auditor):
    root, output, _ = checkpoint_fixture
    plan_path, results_path = root/'reports/secondary/plan.json', root/'reports/secondary/results.json'
    plan = read_json(plan_path)
    plan['jobs'].append({**deepcopy(plan['jobs'][0]), 'arm': 'LATER_SOFTWARE_ARM'})
    write_json(plan_path, plan)
    first = auditor.audit(root)
    assert first['verified_run_count'] == 1 and first['excluded_run_count'] == 1
    established_files = first['verified_runs'][0]['files']
    later_output = output.parent/'LATER_SOFTWARE_ARM-seed17'
    later_output.mkdir()
    for name in ('final.pt', 'checkpoint.pt', 'summary.json'):
        (later_output/name).write_bytes((output/name).read_bytes())
    results = read_json(results_path)
    results.append({**deepcopy(results[0]), 'arm': 'LATER_SOFTWARE_ARM'})
    write_json(results_path, results)
    second = auditor.audit(root)
    assert second['status'] == 'passed' and second['all_planned_runs_verified']
    assert second['verified_run_count'] == 2 and second['previous_verified_run_count'] == 1
    assert second['verified_runs'][0]['files'] == established_files


def test_reaudit_rejects_different_plan_and_preserves_established_inventory(checkpoint_fixture, auditor):
    root, _, _ = checkpoint_fixture
    auditor.audit(root)
    destination = root/'reports/secondary/checkpoints.json'
    established = destination.read_bytes()
    plan_path = root/'reports/secondary/plan.json'
    plan = read_json(plan_path)
    plan['software_fixture_changed_plan'] = True
    write_json(plan_path, plan)
    attempted = auditor.audit(root)
    assert attempted['status'] == 'failed'
    assert 'different frozen plan' in attempted['global_error']['error']
    assert destination.read_bytes() == established


def test_reaudit_rejects_disappearing_previously_verified_result(checkpoint_fixture, auditor):
    root, _, _ = checkpoint_fixture
    auditor.audit(root)
    destination = root/'reports/secondary/checkpoints.json'
    established = destination.read_bytes()
    write_json(root/'reports/secondary/results.json', [])
    attempted = auditor.audit(root)
    assert attempted['status'] == 'failed' and attempted['failed_run_count'] == 1
    assert 'no longer present' in attempted['failed_runs'][0]['error']
    assert destination.read_bytes() == established
