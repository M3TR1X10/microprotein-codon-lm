"""Software-only synthetic fixtures. Nothing here is research training data."""
import hashlib
import json
import math
from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pytest
import torch
import microprotein_lm.secondary_train as secondary_engine

from microprotein_lm.io import read_json, sha256, write_json
from microprotein_lm.secondary_train import (
    SecondaryCorpus, build_secondary_model, fit_training_baselines,
    objective_loss, train_secondary, training_diagnostics,
)


def write_fixture(root, records):
    root.mkdir(exist_ok=True)
    path = root/'cohort.jsonl'
    path.write_text(''.join(json.dumps(r)+'\n' for r in records), encoding='utf-8')
    write_json(root/'manifest.json', {'cohort_sha256': sha256(path), 'stage': 'training_only',
                                    'validation_sequences': 0, 'test_sequences': 0})


@pytest.fixture
def cohort(tmp_path):
    # Family A has three CDS and two taxa, B has one longer CDS. Family B must
    # not manufacture apparent within-family variability from pooled positions.
    records = []
    for rna, family, tax, code in [
        ('AUGAAAUAA', 'FIXTURE_A', 9606, 1),
        ('AUGAACUAA', 'FIXTURE_A', 9606, 1),
        ('AUGAAAUAG', 'FIXTURE_A', 10090, 2),
        ('AUGCCCGGGUAA', 'FIXTURE_B', 10090, 2),
    ]:
        records.append({'rna': rna, 'family': family, 'tax_id': tax, 'translation_table': code,
                        'sequence_sha256': hashlib.sha256(rna.encode()).hexdigest()})
    root = tmp_path/'fixture_cohort'
    write_fixture(root, records)
    return root


def tiny_config():
    return {'stage': 'training_only', 'dataset': 'SOFTWARE_FIXTURE_ONLY',
            'model': {'d_model': 16, 'n_heads': 2, 'n_layers': 1, 'dropout': 0.1},
            'max_nt_context': 12, 'batch_size': 2, 'gradient_accumulation': 2,
            'updates': 4, 'warmup_updates': 1, 'learning_rate': 0.001,
            'min_learning_rate': 0.0001, 'weight_decay': 0.01, 'gradient_clip': 1.,
            'diagnostic_interval': 2, 'checkpoint_interval': 2, 'device': 'cpu',
            'precision': 'float32', 'deterministic': True, 'cpu_threads': 1,
            'loss': 'token', 'position_encoding': True}


def test_matching_targets_variable_lengths_and_exact_vocabularies(cohort):
    base, codon = SecondaryCorpus(cohort, 'base'), SecondaryCorpus(cohort, 'codon')
    assert len(base.tokenizer.vocabulary) == 4
    assert len(codon.tokenizer.vocabulary) == 64
    xb, yb = base.batch([0, 1, 2, 3])
    xc, yc = codon.batch([0, 1, 2, 3])
    assert (yb != -100).sum() == 3*(yc != -100).sum()
    assert (yb[:, :2] == -100).all()
    assert (yb[0, 8:] == -100).all()
    assert (yc[0, 2:] == -100).all()
    for row in range(4):
        assert base.tokenizer.decode(yb[row][yb[row] != -100].tolist()) == codon.tokenizer.decode(yc[row][yc[row] != -100].tolist())
    assert base.identity == codon.identity


def test_variable_codons_are_within_family_and_taxon_and_match_base_region(cohort):
    base, codon = SecondaryCorpus(cohort, 'base'), SecondaryCorpus(cohort, 'codon')
    assert [int(mask.sum()) for mask in codon.variable_masks] == [1, 1, 0, 0]
    assert [int(mask.sum()) for mask in base.variable_masks] == [3, 3, 0, 0]
    assert np.flatnonzero(base.variable_masks[0]).tolist() == [2, 3, 4]
    # Different stop codon in the singleton other taxon does not create variation.
    assert not codon.variable_masks[2].any()


def test_family_macro_has_fixed_normalizer_and_accumulation_equivalence(cohort):
    corpus = SecondaryCorpus(cohort, 'codon')
    _, targets = corpus.batch([0, 1, 2, 3])
    torch.manual_seed(10)
    logits = torch.randn((*targets.shape, 64), requires_grad=True)
    full = objective_loss(logits, targets, [0, 1, 2, 3], corpus, 'family_macro', 4, 9)
    first = objective_loss(logits[:1], targets[:1], [0], corpus, 'family_macro', 4, 9)
    rest = objective_loss(logits[1:], targets[1:], [1, 2, 3], corpus, 'family_macro', 4, 9)
    assert float(full.detach()) == pytest.approx(float((first+rest).detach()))
    per_token = torch.nn.functional.cross_entropy(logits.reshape(-1, 64), targets.reshape(-1),
                                                  ignore_index=-100, reduction='none').reshape_as(targets)
    expected = (per_token[:3].sum()/6 + per_token[3:].sum()/3)/2
    assert torch.allclose(full, expected)
    # Uniform predictions have log(vocabulary) loss under either weighting.
    uniform = torch.zeros_like(logits)
    for objective in ('token', 'family_macro'):
        loss = objective_loss(uniform, targets, [0, 1, 2, 3], corpus, objective, 4, 9)
        assert float(loss) == pytest.approx(math.log(64))
    # Whole-corpus averaged gradients equal summed accumulation contributions.
    grad_full = torch.autograd.grad(full, logits, retain_graph=True)[0]
    grad_parts = torch.autograd.grad(first+rest, logits)[0]
    assert torch.allclose(grad_full, grad_parts)


@pytest.mark.parametrize('position_encoding', [True, False])
def test_position_ablation_keeps_shape_and_future_causality(cohort, position_encoding):
    corpus = SecondaryCorpus(cohort, 'base')
    config = tiny_config()
    config['position_encoding'] = position_encoding
    config['model']['dropout'] = 0
    torch.manual_seed(1)
    model = build_secondary_model(config, corpus).eval()
    assert bool(model.position.count_nonzero()) == position_encoding
    x, _ = corpus.batch([3])
    changed = x.clone()
    changed[:, 6:] = (changed[:, 6:]+1) % 4
    a, b = model(x), model(changed)
    assert a.shape == (*x.shape, 4)
    assert torch.allclose(a[:, :6], b[:, :6], atol=1e-7)


def test_diagnostics_and_baselines_are_training_only_and_cached(cohort):
    corpus = SecondaryCorpus(cohort, 'codon')
    config = tiny_config()
    model = build_secondary_model(config, corpus)
    diagnostics = training_diagnostics(model, corpus, 2, 'cpu')
    assert diagnostics['scored_training_tokens'] == 9
    assert diagnostics['scored_training_bases'] == 27
    assert set(diagnostics['per_family']) == {'FIXTURE_A', 'FIXTURE_B'}
    assert set(diagnostics['per_code']) == {'1', '2'}
    assert set(diagnostics['per_taxon']) == {'9606', '10090'}
    assert diagnostics['variable_positions']['scored_training_tokens'] == 2
    assert diagnostics['variable_positions']['per_family']['FIXTURE_B']['training_bits_per_base'] is None
    assert diagnostics['family_macro_bits_per_base'] == pytest.approx(
        sum(v['training_bits_per_base'] for v in diagnostics['per_family'].values())/2)
    path = cohort/'baselines-codon.json'
    first = fit_training_baselines(corpus)
    assert b'\r\n' not in path.read_bytes()  # Portable tracked/report hashes on Windows and Linux.
    stamp = path.stat().st_mtime_ns
    second = fit_training_baselines(corpus)
    assert first == second and path.stat().st_mtime_ns == stamp
    assert first['prior_total_mass'] == 2
    assert first['prior_per_category'] == 2/64
    assert first['target_bases'] == 27
    assert 'family label' in first['oracle_warning']
    empirical = first['metrics']['empirical']
    assert empirical['family_position']['training_bits_per_base'] <= empirical['position']['training_bits_per_base']


@pytest.mark.parametrize('objective,position', [('token', True), ('family_macro', False)])
def test_resume_exact_sampling_identity_and_retention(cohort, tmp_path, objective, position):
    config = tiny_config()
    config.update(loss=objective, position_encoding=position)
    full, part = tmp_path/'full', tmp_path/'part'
    summary = train_secondary(config, 'codon', 7, cohort, full)
    partial = train_secondary(config, 'codon', 7, cohort, part, stop_after=2)
    assert partial['complete'] is False and not (part/'final.pt').exists()
    resumed = train_secondary(config, 'codon', 7, cohort, part, resume=True)
    assert resumed['complete'] and (part/'final.pt').exists()
    a = torch.load(full/'checkpoint.pt', weights_only=True)
    b = torch.load(part/'checkpoint.pt', weights_only=True)
    assert summary['bases_seen'] == resumed['bases_seen']
    assert summary['sampling_trace_sha256'] == resumed['sampling_trace_sha256']
    assert sum(v['target_bases'] for v in summary['per_family_exposures'].values()) == summary['bases_seen']
    for key in a['model']:
        assert torch.equal(a['model'][key], b['model'][key]), key
    assert summary['final'] == resumed['final']
    assert set(p.name for p in part.glob('*.pt')) == {'checkpoint.pt', 'final.pt'}
    with pytest.raises(FileExistsError):
        train_secondary(config, 'codon', 7, cohort, full)


def test_same_sampling_across_objective_capacity_and_vocabulary(cohort, tmp_path):
    config = tiny_config()
    a = train_secondary(config, 'codon', 11, cohort, tmp_path/'codon')
    modified = deepcopy(config)
    modified['model']['d_model'] = 8
    modified['loss'] = 'family_macro'
    modified['position_encoding'] = False
    b = train_secondary(modified, 'base', 11, cohort, tmp_path/'base')
    assert a['sampling_trace_sha256'] == b['sampling_trace_sha256']
    assert a['sequences_seen'] == b['sequences_seen']
    assert a['bases_seen'] == b['bases_seen']
    assert a['per_family_exposures'] == b['per_family_exposures']


@pytest.mark.parametrize('field,new_value', [('family', 'CHANGED_FAMILY'), ('tax_id', 10116), ('translation_table', 3)])
def test_resume_rejects_changed_biological_metadata_even_with_updated_hash(cohort, tmp_path, field, new_value):
    config = tiny_config()
    output = tmp_path/'partial'
    train_secondary(config, 'codon', 3, cohort, output, stop_after=2)
    records = [json.loads(line) for line in (cohort/'cohort.jsonl').read_text().splitlines()]
    records[0][field] = new_value
    write_fixture(cohort, records)
    with pytest.raises(ValueError, match='identical'):
        train_secondary(config, 'codon', 3, cohort, output, resume=True)


def test_resume_rejects_manifest_only_change_and_cache_refits(cohort, tmp_path):
    config = tiny_config()
    output = tmp_path/'partial'
    train_secondary(config, 'codon', 3, cohort, output, stop_after=2)
    old = read_json(cohort/'baselines-codon.json')['cache_identity']
    manifest = read_json(cohort/'manifest.json')
    manifest['provenance_note'] = 'Changed after freezing'
    write_json(cohort/'manifest.json', manifest)
    changed = fit_training_baselines(SecondaryCorpus(cohort, 'codon'))
    assert changed['cache_identity'] != old
    with pytest.raises(ValueError, match='identical'):
        train_secondary(config, 'codon', 3, cohort, output, resume=True)


def test_reader_rejects_checksum_duplicates_and_nontraining_scope(cohort):
    records = [json.loads(line) for line in (cohort/'cohort.jsonl').read_text().splitlines()]
    records[0]['rna'] = 'AUGGGGUAA'
    write_fixture(cohort, records)
    with pytest.raises(ValueError, match='Sequence checksum'):
        SecondaryCorpus(cohort, 'base')
    records[0]['sequence_sha256'] = hashlib.sha256(records[0]['rna'].encode()).hexdigest()
    write_fixture(cohort, records+[records[0]])
    with pytest.raises(ValueError, match='distinct'):
        SecondaryCorpus(cohort, 'codon')
    manifest = read_json(cohort/'manifest.json')
    manifest['test_sequences'] = 1
    write_json(cohort/'manifest.json', manifest)
    with pytest.raises(ValueError, match='Validation and test'):
        SecondaryCorpus(cohort, 'codon')


def test_runtime_budget_saves_orderly_partial_checkpoint_and_resumes(cohort, tmp_path):
    config = tiny_config()
    full, guarded = tmp_path/'full', tmp_path/'guarded'
    train_secondary(config, 'codon', 23, cohort, full)
    partial = train_secondary(config, 'codon', 23, cohort, guarded, max_seconds=1e-12)
    assert partial['complete'] is False and partial['stop_reason'] == 'max_seconds'
    assert partial['updates'] == 1 and partial['final']['update'] == 1
    assert (guarded/'checkpoint.pt').exists() and not (guarded/'final.pt').exists()
    resumed = train_secondary(config, 'codon', 23, cohort, guarded, resume=True)
    assert resumed['complete'] and resumed['stop_reason'] == 'completed'
    a = torch.load(full/'final.pt', weights_only=True)
    b = torch.load(guarded/'final.pt', weights_only=True)
    for key in a['model']:
        assert torch.equal(a['model'][key], b['model'][key]), key


def test_progress_clock_does_not_add_diagnostics_or_change_training(cohort, tmp_path, monkeypatch, capsys):
    config = tiny_config()
    ordinary, visible = tmp_path/'ordinary', tmp_path/'visible'
    train_secondary(config, 'codon', 31, cohort, ordinary)
    capsys.readouterr()
    ticks = [0.]
    def accelerated_clock():
        ticks[0] += 61
        return ticks[0]
    monkeypatch.setattr(secondary_engine, 'time', SimpleNamespace(monotonic=accelerated_clock))
    train_secondary(config, 'codon', 31, cohort, visible)
    output = capsys.readouterr().out
    assert 'PROGRESS secondary SOFTWARE_FIXTURE_ONLY codon seed=31 update=1/4 elapsed=' in output
    a = torch.load(ordinary/'final.pt', weights_only=True)
    b = torch.load(visible/'final.pt', weights_only=True)
    assert a['history'] == b['history']
    assert [entry['update'] for entry in b['history']] == [0, 2, 4]
    assert a['sampling_trace_sha256'] == b['sampling_trace_sha256']
    assert torch.equal(a['torch_rng'], b['torch_rng'])
    for key in a['model']:
        assert torch.equal(a['model'][key], b['model'][key]), key
