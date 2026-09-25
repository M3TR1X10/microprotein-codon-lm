"""Synthetic reporting values only; no research training or outcome generation."""
from copy import deepcopy
import importlib.util
from pathlib import Path

import pytest

from microprotein_lm.io import sha256
from microprotein_lm.secondary_io import write_json


@pytest.fixture
def report():
    path = Path(__file__).resolve().parents[1]/'scripts/build_secondary_report.py'
    spec = importlib.util.spec_from_file_location('secondary_report_software_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fixture_summary(arm, bits, macro=None, variable=.4, seed=17):
    # Deliberately artificial SOFTWARE TEST values, never saved under real runs.
    stratum = {'training_bits_per_base': bits, 'scored_training_bases': 36}
    variable_stratum = {'training_bits_per_base': variable,
                        'scored_training_bases': 9 if variable is not None else 0}
    metric = {'training_bits_per_base': bits, 'family_macro_bits_per_base': bits if macro is None else macro,
              'scored_training_bases': 36, 'per_family': {'FIXTURE': deepcopy(stratum)},
              'per_taxon': {'9606': deepcopy(stratum)}, 'per_code': {'2': deepcopy(stratum)},
              'variable_positions': {**variable_stratum, 'per_family': {'FIXTURE': deepcopy(variable_stratum)}}}
    return {'arm': arm, 'role': 'matched', 'seed': seed, 'cohort_sha256': 'SOFTWARE_COHORT',
            'sampling_trace_sha256': 'SOFTWARE_DRAWS', 'updates': 4, 'complete': True,
            'sequences_seen': 8, 'bases_seen': 72, 'parameters': 10, 'sequences': 4,
            'sequence_exposures_per_unique_cds': 2, 'elapsed_seconds': 1.,
            'per_family_exposures': {'FIXTURE': {'sequences': 8, 'target_bases': 72}},
            'initial': {**deepcopy(metric), 'training_bits_per_base': 2., 'update': 0},
            'final': {**deepcopy(metric), 'update': 4}}


def fixture_groups(report):
    return {arm: [] for arm in report.PRIMARY+report.CONTROLS}


def fixture_plan(report):
    return {'aggregate_budget_seconds': 100,
            'jobs': [{'arm': arm, 'seed': seed} for arm in report.PRIMARY+report.CONTROLS
                     for seed in (17, 29, 43)]}


def test_paired_directions_use_correct_metric_and_only_common_complete_seeds(report):
    groups = fixture_groups(report)
    for seed, codon in [(17, .2), (29, .6)]:
        groups['human_atp8_base'].append(fixture_summary('human_atp8_base', .5, seed=seed))
        groups['human_atp8_codon'].append(fixture_summary('human_atp8_codon', codon, seed=seed))
    groups['human_atp8_base'].append(fixture_summary('human_atp8_base', .7, seed=43))
    groups['vertebrate_complex_codon'] = [fixture_summary('vertebrate_complex_codon', .3, macro=.2)]
    groups['vertebrate_family_macro'] = [fixture_summary('vertebrate_family_macro', .4, macro=.1)]
    contrasts = report.paired_contrasts(groups, fixture_plan(report))
    token, weighted = contrasts[0], contrasts[2]
    assert token['planned'] == 3 and [p['seed'] for p in token['pairs']] == [17, 29]
    assert [p['overall_delta'] for p in token['pairs']] == pytest.approx([-.3, .1])
    assert weighted['metric'] == 'macro_delta' and weighted['direction'] == -1
    assert weighted['pairs'][0]['overall_delta'] > 0 and weighted['pairs'][0]['macro_delta'] < 0
    assert report.expectation_outcome(-.1, -1, complete=False).startswith('Provisional: met')
    assert report.expectation_outcome(.1, -1) == 'Reversed'
    assert report.expectation_outcome(0, 1) == 'Unchanged'
    assert report.expectation_outcome(None, 1) == 'Pending'


@pytest.mark.parametrize('field', ['cohort_sha256', 'sampling_trace_sha256', 'updates', 'sequences_seen', 'bases_seen'])
def test_paired_scalar_exposure_changes_are_rejected(report, field):
    ref = fixture_summary('human_atp8_base', .5)
    changed = deepcopy(ref)
    changed[field] = 'CHANGED' if isinstance(changed[field], str) else changed[field]+1
    with pytest.raises(ValueError, match='Paired exposure mismatch'):
        report.verify_paired_exposure(changed, ref, 'SOFTWARE FIXTURE')


@pytest.mark.parametrize('category', ['per_family', 'per_taxon', 'per_code', 'variable_positions', 'family_exposure'])
def test_paired_support_and_family_exposure_changes_are_rejected(report, category):
    ref = fixture_summary('vertebrate_complex_codon', .5)
    changed = deepcopy(ref)
    if category == 'family_exposure':
        changed['per_family_exposures']['FIXTURE']['target_bases'] += 3
    elif category == 'variable_positions':
        changed['final'][category]['scored_training_bases'] += 3
    else:
        first = next(iter(changed['final'][category]))
        changed['final'][category][first]['scored_training_bases'] += 3
    with pytest.raises(ValueError, match='Paired'):
        report.verify_paired_exposure(changed, ref, 'SOFTWARE FIXTURE')


def test_zero_variable_support_stays_missing(report):
    groups = fixture_groups(report)
    groups['human_atp8_base'] = [fixture_summary('human_atp8_base', .5, variable=None)]
    groups['human_atp8_codon'] = [fixture_summary('human_atp8_codon', .3, variable=None)]
    pair = report.paired_contrasts(groups, fixture_plan(report))[0]['pairs'][0]
    assert pair['variable_delta'] is None and pair['variable_target_bases'] == 0


def test_execution_identity_and_status_are_checked(report):
    good = {'plan_sha256': 'FROZEN', 'elapsed_seconds': 100., 'status': 'timed_out'}
    report.validate_execution(good, 'FROZEN')
    with pytest.raises(ValueError, match='different frozen plan'):
        report.validate_execution(good, 'OTHER')
    with pytest.raises(ValueError, match='invalid elapsed time'):
        report.validate_execution({**good, 'elapsed_seconds': -1}, 'FROZEN')
    with pytest.raises(ValueError, match='Unknown scheduler'):
        report.validate_execution({**good, 'status': 'invented'}, 'FROZEN')


def test_report_smoke_exposes_partial_pairs_support_and_execution_failure(report, tmp_path):
    root = tmp_path/'SOFTWARE_REPORT_FIXTURE_ONLY'
    plan = fixture_plan(report)
    write_json(root/'reports/secondary/plan.json', plan)
    write_json(root/'configs/secondary/biology.json', {'taxa': {'9606': 'Fixture organism'}})
    write_json(root/'reports/secondary/execution.json', {
        'plan_sha256': sha256(root/'reports/secondary/plan.json'), 'elapsed_seconds': 100.,
        'status': 'timed_out', 'complete': False,
        'recovery_accounting': {'extra_wall_seconds_charged': 30}})
    cohorts = {'views': {'human_atp8': {'n_sequences': 4, 'family_counts': {'FIXTURE': 4},
                  'taxon_counts': {'9606': 4}, 'target_nt': 36}},
               'overlap_distinct_rna': {'human_atp8__human_complex': 3,
                  'human_complex__mammal_complex': 2, 'mammal_complex__vertebrate_complex': 3}}
    diversity = {'human_atp8': {'n_unique_peptides': 2, 'variable_target_bases': 9,
                               'observed_target_codon_classes': 4}}
    runs = [{'summary': fixture_summary('human_atp8_base', .5)},
            {'summary': fixture_summary('human_atp8_codon', .2)}]
    report.build_markdown(root, plan, cohorts, runs, {}, diversity)
    path = root/'reports/secondary/report.md'
    text = path.read_text(encoding='utf-8')
    assert '## Controlled contrasts and prespecified expectations' in text
    assert 'One-factor contrasts' not in text
    assert 'Tokenization: codon minus base' in text and 'Provisional: met' in text
    assert '## Per-taxon training fit' in text and 'Fixture organism (9606)' in text
    assert '## Per-genetic-code training fit' in text
    assert '1/4 CDS (25.00%)' in text and '2/4 (50.00%)' in text
    assert 'Scheduler status: **timed_out**' in text and '0.50 additional minutes' in text
    assert 'Variable support is comparable within each pair only' in text
    assert b'\r\n' not in path.read_bytes()
