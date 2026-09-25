"""Synthetic software fixtures only; these sequences are never research data."""
import copy
import hashlib
import json
from pathlib import Path
import random

import pytest

from microprotein_lm.gap_design import (
    HomologyGraphAuditor, audit_homology_graph, case_regions, load_verified_holdout,
    make_cases, multi_threshold_sequence_clusters, resample_sequence_clusters,
    select_panel, sequence_clusters, sequence_identity, stable_source_accessions, verify_holdout)


def config():
    return json.loads((Path(__file__).parents[1]/'configs/gap-evaluation.json').read_text())


def record(index, family='ATP8', tax_id=9606, sense_codons=68):
    rng = random.Random(index)
    rna = 'AUG' + ''.join(rng.choice(['AAA', 'CCU', 'GGC', 'UUU', 'CAU', 'AAC'])
                        for _ in range(sense_codons-1)) + 'UAA'
    return {'rna': rna, 'sequence_sha256': hashlib.sha256(rna.encode()).hexdigest(),
            'family': family, 'tax_id': tax_id, 'taxa': [tax_id],
            'translation_table': 2 if family == 'ATP8' else 1,
            'provenance': [{'accession': f'SOFTWARE_PROTEIN_{index}.1',
                            'parent_accession': f'SOFTWARE_PARENT_{index}.1'}]}


def test_panel_selection_is_input_order_independent_balanced_and_bounded():
    records = [record(i + t*1000, tax_id=t) for t in range(1, 9) for i in range(20)]
    original = copy.deepcopy(records)
    selected = select_panel(records, config()['panel'])
    reverse = select_panel(list(reversed(records)), config()['panel'])
    assert selected == reverse and records == original
    assert len(selected) == 96
    assert {sum(r['tax_id'] == t for r in selected) for t in range(1, 9)} == {12}
    assert len(select_panel(records[:20], config()['panel'])) == 16


def test_primary_lengths_share_one_anchor_and_exact_checked_coordinates():
    r = record(1)
    design = make_cases([r], config())
    assert [c['gap_codons'] for c in design['primary']] == [1, 2, 4, 8, 16]
    assert len({c['anchor_id'] for c in design['primary']}) == 1
    assert len({c['start_codon'] for c in design['primary']}) == 1
    previous = ''
    for case in design['primary']:
        parts = case_regions(r, case)
        assert parts['truth'].startswith(previous)
        assert len(parts['truth']) == case['target_bases'] == 3*case['gap_codons']
        assert parts['offset_nt'] == 0 and not parts['right']
        assert parts['prefix'] + parts['truth'] == r['rna'][:case['gap_end_nt']]
        assert case['gap_end_nt'] <= len(r['rna'])-3
        previous = parts['truth']
        assert 'rna' not in case and 'truth' not in case


def test_focused_context_is_support_matched_and_rerank_reuses_bank_keys():
    r = record(3)
    design = make_cases([r], config())
    focus, rerank = design['focused_context'], design['right_flank']
    assert len(focus) == 5*7 and len(rerank) == 5*2*5
    assert len({c['anchor_id'] for c in focus+rerank}) == 1
    for case in focus+rerank:
        assert case['start_codon'] >= 32
        parts = case_regions(r, case)
        assert len(parts['prefix']) == 3*case['observed_left_codons']
        assert parts['offset_nt'] == 3*(case['start_codon']-case['observed_left_codons'])
        assert len(parts['right']) == 3*case['right_context_codons']
        if case['left_context_codons'] != 'full':
            assert case['observed_left_codons'] == case['left_context_codons']
    for gap in config()['primary']['gap_codons']:
        for left in (8, 'full'):
            cells = [c for c in rerank if (c['gap_codons'], c['left_context_codons']) == (gap, left)]
            assert len({c['candidate_bank_key'] for c in cells}) == 1
            assert len({c['case_id'] for c in cells}) == 5


def test_short_sequences_are_omitted_not_silently_shortened():
    short = record(5, sense_codons=10)
    middle = record(6, sense_codons=40)
    design = make_cases([short, middle], config())
    assert len(design['primary']) == 5
    assert not design['focused_context'] and not design['right_flank']
    assert len(design['support']['omitted_primary']) == 1
    assert len(design['support']['omitted_focused_context']) == 2


def test_exploratory_gaps_keep_primary_anchor_and_disclose_support():
    records = [record(i, sense_codons=100) for i in range(10, 20)]
    design = make_cases(records, config())
    primary = {c['sequence_sha256']: c for c in design['primary']}
    for case in design['exploratory']:
        original = primary[case['sequence_sha256']]
        assert case['anchor_id'] == original['anchor_id']
        assert case['start_codon'] == original['start_codon']
        assert case['gap_codons'] in (32, 48, 64)
        assert case['gap_end_nt'] <= 300
    assert design['support']['stages']['exploratory']['n_cases'] == len(design['exploratory'])


def test_empty_panel_remains_empty_without_training_substitution():
    design = make_cases([], config())
    assert design['selected_records'] == []
    assert all(design[name] == [] for name in ('primary', 'focused_context', 'right_flank', 'exploratory'))
    assert design['clusters']['membership'] == {}


def test_case_tampering_and_cross_record_use_rejected():
    r = record(25)
    case = make_cases([r], config())['focused_context'][0]
    changed = {**case, 'left_start_nt': 0}
    with pytest.raises(ValueError, match='coordinates'):
        case_regions(r, changed)
    with pytest.raises(ValueError, match='metadata'):
        case_regions(record(26), case)


def test_invalid_record_hash_and_duplicate_cds_rejected():
    r = record(30)
    with pytest.raises(ValueError, match='checksum'):
        make_cases([{**r, 'rna': r['rna'].replace('AAA', 'AAC', 1)}], config())
    with pytest.raises(ValueError, match='Distinct'):
        make_cases([r, r], config())
    with pytest.raises(ValueError, match='genetic code'):
        make_cases([{**r, 'translation_table': 1}], config())


@pytest.mark.parametrize('mutation,match', [
    ('duplicate_rna', 'RNA was present'), ('source_version', 'source/parent'),
    ('foreign_taxon', 'taxonomy'), ('foreign_family', 'family'),
    ('no_provenance', 'provenance'), ('bad_hash', 'checksum')])
def test_strict_holdout_exclusions_and_identity(mutation, match):
    training, test = record(40), record(41)
    if mutation == 'duplicate_rna':
        test['rna'], test['sequence_sha256'] = training['rna'], training['sequence_sha256']
    elif mutation == 'source_version':
        test['provenance'][0]['parent_accession'] = ' software_parent_40.9 '
    elif mutation == 'foreign_taxon':
        test['tax_id'], test['taxa'] = 999, [999]
    elif mutation == 'foreign_family':
        test['family'], test['translation_table'] = 'UNAPPROVED', 1
    elif mutation == 'no_provenance':
        test['provenance'] = []
    elif mutation == 'bad_hash':
        test['sequence_sha256'] = '0'*64
    with pytest.raises(ValueError, match=match):
        verify_holdout([test], [training], allowed_taxa=[9606], allowed_families=['ATP8'])


def test_stable_accessions_normalize_versions_and_do_not_use_study_id():
    r = record(50)
    r['provenance'] = [{'accession': ' A1.2, B2.1 ', 'parent_accession': 'c3.12',
                        'study_accession': 'SAME_PROJECT'}]
    assert stable_source_accessions(r) == {'A1', 'B2', 'C3'}


def test_holdout_audit_reports_unique_unit_and_metadata_identity():
    training, test = record(51), record(52)
    result = verify_holdout([test], [training], allowed_taxa=[9606], allowed_families=['ATP8'])
    assert result['n_sequences'] == result['training_union_sequences'] == 1
    assert len(result['ordered_sequence_metadata_sha256']) == 64
    with pytest.raises(ValueError, match='duplicate RNA'):
        verify_holdout([test, test], [training], allowed_taxa=[9606], allowed_families=['ATP8'])


def test_clusters_are_connected_components_not_only_greedy_representatives():
    # 303 nt: 3 changes are >=99% identity. The endpoints differ by 6 bases.
    base = record(61, sense_codons=100)
    base['rna'] = 'AUG' + 'AAA'*99 + 'UAA'
    rows = []
    for i, positions in enumerate(([], [3, 6, 9], [3, 6, 9, 12, 15, 18])):
        row = copy.deepcopy(base)
        seq = list(row['rna'])
        for pos in positions:
            seq[pos] = 'C'
        row['rna'] = ''.join(seq)
        row['sequence_sha256'] = hashlib.sha256(row['rna'].encode()).hexdigest()
        rows.append(row)
    result = sequence_clusters(rows)
    assert len(result['clusters']) == 1
    assert len(result['clusters'][0]['sequence_sha256s']) == 3
    assert sequence_clusters(list(reversed(rows))) == result


def test_clusters_never_merge_taxa_families_or_lengths():
    records = [record(70), record(71, tax_id=10090), record(72, family='ATP5ME'),
               record(73, sense_codons=60)]
    assert len(sequence_clusters(records, threshold=.01)['clusters']) == 4


def test_inconsistent_context_grid_rejected():
    c = config()
    c['focused_context']['minimum_available_left_codons'] = 1
    with pytest.raises(ValueError, match='every declared left'):
        make_cases([record(80)], c)


def test_file_backed_loader_rejects_unfrozen_stage_without_loading_models(tmp_path):
    target = tmp_path/'data/processed/gap'
    target.mkdir(parents=True)
    (target/'manifest.json').write_text(json.dumps({'stage': 'training_only'}))
    with pytest.raises(ValueError, match='not frozen'):
        load_verified_holdout(tmp_path)


def test_sequence_identity_and_multi_threshold_clustering():
    r1 = record(100, sense_codons=30)
    
    r2 = copy.deepcopy(r1)
    seq2 = list(r2['rna'])
    # Toggle base at pos 3 to guarantee 1 diff out of 96 nt (identity ~0.9895)
    seq2[3] = 'G' if seq2[3] != 'G' else 'C'
    r2['rna'] = ''.join(seq2)
    r2['sequence_sha256'] = hashlib.sha256(r2['rna'].encode()).hexdigest()

    r3 = copy.deepcopy(r1)
    seq3 = list(r3['rna'])
    for pos in (3, 6, 9, 12, 15, 18):
        seq3[pos] = 'G' if seq3[pos] != 'G' else 'C'
    # 6 diffs out of 96 nt (identity ~0.9375)
    r3['rna'] = ''.join(seq3)
    r3['sequence_sha256'] = hashlib.sha256(r3['rna'].encode()).hexdigest()

    records = [r1, r2, r3]
    mt = multi_threshold_sequence_clusters(records, thresholds=(0.99, 0.95, 0.90))

    assert '0.99' in mt['thresholds']
    assert '0.95' in mt['thresholds']
    assert '0.90' in mt['thresholds']

    # At 0.99 threshold: r1 and r2 are separate clusters (identity < 0.99)
    assert mt['summary']['0.99']['n_clusters'] == 3
    # At 0.95 threshold: r1 and r2 are in the same cluster (identity ~0.9895 >= 0.95)
    assert mt['summary']['0.95']['n_clusters'] == 2
    # At 0.90 threshold: all 3 are in 1 cluster
    assert mt['summary']['0.90']['n_clusters'] == 1


def test_resample_sequence_clusters_bootstrap_and_eligibility():
    records = [record(i, sense_codons=30) for i in range(110, 114)]
    res = resample_sequence_clusters(records, thresholds=(0.99, 0.95, 0.90), n_draws=100, seed=42, min_clusters_per_stratum=5)

    assert res['bootstrap_seed'] == 42
    assert res['n_draws'] == 100
    # Stratum has 4 clusters (< 5 min clusters requirement), so bands_available should be False
    assert res['resample_by_threshold']['0.99']['bands_available'] is False
    assert 'draw_summary' in res['resample_by_threshold']['0.99']
    summary = res['resample_by_threshold']['0.99']['draw_summary']
    assert summary['sequence_count_mean'] == 4.0
    assert len(summary['sequence_count_ci_95']) == 2

    # Deterministic seed reproducibility check
    res_repeat = resample_sequence_clusters(records, thresholds=(0.99,), n_draws=100, seed=42, min_clusters_per_stratum=5)
    assert res_repeat['resample_by_threshold']['0.99'] == res['resample_by_threshold']['0.99']


def test_homology_graph_auditor_detects_leakage_and_structure():
    t_rec = record(200, sense_codons=30)
    train_rec = copy.deepcopy(t_rec)  # 100% identical training sequence (leakage)

    t_rec2 = record(201, sense_codons=30)  # Distinct sequence

    report = audit_homology_graph([t_rec, t_rec2], training_records=[train_rec], thresholds=(0.99, 0.95, 0.90))

    assert report['n_test_sequences'] == 2
    assert report['n_training_sequences'] == 1
    assert len(report['leakage_warnings']) > 0
    assert report['leakage_warnings'][0]['identity'] == 1.0

    auditor = HomologyGraphAuditor(thresholds=(0.99, 0.95, 0.90))
    aud_res = auditor.audit([t_rec2], training_records=[])
    assert len(aud_res['leakage_warnings']) == 0
    assert aud_res['thresholds']['0.99']['n_leakage_edges'] == 0
    assert 'graph_audit_sha256' in aud_res

