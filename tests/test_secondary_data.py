"""Synthetic selection invariants, not biological validation examples."""
import hashlib
from collections import Counter

from microprotein_lm.secondary_data import select_views


def test_matched_views_are_deduplicated_balanced_and_preserve_quotas():
    rows = []
    for family, n in [('ATP8', 90), ('small_a', 9), ('small_b', 2)]:
        for tax in [9606, 10090, 10116, 9913, 9031, 7955]:
            if family != 'ATP8' and tax in [9031, 7955]:
                continue
            for i in range(n):
                rows.append({'family': family, 'tax_id': tax,
                             'sequence_sha256': hashlib.sha256(f'{family}:{tax}:{i}'.encode()).hexdigest()})
    config = {'selection_seed': 12, 'mammals': [9606, 10090, 10116, 9913],
              'nonmammalian_vertebrates': [9031, 7955]}
    views, quotas = select_views(rows, config)
    reordered, _ = select_views(list(reversed(rows)), config)
    assert views == reordered
    for name in ['human_complex', 'mammal_complex', 'vertebrate_complex']:
        assert len(views[name]) == 90
        assert len({r['sequence_sha256'] for r in views[name]}) == 90
        assert Counter(r['family'] for r in views[name]) == quotas
    mm = {r['sequence_sha256'] for r in views['mammal_complex'] if r['tax_id'] == 9606}
    vv = {r['sequence_sha256'] for r in views['vertebrate_complex'] if r['tax_id'] == 9606}
    assert mm == vv
    assert {r['tax_id'] for r in views['mammal_complex']} == {9606, 10090, 10116, 9913}
    assert {r['tax_id'] for r in views['vertebrate_complex']} == {9606, 10090, 10116, 9913, 9031, 7955}


def test_prepared_views_refuse_changed_selection_or_provenance(tmp_path):
    import json
    import pytest
    from microprotein_lm.secondary_data import prepare_secondary

    config = {'selection_seed': 12, 'mammals': [9606, 10090, 10116],
              'nonmammalian_vertebrates': [9031]}
    for name, content in {
        'configs/secondary/biology.json': json.dumps(config),
        'docs/secondary-protocol.md': 'Software fixture protocol',
        'data/raw/secondary/genomes/enriched-records.jsonl': 'Software fixture pool',
        'reports/secondary/acquisition.json': '{}',
        'reports/secondary/genome-acquisition.json': '{}',
    }.items():
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    rows = []
    for tax in [9606, 10090, 10116, 9031]:
        for family, count in [('ATP8', 12), ('small', 2)]:
            for i in range(count):
                identifier = f'{family}:{tax}:{i}'
                rows.append({'family': family, 'tax_id': tax, 'rna': 'AUGUUUUAA',
                             'protein': 'MF', 'translation_table': 2,
                             'provenance': [{'source': 'synthetic accession fixture'}],
                             'sequence_sha256': hashlib.sha256(identifier.encode()).hexdigest()})
    # Deliberately synthetic identifiers test selection freezing only; these
    # fixtures never enter the training engine or biological quality pipeline.
    first = prepare_secondary(tmp_path, rows)
    assert prepare_secondary(tmp_path, rows) == first
    frozen = tmp_path / 'data/processed/secondary/human_atp8/cohort.jsonl'
    original = frozen.read_bytes()
    assert b'provenance' not in original
    assert first['views']['human_atp8']['provenance_join_key'] == 'sequence_sha256'
    (tmp_path / 'reports/secondary/genome-acquisition.json').write_text('{"changed":true}')
    with pytest.raises(ValueError, match='Frozen selection metadata differs'):
        prepare_secondary(tmp_path, rows)
    assert frozen.read_bytes() == original
