import hashlib
import json
import pytest

from microprotein_lm.gap_acquire import (novelty_descriptors, outer_query,
    split_untouched, stable_accessions, unreviewed_targets, guard_frozen_manifest)


def row(rna, source='NEW.1', parent='PNEW.1', family='ATP8', protein='MI'):
    return {'rna': rna, 'sequence_sha256': hashlib.sha256(rna.encode()).hexdigest(),
        'family': family, 'tax_id': 9606, 'protein': protein,
        'provenance': [{'accession': source, 'parent_accession': parent}]}


def test_union_exclusion_is_rna_and_stable_source_based():
    train = [row('AUGAUUUAA', source='OLD.1', parent='POLD.1')]
    candidates = [row('AUGAUUUAA'), row('AUGAUCUAA', source='OLD.2'),
                  row('AUGAUCUAG', parent='POLD.2'), row('AUGAUCUAA')]
    kept, excluded = split_untouched(candidates, train)
    assert len(kept) == 1
    assert [r['reason'] for r in excluded] == ['training_rna_duplicate',
        'training_source_accession_overlap', 'training_source_accession_overlap']


def test_parent_multireference_overlap_and_novelty_are_explicit():
    candidate = row('AUGAUCUAA', parent='P1.2,P2.4')
    assert stable_accessions(candidate) == {'NEW', 'P1', 'P2'}
    annotated = novelty_descriptors([candidate], [row('AUGAUUUAA')])[0]
    assert annotated['novelty']['minimum_nt_edit_distance'] == 1
    assert annotated['novelty']['exact_training_peptide']
    assert annotated['novelty']['nearest_nt_normalized_edit_distance'] == 1 / 9
    assert annotated['novelty']['same_family_taxon_length_minimum_nt_hamming'] == 1
    assert annotated['novelty']['same_family_taxon_length_maximum_nt_identity'] == 8 / 9


def test_leads_do_not_create_new_references():
    biology = {'taxa': {'9606': 'human'}, 'families': {'ATP8': {'interpro': 'IPR1'}}}
    entry = {'primaryAccession': 'U1', 'entryType': 'UniProtKB unreviewed (TrEMBL)',
        'organism': {'taxonId': 9606}, 'uniProtKBCrossReferences': [
            {'database': 'InterPro', 'id': 'IPR1'}, {'database': 'EMBL', 'id': 'P1',
                'properties': [{'key': 'ProteinId', 'value': 'C1.2'}]}]}
    targets, rejected = unreviewed_targets([entry], biology, [])
    assert not targets and rejected[0]['reason'] == 'no_frozen_own_species_reviewed_reference'
    refs = [{'tax_id': 9606, 'family': 'ATP8', 'uniprot': 'R1'}]
    targets, rejected = unreviewed_targets([entry], biology, refs)
    assert not rejected and targets['C1.2'][0]['reference']['uniprot'] == 'R1'


def test_outer_query_does_not_reuse_training_parent_interval():
    query = outer_query(9606, {'outer_parent_min_nt': 150, 'outer_parent_max_nt': 100000})
    assert 'tax_eq(9606)' in query and '(base_count<14000 OR base_count>20000)' in query
    assert 'base_count>=150' in query and 'base_count<=100000' in query


def test_frozen_guard_rejects_changed_cohort_or_source_identity(tmp_path):
    path = tmp_path / 'manifest.json'
    identity = {'cohort_sha256': 'original', 'training_pool_sha256': 'blocked-union'}
    path.write_text(json.dumps({'identity': identity}))
    guard_frozen_manifest(path, identity)
    for key in identity:
        with pytest.raises(ValueError, match='frozen test cohort'):
            guard_frozen_manifest(path, {**identity, key: 'different'})
