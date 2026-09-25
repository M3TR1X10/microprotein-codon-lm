"""Offline acquisition checks; all synthetic sequences are software fixtures."""
import hashlib
from copy import deepcopy

import openpyxl
import pytest
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from microprotein_lm.secondary_acquire import (
    fetch_records, match_metadata, merge_records, preserve_identical_report, resolve_impi,
    validate_secondary_record, write_canonical_json, write_frozen_bytes,
)


def fixture_record(tax_id=10090, family='ATP8', dna='ATGTGATAA', protein='MW', table=2):
    ref = {'tax_id': tax_id, 'family': family, 'uniprot': 'FIXTURE',
           'translation_table': table, 'reference_protein': protein}
    rec = SeqRecord(Seq(dna), id='FIXTURE.1')
    # Coordinates intentionally refer to a parent record, not the already sliced CDS.
    rec.features = [SeqFeature(FeatureLocation(100, 100 + len(dna)), type='source',
                    qualifiers={'db_xref': [f'taxon:{tax_id}']}),
                    SeqFeature(FeatureLocation(100, 100 + len(dna)), type='CDS',
                    qualifiers={'protein_id': ['FIXTURE.1'], 'codon_start': ['1'],
                                'transl_table': [str(table)], 'translation': [protein]})]
    return rec, ref


def test_own_species_reference_genetic_code_and_no_second_coordinate_extraction():
    rec, ref = fixture_record()
    result, reason = validate_secondary_record(rec, ref)
    assert reason is None and result['protein'] == 'MW' and result['tax_id'] == 10090
    rec.features[0].qualifiers['db_xref'] = ['taxon:9606']
    assert validate_secondary_record(rec, ref)[1] == 'wrong_taxon'
    rec, ref = fixture_record()
    rec.features[1].qualifiers['transl_table'] = ['1']
    assert validate_secondary_record(rec, ref)[1] == 'unexpected_genetic_code'
    rec, ref = fixture_record(family='ATP5ME', dna='ATGTGGTAA', table=1)
    del rec.features[1].qualifiers['transl_table']
    assert validate_secondary_record(rec, ref)[1] is None


@pytest.mark.parametrize('dna,reason', [('ATGNNNTAA', 'ambiguous_bases'),
    ('ATGTGATA', 'incomplete_codon'), ('ATGTAATAA', 'invalid_complete_cds'),
    ('ATGTGATGG', 'invalid_complete_cds'), ('TTTTGATAA', 'invalid_complete_cds')])
def test_incomplete_ambiguous_stop_start_rejected(dna, reason):
    assert validate_secondary_record(*fixture_record(dna=dna))[1] == reason


def test_partial_translation_and_index_checksum_rejected():
    rec, ref = fixture_record()
    rec.features[1].location = FeatureLocation(BeforePosition(100), 109)
    assert validate_secondary_record(rec, ref)[1] == 'partial_location'
    rec, ref = fixture_record()
    rec.features[1].qualifiers['translation'] = ['ML']
    assert validate_secondary_record(rec, ref)[1] == 'translation_mismatch'
    rec, ref = fixture_record()
    metadata = {'protein_id': rec.id, 'base_count': '9', 'sequence_md5': 'incorrect'}
    assert validate_secondary_record(rec, ref, metadata)[1] == 'sequence_hash_mismatch'
    ref['reference_protein'] = 'MA'
    assert validate_secondary_record(rec, ref)[1] == 'reference_identity_below_threshold'
    ref['reference_protein'] = 'MWA'
    assert validate_secondary_record(rec, ref)[1] == 'reference_length_mismatch'


def test_global_duplicates_keep_taxa_and_human_primary_regardless_of_order():
    result, _ = validate_secondary_record(*fixture_record())
    result['provenance'] = [{'accession': 'MOUSE.1'}]
    human = deepcopy(result)
    human['tax_id'] = 9606
    human['provenance'] = [{'accession': 'HUMAN.1'}]
    merged = merge_records([result, human], [9606, 10090])
    assert len(merged) == 1
    assert merged[0]['tax_id'] == 9606 and merged[0]['taxa'] == [9606, 10090]
    assert len(merged[0]['provenance']) == 2
    assert merge_records(merged, [9606, 10090]) == merged
    human['family'] = 'ATP5ME'
    with pytest.raises(ValueError, match='cross-family'):
        merge_records([result, human], [9606, 10090])


def test_impi_compound_symbols_resolve_without_substring(tmp_path):
    path = tmp_path / 'fixture.xlsx'
    book = openpyxl.Workbook()
    sheet = book.active
    sheet.title = 'IMPI-2021Q4pre'
    sheet.append(['Ensembl', 'Symbol', 'Name', 'IMPI Class'])
    sheet.append(['ENSG_FIXTURE', 'ATP5ME; ATP5I', 'subunit e; old name', 'Verified mitochondrial'])
    book.save(path)
    result = resolve_impi(path, {'families': {'ATP5ME': {'human_gene': 'ATP5ME'}}})
    assert result['ATP5ME']['Ensembl'] == 'ENSG_FIXTURE'
    with pytest.raises(ValueError, match='Expected one'):
        resolve_impi(path, {'families': {'ATP5M': {'human_gene': 'ATP5M'}}})


def test_exact_metadata_mapping_does_not_confuse_subunits():
    refs = [{'genes': ['ATP5ME', 'ATP5I'], 'products': ['ATP synthase membrane subunit e']}]
    row = {'gene': '', 'gene_synonym': '', 'product': 'ATP synthase membrane subunit e', 'protein_id': 'TEST.1'}
    assert match_metadata(row, refs) is refs[0]
    row['product'] = 'ATP synthase membrane subunit epsilon'
    assert match_metadata(row, refs) is None


def test_partially_missing_browser_batch_is_resolved_individually(tmp_path):
    rec, _ = fixture_record()
    rec.annotations['molecule_type'] = 'DNA'
    path = tmp_path / 'fixture.embl'
    SeqIO.write([rec], path, 'embl')

    class FixtureDownload:
        def get(self, name, url):
            if url.endswith('/MISSING.1'):
                response = requests.Response()
                response.status_code = 404
                raise requests.HTTPError('missing fixture', response=response)
            return path

    missing = []
    observed = list(fetch_records(FixtureDownload(), ['FIXTURE.1', 'MISSING.1'], missing))
    assert [r.id for r, _ in observed] == ['FIXTURE.1']
    assert [r['accession'] for r in missing] == ['MISSING.1']


def test_network_failure_is_not_a_biological_zero():
    class FixtureDownload:
        def get(self, name, url):
            response = requests.Response()
            response.status_code = 503
            raise requests.HTTPError('temporary service failure', response=response)

    with pytest.raises(requests.HTTPError):
        list(fetch_records(FixtureDownload(), ['FIXTURE.1'], []))


def test_canonical_newlines_report_identity_and_prewrite_freeze(tmp_path):
    report = tmp_path / 'reports/secondary/acquisition.json'
    original = {'eligible_records_sha256': hashlib.sha256(b'original\n').hexdigest(), 'created_utc': 'original'}
    write_canonical_json(report, original)
    old_bytes = report.read_bytes()
    assert b'\r' not in old_bytes
    regenerated = {**original, 'created_utc': 'later'}
    preserve_identical_report(report, regenerated, ['eligible_records_sha256'])
    assert report.read_bytes() == old_bytes
    (report.parent / 'plan.json').write_bytes(b'{}\n')
    corpus = tmp_path / 'data/cohort.jsonl'
    write_frozen_bytes(tmp_path, corpus, b'original\n', report)
    with pytest.raises(ValueError, match='Frozen acquisition'):
        write_frozen_bytes(tmp_path, corpus, b'changed\n', report)
    assert corpus.read_bytes() == b'original\n'
    with pytest.raises(ValueError, match='after the training plan'):
        preserve_identical_report(report, {'eligible_records_sha256': 'different'}, ['eligible_records_sha256'])
    assert report.read_bytes() == old_bytes
