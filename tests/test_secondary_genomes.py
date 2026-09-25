"""Parent-genome extraction fixtures never enter training acquisition."""
import hashlib

from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from microprotein_lm.secondary_genomes import extract_atp8


def parent_fixture(join=False, reverse=False):
    coding = Seq('ATGAAATAA')
    if join:
        dna = Seq('GGATGCCCAAATAATT')
        location = CompoundLocation([FeatureLocation(2, 5, strand=1), FeatureLocation(8, 14, strand=1)])
    elif reverse:
        dna = Seq('GG') + coding.reverse_complement() + Seq('TT')
        location = FeatureLocation(2, 11, strand=-1)
    else:
        dna = Seq('GG') + coding + Seq('TT')
        location = FeatureLocation(2, 11, strand=1)
    parent = SeqRecord(dna, id='PARENT.2')
    parent.features = [SeqFeature(FeatureLocation(0, len(dna)), type='source',
        qualifiers={'db_xref': ['taxon:10090']}), SeqFeature(location, type='CDS',
        qualifiers={'gene': ['Atp8'], 'protein_id': ['PROTEIN.3'], 'transl_table': ['2'],
                    'codon_start': ['1'], 'translation': ['MK']})]
    reference = {'tax_id': 10090, 'family': 'ATP8', 'reference_protein': 'MK',
                 'translation_table': 2, 'genes': ['Atp8'], 'products': ['ATP synthase 8']}
    metadata = {'accession': 'PARENT', 'sequence_version': '2', 'tax_id': '10090',
                'base_count': str(len(dna)), 'sequence_md5': hashlib.md5(str(dna).encode()).hexdigest()}
    return parent, reference, metadata


def test_parent_feature_extraction_respects_strand_join_and_parent_coordinates():
    for settings in [{}, {'reverse': True}, {'join': True}]:
        parent, ref, meta = parent_fixture(**settings)
        result, reason = extract_atp8(parent, ref, meta)
        assert reason is None and result['rna'] == 'AUGAAAUAA'
        assert result['accession'] == 'PROTEIN.3'
        assert result['original_feature_location'] == str(parent.features[1].location)
        assert result['feature_parts'][0]['strand'] == (-1 if settings.get('reverse') else 1)


def test_parent_identity_checksum_and_gene_gate():
    parent, ref, meta = parent_fixture()
    meta['sequence_version'] = '1'
    assert extract_atp8(parent, ref, meta)[1] == 'parent_accession_mismatch'
    parent, ref, meta = parent_fixture()
    meta['sequence_md5'] = 'wrong'
    assert extract_atp8(parent, ref, meta)[1] == 'parent_hash_mismatch'
    parent, ref, meta = parent_fixture()
    parent.features[1].qualifiers['gene'] = ['ND4L']
    assert extract_atp8(parent, ref, meta)[1] == 'missing_annotated_atp8'
    parent, ref, meta = parent_fixture()
    parent.seq = Seq(None, length=len(parent))
    assert extract_atp8(parent, ref, meta)[1] == 'undefined_parent_sequence'


def test_parent_route_keeps_cds_quality_gates():
    parent, ref, meta = parent_fixture()
    parent.features[1].qualifiers['transl_except'] = ['fixture exception']
    assert extract_atp8(parent, ref, meta)[1] == 'transl_except'
    parent, ref, meta = parent_fixture()
    parent.features[1].qualifiers['translation'] = ['MM']
    assert extract_atp8(parent, ref, meta)[1] == 'translation_mismatch'


def test_parent_coordinates_must_be_exact_bounded_and_unambiguous():
    parent, ref, meta = parent_fixture()
    parent.features[1].location = FeatureLocation(BeforePosition(2), 11, strand=1)
    assert extract_atp8(parent, ref, meta)[1] == 'partial_location'
    parent.features[1].location = FeatureLocation(2, len(parent) + 3, strand=1)
    assert extract_atp8(parent, ref, meta)[1] == 'out_of_bounds_cds_location'
    parent.features[1].location = FeatureLocation(2, 11, strand=None)
    assert extract_atp8(parent, ref, meta)[1] == 'unknown_cds_strand'
    parent, ref, meta = parent_fixture(join=True)
    parent.features[1].location.operator = 'order'
    assert extract_atp8(parent, ref, meta)[1] == 'unsupported_location_operator'
    parent, ref, meta = parent_fixture()
    parent.features[0].qualifiers['db_xref'] = ['taxon:9606']
    assert extract_atp8(parent, ref, meta)[1] == 'wrong_taxon'
