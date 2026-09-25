import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from microprotein_lm.quality import (
    TaxonomyResolver,
    CellularCompartmentResolver,
    ThreeTierQualityCategorizer,
    categorize_quality_tier,
    validate_record,
)


def test_taxonomy_resolver_human():
    resolver = TaxonomyResolver()
    res = resolver.get_lineage(9606)
    assert res['tax_id'] == 9606
    assert 'Homo sapiens' in res['scientific_name']
    assert res['is_eukaryota'] is True
    assert res['is_human'] is True
    assert len(res['lineage']) > 0


def test_taxonomy_resolver_mouse():
    resolver = TaxonomyResolver()
    res = resolver.get_lineage(10090)
    assert res['tax_id'] == 10090
    assert res['is_human'] is False
    assert res['is_eukaryota'] is True


def test_cellular_compartment_resolver():
    resolver = CellularCompartmentResolver()

    # Mitochondrion
    res_mito = resolver.resolve_compartment('mitochondrion')
    assert res_mito['primary_compartment'] == 'mitochondrion'
    assert res_mito['go_id'] == 'GO:0005739'
    assert res_mito['is_organelle'] is True

    # Nucleus
    res_nuc = resolver.resolve_compartment('nucleus')
    assert res_nuc['primary_compartment'] == 'nucleus'
    assert res_nuc['go_id'] in ('GO:0005634', 'GO:0140707')
    assert res_nuc['is_organelle'] is True

    # Cytosol
    res_cyto = resolver.resolve_compartment('cytosol')
    assert res_cyto['primary_compartment'] == 'cytosol'
    assert res_cyto['go_id'] == 'GO:0005829'


def test_three_tier_quality_categorizer():
    categorizer = ThreeTierQualityCategorizer()

    # Tier 1: Gold / Swiss-Prot
    tier1_rec = {
        'protein': 'MPSV',
        'reference_identity': 1.0,
        'initiation_type': 'CANONICAL',
        'initiation_codon': 'AUG',
    }
    ref_tier1 = {'reviewed': True, 'uniprot': 'P12345'}
    res1 = categorizer.categorize(tier1_rec, reference=ref_tier1)
    assert res1['tier'] == 1
    assert 'Tier 1' in res1['tier_name']
    assert len(res1['rationale']) > 0

    # Tier 2: Genomic CDS
    tier2_rec = {
        'protein': 'MPSV',
        'reference_identity': 0.88,
        'initiation_type': 'CANONICAL',
        'initiation_codon': 'AUG',
    }
    ref_tier2 = {'reviewed': False, 'uniprot': 'A0A000'}
    res2 = categorizer.categorize(tier2_rec, reference=ref_tier2)
    assert res2['tier'] == 2
    assert 'Tier 2' in res2['tier_name']

    # Tier 3: sORF lead / Non-canonical initiation
    tier3_rec = {
        'protein': 'MPSV',
        'reference_identity': 0.80,
        'initiation_type': 'NEAR_COGNATE',
        'initiation_codon': 'CUG',
        'is_sorf_lead': True,
    }
    res3 = categorizer.categorize(tier3_rec)
    assert res3['tier'] == 3
    assert 'Tier 3' in res3['tier_name']


def test_validate_record_includes_taxonomy_and_compartment():
    record = SeqRecord(Seq("ATGCCAAGCGTTTAG"), id="TEST_PROT")
    feature = SeqFeature(FeatureLocation(0, 15), type="CDS")
    feature.qualifiers = {
        'protein_id': ['TEST_PROT'],
        'transl_table': ['1'],
        'translation': ['MPSV']
    }
    source_feature = SeqFeature(FeatureLocation(0, 15), type="source")
    source_feature.qualifiers = {'db_xref': ['taxon:9606']}
    record.features = [feature, source_feature]

    metadata = {
        'protein_id': 'TEST_PROT',
        'base_count': 15,
        'parent_accession': 'TEST_PARENT'
    }

    reference = {
        'gene': 'TEST_GENE',
        'uniprot': 'P00000',
        'reference_protein': 'MPSV',
        'reviewed': True,
        'impi_class': 'Verified mitochondrial'
    }

    res, reason = validate_record(record, metadata, reference)
    assert reason is None
    assert 'taxonomy_lineage' in res
    assert res['taxonomy_lineage']['tax_id'] == 9606
    assert 'cellular_compartment' in res
    assert res['cellular_compartment']['primary_compartment'] == 'mitochondrion'
    assert 'quality_tier' in res
    assert res['quality_tier'] == 1
    assert 'Tier 1' in res['quality_tier_name']
