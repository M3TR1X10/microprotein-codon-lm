import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from microprotein_lm.quality import DynamicIdentityGate, validate_record

def test_dynamic_identity_gate():
    references = ["MPSV", "MPSVQQ"] # Includes a canonical and an isoform
    gate = DynamicIdentityGate(references, min_identity=0.7)
    
    # Exact match
    res, reason = gate.evaluate("MPSV")
    assert reason is None
    assert res['variant_class'] == 'EXACT_MATCH'
    assert res['identity'] == 1.0

    # Substitution
    res, reason = gate.evaluate("MQSV")
    assert reason is None
    assert res['variant_class'] == 'SUBSTITUTION_ONLY'
    assert res['identity'] == 0.75 # 3/4
    
    # N-Terminal Extension
    res, reason = gate.evaluate("MMPSV")
    assert reason is None
    assert res['variant_class'] == 'N_TERMINAL_EXTENSION'
    assert res['identity'] == 1.0
    
    # C-Terminal Truncation
    res, reason = gate.evaluate("MPS")
    assert reason is None
    assert res['variant_class'] == 'C_TERMINAL_TRUNCATION'
    assert res['identity'] == 0.75
    
    # Isoform Match
    res, reason = gate.evaluate("MPSVQQ")
    assert reason is None
    assert res['variant_class'] == 'EXACT_MATCH'
    assert res['identity'] == 1.0

def test_validate_record_with_isoforms():
    # Construct a dummy record with a CDS feature
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
    
    # Reference with just one canonical protein
    reference = {
        'gene': 'TEST_GENE',
        'uniprot': 'P00000',
        'reference_protein': 'MPSV'
    }
    
    res, reason = validate_record(record, metadata, reference)
    assert reason is None
    assert res['variant_class'] == 'EXACT_MATCH'
    
    # Reference with isoforms dictionary
    reference_isoforms = {
        'gene': 'TEST_GENE',
        'uniprot': 'P00000',
        'isoforms': {
            'P00000-1': 'MPSVQ',
            'P00000-2': 'MPSV'
        }
    }
    
    res2, reason2 = validate_record(record, metadata, reference_isoforms)
    assert reason2 is None
    assert res2['variant_class'] == 'EXACT_MATCH'
