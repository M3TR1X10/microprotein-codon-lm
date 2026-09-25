import pytest
from microprotein_lm.translation_engine import TranslationEngine


def test_table1_standard_translation():
    engine = TranslationEngine()
    # Table 1: ATG (M), GCT (A), TAA (Stop)
    seq = "ATGGCTTAA"
    protein = engine.translate_cds(seq, table=1, cds=True)
    assert protein == "MA"

    # Test RNA input
    seq_rna = "AUGGCUUAA"
    protein_rna = engine.translate_cds(seq_rna, table=1, cds=True)
    assert protein_rna == "MA"


def test_table2_vertebrate_mitochondrial():
    engine = TranslationEngine()
    # Table 2: AUA is Met (M), UGA is Trp (W), AGA/AGG are Stop (*)
    # Sequence: AUG (M), AUA (M), UGA (W), AGA (Stop)
    seq = "ATGATATGAAGA"
    protein = engine.translate_cds(seq, table=2, cds=True)
    assert protein == "MMW"

    assert engine.is_valid_stop("AGA", table=2) is True
    assert engine.is_valid_stop("AGG", table=2) is True
    assert engine.is_valid_stop("TGA", table=2) is False


def test_initiation_types_and_validation():
    engine = TranslationEngine()
    # Table 1 canonical vs near cognate
    assert engine.get_initiation_type("AUG", table=1) == "canonical"
    assert engine.get_initiation_type("GUG", table=1) == "near_cognate"
    assert engine.get_initiation_type("CCC", table=1) == "invalid"

    assert engine.is_valid_initiation("AUG", table=1) is True
    assert engine.is_valid_initiation("GUG", table=1) is True
    assert engine.is_valid_initiation("CCC", table=1) is False


def test_translation_engine_tables_coverage():
    engine = TranslationEngine()
    # Verify tables 1 through 33 are supported without KeyError
    for t in range(1, 34):
        starts = engine.get_initiation_codons(table=t)
        stops = engine.get_stop_codons(table=t)
        table_map = engine.get_table_map(table=t)
        assert len(starts) > 0
        assert len(stops) > 0
        assert len(table_map) == 64


def test_diagnose_cds():
    engine = TranslationEngine()
    # Valid CDS
    diag = engine.diagnose_cds("AUGGCUUAA", table=1)
    assert diag['is_valid'] is True
    assert diag['errors'] == []
    assert diag['initiation_type'] == "canonical"
    assert diag['translated_protein'] == "MA"

    # Invalid length
    diag_len = engine.diagnose_cds("AUGGCUUA", table=1)
    assert diag_len['is_valid'] is False
    assert "incomplete_codon" in diag_len['errors']

    # Internal stop codon
    diag_stop = engine.diagnose_cds("AUGUAAAGCUAA", table=1)
    assert diag_stop['is_valid'] is False
    assert any("internal_stop" in e for e in diag_stop['errors'])
