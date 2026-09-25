"""Unit tests for Ribo-seq translation evidence integration, tRNA Adaptation Indexing (tAI), and codon usage bias."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition

from microprotein_lm.translation_evidence import (
    calculate_rscu,
    calculate_cai,
    calculate_tai,
    get_per_codon_tai,
    compute_periodicity_score,
    compute_tis_peak_ratio,
    compute_cds_coverage_fraction,
    determine_evidence_tier,
    annotate_record_evidence,
    HUMAN_TRNA_COPY_NUMBERS,
)
from microprotein_lm.quality import validate_record


def test_rscu_calculation():
    # Test sequence with multiple codons
    codons = ["ATG", "GCC", "GCC", "GCU", "TAA"]
    rscu = calculate_rscu(codons)
    assert isinstance(rscu, dict)
    assert "GCC" in rscu
    assert "GCU" in rscu
    # GCC count = 2, GCU count = 1. Total Ala codons = 3, n_syn = 4.
    # Expected per codon = 3/4 = 0.75. RSCU(GCC) = 2 / 0.75 = 2.666...
    assert pytest.approx(rscu["GCC"], 0.01) == 8.0 / 3.0


def test_cai_calculation():
    seq = "ATGCCCGCCGCGGCUGAA"
    cai = calculate_cai(seq)
    assert 0.0 <= cai <= 1.0
    assert isinstance(cai, float)


def test_tai_calculation():
    seq = "AUGGCCGCCGCUGAACUGUAA"
    tai = calculate_tai(seq)
    assert 0.0 < tai <= 1.0
    assert isinstance(tai, float)


def test_per_codon_tai():
    seq = "AUGGCCGCU"
    weights = get_per_codon_tai(seq)
    assert len(weights) == 3
    assert all(0.0 <= w <= 1.0 for w in weights)


def test_ribo_periodicity_score():
    # Perfect frame 0 coverage
    frame0_cov = [10.0, 0.0, 0.0, 15.0, 0.0, 0.0, 20.0, 0.0, 0.0]
    score0 = compute_periodicity_score(frame0_cov)
    assert score0 == 1.0

    # Off-frame flat coverage
    flat_cov = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
    score_flat = compute_periodicity_score(flat_cov)
    assert pytest.approx(score_flat, 0.01) == 1.0 / 3.0


def test_ribo_tis_peak_ratio():
    # Elevated TIS peak at position 0..9 (3 codons)
    cov = [100.0] * 9 + [10.0] * 90
    ratio = compute_tis_peak_ratio(cov, tis_offset=0, window_nt=9)
    assert pytest.approx(ratio, 0.1) == 10.0


def test_ribo_cds_coverage_fraction():
    cov = [1.0, 0.0, 2.0, 0.0, 5.0, 0.0]
    frac = compute_cds_coverage_fraction(cov, min_reads=0.5)
    assert pytest.approx(frac, 0.01) == 0.5


def test_evidence_tier_determination():
    # Tier 1: MS/MS evidence
    t1 = determine_evidence_tier(has_ms_evidence=True)
    assert t1 == 1

    # Tier 2: Active Ribo-seq evidence (high frame 0 periodicity)
    ribo_high_period = [10.0, 0.0, 0.0] * 10
    t2 = determine_evidence_tier(has_ms_evidence=False, ribo_coverage=ribo_high_period)
    assert t2 == 2

    # Tier 3: Conservation evidence
    t3 = determine_evidence_tier(has_ms_evidence=False, ribo_coverage=None, is_conserved=True)
    assert t3 == 3

    # Tier 4: Computational sORF prediction / baseline
    t4 = determine_evidence_tier(has_ms_evidence=False, ribo_coverage=None, is_conserved=False)
    assert t4 == 4


def test_annotate_record_evidence():
    rec = {"rna": "AUGGCCGCUGAACUGUAA", "gene": "TEST"}
    ribo_cov = [10.0, 0.0, 0.0] * 6
    annotated = annotate_record_evidence(rec, ribo_coverage=ribo_cov, has_ms_evidence=False, is_conserved=True)

    assert "tai_score" in annotated
    assert "cai_score" in annotated
    assert "evidence_tier" in annotated
    assert "ribo_periodicity_score" in annotated
    assert "ribo_tis_peak_ratio" in annotated
    assert "ribo_coverage_fraction" in annotated
    assert annotated["evidence_tier"] == 2  # Active Ribo-seq periodicity


def test_quality_validate_record_integration():
    dna = "ATGGCCGCTGAAGAA" + "TAA"  # 18 nt = 6 codons = 5 aa + stop
    protein = "MAAEE"
    seq_obj = Seq(dna)

    cds_feat = SeqFeature(
        FeatureLocation(ExactPosition(0), ExactPosition(18), strand=1),
        type="CDS",
        qualifiers={
            "protein_id": ["TEST.1"],
            "transl_table": ["1"],
            "codon_start": ["1"],
            "translation": [protein],
        }
    )

    src_feat = SeqFeature(
        FeatureLocation(ExactPosition(0), ExactPosition(18), strand=1),
        type="source",
        qualifiers={"db_xref": ["taxon:9606"]}
    )

    record = SeqRecord(seq_obj, id="TEST.1", features=[cds_feat, src_feat])
    metadata = {
        "protein_id": "TEST.1",
        "parent_accession": "PARENT123",
        "base_count": len(dna),
    }
    reference = {
        "gene": "TESTGENE",
        "uniprot": "P12345",
        "reference_protein": protein,
    }

    res, reason = validate_record(record, metadata, reference)
    assert reason is None
    assert res is not None
    assert "tai_score" in res
    assert "cai_score" in res
    assert "evidence_tier" in res
    assert 0.0 < res["tai_score"] <= 1.0
    assert 0.0 <= res["cai_score"] <= 1.0
