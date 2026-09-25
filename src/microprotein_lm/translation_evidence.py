"""Ribo-seq translation evidence integration, tRNA Adaptation Indexing (tAI), and codon usage bias (CAI/RSCU)."""

import math
from typing import Dict, List, Optional, Tuple, Union, Any


# Standard genetic code dictionary (RNA codons -> single-letter amino acid code)
CODON_TO_AA: Dict[str, str] = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAU": "N",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACU": "T",
    "AGA": "R", "AGC": "S", "AGG": "R", "AGU": "S",
    "AUA": "I", "AUC": "I", "AUG": "M", "AUU": "I",
    "CAA": "Q", "CAC": "H", "CAG": "Q", "CAU": "H",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R",
    "CUA": "L", "CUC": "L", "CUG": "L", "CUU": "L",
    "GAA": "E", "GAC": "D", "GAG": "E", "GAU": "D",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCU": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGU": "G",
    "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V",
    "UAA": "*", "UAC": "Y", "UAG": "*", "UAU": "Y",
    "UCA": "S", "UCC": "S", "UCG": "S", "UCU": "S",
    "UGA": "*", "UGC": "C", "UGG": "W", "UGU": "C",
    "UUA": "L", "UUC": "F", "UUG": "L", "UUU": "F",
}

# Amino acid to synonymous RNA codons list
AA_TO_CODONS: Dict[str, List[str]] = {}
for codon, aa in CODON_TO_AA.items():
    if aa != "*":
        AA_TO_CODONS.setdefault(aa, []).append(codon)

# Human nuclear tRNA gene copy numbers (GRCh38 / GtRNAdb) per codon (or wobble-adjusted)
HUMAN_TRNA_COPY_NUMBERS: Dict[str, int] = {
    "AAA": 19, "AAC": 14, "AAG": 22, "AAU": 14,
    "ACA": 10, "ACC": 16, "ACG": 6,  "ACU": 12,
    "AGA": 8,  "AGC": 11, "AGG": 5,  "AGU": 11,
    "AUA": 5,  "AUC": 16, "AUG": 20, "AUU": 16,
    "CAA": 10, "CAC": 17, "CAG": 23, "CAU": 17,
    "CCA": 8,  "CCC": 12, "CCG": 4,  "CCU": 10,
    "CGA": 6,  "CGC": 10, "CGG": 6,  "CGU": 8,
    "CUA": 6,  "CUC": 14, "CUG": 21, "CUU": 10,
    "GAA": 15, "GAC": 17, "GAG": 10, "GAU": 17,
    "GCA": 11, "GCC": 26, "GCG": 5,  "GCU": 26,
    "GGA": 11, "GGC": 16, "GGG": 7,  "GGU": 16,
    "GUA": 5,  "GUC": 12, "GUG": 18, "GUU": 10,
    "UAC": 14, "UAU": 14, "UCA": 10, "UCC": 16,
    "UCG": 5,  "UCU": 12, "UGC": 10, "UGG": 10,
    "UGU": 10, "UUA": 6,  "UUC": 16, "UUG": 12,
    "UUU": 16,
}

# Reference human relative synonymous codon usage (RSCU) values (Kazusa codon usage database)
HUMAN_REFERENCE_RSCU: Dict[str, float] = {
    "AAA": 0.86, "AAC": 1.07, "AAG": 1.14, "AAU": 0.93,
    "ACA": 0.80, "ACC": 1.43, "ACG": 0.44, "ACU": 0.98,
    "AGA": 1.25, "AGC": 1.22, "AGG": 1.21, "AGU": 0.90,
    "AUA": 0.51, "AUC": 1.41, "AUG": 1.00, "AUU": 1.08,
    "CAA": 0.54, "CAC": 1.17, "CAG": 1.46, "CAU": 0.83,
    "CCA": 0.80, "CCC": 1.28, "CCG": 0.44, "CCU": 1.14,
    "CGA": 0.65, "CGC": 1.13, "CGG": 1.20, "CGU": 0.54,
    "CUA": 0.43, "CUC": 1.18, "CUG": 2.37, "CUU": 0.77,
    "GAA": 0.84, "GAC": 1.08, "GAG": 1.16, "GAU": 0.92,
    "GCA": 0.86, "GCC": 1.60, "GCG": 0.44, "GCU": 1.05,
    "GGA": 0.98, "GGC": 1.34, "GGG": 0.97, "GGU": 0.65,
    "GUA": 0.47, "GUC": 1.43, "GUG": 2.81, "GUU": 0.72,
    "UAC": 1.14, "UAU": 0.86, "UCA": 0.90, "UCC": 1.33,
    "UCG": 0.33, "UCU": 1.14, "UGC": 1.09, "UGG": 1.00,
    "UGU": 0.91, "UUA": 0.46, "UUC": 1.09, "UUG": 0.76,
    "UUU": 0.91,
}

# Standard wobble parameters (dos Reis et al. 2004)
# Maps (codon_third_base, anticodon_first_base) -> (1 - s_ij) weight
WOBBLE_EFFICIENCY_WEIGHTS: Dict[str, float] = {
    "WC": 1.0,        # Watson-Crick match
    "I_U": 0.5606,    # Inosine with U
    "I_C": 0.2901,    # Inosine with C
    "I_A": 0.9999,    # Inosine with A
    "G_U": 0.4100,    # G with U
    "U_G": 0.2800,    # U with G
}


def _clean_sequence_to_codons(sequence: Union[str, List[str]]) -> List[str]:
    """Converts input string or list of codons into clean 3-character uppercase RNA codon list."""
    if isinstance(sequence, str):
        rna = sequence.upper().replace('T', 'U').strip()
        # Keep complete codons
        codons = [rna[i:i+3] for i in range(0, len(rna) - len(rna) % 3, 3)]
    else:
        codons = [c.upper().replace('T', 'U') for c in sequence]
    return [c for c in codons if len(c) == 3 and set(c).issubset(set("ACGU"))]


def calculate_rscu(codons: List[str]) -> Dict[str, float]:
    """Calculates Relative Synonymous Codon Usage (RSCU) for a set or sequence of codons.

    RSCU_ij = X_ij / ( (1 / n_i) * sum_{k=1}^{n_i} X_ik )

    Args:
        codons: List of 3-letter RNA codon strings.

    Returns:
        Dict mapping codon to RSCU float value.
    """
    clean_codons = _clean_sequence_to_codons(codons)
    codon_counts: Dict[str, int] = {c: 0 for c in CODON_TO_AA if CODON_TO_AA[c] != "*"}
    for c in clean_codons:
        if c in codon_counts:
            codon_counts[c] += 1

    rscu: Dict[str, float] = {}
    for aa, syn_codons in AA_TO_CODONS.items():
        total_aa_count = sum(codon_counts[c] for c in syn_codons)
        n_syn = len(syn_codons)
        expected_per_codon = total_aa_count / float(n_syn) if n_syn > 0 else 0.0

        for c in syn_codons:
            if expected_per_codon > 0:
                rscu[c] = codon_counts[c] / expected_per_codon
            else:
                rscu[c] = 1.0  # Default neutral RSCU if amino acid never observed

    return rscu


def calculate_cai(
    sequence: Union[str, List[str]],
    reference_rscu: Optional[Dict[str, float]] = None
) -> float:
    """Calculates Codon Adaptation Index (CAI) for a coding sequence.

    CAI = exp( (1 / L) * sum_{i=1}^L ln(w_i) )
    where w_i = RSCU_i / max(RSCU_synonymous).

    Args:
        sequence: RNA nucleotide string or list of 3-letter codons.
        reference_rscu: Optional dictionary of reference RSCU values. Uses human reference if None.

    Returns:
        Float CAI score between 0.0 and 1.0.
    """
    codons = _clean_sequence_to_codons(sequence)
    # Exclude stop codons and ATG/TGG (single-codon amino acids) from CAI calculation
    sense_codons = [c for c in codons if CODON_TO_AA.get(c, "*") not in ("*", "M", "W")]
    if not sense_codons:
        return 1.0

    rscu_ref = reference_rscu if reference_rscu is not None else HUMAN_REFERENCE_RSCU

    # Compute relative adaptiveness w_i for each codon relative to max in its synonymous group
    max_rscu_per_aa: Dict[str, float] = {}
    for aa, syn_codons in AA_TO_CODONS.items():
        max_rscu_per_aa[aa] = max(rscu_ref.get(c, 1.0) for c in syn_codons)

    log_w_sum = 0.0
    valid_count = 0
    for c in sense_codons:
        aa = CODON_TO_AA.get(c, "*")
        if aa in max_rscu_per_aa:
            rscu_val = rscu_ref.get(c, 1.0)
            max_rscu = max_rscu_per_aa[aa]
            w = rscu_val / max_rscu if max_rscu > 0 else 1.0
            w = max(w, 0.01)  # Bound away from 0 to prevent log(0)
            log_w_sum += math.log(w)
            valid_count += 1

    if valid_count == 0:
        return 1.0

    return math.exp(log_w_sum / valid_count)


def get_per_codon_tai(
    sequence: Union[str, List[str]],
    trna_copy_numbers: Optional[Dict[str, int]] = None,
    custom_wobble_weights: Optional[Dict[str, float]] = None
) -> List[float]:
    """Computes per-codon tRNA Adaptation Index (w_k relative adaptiveness) vector for a sequence.

    Args:
        sequence: RNA nucleotide string or list of codons.
        trna_copy_numbers: Optional dict of codon/anticodon tRNA gene copy numbers.
        custom_wobble_weights: Optional dict of wobble weights per codon.

    Returns:
        List of float relative adaptiveness values w_k for each codon in the sequence.
    """
    codons = _clean_sequence_to_codons(sequence)
    if not codons:
        return []

    trna_counts = trna_copy_numbers if trna_copy_numbers is not None else HUMAN_TRNA_COPY_NUMBERS

    # Calculate absolute adaptiveness W_c for all 61 sense codons
    W_scores: Dict[str, float] = {}
    for c in CODON_TO_AA:
        if CODON_TO_AA[c] == "*":
            continue
        if custom_wobble_weights and c in custom_wobble_weights:
            W_scores[c] = custom_wobble_weights[c]
        else:
            # Direct or wobble-adjusted tRNA copy count lookup
            W_scores[c] = float(trna_counts.get(c, 10))

    max_W = max(W_scores.values()) if W_scores else 1.0
    if max_W <= 0:
        max_W = 1.0

    # Compute relative adaptiveness w_c = W_c / max_W
    w_scores: Dict[str, float] = {c: max(W / max_W, 0.01) for c, W in W_scores.items()}

    return [w_scores.get(c, 0.5) for c in codons]


def calculate_tai(
    sequence: Union[str, List[str]],
    trna_copy_numbers: Optional[Dict[str, int]] = None,
    custom_wobble_weights: Optional[Dict[str, float]] = None
) -> float:
    """Calculates overall tRNA Adaptation Index (tAI) for a coding sequence.

    tAI = exp( (1 / L) * sum_{k=1}^L ln(w_k) )

    Args:
        sequence: RNA nucleotide string or list of codons.
        trna_copy_numbers: Optional dict of tRNA gene copy numbers.
        custom_wobble_weights: Optional dict of custom codon weights.

    Returns:
        Float tAI score between 0.0 and 1.0.
    """
    per_codon_w = get_per_codon_tai(sequence, trna_copy_numbers, custom_wobble_weights)
    if not per_codon_w:
        return 1.0

    log_w_sum = sum(math.log(max(w, 0.01)) for w in per_codon_w)
    return math.exp(log_w_sum / len(per_codon_w))


def compute_periodicity_score(coverage: Union[List[float], Any]) -> float:
    """Computes 3-nt translation periodicity score (frame 0 read fraction) from Ribo-seq footprint coverage.

    Frame 0 fraction: P_0 = Frame_0_reads / Total_reads.
    A score of 1.0 means perfect 3-nt frame-0 alignment; 0.33 represents random background.

    Args:
        coverage: List or array of numerical P-site footprint counts per nucleotide.

    Returns:
        Float periodicity score between 0.0 and 1.0.
    """
    if coverage is None or len(coverage) < 3:
        return 0.33

    cov_list = [float(x) for x in coverage]
    f0 = sum(cov_list[i] for i in range(0, len(cov_list), 3))
    f1 = sum(cov_list[i] for i in range(1, len(cov_list), 3))
    f2 = sum(cov_list[i] for i in range(2, len(cov_list), 3))
    total = f0 + f1 + f2

    if total <= 0:
        return 0.33

    return f0 / total


def compute_tis_peak_ratio(
    coverage: Union[List[float], Any],
    tis_offset: int = 0,
    window_nt: int = 9
) -> float:
    """Computes Translation Initiation Site (TIS) Ribo-seq peak ratio.

    Peak ratio = (Mean footprint density in TIS window) / (Mean background footprint density in CDS).

    Args:
        coverage: List or array of numerical footprint counts.
        tis_offset: Start nucleotide position of TIS.
        window_nt: Length of TIS peak window in nucleotides (default 9 nt = 3 codons).

    Returns:
        Float peak ratio >= 0.0.
    """
    if coverage is None or len(coverage) == 0:
        return 1.0

    cov_list = [float(x) for x in coverage]
    total_len = len(cov_list)

    start_idx = max(0, tis_offset)
    end_idx = min(total_len, start_idx + window_nt)

    if start_idx >= total_len or end_idx <= start_idx:
        return 1.0

    tis_reads = cov_list[start_idx:end_idx]
    tis_mean = sum(tis_reads) / float(len(tis_reads))

    background_reads = cov_list[:start_idx] + cov_list[end_idx:]
    if not background_reads:
        return tis_mean

    bg_mean = sum(background_reads) / float(len(background_reads))
    if bg_mean <= 0:
        return tis_mean if tis_mean > 0 else 1.0

    return tis_mean / bg_mean


def compute_cds_coverage_fraction(
    coverage: Union[List[float], Any],
    min_reads: float = 0.1
) -> float:
    """Computes fraction of nucleotide or codon positions covered by Ribo-seq footprints.

    Args:
        coverage: List or array of footprint counts.
        min_reads: Minimum footprint count threshold to consider a position covered.

    Returns:
        Float fraction covered between 0.0 and 1.0.
    """
    if coverage is None or len(coverage) == 0:
        return 0.0

    cov_list = [float(x) for x in coverage]
    covered = sum(1 for x in cov_list if x >= min_reads)
    return covered / float(len(cov_list))


def determine_evidence_tier(
    has_ms_evidence: bool = False,
    ribo_coverage: Optional[List[float]] = None,
    periodicity_score: Optional[float] = None,
    tis_peak_ratio: Optional[float] = None,
    is_conserved: bool = False
) -> int:
    """Determines the evidence tier (1..4) according to the 4-tier evidence hierarchy.

    Hierarchy rules:
        - Tier 1 (1.0 weight): MS/MS peptide evidence verified (+ optional Ribo-seq).
        - Tier 2 (0.75 weight): High-confidence Ribo-seq translation evidence
          (periodicity >= 0.5, TIS peak ratio >= 2.0, or CDS coverage >= 0.5).
        - Tier 3 (0.50 weight): High sequence evolutionary conservation.
        - Tier 4 (0.10 weight): Computational sORF prediction / baseline without direct empirical proof.

    Args:
        has_ms_evidence: True if mass spectrometry peptide match exists.
        ribo_coverage: Optional list of Ribo-seq footprint counts.
        periodicity_score: Optional pre-calculated periodicity score.
        tis_peak_ratio: Optional pre-calculated TIS peak ratio.
        is_conserved: True if sequence is evolutionarily conserved across species.

    Returns:
        Integer tier 1, 2, 3, or 4.
    """
    if has_ms_evidence:
        return 1

    # Check Ribo-seq translation evidence for Tier 2
    if periodicity_score is None and ribo_coverage is not None:
        periodicity_score = compute_periodicity_score(ribo_coverage)

    if tis_peak_ratio is None and ribo_coverage is not None:
        tis_peak_ratio = compute_tis_peak_ratio(ribo_coverage)

    cov_frac = compute_cds_coverage_fraction(ribo_coverage) if ribo_coverage is not None else 0.0

    is_ribo_active = (
        (periodicity_score is not None and periodicity_score >= 0.50) or
        (tis_peak_ratio is not None and tis_peak_ratio >= 2.0) or
        (cov_frac >= 0.50)
    )

    if is_ribo_active:
        return 2

    if is_conserved:
        return 3

    return 4


def annotate_record_evidence(
    record: Dict[str, Any],
    ribo_coverage: Optional[List[float]] = None,
    has_ms_evidence: bool = False,
    is_conserved: bool = False
) -> Dict[str, Any]:
    """Annotates a sequence metadata dict with translation evidence, tAI, CAI, and evidence tier.

    Args:
        record: Sequence record dictionary containing 'rna' or 'dna' string.
        ribo_coverage: Optional list of Ribo-seq P-site coverage values.
        has_ms_evidence: True if mass spectrometry proof exists.
        is_conserved: True if evolutionary conservation proof exists.

    Returns:
        Updated dictionary containing additional evidence and usage index fields.
    """
    updated = dict(record)
    seq = updated.get("rna", updated.get("dna", ""))

    tai_score = calculate_tai(seq)
    cai_score = calculate_cai(seq)

    periodicity = compute_periodicity_score(ribo_coverage) if ribo_coverage is not None else 0.33
    tis_peak = compute_tis_peak_ratio(ribo_coverage) if ribo_coverage is not None else 1.0
    cov_frac = compute_cds_coverage_fraction(ribo_coverage) if ribo_coverage is not None else 0.0

    tier = determine_evidence_tier(
        has_ms_evidence=has_ms_evidence,
        ribo_coverage=ribo_coverage,
        periodicity_score=periodicity,
        tis_peak_ratio=tis_peak,
        is_conserved=is_conserved
    )

    updated["tai_score"] = tai_score
    updated["cai_score"] = cai_score
    updated["ribo_periodicity_score"] = periodicity
    updated["ribo_tis_peak_ratio"] = tis_peak
    updated["ribo_coverage_fraction"] = cov_frac
    updated["evidence_tier"] = tier
    if ribo_coverage is not None:
        updated["ribo_coverage"] = ribo_coverage

    return updated
