"""Reproducible, training-only microprotein language-model experiments."""

from .extended_window import WindowConfig, ExtendedWindowExtractor
from .initiation import DEFAULT_INITIATION_PRIORS, NEAR_COGNATE_START_CODONS, InitiationPrior
from .multitrack import MultiTrackSample, MultiTrackDataset, collate_multitrack_batch
from .evidence_loss import EVIDENCE_TIER_WEIGHTS, get_evidence_weight, EvidenceWeightedLoss
from .translation_evidence import (
    calculate_tai,
    calculate_cai,
    calculate_rscu,
    get_per_codon_tai,
    compute_periodicity_score,
    compute_tis_peak_ratio,
    compute_cds_coverage_fraction,
    determine_evidence_tier,
    annotate_record_evidence,
    HUMAN_TRNA_COPY_NUMBERS,
)
from .tokenization import (
    BioTokenizer,
    Tokenizer,
    PAD_TOKEN,
    UNK_TOKEN,
    BOS_TOKEN,
    EOS_TOKEN,
    MASK_TOKEN,
    FRAME_0_TOKEN,
    FRAME_1_TOKEN,
    FRAME_2_TOKEN,
    CONTROL_TOKENS,
    FRAME_TOKENS,
    IUPAC_DEGENERATE_BASES,
)

__version__ = "0.1.0"

__all__ = [
    "BioTokenizer",
    "Tokenizer",
    "PAD_TOKEN",
    "UNK_TOKEN",
    "BOS_TOKEN",
    "EOS_TOKEN",
    "MASK_TOKEN",
    "FRAME_0_TOKEN",
    "FRAME_1_TOKEN",
    "FRAME_2_TOKEN",
    "CONTROL_TOKENS",
    "FRAME_TOKENS",
    "IUPAC_DEGENERATE_BASES",
    "WindowConfig",
    "ExtendedWindowExtractor",
    "DEFAULT_INITIATION_PRIORS",
    "NEAR_COGNATE_START_CODONS",
    "InitiationPrior",
    "MultiTrackSample",
    "MultiTrackDataset",
    "collate_multitrack_batch",
    "EVIDENCE_TIER_WEIGHTS",
    "get_evidence_weight",
    "EvidenceWeightedLoss",
    "calculate_tai",
    "calculate_cai",
    "calculate_rscu",
    "get_per_codon_tai",
    "compute_periodicity_score",
    "compute_tis_peak_ratio",
    "compute_cds_coverage_fraction",
    "determine_evidence_tier",
    "annotate_record_evidence",
    "HUMAN_TRNA_COPY_NUMBERS",
]

