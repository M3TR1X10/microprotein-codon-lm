"""Reproducible, training-only microprotein language-model experiments."""

from .extended_window import WindowConfig, ExtendedWindowExtractor
from .initiation import DEFAULT_INITIATION_PRIORS, NEAR_COGNATE_START_CODONS, InitiationPrior
from .multitrack import MultiTrackSample, MultiTrackDataset, collate_multitrack_batch
from .evidence_loss import EVIDENCE_TIER_WEIGHTS, get_evidence_weight, EvidenceWeightedLoss

__version__ = "0.1.0"

__all__ = [
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
]
