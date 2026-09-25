"""Multi-tier experimental validation evidence hierarchy and weighted loss functions."""

from typing import Dict, Optional, Union
import torch
import torch.nn as nn
import torch.nn.functional as F


# 4-Tier evidence hierarchy with defined weights
EVIDENCE_TIER_WEIGHTS: Dict[Union[int, str], float] = {
    1: 1.0,   # Tier 1: MS/MS + Ribo-seq
    2: 0.75,  # Tier 2: Ribo-seq TIS
    3: 0.50,  # Tier 3: Conservation
    4: 0.10,  # Tier 4: Computational / sORF prediction
    "tier1": 1.0,
    "tier2": 0.75,
    "tier3": 0.50,
    "tier4": 0.10,
}


def get_evidence_weight(tier: Union[int, str], custom_weights: Optional[Dict] = None) -> float:
    """Returns the evidence weight for a given tier (1..4 or string name)."""
    weights = dict(EVIDENCE_TIER_WEIGHTS)
    if custom_weights:
        weights.update(custom_weights)
    if tier not in weights:
        raise ValueError(f"Unknown evidence tier: {tier}. Must be 1, 2, 3, 4 or tier1..tier4.")
    return weights[tier]


class EvidenceWeightedLoss(nn.Module):
    """Computes evidence-weighted cross-entropy loss L_weighted across dataset samples.

    Loss formula:
        L_weighted = sum(w_i * L_i) / sum(w_i)
    where w_i is the evidence tier weight for sample i (optionally scaled by initiation prior E_init),
    and L_i is the per-sample cross-entropy token loss.
    """

    def __init__(
        self,
        tier_weights: Optional[Dict] = None,
        ignore_index: int = -100,
        use_initiation_priors: bool = False
    ):
        super().__init__()
        self.tier_weights = dict(EVIDENCE_TIER_WEIGHTS)
        if tier_weights:
            self.tier_weights.update(tier_weights)
        self.ignore_index = ignore_index
        self.use_initiation_priors = use_initiation_priors

    def forward(
        self,
        logits: torch.Tensor,
        targets: torch.Tensor,
        evidence_tiers: torch.Tensor,
        initiation_weights: Optional[torch.Tensor] = None
    ) -> torch.Tensor:
        """
        Args:
            logits: Predicted logits tensor of shape (B, T, V) or (B*T, V).
            targets: Target class indices tensor of shape (B, T) or (B*T,).
            evidence_tiers: Tensor of shape (B,) containing integer tier numbers (1, 2, 3, 4).
            initiation_weights: Optional tensor of shape (B,) containing initiation prior weights E_init.

        Returns:
            Scalar PyTorch Tensor representing the evidence-weighted loss L_weighted.
        """
        if logits.dim() == 3:
            b, t, v = logits.shape
            logits_flat = logits.reshape(-1, v)
            targets_flat = targets.reshape(-1)
        else:
            b = evidence_tiers.shape[0]
            v = logits.shape[-1]
            t = logits.shape[0] // b
            logits_flat = logits
            targets_flat = targets

        # Per-token cross entropy without reduction (unreduced loss)
        per_token_loss = F.cross_entropy(logits_flat, targets_flat, ignore_index=self.ignore_index, reduction='none')
        per_token_loss = per_token_loss.view(b, t)

        # Compute valid token mask and per-sample unweighted mean loss
        valid_mask = (targets.view(b, t) != self.ignore_index).float()
        sample_valid_tokens = valid_mask.sum(dim=1)
        sample_token_loss_sum = (per_token_loss * valid_mask).sum(dim=1)

        # Per-sample loss (avoid division by zero)
        sample_loss = torch.where(
            sample_valid_tokens > 0,
            sample_token_loss_sum / (sample_valid_tokens + 1e-8),
            torch.zeros_like(sample_token_loss_sum)
        )

        # Get evidence weights w_i for each sample in batch
        sample_weights = torch.tensor(
            [self.tier_weights.get(int(t_idx.item()), 0.10) for t_idx in evidence_tiers],
            dtype=logits.dtype,
            device=logits.device
        )

        # Optionally scale sample weight by E_init prior weight
        if self.use_initiation_priors and initiation_weights is not None:
            sample_weights = sample_weights * initiation_weights.to(logits.device)

        # Weighted mean loss across batch
        total_weight = sample_weights.sum()
        if total_weight > 0:
            weighted_loss = (sample_weights * sample_loss).sum() / total_weight
        else:
            weighted_loss = sample_loss.mean()

        return weighted_loss
