import pytest
import torch
from microprotein_lm.evidence_loss import (
    EVIDENCE_TIER_WEIGHTS,
    get_evidence_weight,
    EvidenceWeightedLoss,
)


def test_evidence_tier_weights():
    assert get_evidence_weight(1) == 1.0
    assert get_evidence_weight(2) == 0.75
    assert get_evidence_weight(3) == 0.50
    assert get_evidence_weight(4) == 0.10
    assert get_evidence_weight("tier1") == 1.0
    assert get_evidence_weight("tier4") == 0.10

    with pytest.raises(ValueError):
        get_evidence_weight(5)


def test_evidence_weighted_loss_computation():
    criterion = EvidenceWeightedLoss()

    # Batch of 2 samples, sequence length 4, vocab size 5
    torch.manual_seed(42)
    logits = torch.randn(2, 4, 5)
    targets = torch.tensor([
        [1, 2, 3, 4],
        [0, 1, 2, 3]
    ])

    # Sample 0 is Tier 1 (weight 1.0), Sample 1 is Tier 4 (weight 0.10)
    tiers = torch.tensor([1, 4])

    loss_weighted = criterion(logits, targets, tiers)
    assert isinstance(loss_weighted, torch.Tensor)
    assert loss_weighted.dim() == 0  # Scalar

    # Calculate unweighted loss per sample manually
    loss_s0 = torch.nn.functional.cross_entropy(logits[0], targets[0]).item()
    loss_s1 = torch.nn.functional.cross_entropy(logits[1], targets[1]).item()

    # Expected weighted mean: (1.0 * loss_s0 + 0.10 * loss_s1) / (1.0 + 0.10)
    expected = (1.0 * loss_s0 + 0.10 * loss_s1) / 1.10
    assert torch.isclose(loss_weighted, torch.tensor(expected, dtype=torch.float32), atol=1e-4)


def test_evidence_weighted_loss_with_initiation_priors():
    criterion = EvidenceWeightedLoss(use_initiation_priors=True)

    logits = torch.randn(2, 4, 5)
    targets = torch.tensor([
        [1, 2, 3, 4],
        [0, 1, 2, 3]
    ])
    tiers = torch.tensor([1, 1])  # Both Tier 1 (weight 1.0)
    # Sample 0 has canonical AUG (init_weight 1.0), Sample 1 has non-canonical ACG (init_weight 0.3)
    init_weights = torch.tensor([1.0, 0.3])

    loss_weighted = criterion(logits, targets, tiers, initiation_weights=init_weights)
    loss_s0 = torch.nn.functional.cross_entropy(logits[0], targets[0]).item()
    loss_s1 = torch.nn.functional.cross_entropy(logits[1], targets[1]).item()

    expected = (1.0 * 1.0 * loss_s0 + 1.0 * 0.3 * loss_s1) / (1.0 + 0.3)
    assert torch.isclose(loss_weighted, torch.tensor(expected, dtype=torch.float32), atol=1e-4)


def test_evidence_weighted_loss_ignore_index():
    criterion = EvidenceWeightedLoss(ignore_index=-100)

    logits = torch.randn(1, 4, 5)
    targets = torch.tensor([[1, 2, -100, -100]])
    tiers = torch.tensor([1])

    loss = criterion(logits, targets, tiers)
    assert not torch.isnan(loss)
    # Should equal cross entropy over first 2 valid positions
    loss_valid = torch.nn.functional.cross_entropy(logits[0, :2], targets[0, :2])
    assert torch.isclose(loss, loss_valid, atol=1e-4)
