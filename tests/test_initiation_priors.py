import pytest
import torch
from microprotein_lm.initiation import (
    DEFAULT_INITIATION_PRIORS,
    NEAR_COGNATE_START_CODONS,
    InitiationPrior,
)
from microprotein_lm.tokenization import Tokenizer


def test_default_priors():
    prior = InitiationPrior()
    assert prior.get_weight("AUG") == 1.0
    assert prior.get_weight("CUG") == 0.6
    assert prior.get_weight("GUG") == 0.5
    assert prior.get_weight("UUG") == 0.4
    assert prior.get_weight("ACG") == 0.3
    assert prior.get_weight("AUU") == 0.2
    # DNA format string handling
    assert prior.get_weight("ATG") == 1.0
    # Unknown codon default
    assert prior.get_weight("CCC") == 0.0


def test_custom_priors():
    custom = {"AUG": 1.0, "CUG": 0.9, "XYZ": 0.8}
    prior = InitiationPrior(prior_table=custom, default_weight=0.05)
    assert prior.get_weight("CUG") == 0.9
    assert prior.get_weight("XYZ") == 0.8
    assert prior.get_weight("AAA") == 0.05


def test_prior_tensor():
    tokenizer = Tokenizer("codon")
    prior = InitiationPrior()
    tensor = prior.get_prior_tensor(tokenizer.vocabulary)
    assert isinstance(tensor, torch.Tensor)
    assert tensor.shape == (64,)
    aug_idx = tokenizer.token_to_id["AUG"]
    cug_idx = tokenizer.token_to_id["CUG"]
    assert tensor[aug_idx].item() == pytest.approx(1.0)
    assert tensor[cug_idx].item() == pytest.approx(0.6)


def test_calculate_frame_offsets_base():
    prior = InitiationPrior()
    # 5' UTR of len 5, TIS at idx 5, CDS of 6 bases
    # total len = 11 bases, tis_offset = 5
    offsets = prior.calculate_frame_offsets(seq_len=11, tis_offset=5, token_width=1)
    assert len(offsets) == 11
    # 5' UTR positions (0..4) should be -1
    assert offsets[:5] == [-1, -1, -1, -1, -1]
    # CDS positions (5..10) should be 0, 1, 2, 0, 1, 2
    assert offsets[5:] == [0, 1, 2, 0, 1, 2]


def test_calculate_frame_offsets_codon():
    prior = InitiationPrior()
    # 2 UTR codons, TIS at codon idx 2, 3 CDS codons
    # total len = 5 codons, tis_offset = 2
    offsets = prior.calculate_frame_offsets(seq_len=5, tis_offset=2, token_width=3)
    assert len(offsets) == 5
    assert offsets[:2] == [-1, -1]
    assert offsets[2:] == [0, 0, 0]
