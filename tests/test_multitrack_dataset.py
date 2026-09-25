import pytest
import torch
from microprotein_lm.tokenization import Tokenizer
from microprotein_lm.extended_window import WindowConfig
from microprotein_lm.multitrack import (
    MultiTrackSample,
    MultiTrackDataset,
    collate_multitrack_batch,
)


def test_multitrack_sample_dataclass():
    sample = MultiTrackSample(
        sequence_tokens=torch.tensor([0, 1, 2]),
        segment_ids=torch.tensor([0, 1, 2]),
        frame_offsets=torch.tensor([-1, 0, 1]),
        ribo_coverage=torch.tensor([0.0, 1.5, 3.2]),
        evidence_tier=1,
        initiation_weight=1.0,
    )
    assert sample.sequence_tokens.shape == (3,)
    assert sample.evidence_tier == 1
    assert sample.initiation_weight == 1.0


def test_multitrack_dataset_base_mode():
    records = [
        {
            'rna': "A" * 50 + "AUGGCGUAA" + "C" * 30,
            'tis_start': 50,
            'cds_len': 9,
            'evidence_tier': 1,
            'ribo_coverage': [0.1] * 89,
        },
        {
            'rna': "G" * 50 + "CUGGCGUAA" + "U" * 30,
            'tis_start': 50,
            'cds_len': 9,
            'evidence_tier': 2,
            'ribo_coverage': [0.5] * 89,
        }
    ]
    tokenizer = Tokenizer('base')
    window_cfg = WindowConfig(utr5_len=50, utr3_len=30)
    dataset = MultiTrackDataset(records, tokenizer, window_config=window_cfg)

    assert len(dataset) == 2
    s0 = dataset[0]
    assert s0.sequence_tokens.shape[0] == 89
    assert s0.segment_ids.shape[0] == 89
    assert s0.frame_offsets.shape[0] == 89
    assert s0.ribo_coverage.shape[0] == 89
    assert s0.evidence_tier == 1
    assert s0.initiation_weight == 1.0  # AUG

    s1 = dataset[1]
    assert s1.evidence_tier == 2
    assert s1.initiation_weight == 0.6  # CUG near-cognate start


def test_multitrack_dataset_codon_mode():
    records = [
        {
            'rna': "A" * 48 + "AUGGCGUAA" + "C" * 30,  # 5' UTR = 48 nt (16 codons)
            'tis_start': 48,
            'cds_len': 9,
            'evidence_tier': 3,
        }
    ]
    tokenizer = Tokenizer('codon')
    window_cfg = WindowConfig(utr5_len=48, utr3_len=30)
    dataset = MultiTrackDataset(records, tokenizer, window_config=window_cfg)

    s = dataset[0]
    # Total nucleotides = 48 + 9 + 30 = 87 nt = 29 codons
    assert s.sequence_tokens.shape[0] == 29
    assert s.segment_ids.shape[0] == 29
    assert s.frame_offsets.shape[0] == 29
    assert s.ribo_coverage.shape[0] == 29
    assert s.evidence_tier == 3


def test_collate_multitrack_batch():
    s1 = MultiTrackSample(
        sequence_tokens=torch.tensor([1, 2, 3]),
        segment_ids=torch.tensor([0, 1, 2]),
        frame_offsets=torch.tensor([-1, 0, 1]),
        ribo_coverage=torch.tensor([0.0, 1.0, 2.0]),
        evidence_tier=1,
        initiation_weight=1.0,
    )
    s2 = MultiTrackSample(
        sequence_tokens=torch.tensor([4, 5]),
        segment_ids=torch.tensor([0, 2]),
        frame_offsets=torch.tensor([-1, 0]),
        ribo_coverage=torch.tensor([0.5, 0.5]),
        evidence_tier=4,
        initiation_weight=0.3,
    )

    batch = collate_multitrack_batch([s1, s2], pad_token_id=0)

    assert batch['sequence_tokens'].shape == (2, 3)
    assert batch['segment_ids'].shape == (2, 3)
    assert batch['frame_offsets'].shape == (2, 3)
    assert batch['ribo_coverage'].shape == (2, 3)
    assert batch['evidence_tiers'].shape == (2,)
    assert batch['initiation_weights'].shape == (2,)

    # Check padding values for second sample (length 2 padded to 3)
    assert batch['sequence_tokens'][1, 2].item() == 0
    assert batch['segment_ids'][1, 2].item() == 0
    assert batch['frame_offsets'][1, 2].item() == -1
    assert batch['ribo_coverage'][1, 2].item() == 0.0
    assert batch['padding_mask'][1, 2].item() is False
    assert batch['padding_mask'][0, 2].item() is True
