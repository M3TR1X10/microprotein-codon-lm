"""Multi-track dataset representation and batch collation."""

from dataclasses import dataclass
from typing import List, Dict, Any, Optional
import torch
from torch.utils.data import Dataset

from .tokenization import Tokenizer
from .extended_window import ExtendedWindowExtractor, WindowConfig
from .initiation import InitiationPrior


@dataclass
class MultiTrackSample:
    """Container for a multi-track sequence sample."""
    sequence_tokens: torch.Tensor  # Track 1: LongTensor (L,)
    segment_ids: torch.Tensor      # Track 2: LongTensor (L,) - 0=5' UTR, 1=Kozak, 2=CDS, 3=3' UTR
    frame_offsets: torch.Tensor    # Track 3: LongTensor (L,) - -1, 0, 1, 2
    ribo_coverage: torch.Tensor    # Track 4: FloatTensor (L,) - Ribo-seq P-site coverage vector
    evidence_tier: int = 1         # Tier 1..4 (for evidence-weighted loss)
    initiation_weight: float = 1.0 # Empirical E_init prior weight for active TIS


class MultiTrackDataset(Dataset):
    """Dataset producing multi-track tensor representations for microprotein sequences."""

    def __init__(
        self,
        records: List[Dict[str, Any]],
        tokenizer: Tokenizer,
        window_config: Optional[WindowConfig] = None,
        prior_table: Optional[Dict[str, float]] = None
    ):
        """
        Args:
            records: List of dicts, each containing:
                - 'rna': Full or transcript RNA sequence
                - 'tis_start': Start offset of TIS in rna (default 0 if window pre-sliced)
                - 'cds_len': CDS length in nucleotides
                - 'ribo_coverage': List/array of floats for Ribo-seq coverage (optional)
                - 'evidence_tier': Tier int 1..4 (optional, default 1)
            tokenizer: Tokenizer instance (mode='base' or 'codon')
            window_config: Optional WindowConfig instance
            prior_table: Optional dict of custom TIS prior weights
        """
        self.records = records
        self.tokenizer = tokenizer
        self.window_extractor = ExtendedWindowExtractor(window_config)
        self.initiation_prior = InitiationPrior(prior_table)
        self.samples = [self._process_record(r) for r in records]

    def _process_record(self, record: Dict[str, Any]) -> MultiTrackSample:
        rna = record['rna'].upper().replace('T', 'U')
        tis_start = record.get('tis_start', 0)
        cds_len = record.get('cds_len', len(rna) - tis_start)

        # Extract extended window
        extracted = self.window_extractor.extract_window(rna, tis_start, cds_len)
        ext_seq = extracted['extended_seq']
        raw_seg_ids = extracted['segment_ids']
        tis_offset = extracted['tis_offset_in_extended']

        # Determine start codon & E_init weight
        start_codon = ext_seq[tis_offset:tis_offset + 3] if len(ext_seq) >= tis_offset + 3 else "AUG"
        init_weight = self.initiation_prior.get_weight(start_codon)

        # Tokenization & track construction based on mode
        if self.tokenizer.mode == 'base':
            tokens = self.tokenizer.encode(ext_seq)
            seq_tokens = torch.tensor(tokens, dtype=torch.long)
            seg_ids = torch.tensor(raw_seg_ids, dtype=torch.long)
            frame_offsets = torch.tensor(
                self.initiation_prior.calculate_frame_offsets(len(tokens), tis_offset, token_width=1),
                dtype=torch.long
            )

            # Ribo-seq coverage
            raw_ribo = record.get('ribo_coverage', None)
            if raw_ribo is None or len(raw_ribo) != len(ext_seq):
                ribo = torch.zeros(len(ext_seq), dtype=torch.float32)
            else:
                ribo = torch.tensor(raw_ribo, dtype=torch.float32)

        else:  # codon mode
            # Ensure length is divisible by 3
            rem = len(ext_seq) % 3
            if rem != 0:
                ext_seq = ext_seq[:-rem]
                raw_seg_ids = raw_seg_ids[:-rem]

            tokens = self.tokenizer.encode(ext_seq)
            seq_tokens = torch.tensor(tokens, dtype=torch.long)

            # Downsample segment IDs to codon level (take segment ID of first nucleotide in codon)
            codon_seg_ids = [raw_seg_ids[i] for i in range(0, len(raw_seg_ids), 3)]
            seg_ids = torch.tensor(codon_seg_ids, dtype=torch.long)

            codon_tis_offset = tis_offset // 3
            frame_offsets = torch.tensor(
                self.initiation_prior.calculate_frame_offsets(len(tokens), codon_tis_offset, token_width=3),
                dtype=torch.long
            )

            raw_ribo = record.get('ribo_coverage', None)
            if raw_ribo is None or len(raw_ribo) != len(ext_seq):
                ribo = torch.zeros(len(tokens), dtype=torch.float32)
            else:
                # Average or sum Ribo coverage per codon
                codon_ribo = [sum(raw_ribo[i:i+3]) / 3.0 for i in range(0, len(raw_ribo) - rem, 3)]
                ribo = torch.tensor(codon_ribo, dtype=torch.float32)

        evidence_tier = record.get('evidence_tier', 1)

        return MultiTrackSample(
            sequence_tokens=seq_tokens,
            segment_ids=seg_ids,
            frame_offsets=frame_offsets,
            ribo_coverage=ribo,
            evidence_tier=evidence_tier,
            initiation_weight=init_weight
        )

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> MultiTrackSample:
        return self.samples[idx]


def collate_multitrack_batch(samples: List[MultiTrackSample], pad_token_id: int = 0) -> Dict[str, torch.Tensor]:
    """Collates a list of MultiTrackSample objects into batched 2D multi-track tensors.

    Args:
        samples: List of MultiTrackSample instances.
        pad_token_id: Token ID used for padding sequence tokens.

    Returns:
        Dict with batched tensors:
            - 'sequence_tokens': LongTensor of shape (B, Max_L)
            - 'segment_ids': LongTensor of shape (B, Max_L), padded with 0
            - 'frame_offsets': LongTensor of shape (B, Max_L), padded with -1
            - 'ribo_coverage': FloatTensor of shape (B, Max_L), padded with 0.0
            - 'evidence_tiers': LongTensor of shape (B,)
            - 'initiation_weights': FloatTensor of shape (B,)
            - 'padding_mask': BoolTensor of shape (B, Max_L) (True where valid, False where padded)
    """
    batch_size = len(samples)
    max_len = max(s.sequence_tokens.shape[0] for s in samples)

    seq_batch = torch.full((batch_size, max_len), pad_token_id, dtype=torch.long)
    seg_batch = torch.zeros((batch_size, max_len), dtype=torch.long)
    frame_batch = torch.full((batch_size, max_len), -1, dtype=torch.long)
    ribo_batch = torch.zeros((batch_size, max_len), dtype=torch.float32)
    mask_batch = torch.zeros((batch_size, max_len), dtype=torch.bool)

    tiers = torch.tensor([s.evidence_tier for s in samples], dtype=torch.long)
    init_weights = torch.tensor([s.initiation_weight for s in samples], dtype=torch.float32)

    for i, sample in enumerate(samples):
        l = sample.sequence_tokens.shape[0]
        seq_batch[i, :l] = sample.sequence_tokens
        seg_batch[i, :l] = sample.segment_ids
        frame_batch[i, :l] = sample.frame_offsets
        ribo_batch[i, :l] = sample.ribo_coverage
        mask_batch[i, :l] = True

    return {
        'sequence_tokens': seq_batch,
        'segment_ids': seg_batch,
        'frame_offsets': frame_batch,
        'ribo_coverage': ribo_batch,
        'evidence_tiers': tiers,
        'initiation_weights': init_weights,
        'padding_mask': mask_batch,
    }
