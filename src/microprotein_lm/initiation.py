"""Non-canonical initiation codon (TIS) dynamics and empirical prior tensors."""

from typing import Dict, List, Optional
import torch

# Empirical initiation efficiency priors mapping start codons to initiation weights
DEFAULT_INITIATION_PRIORS: Dict[str, float] = {
    "AUG": 1.0,
    "CUG": 0.6,
    "GUG": 0.5,
    "UUG": 0.4,
    "ACG": 0.3,
    "AUU": 0.2,
}

NEAR_COGNATE_START_CODONS: List[str] = list(DEFAULT_INITIATION_PRIORS.keys())


class InitiationPrior:
    """Manages empirical initiation efficiency priors (E_init) for canonical and near-cognate start codons."""

    def __init__(self, prior_table: Optional[Dict[str, float]] = None, default_weight: float = 0.0):
        self.priors = dict(DEFAULT_INITIATION_PRIORS)
        if prior_table:
            self.priors.update(prior_table)
        self.default_weight = default_weight

    def get_weight(self, codon: str) -> float:
        """Returns empirical initiation weight for a given codon string."""
        codon = codon.upper().replace('T', 'U')
        return self.priors.get(codon, self.default_weight)

    def get_prior_tensor(self, vocabulary: List[str]) -> torch.Tensor:
        """Creates E_init prior tensor aligned with a tokenizer vocabulary list.

        Args:
            vocabulary: List of token strings (e.g. 64 codons or 4 bases).

        Returns:
            1D PyTorch FloatTensor of shape (len(vocabulary),) containing initiation weights.
        """
        weights = [self.get_weight(tok) for tok in vocabulary]
        return torch.tensor(weights, dtype=torch.float32)

    def calculate_frame_offsets(
        self,
        seq_len: int,
        tis_offset: int,
        token_width: int = 1
    ) -> List[int]:
        """Calculates reading-phase coordinate tracks relative to active TIS.

        Args:
            seq_len: Total length of token sequence.
            tis_offset: Token index where active TIS starts.
            token_width: 1 for nucleotide/base level, 3 for codon level.

        Returns:
            List of frame offset integers (-1 for UTR/out-of-frame, 0, 1, 2 for reading frames).
        """
        frame_offsets = []
        for i in range(seq_len):
            if token_width == 3:
                # Codon mode: i is codon index
                rel_codon = i - tis_offset
                if rel_codon < 0:
                    frame_offsets.append(-1)  # 5' UTR / untranslated
                else:
                    frame_offsets.append(0)   # In-frame codon position 0
            else:
                # Base mode: i is nucleotide index
                rel_base = i - tis_offset
                if rel_base < 0:
                    frame_offsets.append(-1)
                else:
                    frame_offsets.append(rel_base % 3)
        return frame_offsets
