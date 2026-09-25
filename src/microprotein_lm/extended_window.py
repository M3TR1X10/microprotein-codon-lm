"""Extended sequence window utilities (5' UTR + Kozak + CDS + 3' UTR)."""

from dataclasses import dataclass
from typing import Dict, List, Any, Optional
import torch

from .tokenization import Tokenizer


@dataclass
class WindowConfig:
    utr5_len: int = 50
    kozak_upstream: int = 6
    kozak_downstream: int = 4  # -6 to +4 nt around TIS (length 10 nt)
    utr3_len: int = 30
    pad_char: str = 'A'


class ExtendedWindowExtractor:
    """Extracts and pads extended sequence windows including 5' UTR, Kozak, CDS, and 3' UTR."""

    def __init__(self, config: Optional[WindowConfig] = None):
        self.config = config or WindowConfig()

    def extract_window(
        self,
        transcript: str,
        tis_start: int,
        cds_len: int
    ) -> Dict[str, Any]:
        """Extract region slices from a full transcript sequence given TIS start index and CDS length.

        Args:
            transcript: Complete RNA nucleotide sequence (A, U, C, G).
            tis_start: 0-indexed position where the Translation Initiation Site (TIS) begins.
            cds_len: Length of the CDS in nucleotides.

        Returns:
            Dictionary containing:
                'utr5': 5' UTR sequence string (length config.utr5_len)
                'kozak': Kozak window sequence string (-6 to +4 around TIS)
                'cds': Coding sequence string
                'utr3': 3' UTR sequence string (length config.utr3_len)
                'extended_seq': Combined extended sequence string (5' UTR + CDS + 3' UTR)
                'segment_ids': List of integer segment IDs corresponding to each nucleotide:
                                0 = 5' UTR, 1 = Kozak, 2 = CDS, 3 = 3' UTR
                'tis_offset_in_extended': Index of TIS start in extended_seq
        """
        transcript = transcript.upper().replace('T', 'U')
        cfg = self.config

        # 5' UTR extraction with padding if needed
        utr5_start = tis_start - cfg.utr5_len
        if utr5_start < 0:
            pad_left = cfg.pad_char * (-utr5_start)
            utr5_seq = pad_left + transcript[0:tis_start]
        else:
            utr5_seq = transcript[utr5_start:tis_start]

        # Kozak extraction (-6 to +4 nt around tis_start)
        k_start = tis_start - cfg.kozak_upstream
        k_end = tis_start + cfg.kozak_downstream
        if k_start < 0:
            kozak_seq = (cfg.pad_char * (-k_start)) + transcript[0:min(k_end, len(transcript))]
        else:
            kozak_seq = transcript[k_start:min(k_end, len(transcript))]
            if len(kozak_seq) < (cfg.kozak_upstream + cfg.kozak_downstream):
                kozak_seq += cfg.pad_char * ((cfg.kozak_upstream + cfg.kozak_downstream) - len(kozak_seq))

        # CDS extraction
        cds_end = tis_start + cds_len
        cds_seq = transcript[tis_start:cds_end]

        # 3' UTR extraction with padding if needed
        utr3_end = cds_end + cfg.utr3_len
        if utr3_end > len(transcript):
            utr3_seq = transcript[cds_end:len(transcript)]
            utr3_seq += cfg.pad_char * (cfg.utr3_len - len(utr3_seq))
        else:
            utr3_seq = transcript[cds_end:utr3_end]

        # Extended sequence assembly
        extended_seq = utr5_seq + cds_seq + utr3_seq
        tis_offset = len(utr5_seq)

        # Build nucleotide-level segment IDs
        # 0 = 5' UTR, 1 = Kozak, 2 = CDS, 3 = 3' UTR
        segment_ids = []

        # 5' UTR positions
        for i in range(len(utr5_seq)):
            rel_pos = i - tis_offset
            if -cfg.kozak_upstream <= rel_pos < 0:
                segment_ids.append(1)  # Kozak
            else:
                segment_ids.append(0)  # 5' UTR

        # CDS positions
        for i in range(len(cds_seq)):
            rel_pos = i  # relative to TIS start in CDS
            if rel_pos < cfg.kozak_downstream:
                segment_ids.append(1)  # Kozak
            else:
                segment_ids.append(2)  # CDS

        # 3' UTR positions
        for _ in range(len(utr3_seq)):
            segment_ids.append(3)  # 3' UTR

        return {
            'utr5': utr5_seq,
            'kozak': kozak_seq,
            'cds': cds_seq,
            'utr3': utr3_seq,
            'extended_seq': extended_seq,
            'segment_ids': segment_ids,
            'tis_offset_in_extended': tis_offset,
        }
