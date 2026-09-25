"""BioTokenizer Special Token Framework with control tokens, IUPAC degenerate bases, and legacy compatibility."""

import re
from itertools import product

# Explicit Control Tokens
PAD_TOKEN = "[PAD]"
UNK_TOKEN = "[UNK]"
BOS_TOKEN = "[BOS]"
EOS_TOKEN = "[EOS]"
MASK_TOKEN = "[MASK]"
FRAME_0_TOKEN = "[FRAME_0]"
FRAME_1_TOKEN = "[FRAME_1]"
FRAME_2_TOKEN = "[FRAME_2]"

CONTROL_TOKENS = [
    PAD_TOKEN,
    UNK_TOKEN,
    BOS_TOKEN,
    EOS_TOKEN,
    MASK_TOKEN,
    FRAME_0_TOKEN,
    FRAME_1_TOKEN,
    FRAME_2_TOKEN,
]

FRAME_TOKENS = {
    0: FRAME_0_TOKEN,
    1: FRAME_1_TOKEN,
    2: FRAME_2_TOKEN,
}

# IUPAC Degenerate Base Definitions
IUPAC_DEGENERATE_BASES = ["N", "R", "Y"]
STANDARD_BASES = ["A", "U", "C", "G"]

IUPAC_MAP = {
    "A": {"A"},
    "U": {"U"},
    "C": {"C"},
    "G": {"G"},
    "T": {"U"},
    "N": {"A", "U", "C", "G"},
    "R": {"A", "G"},
    "Y": {"C", "U"},
}


class BioTokenizer:
    """BioTokenizer supporting control tokens, IUPAC degenerate bases, and flexible encoding/decoding."""

    def __init__(
        self,
        mode="base",
        include_special_tokens=True,
        include_degenerate=True,
        handle_dna=True,
        strict_legacy=False,
    ):
        if mode not in ("base", "codon"):
            raise ValueError("mode must be base or codon")
        self.mode = mode
        self.width = 1 if mode == "base" else 3
        self.include_special_tokens = include_special_tokens
        self.include_degenerate = include_degenerate
        self.handle_dna = handle_dna
        self.strict_legacy = strict_legacy

        self.vocabulary = []

        # 1. Special Control Tokens
        if self.include_special_tokens:
            self.vocabulary.extend(CONTROL_TOKENS)

        # 2. Standard Vocabularies
        self.standard_bases = list(STANDARD_BASES)
        self.standard_codons = ["".join(c) for c in product(STANDARD_BASES, repeat=3)]

        if mode == "base":
            self.vocabulary.extend(self.standard_bases)
            if self.include_degenerate:
                self.vocabulary.extend(IUPAC_DEGENERATE_BASES)
        else:  # codon mode
            self.vocabulary.extend(self.standard_codons)
            if self.include_degenerate:
                all_iupac_bases = STANDARD_BASES + IUPAC_DEGENERATE_BASES
                all_triplets = ["".join(c) for c in product(all_iupac_bases, repeat=3)]
                degenerate_codons = [c for c in all_triplets if c not in set(self.standard_codons)]
                self.vocabulary.extend(degenerate_codons)

        self.token_to_id = {token: i for i, token in enumerate(self.vocabulary)}
        self.id_to_token = {i: token for i, token in enumerate(self.vocabulary)}

        # Helper ID lookups
        self.pad_id = self.token_to_id.get(PAD_TOKEN)
        self.unk_id = self.token_to_id.get(UNK_TOKEN)
        self.bos_id = self.token_to_id.get(BOS_TOKEN)
        self.eos_id = self.token_to_id.get(EOS_TOKEN)
        self.mask_id = self.token_to_id.get(MASK_TOKEN)
        self.frame_0_id = self.token_to_id.get(FRAME_0_TOKEN)
        self.frame_1_id = self.token_to_id.get(FRAME_1_TOKEN)
        self.frame_2_id = self.token_to_id.get(FRAME_2_TOKEN)

    def is_control_token(self, token_or_id):
        if isinstance(token_or_id, int):
            if 0 <= token_or_id < len(self.vocabulary):
                token = self.id_to_token[token_or_id]
            else:
                return False
        else:
            token = str(token_or_id)
        return token in CONTROL_TOKENS

    def is_degenerate(self, token_or_id):
        if isinstance(token_or_id, int):
            if 0 <= token_or_id < len(self.vocabulary):
                token = self.id_to_token[token_or_id]
            else:
                return False
        else:
            token = str(token_or_id)

        if token in CONTROL_TOKENS:
            return False
        return any(b in IUPAC_DEGENERATE_BASES for b in token)

    def expand_degenerate(self, token_or_sequence):
        """Expands a degenerate token or sequence into all standard base or codon permutations."""
        if isinstance(token_or_sequence, list):
            expansions = [self.expand_degenerate(item) for item in token_or_sequence]
            return ["".join(p) for p in product(*expansions)]

        token = token_or_sequence
        if isinstance(token, int):
            token = self.id_to_token[token]

        if token in CONTROL_TOKENS:
            return [token]

        char_sets = [IUPAC_MAP.get(c, {c}) for c in token]
        expanded = ["".join(p) for p in product(*char_sets)]
        return expanded

    def encode(self, rna, add_bos=False, add_eos=False, frame=None):
        if self.strict_legacy:
            if not rna or not isinstance(rna, str) or set(rna) - set('AUCG'):
                raise ValueError('Expected nonempty uppercase RNA containing only A,U,C,G')
            if len(rna) % self.width:
                raise ValueError('Incomplete codon')
            return [self.token_to_id[rna[i:i+self.width]] for i in range(0, len(rna), self.width)]

        if rna is None:
            raise ValueError('Input sequence cannot be None')

        tokens = []
        if frame is not None:
            if frame not in (0, 1, 2):
                raise ValueError('frame must be 0, 1, or 2')
            tokens.append(FRAME_TOKENS[frame])

        if add_bos:
            tokens.append(BOS_TOKEN)

        if isinstance(rna, list):
            tokens.extend(rna)
        elif isinstance(rna, str):
            pattern = r"(\[[A-Z0-9_]+\])"
            parts = re.split(pattern, rna)
            for part in parts:
                if not part:
                    continue
                if part in CONTROL_TOKENS:
                    tokens.append(part)
                elif part.startswith("[") and part.endswith("]"):
                    if self.unk_id is not None:
                        tokens.append(UNK_TOKEN)
                    else:
                        raise ValueError(f"Unknown control token: {part}")
                else:
                    seq = part
                    if self.handle_dna:
                        seq = seq.replace("T", "U")
                    if len(seq) % self.width != 0:
                        raise ValueError(f"Incomplete sequence chunk: length {len(seq)} not divisible by width {self.width}")
                    for i in range(0, len(seq), self.width):
                        tokens.append(seq[i:i+self.width])

        if add_eos:
            tokens.append(EOS_TOKEN)

        encoded_ids = []
        for tok in tokens:
            if tok in self.token_to_id:
                encoded_ids.append(self.token_to_id[tok])
            elif self.unk_id is not None:
                encoded_ids.append(self.unk_id)
            else:
                raise ValueError(f"Token '{tok}' not in vocabulary")

        return encoded_ids

    def decode(self, ids, skip_special_tokens=False):
        if self.strict_legacy:
            if any(not isinstance(i, int) or i < 0 or i >= len(self.vocabulary) for i in ids):
                raise ValueError('Token ID outside vocabulary')
            return ''.join(self.vocabulary[i] for i in ids)

        res = []
        for i in ids:
            if not isinstance(i, int) or i < 0 or i >= len(self.vocabulary):
                raise ValueError(f'Token ID {i} outside vocabulary range (0-{len(self.vocabulary)-1})')
            tok = self.vocabulary[i]
            if skip_special_tokens and tok in CONTROL_TOKENS:
                continue
            res.append(tok)
        return ''.join(res)


class Tokenizer(BioTokenizer):
    """Legacy tokenizer interface maintaining exact 4-base / 64-codon fixed vocabularies."""

    def __init__(self, mode):
        super().__init__(
            mode=mode,
            include_special_tokens=False,
            include_degenerate=False,
            handle_dna=False,
            strict_legacy=True,
        )


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
    "STANDARD_BASES",
    "IUPAC_MAP",
]
