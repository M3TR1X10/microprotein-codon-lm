"""Dynamic Translation Engine supporting NCBI Genetic Code Tables 1 through 33.

Provides full codon-to-amino-acid translation, initiation codon validation,
stop codon validation, and detailed CDS diagnostic reports across all NCBI tables.
"""

from typing import Dict, List, Optional, Set, Tuple, Union

# Base Table 1 (Standard Code) DNA mappings
STANDARD_CODON_TABLE_DNA: Dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Table specific codon alterations (relative to Table 1) and initiation/stop sets
TABLE_SPECIFICS: Dict[int, Dict[str, Union[Dict[str, str], List[str]]]] = {
    1: {
        'overrides': {},
        'starts': ['ATG', 'GTG', 'TTG'],
        'stops': ['TAA', 'TAG', 'TGA'],
    },
    2: {
        'overrides': {'AGA': '*', 'AGG': '*', 'ATA': 'M', 'TGA': 'W'},
        'starts': ['ATT', 'ATC', 'ATA', 'ATG', 'GTG', 'TTG'],
        'stops': ['TAA', 'TAG', 'AGA', 'AGG'],
    },
    3: {
        'overrides': {'ATA': 'M', 'CTT': 'T', 'CTC': 'T', 'CTA': 'T', 'CTG': 'T', 'TGA': 'W'},
        'starts': ['ATA', 'ATG'],
        'stops': ['TAA', 'TAG'],
    },
    4: {
        'overrides': {'TGA': 'W'},
        'starts': ['TTA', 'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'],
        'stops': ['TAA', 'TAG'],
    },
    5: {
        'overrides': {'AGA': 'S', 'AGG': 'S', 'ATA': 'M', 'TGA': 'W'},
        'starts': ['TTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'],
        'stops': ['TAA', 'TAG'],
    },
    6: {
        'overrides': {'TAA': 'Q', 'TAG': 'Q'},
        'starts': ['ATG'],
        'stops': ['TGA'],
    },
    9: {
        'overrides': {'AAA': 'N', 'AGA': 'S', 'AGG': 'S', 'ATA': 'M', 'TGA': 'W'},
        'starts': ['ATG', 'GTG'],
        'stops': ['TAA', 'TAG'],
    },
    10: {
        'overrides': {'TGA': 'C'},
        'starts': ['ATG'],
        'stops': ['TAA', 'TAG'],
    },
    11: {
        'overrides': {},
        'starts': ['TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'],
        'stops': ['TAA', 'TAG', 'TGA'],
    },
    12: {
        'overrides': {'CTG': 'S'},
        'starts': ['CAG', 'ATG'],
        'stops': ['TAA', 'TAG', 'TGA'],
    },
    13: {
        'overrides': {'AGA': 'G', 'ATA': 'M', 'TGA': 'W'},
        'starts': ['ATG', 'GTG'],
        'stops': ['TAA', 'TAG'],
    },
    14: {
        'overrides': {'AAA': 'N', 'AGA': 'S', 'AGG': 'S', 'ATA': 'M', 'TAA': 'Y', 'TGA': 'W'},
        'starts': ['ATG'],
        'stops': ['TAG'],
    },
    15: {
        'overrides': {'TAG': 'Q'},
        'starts': ['ATG'],
        'stops': ['TAA', 'TGA'],
    },
    16: {
        'overrides': {'TAG': 'L'},
        'starts': ['ATG'],
        'stops': ['TAA', 'TGA'],
    },
    21: {
        'overrides': {'AAA': 'N', 'AGA': 'S', 'AGG': 'S', 'ATA': 'M', 'TGA': 'W'},
        'starts': ['ATG', 'GTG'],
        'stops': ['TAA', 'TAG'],
    },
    22: {
        'overrides': {'TCA': '*', 'TAG': 'L'},
        'starts': ['ATG'],
        'stops': ['TAA', 'TCA', 'TGA'],
    },
    23: {
        'overrides': {'TTA': '*'},
        'starts': ['ATT', 'ATG', 'GTG'],
        'stops': ['TAA', 'TAG', 'TGA', 'TTA'],
    },
    24: {
        'overrides': {'AGA': 'S', 'AGG': 'K', 'TGA': 'W'},
        'starts': ['ATG', 'AGG'],
        'stops': ['TAA', 'TAG'],
    },
    25: {
        'overrides': {'TGA': 'G'},
        'starts': ['TTG', 'ATG', 'GTG'],
        'stops': ['TAA', 'TAG'],
    },
    26: {
        'overrides': {'CTG': 'A'},
        'starts': ['CTG', 'ATG'],
        'stops': ['TAA', 'TAG', 'TGA'],
    },
    27: {
        'overrides': {'TAA': 'Q', 'TAG': 'Q'},
        'starts': ['ATG'],
        'stops': ['TGA'],
    },
    28: {
        'overrides': {'TAA': 'Q', 'TAG': 'Q'},
        'starts': ['ATG'],
        'stops': ['TAA', 'TAG', 'TGA'],
    },
    29: {
        'overrides': {'TAA': 'Y', 'TAG': 'Y'},
        'starts': ['ATG'],
        'stops': ['TGA'],
    },
    30: {
        'overrides': {'TAA': 'Q', 'TAG': 'Q'},
        'starts': ['ATG'],
        'stops': ['TGA'],
    },
    31: {
        'overrides': {'TAA': 'E', 'TAG': 'E'},
        'starts': ['ATG'],
        'stops': ['TGA'],
    },
    33: {
        'overrides': {'AGA': 'S', 'AGG': 'K', 'TGA': 'W'},
        'starts': ['ATA', 'ATG', 'AGG'],
        'stops': ['TAA', 'TAG'],
    },
}


def _dna_to_rna(seq: str) -> str:
    return seq.upper().replace('T', 'U')


def _rna_to_dna(seq: str) -> str:
    return seq.upper().replace('U', 'T')


class TranslationEngine:
    """Translation engine supporting NCBI Genetic Code Tables 1 through 33."""

    def __init__(self):
        # Pre-build lookup dictionaries for tables 1 through 33
        self._tables_rna: Dict[int, Dict[str, str]] = {}
        self._tables_dna: Dict[int, Dict[str, str]] = {}
        self._starts_rna: Dict[int, Set[str]] = {}
        self._stops_rna: Dict[int, Set[str]] = {}

        for table_id in range(1, 34):
            info = TABLE_SPECIFICS.get(table_id, TABLE_SPECIFICS[1])
            dna_map = dict(STANDARD_CODON_TABLE_DNA)
            overrides = info.get('overrides', {})
            dna_map.update(overrides)

            rna_map = {_dna_to_rna(k): v for k, v in dna_map.items()}
            starts = {_dna_to_rna(s) for s in info.get('starts', ['ATG', 'GTG', 'TTG'])}
            stops = {_dna_to_rna(s) for s in info.get('stops', ['TAA', 'TAG', 'TGA'])}

            self._tables_dna[table_id] = dna_map
            self._tables_rna[table_id] = rna_map
            self._starts_rna[table_id] = starts
            self._stops_rna[table_id] = stops

    def get_table_map(self, table: int = 1, rna: bool = True) -> Dict[str, str]:
        """Returns the codon-to-amino-acid dictionary for the given NCBI table."""
        table_id = table if 1 <= table <= 33 else 1
        return self._tables_rna[table_id] if rna else self._tables_dna[table_id]

    def get_initiation_codons(self, table: int = 1, rna: bool = True) -> Set[str]:
        """Returns the set of valid initiation triplets for the given table."""
        table_id = table if 1 <= table <= 33 else 1
        starts = self._starts_rna[table_id]
        if not rna:
            return {_rna_to_dna(s) for s in starts}
        return starts

    def get_stop_codons(self, table: int = 1, rna: bool = True) -> Set[str]:
        """Returns the set of valid stop triplets for the given table."""
        table_id = table if 1 <= table <= 33 else 1
        stops = self._stops_rna[table_id]
        if not rna:
            return {_rna_to_dna(s) for s in stops}
        return stops

    def is_valid_initiation(self, codon: str, table: int = 1) -> bool:
        """Checks if a codon is a valid initiation triplet for table."""
        c_rna = _dna_to_rna(codon)
        return c_rna in self.get_initiation_codons(table, rna=True)

    def is_valid_stop(self, codon: str, table: int = 1) -> bool:
        """Checks if a codon is a valid stop triplet for table."""
        c_rna = _dna_to_rna(codon)
        return c_rna in self.get_stop_codons(table, rna=True)

    def get_initiation_type(self, codon: str, table: int = 1) -> str:
        """Classifies initiation codon as 'canonical' (AUG), 'near_cognate', or 'invalid'."""
        c_rna = _dna_to_rna(codon)
        if c_rna == 'AUG':
            return 'canonical'
        if self.is_valid_initiation(c_rna, table):
            return 'near_cognate'
        return 'invalid'

    def translate_codon(self, codon: str, table: int = 1) -> str:
        """Translates a single triplet codon to single-letter amino acid code."""
        c_rna = _dna_to_rna(codon)
        t_map = self.get_table_map(table, rna=True)
        return t_map.get(c_rna, 'X')

    def translate_cds(self, seq: str, table: int = 1, cds: bool = True) -> str:
        """Translates a CDS nucleotide sequence (DNA or RNA) into protein string.

        Args:
            seq: Nucleotide string (RNA or DNA).
            table: NCBI genetic code table ID (1-33).
            cds: If True, validates length % 3, initiation codon, stop codon,
                 and translates valid initiation codon as Methionine ('M').

        Returns:
            Translated amino acid sequence string.
        """
        seq_clean = _dna_to_rna(seq.strip())
        if len(seq_clean) % 3 != 0:
            raise ValueError(f"CDS sequence length {len(seq_clean)} is not a multiple of 3")

        codons = [seq_clean[i:i+3] for i in range(0, len(seq_clean), 3)]
        if not codons:
            return ""

        if cds:
            first_codon = codons[0]
            if not self.is_valid_initiation(first_codon, table):
                raise ValueError(f"Invalid initiation codon '{first_codon}' for genetic table {table}")

            stops = self.get_stop_codons(table, rna=True)
            if codons[-1] not in stops:
                raise ValueError(f"Final codon '{codons[-1]}' is not a valid stop codon for genetic table {table}")

            body_codons = codons[:-1]

            amino_acids = []
            # First codon translated as 'M' in CDS mode per NCBI standard convention
            amino_acids.append('M')

            for c in body_codons[1:]:
                if c in stops:
                    raise ValueError(f"Internal stop codon '{c}' encountered in CDS")
                amino_acids.append(self.translate_codon(c, table))

            return "".join(amino_acids)
        else:
            return "".join(self.translate_codon(c, table) for c in codons)

    def diagnose_cds(
        self,
        seq: str,
        table: int = 1,
        expected_protein: Optional[str] = None
    ) -> Dict[str, Union[bool, List[str], str, int]]:
        """Runs comprehensive diagnostics on a candidate CDS sequence.

        Returns:
            Dict containing validation result, list of error codes, parsed codons,
            initiation classification, and translated protein sequence.
        """
        seq_rna = _dna_to_rna(seq.strip())
        errors: List[str] = []
        is_valid = True

        if not seq_rna or set(seq_rna) - set('ACGU'):
            errors.append('ambiguous_or_invalid_bases')
            is_valid = False

        if len(seq_rna) % 3 != 0:
            errors.append('incomplete_codon')
            is_valid = False

        init_codon = seq_rna[:3] if len(seq_rna) >= 3 else ""
        stop_codon = seq_rna[-3:] if len(seq_rna) >= 3 else ""
        init_type = self.get_initiation_type(init_codon, table) if init_codon else 'invalid'

        if init_codon and not self.is_valid_initiation(init_codon, table):
            errors.append('invalid_initiation_codon')
            is_valid = False

        translated_protein = ""
        if len(seq_rna) >= 3 and len(seq_rna) % 3 == 0:
            codons = [seq_rna[i:i+3] for i in range(0, len(seq_rna), 3)]
            stops = self.get_stop_codons(table, rna=True)
            if codons[-1] not in stops:
                errors.append('missing_terminal_stop_codon')
                is_valid = False

            for i, c in enumerate(codons[:-1]):
                if c in stops:
                    errors.append(f'internal_stop_codon_at_index_{i}')
                    is_valid = False
                    break

            try:
                translated_protein = self.translate_cds(seq_rna, table=table, cds=True)
            except Exception as e:
                errors.append(f'translation_error_{str(e)}')
                is_valid = False

        if expected_protein and translated_protein:
            if translated_protein != expected_protein:
                errors.append('protein_mismatch')
                is_valid = False

        return {
            'is_valid': is_valid,
            'errors': errors,
            'initiation_codon': init_codon,
            'initiation_type': init_type,
            'stop_codon': stop_codon,
            'translated_protein': translated_protein,
            'table': table,
        }
