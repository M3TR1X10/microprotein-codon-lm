"""Fixed vocabularies; no learned merges, special tokens, or unknown symbols."""
from itertools import product


class Tokenizer:
    def __init__(self, mode):
        if mode not in ('base', 'codon'):
            raise ValueError('mode must be base or codon')
        self.mode = mode
        self.width = 1 if mode == 'base' else 3
        self.vocabulary = list('AUCG') if mode == 'base' else [''.join(c) for c in product('AUCG', repeat=3)]
        self.token_to_id = {token:i for i,token in enumerate(self.vocabulary)}

    def encode(self, rna):
        if not rna or set(rna) - set('AUCG'):
            raise ValueError('Expected nonempty uppercase RNA containing only A,U,C,G')
        if len(rna) % self.width:
            raise ValueError('Incomplete codon')
        return [self.token_to_id[rna[i:i+self.width]] for i in range(0,len(rna),self.width)]

    def decode(self, ids):
        if any(not isinstance(i, int) or i < 0 or i >= len(self.vocabulary) for i in ids):
            raise ValueError('Token ID outside vocabulary')
        return ''.join(self.vocabulary[i] for i in ids)

