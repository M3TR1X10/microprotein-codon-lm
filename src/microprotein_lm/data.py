import json
from pathlib import Path

import numpy as np
import torch

from .io import read_json, write_json, sha256
from .tokenization import Tokenizer


def prepare(processed='data/processed'):
    root = Path(processed)
    manifest = read_json(root/'manifest.json')
    if sha256(root/'cohort.jsonl') != manifest['cohort_sha256']:
        raise ValueError('Cohort changed after selection')
    records = [json.loads(line) for line in (root/'cohort.jsonl').read_text(encoding='utf-8').splitlines()]
    if not records or len({r['sequence_sha256'] for r in records}) != len(records):
        raise ValueError('Cohort must be nonempty and contain unique sequences')
    for mode in ('base','codon'):
        tokenizer = Tokenizer(mode)
        arrays = [np.asarray(tokenizer.encode(r['rna']),dtype=np.uint8) for r in records]
        lengths = [len(x) for x in arrays]
        if min(lengths) < (4 if mode == 'base' else 2):
            raise ValueError('Need a first codon and at least one target codon')
        path = root/mode
        path.mkdir(exist_ok=True)
        np.concatenate(arrays).tofile(path/'tokens.bin')
        np.save(path/'offsets.npy',np.asarray([0,*np.cumsum(lengths)],dtype=np.int64))
        write_json(path/'metadata.json',{'mode':mode,'vocabulary':tokenizer.vocabulary,
            'vocab_size':len(tokenizer.vocabulary),'sequences':len(records),
            'cohort_sha256':manifest['cohort_sha256'], 'split':'train',
            'tokens_sha256':sha256(path/'tokens.bin'), 'offsets_sha256':sha256(path/'offsets.npy'),
            'tokens':sum(lengths), 'max_sequence_tokens':max(lengths),
            'predicted_region':'All nucleotides after first complete codon, including terminal stop',
            'sequence_ids':[r['sequence_sha256'] for r in records]})


class Corpus:
    def __init__(self, root, mode):
        self.mode = mode
        self.path = Path(root)/mode
        self.metadata = read_json(self.path/'metadata.json')
        self.tokenizer = Tokenizer(mode)
        if self.metadata['split'] != 'train' or self.metadata['vocabulary'] != self.tokenizer.vocabulary:
            raise ValueError('Unexpected split or vocabulary')
        for name,key in [('tokens.bin','tokens_sha256'),('offsets.npy','offsets_sha256')]:
            if sha256(self.path/name) != self.metadata[key]:
                raise ValueError('Encoded data checksum mismatch')
        self.tokens = np.memmap(self.path/'tokens.bin',dtype=np.uint8,mode='r')
        self.offsets = np.load(self.path/'offsets.npy',allow_pickle=False)
        if self.offsets[0] != 0 or self.offsets[-1] != len(self.tokens) or np.any(np.diff(self.offsets) <= 1):
            raise ValueError('Invalid sequence boundaries')

    def __len__(self):
        return len(self.offsets)-1

    def sequence(self, index):
        return self.tokens[self.offsets[index]:self.offsets[index+1]].astype(np.int64)

    def batch(self, indices, device='cpu'):
        sequences = [self.sequence(int(i)) for i in indices]
        length = max(len(s) for s in sequences)-1
        # Right-fill with existing token 0; it is never a vocabulary PAD token.
        # Causality keeps this fill from influencing earlier valid positions.
        x = torch.zeros((len(indices),length),dtype=torch.long)
        y = torch.full_like(x,-100)
        for i,seq in enumerate(sequences):
            x[i,:len(seq)-1] = torch.from_numpy(seq[:-1])
            y[i,:len(seq)-1] = torch.from_numpy(seq[1:])
            if self.mode == 'base':
                y[i,:2] = -100  # both models condition on the same first codon
        return x.to(device),y.to(device)

