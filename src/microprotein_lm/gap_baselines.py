"""Codon baselines fit only to each frozen model's training cohort."""
from collections import Counter, defaultdict
import math

import numpy as np

from .gap_metrics import codons
from .tokenization import Tokenizer


METHODS = ('uniform', 'unigram', 'bigram', 'position', 'privileged_family_taxon_position')


class GapBaselines:
    def __init__(self, records, prior_total_mass=2.0):
        if prior_total_mass <= 0 or not math.isfinite(prior_total_mass):
            raise ValueError('Prior mass must be positive and finite')
        self.vocabulary = Tokenizer('codon').vocabulary
        self.index = {c: i for i, c in enumerate(self.vocabulary)}
        self.alpha = prior_total_mass / 64
        self.unigram = Counter()
        self.bigram = defaultdict(Counter)
        self.position = defaultdict(Counter)
        self.privileged = defaultdict(Counter)
        if not records:
            raise ValueError('Baselines require training records')
        for record in records:
            seq = codons(record['rna'])
            for pos in range(1, len(seq)):
                c = seq[pos]
                self.unigram[c] += 1
                self.bigram[seq[pos - 1]][c] += 1
                self.position[pos][c] += 1
                self.privileged[(record['family'], record['tax_id'], pos)][c] += 1

    def distribution(self, method, previous, position, record):
        fallback = None
        if method not in METHODS:
            raise ValueError('Unknown baseline')
        if method == 'uniform':
            return np.full(64, 1 / 64), fallback
        if method == 'unigram':
            counts = self.unigram
        elif method == 'bigram':
            counts = self.bigram.get(previous)
            if not counts:
                counts, fallback = self.unigram, 'pooled_unigram'
        elif method == 'position':
            counts = self.position.get(position)
            if not counts:
                counts, fallback = self.unigram, 'pooled_unigram'
        else:
            counts = self.privileged.get((record['family'], record['tax_id'], position))
            if not counts:
                counts, fallback = self.position.get(position), 'pooled_position'
            if not counts:
                counts, fallback = self.unigram, 'pooled_unigram'
        values = np.asarray([counts[c] for c in self.vocabulary], dtype=np.float64) + self.alpha
        return values / values.sum(), fallback

    def score_and_generate(self, method, prefix, truth, start_codon, record):
        actual = codons(truth)
        previous_true = previous_pred = codons(prefix)[-1]
        prediction, logp, scoring_fallbacks, generation_fallbacks = [], 0., Counter(), Counter()
        for i, c in enumerate(actual):
            values, fallback = self.distribution(method, previous_true, start_codon + i, record)
            logp += math.log(values[self.index[c]])
            if fallback:
                scoring_fallbacks[fallback] += 1
            values, fallback = self.distribution(method, previous_pred, start_codon + i, record)
            chosen = self.vocabulary[int(values.argmax())]
            if fallback:
                generation_fallbacks[fallback] += 1
            prediction.append(chosen)
            previous_true, previous_pred = c, chosen
        return {'prediction': ''.join(prediction), 'gap_log_probability': logp,
                'gap_bits_per_base': -logp / (len(truth) * math.log(2)),
                'scoring_fallbacks': dict(scoring_fallbacks),
                'generation_fallbacks': dict(generation_fallbacks),
                'privileged_metadata': method == 'privileged_family_taxon_position'}
