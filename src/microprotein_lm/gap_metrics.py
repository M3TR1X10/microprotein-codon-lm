"""Gap reconstruction outcomes with explicit denominators and genetic codes."""
from collections import Counter
from Bio.Data import CodonTable


def codons(rna):
    if not rna or len(rna) % 3 or set(rna) - set('AUCG'):
        raise ValueError('Expected complete unambiguous RNA codons')
    return [rna[i:i + 3] for i in range(0, len(rna), 3)]


def translate_gap(rna, translation_table):
    """Internal residues: alternative initiation rules do not apply to a gap."""
    table = CodonTable.unambiguous_rna_by_id[translation_table]
    return ''.join('*' if c in table.stop_codons else table.forward_table[c]
                   for c in codons(rna))


def reconstruction(truth, prediction, translation_table):
    if len(truth) != len(prediction):
        raise ValueError('A fixed-length gap must have the same target and output length')
    true_codons, pred_codons = codons(truth), codons(prediction)
    true_aa = translate_gap(truth, translation_table)
    pred_aa = translate_gap(prediction, translation_table)
    codon_matches = [int(a == b) for a, b in zip(true_codons, pred_codons)]
    aa_matches = [int(a == b) for a, b in zip(true_aa, pred_aa)]
    return {'exact_span': int(truth == prediction),
            'correct_bases': sum(a == b for a, b in zip(truth, prediction)),
            'target_bases': len(truth), 'correct_codons': sum(codon_matches),
            'target_codons': len(true_codons), 'codon_matches': codon_matches,
            'correct_amino_acids': sum(aa_matches), 'amino_acid_matches': aa_matches,
            'exact_peptide': int(true_aa == pred_aa),
            'predicted_stop_codons': pred_aa.count('*'),
            'has_premature_stop': int('*' in pred_aa),
            'first_error_codon': next((i + 1 for i, match in enumerate(codon_matches)
                                      if not match), None),
            'synonymous_errors': sum(c == 0 and a == 1 for c, a in zip(codon_matches, aa_matches))}


class TrainingPositionSupport:
    """Descriptive support derived solely from the frozen per-model training CDS."""

    def __init__(self, records):
        self.counts = {}
        self.family_taxa = set()
        self.families = set()
        self.taxa = set()
        for r in records:
            key = (r['family'], r['tax_id'])
            self.family_taxa.add(key)
            self.families.add(r['family'])
            self.taxa.add(r['tax_id'])
            for i, c in enumerate(codons(r['rna'])):
                self.counts.setdefault((*key, i), Counter())[c] += 1

    def describe(self, record, start, truth, prediction):
        true, pred = codons(truth), codons(prediction)
        buckets = {name: {'correct_codons': 0, 'target_codons': 0} for name in
                   ('majority', 'observed_minority', 'unobserved_codon', 'unsupported_position',
                    'observed_rare', 'training_variable', 'training_conserved')}
        labels, frequencies = [], []
        for i, (actual, guess) in enumerate(zip(true, pred)):
            counts = self.counts.get((record['family'], record['tax_id'], start + i))
            if not counts:
                label = 'unsupported_position'
            elif actual not in counts:
                label = 'unobserved_codon'
            elif counts[actual] < max(counts.values()):
                label = 'observed_minority'
            else:
                label = 'majority'
            labels.append(label)
            frequency = counts.get(actual, 0) / sum(counts.values()) if counts else None
            frequencies.append(frequency)
            buckets[label]['target_codons'] += 1
            buckets[label]['correct_codons'] += actual == guess
            if counts:
                conservation = 'training_variable' if len(counts) > 1 else 'training_conserved'
                buckets[conservation]['target_codons'] += 1
                buckets[conservation]['correct_codons'] += actual == guess
                if 0 < frequency <= .01:
                    buckets['observed_rare']['target_codons'] += 1
                    buckets['observed_rare']['correct_codons'] += actual == guess
        return {'family_seen': record['family'] in self.families,
                'taxon_seen': record['tax_id'] in self.taxa,
                'family_taxon_seen': (record['family'], record['tax_id']) in self.family_taxa,
                'codon_support_labels': labels, 'training_allele_frequencies': frequencies,
                'training_support': buckets}
