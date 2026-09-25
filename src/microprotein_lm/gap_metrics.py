"""Gap reconstruction outcomes with explicit denominators and genetic codes."""
import math
from collections import Counter
from Bio.Align import PairwiseAligner
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
        self.msa_support = MSAConservationSupport(records)
        self.domain_support = DomainBoundarySupport()

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
        base_out = {'family_seen': record['family'] in self.families,
                    'taxon_seen': record['tax_id'] in self.taxa,
                    'family_taxon_seen': (record['family'], record['tax_id']) in self.family_taxa,
                    'codon_support_labels': labels, 'training_allele_frequencies': frequencies,
                    'training_support': buckets}
        msa_out = self.msa_support.describe_msa(record, start, truth, prediction)
        domain_out = self.domain_support.describe_domain_boundaries(record, start, len(true))
        return {**base_out, **msa_out, **domain_out}


def calculate_shannon_entropy(counts):
    """Calculate Shannon entropy in bits (base 2) for a frequency distribution."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = 0.0
    for cnt in counts.values():
        if cnt > 0:
            p = cnt / total
            entropy -= p * math.log2(p)
    return entropy


def calculate_conservation_score(counts, max_classes=64):
    """Calculate normalized conservation score [0.0, 1.0] from Shannon entropy."""
    total = sum(counts.values())
    if total <= 1 or max_classes <= 1:
        return 1.0
    entropy = calculate_shannon_entropy(counts)
    max_entropy = math.log2(max_classes)
    return max(0.0, min(1.0, 1.0 - (entropy / max_entropy)))


class MSAConservationSupport:
    """MSA-aligned conservation metrics across family sequences."""

    def __init__(self, records):
        self.family_records = {}
        self.family_msa = {}
        for r in records:
            fam = r.get('family')
            if fam:
                self.family_records.setdefault(fam, []).append(r)
        self._build_msa_alignments()

    def _build_msa_alignments(self):
        aligner = PairwiseAligner(mode='global', match_score=1.0, mismatch_score=-1.0,
                                  open_gap_score=-2.0, extend_gap_score=-0.5)
        for fam, recs in self.family_records.items():
            if not recs:
                continue
            ref_rec = next((r for r in recs if r.get('tax_id') == 9606), recs[0])
            ref_codons = codons(ref_rec['rna'])
            ref_aa = translate_gap(ref_rec['rna'], ref_rec.get('translation_table', 1))

            msa_columns = [Counter() for _ in ref_codons]
            msa_aa_columns = [Counter() for _ in ref_aa]

            for r in recs:
                r_codons = codons(r['rna'])
                r_aa = translate_gap(r['rna'], r.get('translation_table', 1))
                if len(r_codons) == len(ref_codons):
                    for idx, (c, a) in enumerate(zip(r_codons, r_aa)):
                        msa_columns[idx][c] += 1
                        msa_aa_columns[idx][a] += 1
                else:
                    try:
                        alignments = aligner.align(r_aa, ref_aa)
                        if alignments:
                            al = alignments[0]
                            for q_seg, r_seg in zip(al.aligned[0], al.aligned[1]):
                                q_s, q_e = q_seg
                                r_s, r_e = r_seg
                                length = min(q_e - q_s, r_e - r_s)
                                for k in range(length):
                                    if r_s + k < len(ref_codons) and q_s + k < len(r_codons):
                                        msa_columns[r_s + k][r_codons[q_s + k]] += 1
                                        msa_aa_columns[r_s + k][r_aa[q_s + k]] += 1
                    except Exception:
                        for idx, (c, a) in enumerate(zip(r_codons[:len(ref_codons)], r_aa[:len(ref_aa)])):
                            msa_columns[idx][c] += 1
                            msa_aa_columns[idx][a] += 1

            self.family_msa[fam] = {
                'ref_rec': ref_rec,
                'ref_codons': ref_codons,
                'ref_aa': ref_aa,
                'codon_columns': msa_columns,
                'aa_columns': msa_aa_columns,
            }

    def describe_msa(self, record, start, truth, prediction):
        true_c, pred_c = codons(truth), codons(prediction)
        fam = record.get('family')
        msa_data = self.family_msa.get(fam)

        codon_entropies = []
        aa_entropies = []
        codon_conservations = []
        aa_conservations = []
        consensus_codons = []
        consensus_aas = []
        is_consensus_matches = []
        allele_frequencies = []
        msa_labels = []

        for i, (act_c, pred_c_val) in enumerate(zip(true_c, pred_c)):
            pos = start + i
            if not msa_data or pos >= len(msa_data['codon_columns']):
                codon_entropies.append(None)
                aa_entropies.append(None)
                codon_conservations.append(None)
                aa_conservations.append(None)
                consensus_codons.append(None)
                consensus_aas.append(None)
                is_consensus_matches.append(False)
                allele_frequencies.append(0.0)
                msa_labels.append('unsupported_position')
                continue

            c_counts = msa_data['codon_columns'][pos]
            a_counts = msa_data['aa_columns'][pos]
            total_seqs = sum(c_counts.values())

            if total_seqs == 0:
                codon_entropies.append(None)
                aa_entropies.append(None)
                codon_conservations.append(None)
                aa_conservations.append(None)
                consensus_codons.append(None)
                consensus_aas.append(None)
                is_consensus_matches.append(False)
                allele_frequencies.append(0.0)
                msa_labels.append('unsupported_position')
                continue

            c_entropy = calculate_shannon_entropy(c_counts)
            a_entropy = calculate_shannon_entropy(a_counts)
            c_cons = calculate_conservation_score(c_counts, max_classes=64)
            a_cons = calculate_conservation_score(a_counts, max_classes=20)

            top_c = max(c_counts.keys(), key=lambda k: c_counts[k])
            top_a = max(a_counts.keys(), key=lambda k: a_counts[k])

            freq = c_counts.get(act_c, 0) / total_seqs

            if act_c not in c_counts:
                label = 'unobserved_codon'
            elif c_counts[act_c] < max(c_counts.values()):
                label = 'observed_minority'
            else:
                label = 'majority'

            codon_entropies.append(c_entropy)
            aa_entropies.append(a_entropy)
            codon_conservations.append(c_cons)
            aa_conservations.append(a_cons)
            consensus_codons.append(top_c)
            consensus_aas.append(top_a)
            is_consensus_matches.append(act_c == top_c)
            allele_frequencies.append(freq)
            msa_labels.append(label)

        valid_entropies = [e for e in codon_entropies if e is not None]
        valid_cons = [c for c in codon_conservations if c is not None]
        mean_entropy = sum(valid_entropies) / len(valid_entropies) if valid_entropies else 0.0
        mean_conservation = sum(valid_cons) / len(valid_cons) if valid_cons else 1.0

        return {
            'msa_family': fam,
            'msa_codon_entropies': codon_entropies,
            'msa_aa_entropies': aa_entropies,
            'msa_codon_conservation_scores': codon_conservations,
            'msa_aa_conservation_scores': aa_conservations,
            'msa_consensus_codons': consensus_codons,
            'msa_consensus_amino_acids': consensus_aas,
            'msa_is_consensus_match': is_consensus_matches,
            'msa_allele_frequencies': allele_frequencies,
            'msa_labels': msa_labels,
            'mean_msa_codon_entropy': mean_entropy,
            'mean_msa_codon_conservation': mean_conservation,
            'msa_consensus_match_rate': sum(is_consensus_matches) / len(is_consensus_matches) if is_consensus_matches else 0.0,
        }


DEFAULT_FAMILY_DOMAINS = {
    'ATP8': [{'domain_id': 'IPR001421', 'pfam_id': 'PF00137', 'name': 'ATP synthase subunit 8', 'database': 'InterPro', 'start': 0, 'end': 67}],
    'ATP5ME': [{'domain_id': 'IPR008386', 'pfam_id': 'PF05490', 'name': 'ATP synthase subunit e', 'database': 'InterPro', 'start': 0, 'end': 70}],
    'ATP5MF': [{'domain_id': 'IPR019344', 'pfam_id': 'PF10173', 'name': 'ATP synthase subunit f', 'database': 'InterPro', 'start': 0, 'end': 90}],
    'ATP5MJ': [{'domain_id': 'IPR012574', 'pfam_id': 'PF07519', 'name': 'ATP synthase subunit j', 'database': 'InterPro', 'start': 0, 'end': 80}],
    'ATP5MK': [{'domain_id': 'IPR009125', 'pfam_id': 'PF05943', 'name': 'ATP synthase subunit k', 'database': 'InterPro', 'start': 0, 'end': 58}],
    'ATP5F1E': [{'domain_id': 'IPR006721', 'pfam_id': 'PF04604', 'name': 'ATP synthase subunit epsilon', 'database': 'InterPro', 'start': 0, 'end': 50}],
    'ATP5MGL': [{'domain_id': 'IPR006808', 'pfam_id': 'PF04675', 'name': 'ATP synthase subunit g-like', 'database': 'InterPro', 'start': 0, 'end': 95}],
}


class DomainBoundarySupport:
    """InterPro/Pfam Domain Boundary coordinate feature extraction."""

    def __init__(self, family_domains=None):
        self.family_domains = family_domains or DEFAULT_FAMILY_DOMAINS

    def extract_domains(self, record):
        """Extract explicit or family-level default domain annotations."""
        if 'domains' in record and record['domains']:
            return record['domains']
        if 'interpro_domains' in record and record['interpro_domains']:
            return record['interpro_domains']
        if 'pfam_domains' in record and record['pfam_domains']:
            return record['pfam_domains']

        fam = record.get('family')
        if fam and fam in self.family_domains:
            cds_codons = len(record['rna']) // 3 if 'rna' in record else None
            defs = []
            for d in self.family_domains[fam]:
                d_copy = dict(d)
                if cds_codons is not None:
                    d_copy['end'] = min(d_copy['end'], cds_codons - 1)
                defs.append(d_copy)
            return defs
        return []

    def describe_domain_boundaries(self, record, start, gap_codons):
        domains = self.extract_domains(record)
        cds_total_codons = len(record['rna']) // 3 if 'rna' in record else None
        gap_end = start + gap_codons - 1
        gap_center = start + (gap_codons - 1) / 2.0

        if not domains:
            return {
                'has_domain_annotation': False,
                'annotated_domains': [],
                'domain_ids': [],
                'pfam_ids': [],
                'gap_start_codon': start,
                'gap_end_codon': gap_end,
                'gap_center_codon': gap_center,
                'gap_relative_cds_position': gap_center / cds_total_codons if cds_total_codons else None,
                'is_inside_domain': False,
                'is_overlapping_domain': False,
                'straddles_n_terminal_boundary': False,
                'straddles_c_terminal_boundary': False,
                'distance_to_domain_start': None,
                'distance_to_domain_end': None,
                'min_distance_to_boundary': None,
                'relative_domain_position': None,
                'domain_overlap_type': 'no_domain_annotation',
                'per_codon_in_domain': [False] * gap_codons,
                'per_codon_domain_relative_positions': [None] * gap_codons,
            }

        primary_domain = domains[0]
        d_start = primary_domain.get('start', 0)
        d_end = primary_domain.get('end', (cds_total_codons - 1) if cds_total_codons else 100)
        d_len = max(1, d_end - d_start + 1)

        is_inside = (start >= d_start) and (gap_end <= d_end)
        is_overlapping = max(start, d_start) <= min(gap_end, d_end)
        straddles_n = (start < d_start <= gap_end)
        straddles_c = (start <= d_end < gap_end)

        dist_to_start = start - d_start
        dist_to_end = d_end - gap_end
        min_dist_boundary = min(
            abs(start - d_start),
            abs(start - d_end),
            abs(gap_end - d_start),
            abs(gap_end - d_end)
        )

        rel_pos = (gap_center - d_start) / d_len

        if is_inside:
            overlap_type = 'fully_inside_domain'
        elif straddles_n and straddles_c:
            overlap_type = 'straddles_both_boundaries'
        elif straddles_n:
            overlap_type = 'straddles_n_boundary'
        elif straddles_c:
            overlap_type = 'straddles_c_boundary'
        elif gap_end < d_start:
            overlap_type = 'outside_n_terminal'
        else:
            overlap_type = 'outside_c_terminal'

        per_codon_in_domain = []
        per_codon_rel_pos = []
        for k in range(gap_codons):
            pos = start + k
            in_d = (d_start <= pos <= d_end)
            per_codon_in_domain.append(in_d)
            per_codon_rel_pos.append((pos - d_start) / d_len if in_d else None)

        domain_ids = [d.get('domain_id') for d in domains if 'domain_id' in d]
        pfam_ids = [d.get('pfam_id') for d in domains if 'pfam_id' in d]

        return {
            'has_domain_annotation': True,
            'annotated_domains': domains,
            'domain_ids': domain_ids,
            'pfam_ids': pfam_ids,
            'gap_start_codon': start,
            'gap_end_codon': gap_end,
            'gap_center_codon': gap_center,
            'gap_relative_cds_position': gap_center / cds_total_codons if cds_total_codons else None,
            'is_inside_domain': is_inside,
            'is_overlapping_domain': is_overlapping,
            'straddles_n_terminal_boundary': straddles_n,
            'straddles_c_terminal_boundary': straddles_c,
            'distance_to_domain_start': dist_to_start,
            'distance_to_domain_end': dist_to_end,
            'min_distance_to_boundary': min_dist_boundary,
            'relative_domain_position': rel_pos,
            'domain_overlap_type': overlap_type,
            'per_codon_in_domain': per_codon_in_domain,
            'per_codon_domain_relative_positions': per_codon_rel_pos,
        }

