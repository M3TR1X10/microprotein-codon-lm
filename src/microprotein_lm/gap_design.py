"""Outcome-blind, frozen-cohort gap cases and strict holdout identity checks.

This module does not load models, score sequences, acquire biological records or
write files. It never substitutes training records for a missing held-out panel.
"""
from collections import Counter, defaultdict
import copy
import hashlib
import json
from pathlib import Path
import random
import re

from Bio.Data import CodonTable

from .io import read_json, sha256


SOURCE_KEYS = ('accession', 'retrieved_accession', 'parent_accession',
               'retrieved_parent_accession', 'index_parent_accession')
IDENTITY_KEYS = ('sequence_sha256', 'family', 'tax_id', 'translation_table')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True,
                                    separators=(',', ':')).encode()).hexdigest()


def _rank(label, seed, *values):
    return digest([label, seed, *values])


def _positive_integer(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f'{name} must be a positive integer')
    return value


def stable_source_accessions(record):
    """Normalize accession versions, not studies/projects into individuals."""
    found = set()
    provenance = record.get('provenance', [])
    if not isinstance(provenance, list):
        raise ValueError('Provenance must be a list')
    for source in provenance:
        if not isinstance(source, dict):
            raise ValueError('Every provenance entry must be a mapping')
        for key in SOURCE_KEYS:
            value = source.get(key) or ''
            if not isinstance(value, str):
                raise ValueError('Source accession fields must be strings')
            for accession in value.split(','):
                accession = accession.strip().upper()
                if accession:
                    found.add(re.sub(r'\.\d+$', '', accession))
    return found


def validate_record_identity(record):
    rna = record.get('rna')
    if not isinstance(rna, str) or len(rna) < 6 or len(rna) % 3 or set(rna) - set('AUCG'):
        raise ValueError('Expected complete unambiguous RNA codons')
    if hashlib.sha256(rna.encode()).hexdigest() != record.get('sequence_sha256'):
        raise ValueError('RNA checksum does not match record identity')
    if not isinstance(record.get('family'), str) or not record['family'].strip():
        if record.get('gene') in ('MT-ATP8', 'ATP8'):
            record['family'] = 'ATP8'
        else:
            raise ValueError('Family label must be nonempty')
    _positive_integer(record.get('tax_id'), 'tax_id')
    code = record.get('translation_table')
    if type(code) is not int or code not in (1, 2):
        raise ValueError('Only the declared genetic codes 1 and 2 are supported')
    expected = 2 if record['family'] == 'ATP8' else 1
    if code != expected:
        raise ValueError('Family and declared genetic code disagree')
    table = CodonTable.unambiguous_rna_by_id[code]
    triplets = [rna[i:i+3] for i in range(0, len(rna), 3)]
    if triplets[0] not in table.start_codons or triplets[-1] not in table.stop_codons:
        raise ValueError('Expected complete annotated CDS start and terminal stop')
    if any(c in table.stop_codons for c in triplets[1:-1]):
        raise ValueError('Internal stop in an alleged eligible observed CDS')
    if len(triplets) - 1 > 100:
        raise ValueError('CDS exceeds the existing microprotein gate')
    if 'taxa' in record:
        taxa = record['taxa']
        if not isinstance(taxa, list) or any(type(t) is not int or t < 1 for t in taxa):
            raise ValueError('Invalid recorded taxonomic provenance')
        if record['tax_id'] not in taxa:
            raise ValueError('Primary taxon is absent from taxonomic provenance')
    return {k: record[k] for k in IDENTITY_KEYS}


def verify_holdout(records, training_records, *, allowed_taxa, allowed_families):
    """Fail on leakage or invalid identity; selection never silently fixes it.

    Pass provenance-rich training records, not the compact training-only cohort.
    The file-backed loader below additionally verifies that this pool includes
    the original pilot and all five frozen secondary training views.
    """
    allowed_taxa = {int(t) for t in allowed_taxa}
    allowed_families = set(allowed_families)
    blocked, sources = set(), set()
    for record in training_records:
        validate_record_identity(record)
        blocked.add(record['sequence_sha256'])
        sources.update(stable_source_accessions(record))
    if not blocked or not sources:
        raise ValueError('Training RNA and provenance are required for a strict holdout audit')
    identities, seen = [], set()
    for record in records:
        identity = validate_record_identity(record)
        sequence_id = record['sequence_sha256']
        if sequence_id in seen:
            raise ValueError('Test cohort contains duplicate RNA')
        seen.add(sequence_id)
        if record['tax_id'] not in allowed_taxa or set(record.get('taxa', [record['tax_id']])) - allowed_taxa:
            raise ValueError('Test taxonomy exceeds the original species panel')
        if record['family'] not in allowed_families:
            raise ValueError('Test family exceeds the original family panel')
        if sequence_id in blocked:
            raise ValueError('Test RNA was present in the frozen training union')
        accessions = stable_source_accessions(record)
        if not accessions:
            raise ValueError('Test source provenance is required')
        if accessions & sources:
            raise ValueError('Test source/parent accession overlaps frozen training provenance')
        identities.append(identity)
    return {'n_sequences': len(records), 'training_union_sequences': len(blocked),
            'training_union_sha256': digest(sorted(blocked)),
            'training_source_accessions_sha256': digest(sorted(sources)),
            'ordered_sequence_metadata_sha256': digest(identities),
            'individual_independence': 'Not established by accession disjointness'}


def _within(root, name):
    path = (root / name).resolve()
    if path != root and root not in path.parents:
        raise ValueError('Frozen file path leaves the repository')
    return path


def _jsonl(path):
    with path.open(encoding='utf-8') as stream:
        return [json.loads(line) for line in stream if line.strip()]


def load_verified_holdout(root, manifest_path='data/processed/gap/manifest.json',
                          acquisition_config='configs/gap-biology.json',
                          evaluation_config='configs/gap-evaluation.json'):
    """Read the acquisition schema and verify its frozen provenance/identities.

    This checks hashes only; no research checkpoint is loaded. Raw retrieval
    files are checked once per unique path, never inferred from accession names.
    """
    root = Path(root).resolve()
    manifest_file = _within(root, manifest_path)
    manifest = read_json(manifest_file)
    if manifest.get('stage') != 'untouched_test':
        raise ValueError('Gap acquisition is not frozen as an untouched test cohort')
    identity = manifest['identity']
    config_file = _within(root, acquisition_config)
    config = read_json(config_file)
    biology_file = _within(root, config['training_biology_file'])
    training_file = _within(root, config['training_pool_file'])
    references_file = _within(root, config['reference_file'])
    cohort_file = manifest_file.parent / 'cohort.jsonl'
    report_file = root / 'reports/gap/acquisition.json'
    required = [(config_file, identity['config_sha256']),
                (biology_file, identity['training_biology_sha256']),
                (references_file, identity['references_sha256']),
                (training_file, identity['training_pool_sha256']),
                (cohort_file, identity['cohort_sha256']),
                (report_file, manifest['acquisition_report_sha256'])]
    required.extend((_within(root, p), h) for p, h in identity['source_hashes'].items())
    for path, expected in required:
        if sha256(path) != expected:
            raise ValueError(f'Frozen gap source/manifest checksum mismatch: {path.name}')
    if sha256(training_file) != config['training_pool_sha256']:
        raise ValueError('Training pool configuration identity disagrees')
    report = read_json(report_file)
    if report['identity'] != identity:
        raise ValueError('Acquisition report and frozen test identity disagree')
    records, training = _jsonl(cohort_file), _jsonl(training_file)
    if len(records) != manifest['n_sequences'] or len(training) != config['training_pool_sequences']:
        raise ValueError('Frozen sequence counts disagree')
    evaluation = read_json(_within(root, evaluation_config))
    blocked = {r['sequence_sha256'] for r in training}
    for relative in evaluation['exclude_training_cohorts']:
        for record in _jsonl(_within(root, relative)):
            validate_record_identity(record)
            if record['sequence_sha256'] not in blocked:
                raise ValueError('Provenance-rich exclusion pool omits a frozen training sequence')
    biology = read_json(biology_file)
    verification = verify_holdout(records, training, allowed_taxa=biology['taxa'],
                                  allowed_families=biology['families'])
    raw_hashes = {}
    for record in records:
        for source in record['provenance']:
            for name_key, hash_key in (('source_file', 'source_file_sha256'),
                                       ('index_file', 'index_file_sha256')):
                name, expected = source.get(name_key), source.get(hash_key)
                if not name or not expected:
                    raise ValueError('Test provenance lacks raw source/index identity')
                if name in raw_hashes and raw_hashes[name] != expected:
                    raise ValueError('Conflicting raw source hashes in test provenance')
                raw_hashes[name] = expected
    for name, expected in raw_hashes.items():
        if sha256(_within(root, name)) != expected:
            raise ValueError('Test raw source checksum mismatch')
    return {'records': records, 'training_records': training, 'verification': verification,
            'identity': {'manifest_sha256': sha256(manifest_file),
                         'cohort_sha256': sha256(cohort_file),
                         'evaluation_config_sha256': sha256(_within(root, evaluation_config)),
                         **verification}}


def select_panel(records, panel):
    maximum = _positive_integer(panel['max_sequences'], 'max_sequences')
    quota = _positive_integer(panel['max_per_family_taxon'], 'max_per_family_taxon')
    seed = panel['selection_seed']
    groups, seen = defaultdict(list), set()
    for record in records:
        validate_record_identity(record)
        if record['sequence_sha256'] in seen:
            raise ValueError('Distinct test CDS are required before panel selection')
        seen.add(record['sequence_sha256'])
        groups[(record['family'], record['tax_id'])].append(record)
    strata = sorted(groups, key=lambda key: _rank('stratum', seed, *key))
    for key in strata:
        groups[key].sort(key=lambda r: _rank('sequence', seed, r['sequence_sha256']))
    selected = []
    for level in range(quota):
        for key in strata:
            if level < len(groups[key]):
                selected.append(copy.deepcopy(groups[key][level]))
                if len(selected) == maximum:
                    return selected
    return selected


def _anchor(record, minimum_left, maximum_gap, minimum_right, seed, label):
    sense_codons = len(record['rna']) // 3 - 1
    last = sense_codons - maximum_gap - minimum_right
    starts = range(max(1, minimum_left), last + 1)
    return min(starts, key=lambda s: _rank(label, seed, record['sequence_sha256'], s),
               default=None)


def _case(record, stage, start, gap, left, right):
    if type(start) is not int or start < 1:
        raise ValueError('A gap cannot include the initiation codon')
    _positive_integer(gap, 'gap_codons')
    if left == 'full':
        supplied_left = start
    else:
        supplied_left = _positive_integer(left, 'left_context_codons')
        if supplied_left > start:
            raise ValueError('Requested left context is unavailable')
    if type(right) is not int or right < 0:
        raise ValueError('Right context must be a nonnegative integer')
    sense_codons = len(record['rna']) // 3 - 1
    if start + gap + right > sense_codons:
        raise ValueError('Gap/right context crosses the terminal stop boundary')
    anchor_id = digest(['gap_anchor_v1', *[record[k] for k in IDENTITY_KEYS], start])
    row = {k: record[k] for k in IDENTITY_KEYS}
    row.update({'anchor_id': anchor_id, 'start_codon': start, 'gap_codons': gap,
                'left_context_codons': left, 'right_context_codons': right, 'stage': stage,
                'left_start_codon': start - supplied_left,
                'left_start_nt': 3 * (start - supplied_left), 'gap_start_nt': 3 * start,
                'gap_end_nt': 3 * (start + gap), 'right_end_nt': 3 * (start + gap + right),
                'observed_left_codons': supplied_left, 'target_bases': 3 * gap,
                'available_left_codons': start, 'available_right_codons': sense_codons - start - gap,
                'support_key': f'{record["family"]}:{record["tax_id"]}:{record["translation_table"]}'})
    # Right length is deliberately absent from this key: the same beam is reused.
    row['candidate_bank_key'] = digest(['candidate_bank_v1', stage, anchor_id, gap, left])
    row['case_id'] = digest(row)
    return row


def case_regions(record, case):
    """Checked coordinates for the runner; keep truth out of proposal APIs."""
    validate_record_identity(record)
    if any(case[k] != record[k] for k in IDENTITY_KEYS):
        raise ValueError('Case and sequence metadata disagree')
    expected = _case(record, case['stage'], case['start_codon'], case['gap_codons'],
                     case['left_context_codons'], case['right_context_codons'])
    if case != expected:
        raise ValueError('Case coordinates or identity changed after design')
    rna = record['rna']
    return {'prefix': rna[case['left_start_nt']:case['gap_start_nt']],
            'truth': rna[case['gap_start_nt']:case['gap_end_nt']],
            'right': rna[case['gap_end_nt']:case['right_end_nt']],
            'offset_nt': case['left_start_nt']}


def sequence_clusters(records, threshold=0.99):
    if type(threshold) not in (int, float) or not 0 < threshold <= 1:
        raise ValueError('Cluster identity threshold must lie in (0,1]')
    groups = defaultdict(list)
    seen = set()
    for record in records:
        validate_record_identity(record)
        if record['sequence_sha256'] in seen:
            raise ValueError('Clustering requires distinct RNA')
        seen.add(record['sequence_sha256'])
        groups[(record['family'], record['tax_id'], len(record['rna']))].append(record)
    membership, clusters, strata = {}, [], {}
    for key in sorted(groups):
        rows = sorted(groups[key], key=lambda r: r['sequence_sha256'])
        parents = list(range(len(rows)))

        def find(i):
            while parents[i] != i:
                parents[i] = parents[parents[i]]
                i = parents[i]
            return i

        for i, row in enumerate(rows):
            for j in range(i):
                differences = sum(a != b for a, b in zip(row['rna'], rows[j]['rna']))
                if differences <= (1 - threshold) * key[2] + 1e-10:
                    parents[find(i)] = find(j)
        components = defaultdict(list)
        for i, row in enumerate(rows):
            components[find(i)].append(row['sequence_sha256'])
        for members in sorted(components.values()):
            members.sort()
            cluster_id = digest(['gap_cluster_v1', threshold, key, members])
            clusters.append({'cluster_id': cluster_id, 'family': key[0], 'tax_id': key[1],
                             'length_nt': key[2], 'sequence_sha256s': members})
            membership.update({member: cluster_id for member in members})
        strata[f'{key[0]}:{key[1]}:{key[2]}'] = {'n_sequences': len(rows),
                                               'n_clusters': len(components)}
    return {'identity_threshold': threshold, 'membership': membership,
            'clusters': clusters, 'strata': strata,
            'interpretation': 'Similarity components, not verified independent individuals'}


def _support(cases):
    conditions = defaultdict(list)
    for case in cases:
        conditions[(case['gap_codons'], str(case['left_context_codons']),
                    case['right_context_codons'])].append(case)
    cells = []
    for key in sorted(conditions):
        rows = conditions[key]
        cells.append({'gap_codons': rows[0]['gap_codons'],
                      'left_context_codons': rows[0]['left_context_codons'],
                      'right_context_codons': rows[0]['right_context_codons'],
                      'n_cases': len(rows), 'n_sequences': len({r['sequence_sha256'] for r in rows}),
                      'target_bases': sum(r['target_bases'] for r in rows),
                      'family_taxon_sequence_counts': dict(sorted(Counter(
                          f'{r["family"]}:{r["tax_id"]}' for r in rows).items()))})
    return {'n_cases': len(cases), 'n_sequences': len({r['sequence_sha256'] for r in cases}),
            'n_anchors': len({r['anchor_id'] for r in cases}), 'cells': cells}


def make_cases(records, config):
    """Deterministic, support-matched plan; RNA is retained only in selected_records.

    Run verify_holdout/load_verified_holdout before this function on research
    data. This pure case constructor independently rechecks RNA/metadata identity.
    """
    selected = select_panel(records, config['panel'])
    primary, focused, rerank, exploratory = [], [], [], []
    omitted, focused_omitted = [], []
    lengths = config['primary']['gap_codons']
    if not lengths or sorted(set(lengths)) != lengths:
        raise ValueError('Primary gap lengths must be nonempty, unique and increasing')
    for length in lengths:
        _positive_integer(length, 'gap_codons')
    if config['focused_context']['gap_codons'] != lengths or config['right_flank']['gap_codons'] != lengths:
        raise ValueError('Focused comparisons must retain primary gap lengths')
    focus_config, right_config = config['focused_context'], config['right_flank']
    focus_max = _positive_integer(focus_config['max_sequences'], 'focused max_sequences')
    if focus_config['minimum_available_left_codons'] < max(
            x for x in focus_config['left_context_codons'] if x != 'full'):
        raise ValueError('Focused anchor cannot support every declared left context')
    if focus_config['minimum_available_right_after_gap16_codons'] < max(right_config['right_context_codons']):
        raise ValueError('Focused anchor cannot support every declared right context')
    for record in selected:
        start = _anchor(record, 1, max(lengths), 0, config['primary']['gap_seed'], 'primary_anchor')
        if start is None:
            omitted.append({'sequence_sha256': record['sequence_sha256'],
                            'reason': 'insufficient_sense_codons_for_primary_nested_gaps'})
            continue
        primary.extend(_case(record, 'primary', start, length, 'full', 0) for length in lengths)
        for length in config['exploratory']['gap_codons']:
            _positive_integer(length, 'exploratory gap_codons')
            if start + length <= len(record['rna']) // 3 - 1:
                exploratory.append(_case(record, 'exploratory', start, length, 'full', 0))
    eligible_focus = []
    for record in selected:
        if (record['family'], record['tax_id']) != (focus_config['family'], focus_config['tax_id']):
            continue
        start = _anchor(record, focus_config['minimum_available_left_codons'], max(lengths),
                        focus_config['minimum_available_right_after_gap16_codons'],
                        focus_config['gap_seed'], 'focused_anchor')
        if start is None:
            focused_omitted.append({'sequence_sha256': record['sequence_sha256'],
                                   'reason': 'insufficient_common_left_gap_right_support'})
        else:
            eligible_focus.append((record, start))
    eligible_focus.sort(key=lambda pair: _rank('focused_sequence', config['panel']['selection_seed'],
                                               pair[0]['sequence_sha256']))
    for record, start in eligible_focus[:focus_max]:
        for length in lengths:
            for left in focus_config['left_context_codons']:
                focused.append(_case(record, 'focused_context', start, length, left, 0))
            for left in right_config['left_context_codons']:
                for right in right_config['right_context_codons']:
                    rerank.append(_case(record, 'right_flank', start, length, left, right))
    stages = {'primary': primary, 'focused_context': focused,
              'right_flank': rerank, 'exploratory': exploratory}
    all_cases = [case for rows in stages.values() for case in rows]
    if len({case['case_id'] for case in all_cases}) != len(all_cases):
        raise ValueError('Duplicate cases from an invalid experimental grid')
    return {'selected_records': selected, **stages,
            'support': {'n_eligible_records': len(records), 'n_selected_records': len(selected),
                        'selected_family_taxon_counts': dict(sorted(Counter(
                            f'{r["family"]}:{r["tax_id"]}' for r in selected).items())),
                        'omitted_primary': omitted, 'omitted_focused_context': focused_omitted,
                        'stages': {name: _support(rows) for name, rows in stages.items()},
                        'ordered_sequence_metadata_sha256': digest([
                            {k: r[k] for k in IDENTITY_KEYS} for r in selected])},
            'clusters': sequence_clusters(selected, config['uncertainty']['cluster_identity_threshold']),
            'multi_threshold_clusters': multi_threshold_sequence_clusters(selected, (0.99, 0.95, 0.90)),
            'cluster_resampling': resample_sequence_clusters(
                selected, (0.99, 0.95, 0.90),
                n_draws=config.get('uncertainty', {}).get('bootstrap_iterations', 2000),
                seed=config.get('uncertainty', {}).get('bootstrap_seed', 20260929),
                min_clusters_per_stratum=config.get('uncertainty', {}).get('minimum_clusters_per_stratum_for_band', 5))}


def sequence_identity(rna1, rna2):
    """Calculate normalized sequence identity in [0.0, 1.0]."""
    if rna1 == rna2:
        return 1.0
    len1, len2 = len(rna1), len(rna2)
    if len1 == len2:
        mismatches = sum(a != b for a, b in zip(rna1, rna2))
        return 1.0 - (mismatches / len1)
    try:
        from Bio.Align import PairwiseAligner
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 1.0
        aligner.mismatch_score = 0.0
        aligner.open_gap_score = -1.0
        aligner.extend_gap_score = -0.5
        alignments = aligner.align(rna1, rna2)
        if alignments:
            aln = alignments[0]
            c = aln.counts()
            return c.identities / max(len1, len2)
    except Exception:
        pass
    return 1.0 - (_levenshtein_distance(rna1, rna2) / max(len1, len2))


def _levenshtein_distance(s1, s2):
    if len(s1) < len(s2):
        return _levenshtein_distance(s2, s1)
    if len(s2) == 0:
        return len(s1)
    previous_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]


def multi_threshold_sequence_clusters(records, thresholds=(0.99, 0.95, 0.90)):
    """Compute sequence clusters across multiple identity thresholds (e.g. 99%, 95%, 90%)."""
    results = {}
    for threshold in thresholds:
        if type(threshold) not in (int, float) or not 0 < threshold <= 1:
            raise ValueError('Cluster identity threshold must lie in (0,1]')
        key_str = f"{threshold:.2f}"
        results[key_str] = sequence_clusters(records, threshold=threshold)
    return {
        'thresholds': results,
        'summary': {
            f"{t:.2f}": {
                'n_clusters': len(results[f"{t:.2f}"]['clusters']),
                'n_sequences': len(results[f"{t:.2f}"]['membership']),
                'n_strata': len(results[f"{t:.2f}"]['strata']),
                'max_cluster_size': max((len(c['sequence_sha256s']) for c in results[f"{t:.2f}"]['clusters']), default=0)
            }
            for t in thresholds
        }
    }


def resample_sequence_clusters(records, thresholds=(0.99, 0.95, 0.90), n_draws=2000,
                                seed=20260929, min_clusters_per_stratum=5):
    """Stratified cluster bootstrap resampling across specified identity thresholds (99%, 95%, 90%)."""
    _positive_integer(n_draws, 'n_draws')
    _positive_integer(min_clusters_per_stratum, 'min_clusters_per_stratum')

    resample_report = {}
    for threshold in thresholds:
        if type(threshold) not in (int, float) or not 0 < threshold <= 1:
            raise ValueError('Cluster identity threshold must lie in (0,1]')
        key_str = f"{threshold:.2f}"
        clust = sequence_clusters(records, threshold=threshold)
        clusters = clust['clusters']

        strata_clusters = defaultdict(list)
        for cluster in clusters:
            stratum_key = f"{cluster['family']}:{cluster['tax_id']}:{cluster['length_nt']}"
            strata_clusters[stratum_key].append(cluster)

        strata_status = {}
        all_strata_eligible = True
        for s_key, c_list in sorted(strata_clusters.items()):
            n_c = len(c_list)
            eligible = n_c >= min_clusters_per_stratum
            if not eligible:
                all_strata_eligible = False
            strata_status[s_key] = {
                'n_clusters': n_c,
                'n_sequences': sum(len(c['sequence_sha256s']) for c in c_list),
                'bands_eligible': eligible
            }

        draw_sequence_counts = []
        draw_cluster_counts = []

        for d in range(n_draws):
            rng = random.Random(seed + d + int(threshold * 10000))
            draw_seq_count = 0
            draw_cluster_count = 0

            for s_key in sorted(strata_clusters):
                c_list = strata_clusters[s_key]
                if not c_list:
                    continue
                sampled_clusters = rng.choices(c_list, k=len(c_list))
                draw_cluster_count += len(sampled_clusters)
                draw_seq_count += sum(len(c['sequence_sha256s']) for c in sampled_clusters)

            draw_sequence_counts.append(draw_seq_count)
            draw_cluster_counts.append(draw_cluster_count)

        draw_sequence_counts.sort()
        draw_cluster_counts.sort()

        lower_idx = int(0.025 * n_draws)
        median_idx = int(0.50 * n_draws)
        upper_idx = int(0.975 * n_draws)

        mean_seqs = sum(draw_sequence_counts) / max(1, n_draws)
        var_seqs = sum((x - mean_seqs) ** 2 for x in draw_sequence_counts) / max(1, n_draws - 1)
        std_seqs = var_seqs ** 0.5

        resample_report[key_str] = {
            'identity_threshold': threshold,
            'bands_available': all_strata_eligible and len(records) > 0,
            'min_clusters_per_stratum': min_clusters_per_stratum,
            'n_draws': n_draws,
            'strata_status': strata_status,
            'draw_summary': {
                'sequence_count_mean': mean_seqs,
                'sequence_count_std': std_seqs,
                'sequence_count_ci_95': [draw_sequence_counts[lower_idx], draw_sequence_counts[upper_idx]] if draw_sequence_counts else [0, 0],
                'sequence_count_median': draw_sequence_counts[median_idx] if draw_sequence_counts else 0,
                'cluster_count_mean': sum(draw_cluster_counts) / max(1, n_draws),
                'cluster_count_median': draw_cluster_counts[median_idx] if draw_cluster_counts else 0
            }
        }

    return {
        'resample_by_threshold': resample_report,
        'bootstrap_seed': seed,
        'n_draws': n_draws,
        'n_input_records': len(records)
    }


class HomologyGraphAuditor:
    """Automated MMseqs2 / BLAST Homology Graph Auditor.

    Audits sequence similarity graph topology, cross-cohort homology edges, and potential
    leakage between test records and training records across multiple identity thresholds.
    """
    def __init__(self, thresholds=(0.99, 0.95, 0.90), min_identity_report=0.80):
        self.thresholds = sorted(thresholds, reverse=True)
        self.min_identity_report = min_identity_report

    def audit(self, records, training_records=None):
        test_records = records or []
        train_records = training_records or []

        for r in test_records:
            validate_record_identity(r)
        for r in train_records:
            validate_record_identity(r)

        nodes = {}
        for idx, r in enumerate(test_records):
            node_id = f"test_{r['sequence_sha256']}_{idx}"
            nodes[node_id] = {
                'node_id': node_id,
                'sequence_sha256': r['sequence_sha256'],
                'cohort': 'test',
                'family': r['family'],
                'tax_id': r['tax_id'],
                'length_nt': len(r['rna']),
                'rna': r['rna']
            }
        for idx, r in enumerate(train_records):
            node_id = f"train_{r['sequence_sha256']}_{idx}"
            nodes[node_id] = {
                'node_id': node_id,
                'sequence_sha256': r['sequence_sha256'],
                'cohort': 'training',
                'family': r['family'],
                'tax_id': r['tax_id'],
                'length_nt': len(r['rna']),
                'rna': r['rna']
            }

        node_ids = sorted(nodes.keys())
        edges = []

        for i, id1 in enumerate(node_ids):
            n1 = nodes[id1]
            for j in range(i + 1, len(node_ids)):
                id2 = node_ids[j]
                n2 = nodes[id2]

                ident = sequence_identity(n1['rna'], n2['rna'])
                if ident >= self.min_identity_report:
                    is_cross_cohort = (n1['cohort'] != n2['cohort'])
                    is_cross_stratum = (n1['family'] != n2['family']) or (n1['tax_id'] != n2['tax_id'])
                    edges.append({
                        'node1': id1,
                        'node2': id2,
                        'identity': ident,
                        'is_cross_cohort': is_cross_cohort,
                        'is_cross_stratum': is_cross_stratum,
                        'cohorts': tuple(sorted([n1['cohort'], n2['cohort']])),
                        'strata': tuple(sorted([f"{n1['family']}:{n1['tax_id']}", f"{n2['family']}:{n2['tax_id']}"]))
                    })

        threshold_audits = {}
        leakage_warnings = []

        for t in self.thresholds:
            t_key = f"{t:.2f}"
            t_edges = [e for e in edges if e['identity'] >= t - 1e-10]

            parents = {nid: nid for nid in node_ids}
            def find(nid):
                while parents[nid] != nid:
                    parents[nid] = parents[parents[nid]]
                    nid = parents[nid]
                return nid

            for e in t_edges:
                p1, p2 = find(e['node1']), find(e['node2'])
                if p1 != p2:
                    parents[p1] = p2

            components = defaultdict(list)
            for nid in node_ids:
                components[find(nid)].append(nid)

            comp_list = sorted(components.values(), key=lambda c: (-len(c), sorted(c)))

            leakage_edges = [e for e in t_edges if e['is_cross_cohort']]
            cross_stratum_edges = [e for e in t_edges if e['is_cross_stratum']]

            degrees = defaultdict(int)
            for e in t_edges:
                degrees[e['node1']] += 1
                degrees[e['node2']] += 1

            max_degree = max(degrees.values(), default=0)
            isolated_count = sum(1 for nid in node_ids if degrees[nid] == 0)

            threshold_audits[t_key] = {
                'identity_threshold': t,
                'n_nodes': len(node_ids),
                'n_edges': len(t_edges),
                'n_components': len(comp_list),
                'n_leakage_edges': len(leakage_edges),
                'n_cross_stratum_edges': len(cross_stratum_edges),
                'max_component_size': max((len(c) for c in comp_list), default=0),
                'mean_component_size': (len(node_ids) / len(comp_list)) if comp_list else 0.0,
                'max_node_degree': max_degree,
                'isolated_node_count': isolated_count,
                'max_cross_cohort_identity': max((e['identity'] for e in leakage_edges), default=0.0)
            }

            if leakage_edges:
                for le in leakage_edges:
                    leakage_warnings.append({
                        'threshold': t,
                        'test_sha256': le['node1'] if nodes[le['node1']]['cohort'] == 'test' else le['node2'],
                        'training_sha256': le['node2'] if nodes[le['node1']]['cohort'] == 'test' else le['node1'],
                        'identity': le['identity']
                    })

        report = {
            'n_test_sequences': len(test_records),
            'n_training_sequences': len(train_records),
            'n_total_nodes': len(nodes),
            'thresholds': threshold_audits,
            'leakage_warnings': leakage_warnings,
            'auditor': 'MMseqs2/BLAST Homology Graph Auditor v1'
        }
        report['graph_audit_sha256'] = digest(report)
        return report


def audit_homology_graph(records, training_records=None, thresholds=(0.99, 0.95, 0.90),
                         min_identity_report=0.80):
    """Convenience wrapper for HomologyGraphAuditor."""
    auditor = HomologyGraphAuditor(thresholds=thresholds, min_identity_report=min_identity_report)
    return auditor.audit(records, training_records=training_records)

