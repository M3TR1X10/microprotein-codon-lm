"""Outcome-blind, frozen-cohort gap cases and strict holdout identity checks.

This module does not load models, score sequences, acquire biological records or
write files. It never substitutes training records for a missing held-out panel.
"""
from collections import Counter, defaultdict
import copy
import hashlib
import json
from pathlib import Path
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
            'clusters': sequence_clusters(selected, config['uncertainty']['cluster_identity_threshold'])}
