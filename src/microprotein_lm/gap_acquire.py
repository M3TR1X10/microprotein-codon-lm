"""New observed CDS for an untouched gap benchmark; never changes training data.

Unreviewed UniProt records are retrieval leads, not new biological references.
The frozen reviewed own-species references and complete-CDS quality gates apply.
"""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from Bio.Align import PairwiseAligner

from .acquire import Download, ENA_SEARCH, UNIPROT_URL
from .io import now, read_json, sha256
from .secondary_acquire import (fetch_records, merge_records,
    validate_secondary_record, write_canonical_json)
from .secondary_genomes import GENOME_FIELDS, extract_atp8, versioned_accession


def unreviewed_query(biology):
    taxa = '(' + ' OR '.join('organism_id:' + t for t in biology['taxa']) + ')'
    families = '(' + ' OR '.join('xref:interpro-' + f['interpro']
                                 for f in biology['families'].values()) + ')'
    return f'reviewed:false AND length:[1 TO 100] AND {taxa} AND {families}'


def outer_query(tax_id, config):
    return (f'tax_eq({tax_id}) AND organelle="mitochondrion" AND '
            f'base_count>={config["outer_parent_min_nt"]} AND '
            f'base_count<={config["outer_parent_max_nt"]} AND '
            '(base_count<14000 OR base_count>20000)')


def stable_accessions(record):
    """Accession-version changes and database mirrors do not create new samples."""
    found = set()
    for p in record.get('provenance', []):
        for key in ('accession', 'retrieved_accession', 'parent_accession',
                    'retrieved_parent_accession', 'index_parent_accession'):
            for value in (p.get(key) or '').split(','):
                if value:
                    found.add(value.strip().rsplit('.', 1)[0])
    return found


def source_qualifiers(record):
    return [dict(f.qualifiers) for f in record.features if f.type == 'source']


def guard_frozen_manifest(path, identity):
    if path.exists() and read_json(path)['identity'] != identity:
        raise ValueError('Refusing to change a frozen test cohort; create a new benchmark version')


def split_untouched(records, training):
    blocked = {r['sequence_sha256'] for r in training}
    sources = set().union(*(stable_accessions(r) for r in training))
    kept, exclusions = [], []
    for record in records:
        overlaps = sorted(stable_accessions(record) & sources)
        reason = ('training_rna_duplicate' if record['sequence_sha256'] in blocked else
                  'training_source_accession_overlap' if overlaps else None)
        if reason:
            exclusions.append({'sequence_sha256': record['sequence_sha256'],
                'family': record['family'], 'tax_id': record['tax_id'],
                'reason': reason, 'shared_stable_accessions': overlaps})
        else:
            kept.append(record)
    return kept, exclusions


def novelty_descriptors(records, training):
    """Distance annotations never select records or inspect model performance."""
    aligner = PairwiseAligner(mode='global', match_score=0, mismatch_score=-1,
                              open_gap_score=-1, extend_gap_score=-1)
    for row in records:
        same = [r for r in training if r['family'] == row['family']]
        distances = [(int(round(-aligner.score(row['rna'], r['rna']))),
                      r['sequence_sha256'], len(r['rna'])) for r in same]
        distance, nearest, length = min(distances)
        peptide_distance = min(int(round(-aligner.score(row['protein'], r['protein']))) for r in same)
        matched = [r for r in same if r['tax_id'] == row['tax_id'] and len(r['rna']) == len(row['rna'])]
        nt_hamming = min((sum(a != b for a, b in zip(row['rna'], r['rna'])) for r in matched), default=None)
        peptide_hamming = min((sum(a != b for a, b in zip(row['protein'], r['protein']))
                              for r in matched if len(r['protein']) == len(row['protein'])), default=None)
        row['novelty'] = {'comparison': 'all frozen training-pool CDS in the same annotated family',
            'nearest_training_sequence_sha256': nearest,
            'minimum_nt_edit_distance': distance,
            'nearest_nt_normalized_edit_distance': distance / max(length, len(row['rna'])),
            'minimum_peptide_edit_distance': peptide_distance,
            'exact_training_peptide': peptide_distance == 0,
            'same_family_taxon_length_training_sequences': len(matched),
            'same_family_taxon_length_minimum_nt_hamming': nt_hamming,
            'same_family_taxon_length_maximum_nt_identity': None if nt_hamming is None else 1 - nt_hamming / len(row['rna']),
            'same_family_taxon_length_minimum_peptide_hamming': peptide_hamming,
            'test_peptide_group_sha256': hashlib.sha256(row['protein'].encode()).hexdigest()}
    return records


def unreviewed_targets(entries, biology, references):
    reference_by_pair = {(r['tax_id'], r['family']): r for r in references}
    targets, rejected = defaultdict(list), []
    for entry in entries:
        accession = entry['primaryAccession']
        tax = int(entry['organism']['taxonId'])
        if str(tax) not in biology['taxa'] or entry['entryType'] != 'UniProtKB unreviewed (TrEMBL)':
            raise ValueError('Unexpected organism/review status in retrieval query')
        crossrefs = entry.get('uniProtKBCrossReferences', [])
        interpro = {x['id'] for x in crossrefs if x['database'] == 'InterPro'}
        families = [f for f, spec in biology['families'].items() if spec['interpro'] in interpro]
        if len(families) != 1:
            raise ValueError(f'Ambiguous family retrieval lead: {accession}')
        ref = reference_by_pair.get((tax, families[0]))
        if ref is None:
            rejected.append({'uniprot': accession, 'tax_id': tax, 'family': families[0],
                             'reason': 'no_frozen_own_species_reviewed_reference'})
            continue
        found = 0
        for crossref in crossrefs:
            if crossref['database'] != 'EMBL':
                continue
            props = {p['key']: p['value'] for p in crossref.get('properties', [])}
            pid = props.get('ProteinId', '')
            if pid in ('', '-'):
                continue
            found += 1
            targets[pid].append({'reference': ref, 'retrieval_uniprot': accession,
                'parent_accession': crossref['id'], 'cross_reference_status': props.get('Status', ''),
                'lead_fragment_flag': entry.get('proteinDescription', {}).get('flag')})
        if not found:
            rejected.append({'uniprot': accession, 'tax_id': tax, 'family': families[0],
                             'reason': 'no_embl_protein_link'})
    return dict(targets), rejected


def acquire_gap(root, include_outer_parents=False, freeze=False):
    root = Path(root).resolve()
    config_path = root / 'configs/gap-biology.json'
    config = read_json(config_path)
    biology_path = root / config['training_biology_file']
    biology = read_json(biology_path)
    training_path = root / config['training_pool_file']
    if sha256(training_path) != config['training_pool_sha256']:
        raise ValueError('Frozen training union changed')
    training = [json.loads(line) for line in training_path.read_text().splitlines()]
    if len(training) != config['training_pool_sequences']:
        raise ValueError('Unexpected frozen training union count')
    refs_path = root / config['reference_file']
    references = read_json(refs_path)
    if sha256(refs_path) != read_json(root / 'reports/secondary/genome-acquisition.json')['reference_inventory_sha256']:
        raise ValueError('Frozen reviewed references changed')
    panel = biology['mammals'] + biology['nonmammalian_vertebrates']
    raw, output = root / 'data/raw/gap', root / 'data/processed/gap'
    reports = root / 'reports/gap'
    dl = Download(raw)
    source = dl.get('unreviewed-leads.json', UNIPROT_URL,
                    {'query': unreviewed_query(biology), 'format': 'json'})
    entries = read_json(source)['results']
    targets, lead_rejections = unreviewed_targets(entries, biology, references)
    write_canonical_json(raw / 'unreviewed-targets.json', targets)
    accepted, rejected, missing = [], [], []
    source_hash = sha256(source)
    ids = sorted(targets)
    for start in range(0, len(ids), 50):
        for record, path in fetch_records(dl, ids[start:start + 50], missing):
            path_hash = sha256(path)
            for context in targets[record.id]:
                ref = context['reference']
                result, reason = validate_secondary_record(record, ref,
                    min_identity=biology['min_reference_identity'], max_aa=biology['max_amino_acids'])
                if reason:
                    rejected.append({'accession': record.id, 'tax_id': ref['tax_id'],
                        'family': ref['family'], 'reason': reason,
                        'source_file': path.relative_to(root).as_posix(), 'source_file_sha256': path_hash})
                    continue
                parents = sorted({p.ref for f in record.features if f.type == 'CDS'
                                  for p in f.location.parts if p.ref})
                result['provenance'] = [{'accession': record.id, 'retrieved_accession': record.id,
                    'parent_accession': ','.join(parents) or context['parent_accession'],
                    'index_parent_accession': context['parent_accession'],
                    'tax_id': ref['tax_id'], 'reference_uniprot': ref['uniprot'],
                    'reference_identity': result['reference_identity'],
                    'retrieval_uniprot': context['retrieval_uniprot'],
                    'cross_reference_status': context['cross_reference_status'],
                    'source_qualifiers': source_qualifiers(record),
                    'route': 'unreviewed_uniprot_embl_source_cds',
                    'source_file': path.relative_to(root).as_posix(), 'source_file_sha256': path_hash,
                    'index_file': source.relative_to(root).as_posix(), 'index_file_sha256': source_hash}]
                accepted.append(result)
        print(f'Gap CDS leads {min(start+50,len(ids))}/{len(ids)}; accepted routes {len(accepted)}', flush=True)
    inventories, outer_targets = [], {}
    if include_outer_parents:
        atp8_refs = {r['tax_id']: r for r in references if r['family'] == 'ATP8'}
        for tax in panel:
            query = outer_query(tax, config)
            index = dl.get(f'outer-index/{tax}.tsv', ENA_SEARCH,
                {'result': 'sequence', 'query': query, 'fields': GENOME_FIELDS, 'format': 'tsv', 'limit': '0'})
            count = dl.get(f'outer-index/{tax}-count.json', 'https://www.ebi.ac.uk/ena/portal/api/count',
                {'result': 'sequence', 'query': query, 'format': 'json'})
            rows = list(csv.DictReader(index.open(encoding='utf-8'), delimiter='\t'))
            if len(rows) != int(read_json(count)['count']):
                raise ValueError('Outside-range parent inventory count mismatch')
            selected = [row for row in rows if any(term in row['description'].lower()
                        for term in config['outer_parent_description_terms'])]
            inventories.append({'tax_id': tax, 'query': query, 'inventory_records': len(rows),
                'selected_description_records': len(selected),
                'selected_total_nt': sum(int(row['base_count']) for row in selected),
                'index_file': index.relative_to(root).as_posix(), 'index_file_sha256': sha256(index)})
            for row in selected:
                if int(row['tax_id']) != tax or 14000 <= int(row['base_count']) <= 20000:
                    raise ValueError('Unexpected outside-range parent identity')
                acc = versioned_accession(row)
                if acc in outer_targets:
                    raise ValueError('Duplicate parent accession in inventory')
                outer_targets[acc] = (row, inventories[-1])
        write_canonical_json(reports / 'outer-inventory.json', {'inventories': inventories,
            'selected_parents': len(outer_targets),
            'selected_total_nt': sum(int(v[0]['base_count']) for v in outer_targets.values())})
        print('Outside-range inventory: ' + json.dumps(inventories), flush=True)
        # Metadata is inexpensive; stop before unbounded payload retrieval.
        if len(outer_targets) > 10000 or sum(int(v[0]['base_count']) for v in outer_targets.values()) > 100_000_000:
            raise ValueError('Outside-range acquisition exceeds the declared bounded payload budget')
        parent_ids = sorted(outer_targets)
        for start in range(0, len(parent_ids), 100):
            for parent, path in fetch_records(dl, parent_ids[start:start + 100], missing):
                row, inv = outer_targets[parent.id]
                ref = atp8_refs[int(row['tax_id'])]
                result, reason = extract_atp8(parent, ref, row, biology['min_reference_identity'], biology['max_amino_acids'])
                if reason:
                    rejected.append({'parent_accession': parent.id, 'tax_id': ref['tax_id'],
                        'family': 'ATP8', 'reason': reason,
                        'source_file': path.relative_to(root).as_posix(), 'source_file_sha256': sha256(path)})
                    continue
                result['provenance'] = [{'accession': result.pop('accession'), 'parent_accession': parent.id,
                    'tax_id': ref['tax_id'], 'reference_uniprot': ref['uniprot'],
                    'reference_identity': result['reference_identity'],
                    'original_feature_location': result.pop('original_feature_location'),
                    'feature_parts': result.pop('feature_parts'), 'parent_sequence_md5': row['sequence_md5'],
                    'source_qualifiers': source_qualifiers(parent),
                    'route': 'ena_outer_length_parent_feature',
                    'source_file': path.relative_to(root).as_posix(), 'source_file_sha256': sha256(path),
                    'index_file': inv['index_file'], 'index_file_sha256': inv['index_file_sha256']}]
                accepted.append(result)
            print(f'Gap outside parents {min(start+100,len(parent_ids))}/{len(parent_ids)}', flush=True)
    pool = merge_records(accepted, panel)
    untouched, blocked = split_untouched(pool, training)
    novelty_descriptors(untouched, training)
    payload = ''.join(json.dumps(r, sort_keys=True) + '\n' for r in untouched).encode()
    digest = hashlib.sha256(payload).hexdigest()
    identity = {'config_sha256': sha256(config_path), 'training_biology_sha256': sha256(biology_path),
        'references_sha256': sha256(refs_path), 'training_pool_sha256': sha256(training_path),
        'data_card_sha256': sha256(root / 'docs/gap-data-card.md'),
        'cohort_sha256': digest, 'include_outer_parents': include_outer_parents,
        'raw_source_hashes': {p.relative_to(root).as_posix().removesuffix('.source.json'):
                             read_json(p)['sha256'] for p in sorted(raw.rglob('*.source.json'))},
        'source_hashes': {p.relative_to(root).as_posix(): sha256(p) for p in [
            Path(__file__), root / 'src/microprotein_lm/secondary_acquire.py',
            root / 'src/microprotein_lm/secondary_genomes.py', root / 'src/microprotein_lm/acquire.py']}}
    frozen = output / 'manifest.json'
    guard_frozen_manifest(frozen, identity)
    if freeze and not untouched:
        raise ValueError('No genuinely untouched eligible CDS; cannot freeze an empty test cohort')
    output.mkdir(parents=True, exist_ok=True)
    cohort_path = output / 'cohort.jsonl'
    if not cohort_path.exists() or sha256(cohort_path) != digest:
        cohort_path.write_bytes(payload)
    for row in missing:
        if 'source_file' in row:
            row['source_file'] = Path(row['source_file']).relative_to(root).as_posix()
    write_canonical_json(raw / 'rejections.json', rejected)
    write_canonical_json(raw / 'training-overlap-exclusions.json', blocked)
    write_canonical_json(raw / 'unavailable.json', missing)
    report = {'created_utc': now(), 'stage': 'untouched_test_acquisition', 'identity': identity,
        'frozen': frozen.exists() or freeze, 'unreviewed_leads': len(entries),
        'unreviewed_accessions_requested': len(ids), 'lead_exclusions': lead_rejections,
        'outer_parent_inventory': inventories, 'outer_parents_requested': len(outer_targets),
        'source_qc_rejections': dict(Counter(r['reason'] for r in rejected)),
        'unavailable_count': len(missing), 'eligible_before_training_exclusion': len(pool),
        'training_overlap_exclusions': dict(Counter(r['reason'] for r in blocked)),
        'untouched_cds': len(untouched), 'unique_peptides': len({r['protein'] for r in untouched}),
        'same_peptide_as_training_cds': sum(r['novelty']['exact_training_peptide'] for r in untouched),
        'family_taxon_counts': dict(Counter(f'{r["family"]}:{r["tax_id"]}' for r in untouched)),
        'cohort_file': cohort_path.relative_to(root).as_posix(),
        'sources': [{**read_json(p), 'file': p.relative_to(root).as_posix().removesuffix('.source.json')}
                    for p in sorted(raw.rglob('*.source.json'))],
        'limits': 'Exact RNA and stable source accessions are disjoint from all 533 training-pool CDS. '
            'Close homologs and identical peptides can remain; novelty distances are descriptive, not performance-driven filters. '
            'Same species/family does not imply independent individuals or experimentally equivalent function. '
            'No test prediction has been run by this acquisition process.'}
    report_path = reports / 'acquisition.json'
    if frozen.exists():
        existing = read_json(report_path)
        if existing['identity'] != identity:
            raise ValueError('Frozen test report identity mismatch')
    else:
        write_canonical_json(report_path, report)
    if freeze and not frozen.exists():
        write_canonical_json(frozen, {'created_utc': now(), 'stage': 'untouched_test',
            'identity': identity, 'n_sequences': len(untouched), 'acquisition_report_sha256': sha256(report_path),
            'grouping': 'RNA SHA256 for repeated gap scenarios; peptide SHA256 for conservative related-sequence summaries'})
    print(json.dumps({k: report[k] for k in ['unreviewed_leads', 'unreviewed_accessions_requested',
        'outer_parents_requested', 'eligible_before_training_exclusion', 'training_overlap_exclusions',
        'untouched_cds', 'same_peptide_as_training_cds', 'family_taxon_counts']}, indent=2), flush=True)
    return untouched
