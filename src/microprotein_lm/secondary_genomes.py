"""Recover annotated ATP8 CDS from the declared whole-mitochondrial-genome panel.

This third ascertainment route extracts annotated features from parent genomes;
it never reverse translates, guesses reading frames, or repairs missing stops.
"""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path

from Bio.SeqFeature import ExactPosition
from Bio.SeqRecord import SeqRecord

from .acquire import Download, ENA_SEARCH
from .io import now, read_json, sha256
from .secondary_acquire import (
    fetch_records, match_metadata, merge_records, preserve_identical_report,
    validate_secondary_record, write_canonical_json, write_frozen_bytes,
)


GENOME_FIELDS = 'accession,sequence_version,base_count,sequence_md5,tax_id,description,organelle'


def genome_query(tax_id):
    return f'tax_eq({tax_id}) AND organelle="mitochondrion" AND base_count>=14000 AND base_count<=20000'


def versioned_accession(row):
    accession = row['accession']
    if '.' in accession:
        if accession.rsplit('.', 1)[1] != row['sequence_version']:
            raise ValueError('Parent accession and sequence version disagree')
        return accession
    return accession + '.' + row['sequence_version']


def extract_atp8(parent, reference, metadata, min_identity=0.95, max_aa=100):
    """Validate parent identity, then use BioPython's annotated strand/join logic."""
    if parent.id != versioned_accession(metadata):
        return None, 'parent_accession_mismatch'
    if str(reference['tax_id']) != str(metadata['tax_id']):
        return None, 'parent_metadata_taxon_mismatch'
    if not parent.seq.defined:
        return None, 'undefined_parent_sequence'
    dna = str(parent.seq).upper()
    if len(dna) != int(metadata['base_count']):
        return None, 'parent_length_mismatch'
    if metadata.get('sequence_md5') and hashlib.md5(dna.encode()).hexdigest() != metadata['sequence_md5']:
        return None, 'parent_hash_mismatch'
    sources = [f for f in parent.features if f.type == 'source']
    source_taxa = {v for f in sources for v in f.qualifiers.get('db_xref', []) if v.startswith('taxon:')}
    if source_taxa != {f'taxon:{reference["tax_id"]}'}:
        return None, 'wrong_taxon'
    matches = []
    for feature in parent.features:
        if feature.type != 'CDS':
            continue
        q = feature.qualifiers
        row = {'gene': ';'.join(q.get('gene', [])), 'gene_synonym': ';'.join(q.get('gene_synonym', [])),
               'product': ';'.join(q.get('product', [])), 'protein_id': q.get('protein_id', [''])[0]}
        if match_metadata(row, [reference]):
            matches.append(feature)
    if not matches:
        return None, 'missing_annotated_atp8'
    if len(matches) != 1:
        return None, 'multiple_annotated_atp8'
    feature = matches[0]
    accession = feature.qualifiers.get('protein_id', [''])[0]
    if not accession:
        return None, 'missing_protein_accession'
    if feature.location is None or any(p.ref for p in feature.location.parts):
        return None, 'remote_or_missing_cds_location'
    if getattr(feature.location, 'operator', 'join') != 'join':
        return None, 'unsupported_location_operator'
    for part in feature.location.parts:
        if not isinstance(part.start, ExactPosition) or not isinstance(part.end, ExactPosition):
            return None, 'partial_location'
        if not 0 <= int(part.start) < int(part.end) <= len(parent.seq):
            return None, 'out_of_bounds_cds_location'
        if part.strand not in (-1, 1):
            return None, 'unknown_cds_strand'
    try:
        cds = feature.extract(parent.seq)
    except (ValueError, TypeError):
        return None, 'invalid_cds_location'
    derived = SeqRecord(cds, id=accession)
    # Keep original locations so incomplete endpoints remain visible to the same
    # complete-CDS gate used by the first two routes; that gate does not extract.
    derived.features = sources + [feature]
    result, reason = validate_secondary_record(derived, reference, min_identity=min_identity, max_aa=max_aa)
    if reason:
        return None, reason
    result['accession'] = accession
    result['original_feature_location'] = str(feature.location)
    result['feature_parts'] = [{'start_0based': int(p.start), 'end_0based_exclusive': int(p.end),
                               'strand': p.strand} for p in feature.location.parts]
    return result, None


def _fetch_batch(raw, ids):
    missing = []
    records = list(fetch_records(Download(raw), ids, missing))
    return records, missing


def enrich_genomes(root, base_records=None, batch_size=200, workers=2):
    """Return all three routes merged; preserve the original two-route snapshot."""
    if not 1 <= batch_size <= 500 or not 1 <= workers <= 2:
        raise ValueError('Bounded genome retrieval requires batch<=500 and at most two workers')
    root = Path(root).resolve()
    biology_path = root / 'configs/secondary/biology.json'
    biology = read_json(biology_path)
    panel = biology['mammals'] + biology['nonmammalian_vertebrates']
    secondary = root / 'data/raw/secondary'
    raw = secondary / 'genomes'
    report_dir = root / 'reports/secondary'
    base_path = secondary / 'eligible-records.jsonl'
    base_manifest = read_json(report_dir / 'acquisition.json')
    if base_manifest['eligible_records_sha256'] != sha256(base_path):
        raise ValueError('Two-route base corpus changed since its acquisition report')
    if base_records is None:
        base_records = [json.loads(line) for line in base_path.read_text(encoding='utf-8').splitlines()]
    expected_base = [json.loads(line) for line in base_path.read_text(encoding='utf-8').splitlines()]
    if base_records != expected_base:
        raise ValueError('Provided base records differ from the preserved two-route snapshot')
    all_refs = read_json(secondary / 'references.json')
    references = {r['tax_id']: r for r in all_refs if r['family'] == 'ATP8'}
    if set(references) != set(panel):
        raise ValueError('Missing an own-species ATP8 reference in the declared panel')
    download = Download(raw)
    metadata = {}
    alias_groups = {}
    inventories = []
    for tax_id in panel:
        query = genome_query(tax_id)
        path = download.get(f'index/{tax_id}.tsv', ENA_SEARCH,
            {'result': 'sequence', 'query': query, 'fields': GENOME_FIELDS, 'format': 'tsv', 'limit': '0'})
        count_path = download.get(f'index/{tax_id}-count.json', 'https://www.ebi.ac.uk/ena/portal/api/count',
            {'result': 'sequence', 'query': query, 'format': 'json'})
        expected_count = int(read_json(count_path)['count'])
        groups = defaultdict(list)
        seen = set()
        with path.open(encoding='utf-8') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                accession = versioned_accession(row)
                if accession in seen:
                    raise ValueError(f'Duplicate parent in inventory: {accession}')
                seen.add(accession)
                if str(row['tax_id']) != str(tax_id) or not 14000 <= int(row['base_count']) <= 20000:
                    raise ValueError(f'Unexpected parent metadata: {accession}')
                groups[row['sequence_md5'] or accession].append(row)
        if len(seen) != expected_count:
            raise ValueError(f'Parent inventory count mismatch for {tax_id}: {len(seen)} != {expected_count}')
        index_hash = sha256(path)
        for rows in groups.values():
            for row in rows:
                accession = versioned_accession(row)
                if accession in metadata:
                    raise ValueError(f'Parent accession assigned to more than one taxon: {accession}')
                metadata[accession] = row
                # Identical parent nucleotides can have different annotations.
                # Validate every accession; deduplicate only the resulting CDS.
                alias_groups[accession] = {'rows': [row], 'index_file': path.relative_to(root).as_posix(),
                                           'index_file_sha256': index_hash}
        inventories.append({'tax_id': tax_id, 'query': query, 'rows': len(seen),
                            'count_verified': True, 'unique_parent_hashes': len(groups),
                            'index_file': path.relative_to(root).as_posix(), 'index_file_sha256': index_hash})
        print(f'genome index {tax_id}: {len(seen)} parents / {len(groups)} whole-genome hashes', flush=True)
    ids = sorted(metadata)
    write_canonical_json(raw / 'parent-targets.json', {'targets': metadata, 'aliases': alias_groups})
    chunks = [ids[i:i + batch_size] for i in range(0, len(ids), batch_size)]
    accepted, rejected, missing = [], [], []
    source_hashes = {}
    complete = 0
    unique_rna = set()
    # Only workers batches are outstanding, preventing completed BioPython record
    # objects from accumulating in an unbounded future queue.
    with ThreadPoolExecutor(max_workers=workers) as pool:
        iterator = iter(chunks)
        pending = {}
        for _ in range(workers):
            chunk = next(iterator, None)
            if chunk:
                pending[pool.submit(_fetch_batch, raw, chunk)] = chunk
        while pending:
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                chunk = pending.pop(future)
                fetched, unavailable = future.result()
                for unavailable_record in unavailable:
                    if 'source_file' in unavailable_record:
                        unavailable_record['source_file'] = Path(unavailable_record['source_file']).relative_to(root).as_posix()
                missing.extend(unavailable)
                for parent, path in fetched:
                    if path not in source_hashes:
                        source_hashes[path] = sha256(path)
                    row = metadata[parent.id]
                    ref = references[int(row['tax_id'])]
                    result, reason = extract_atp8(parent, ref, row, biology['min_reference_identity'],
                                                 biology['max_amino_acids'])
                    if reason:
                        rejected.append({'parent_accession': parent.id, 'tax_id': ref['tax_id'],
                            'reason': reason, 'source_file': path.relative_to(root).as_posix(),
                            'source_file_sha256': source_hashes[path],
                            'metadata_aliases': len(alias_groups[parent.id]['rows'])})
                        continue
                    result['provenance'] = []
                    group = alias_groups[parent.id]
                    for alias in group['rows']:
                        alias_id = versioned_accession(alias)
                        result['provenance'].append({
                            'accession': result['accession'] if alias_id == parent.id else None,
                            'retrieved_accession': result['accession'], 'parent_accession': alias_id,
                            'retrieved_parent_accession': parent.id, 'tax_id': ref['tax_id'],
                            'reference_uniprot': ref['uniprot'], 'reference_identity': result['reference_identity'],
                            'route': 'ena_parent_genome_feature',
                            'original_feature_location': result['original_feature_location'],
                            'feature_parts': result['feature_parts'],
                            'source_file': path.relative_to(root).as_posix(), 'source_file_sha256': source_hashes[path],
                            'index_file': group['index_file'], 'index_file_sha256': group['index_file_sha256'],
                            'parent_sequence_md5': alias['sequence_md5'],
                            'sequence_equivalent_metadata_alias': alias_id != parent.id})
                    for temporary_key in ['accession', 'original_feature_location', 'feature_parts']:
                        del result[temporary_key]
                    accepted.append(result)
                    unique_rna.add(result['sequence_sha256'])
                complete += len(chunk)
                print(f'genomes {complete}/{len(ids)}; unique eligible ATP8 {len(unique_rna)}; '
                      f'rejected parents {len(rejected)}; unavailable {len(missing)}', flush=True)
                next_chunk = next(iterator, None)
                if next_chunk:
                    pending[pool.submit(_fetch_batch, raw, next_chunk)] = next_chunk
    records = merge_records(base_records + accepted, panel)
    output = raw / 'enriched-records.jsonl'
    payload = ''.join(json.dumps(r, sort_keys=True) + '\n' for r in records).encode('utf-8')
    rejected.sort(key=lambda r: r['parent_accession'])
    missing.sort(key=lambda r: r['accession'])
    write_canonical_json(raw / 'rejections.json', rejected)
    write_canonical_json(raw / 'unavailable-accessions.json', missing)
    manifest = {'created_utc': now(), 'stage': 'training_only', 'route': 'ena_parent_genome_feature',
        'biology_sha256': sha256(biology_path), 'base_eligible_records_sha256': sha256(base_path),
        'reference_inventory_sha256': sha256(secondary / 'references.json'),
        'metadata_inventory': inventories, 'parent_inventory_rows': sum(i['rows'] for i in inventories),
        'requested_parent_accessions': len(ids), 'accepted_parent_records': len(accepted),
        'unavailable_accessions': missing, 'record_rejections': dict(Counter(r['reason'] for r in rejected)),
        'eligible_records_file': output.relative_to(root).as_posix(), 'eligible_records_sha256': hashlib.sha256(payload).hexdigest(),
        'base_unique_cds': len(base_records), 'globally_unique_cds': len(records),
        'additional_unique_cds': len(records) - len(base_records),
        'by_primary_taxon_family': dict(Counter(f'{r["tax_id"]}/{r["family"]}' for r in records)),
        'by_observed_taxon_family': dict(Counter(f'{t}/{r["family"]}' for r in records for t in r['taxa'])),
        'retained_raw_bytes': sum(p.stat().st_size for p in source_hashes),
        'source_registry': [{**read_json(p), 'file': p.relative_to(raw).as_posix().removesuffix('.source.json')}
                            for p in sorted(raw.rglob('*.source.json'))],
        'representative_policy': 'Every indexed parent accession version is fetched and feature-validated independently, '
            'including parents sharing identical DNA; deduplicate only the final eligible RNA and retain all provenance.',
        'scope_limit': 'Annotated ATP8 from complete-length mitochondrial parent records in the declared 14,000–20,000nt index, '
            'plus the two prior routes. This is a declared-route census, not proof of all possible database alleles.'}
    report_path = report_dir / 'genome-acquisition.json'
    identity_keys = ['biology_sha256', 'base_eligible_records_sha256', 'reference_inventory_sha256',
                     'eligible_records_sha256', 'metadata_inventory', 'globally_unique_cds',
                     'by_primary_taxon_family', 'record_rejections']
    preserve_identical_report(report_path, manifest, identity_keys, check_only=True)
    write_frozen_bytes(root, output, payload, report_path)
    preserve_identical_report(report_path, manifest, identity_keys)
    print(json.dumps({k: manifest[k] for k in ['parent_inventory_rows', 'requested_parent_accessions',
        'globally_unique_cds', 'additional_unique_cds', 'by_primary_taxon_family', 'record_rejections']}, indent=2), flush=True)
    return records


if __name__ == '__main__':
    enrich_genomes(Path.cwd())
