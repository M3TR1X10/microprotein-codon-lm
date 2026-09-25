"""Observed CDS acquisition for the preregistered ATP-synthase secondary pilot.

The original pilot remains untouched. Search-index records and UniProt's complete
EMBL cross-reference lists are complementary ascertainment routes, not samples.
"""
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import openpyxl
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition

from .acquire import Download, ENA_SEARCH, FIELDS, UNIPROT_URL, normalize_name
from .io import now, read_json, sha256


def write_canonical_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes((json.dumps(value, indent=2, sort_keys=True) + '\n').encode('utf-8'))


def write_frozen_bytes(root, path, payload, report_path):
    if (root / 'reports/secondary/plan.json').exists():
        digest = hashlib.sha256(payload).hexdigest()
        if report_path.exists() and read_json(report_path)['eligible_records_sha256'] != digest:
            raise ValueError('Frozen acquisition corpus differs from its archived report')
        if path.exists() and sha256(path) != digest:
            raise ValueError('Refusing to overwrite a frozen acquisition corpus with changed records')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def preserve_identical_report(path, manifest, identity_keys, check_only=False):
    """A reproduced corpus does not create a new experimental retrieval receipt.

    Keep the archived report byte-for-byte when scientific identity is unchanged;
    fresh download receipts remain in the new raw cache for separate inspection.
    """
    if path.exists():
        archived = read_json(path)
        if all(archived.get(k) == manifest.get(k) for k in identity_keys):
            return archived
        if (path.parent / 'plan.json').exists():
            raise ValueError('Refusing to change an acquisition report after the training plan was frozen')
    if not check_only:
        write_canonical_json(path, manifest)
    return manifest


def split_names(value):
    return [v.strip() for v in re.split(r'[;,]', str(value or '')) if v.strip()]


def verify_receipt(path):
    source = read_json(path.with_name(path.name + '.source.json'))
    if source['sha256'] != sha256(path):
        raise ValueError(f'Source hash mismatch: {path}')
    return source


def resolve_impi(path, biology):
    """Resolve compound symbols without collapsing distinct genes by substring."""
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    rows = iter(workbook['IMPI-2021Q4pre'].values)
    header = next(rows)
    symbols = defaultdict(list)
    for values in rows:
        row = dict(zip(header, values))
        for name in split_names(row['Symbol']):
            symbols[normalize_name(name)].append(row)
    result = {}
    for family, specification in biology['families'].items():
        matches = [r for r in symbols[normalize_name(specification['human_gene'])]
                   if str(r['IMPI Class']).startswith('Verified mitochondrial')]
        if len(matches) != 1:
            raise ValueError(f'Expected one verified IMPI row for {family}: {len(matches)}')
        row = matches[0]
        result[family] = {key: row[key] for key in ['Symbol', 'Name', 'Ensembl', 'IMPI Class']}
    workbook.close()
    return result


def reference_inventory(entries, biology, impi):
    selected, excluded = [], []
    for entry in entries:
        accession = entry['primaryAccession']
        if accession in biology.get('excluded_uniprot', []):
            excluded.append({'accession': accession, 'reason': 'declared_exclusion'})
            continue
        if entry.get('entryType') != 'UniProtKB reviewed (Swiss-Prot)':
            raise ValueError(f'Unreviewed reference returned: {accession}')
        if entry.get('proteinDescription', {}).get('flag') == 'Fragment':
            excluded.append({'accession': accession, 'reason': 'fragment'})
            continue
        tax_id = int(entry['organism']['taxonId'])
        if str(tax_id) not in biology['taxa']:
            raise ValueError(f'Unexpected reference taxon: {tax_id}')
        crossrefs = entry.get('uniProtKBCrossReferences', [])
        interpro = {x['id'] for x in crossrefs if x['database'] == 'InterPro'}
        families = [name for name, spec in biology['families'].items()
                    if spec['interpro'] in interpro]
        if not families:
            excluded.append({'accession': accession, 'reason': 'outside_declared_families'})
            continue
        if len(families) != 1:
            raise ValueError(f'Ambiguous InterPro family: {accession}')
        family = families[0]
        if tax_id == 9606 and accession != biology['families'][family]['human_uniprot']:
            excluded.append({'accession': accession, 'reason': 'not_declared_human_reference'})
            continue
        sequence = entry['sequence']['value']
        if not 1 <= len(sequence) <= biology['max_amino_acids']:
            raise ValueError(f'Reference outside length query: {accession}')
        genes = sorted({v['value'] for g in entry.get('genes', [])
                        for key in ['geneName', 'synonyms', 'orderedLocusNames', 'orfNames']
                        for v in ([g[key]] if key == 'geneName' and key in g else
                                  g.get(key, []) if key != 'geneName' else [])})
        description = entry.get('proteinDescription', {})
        names = [n['fullName']['value'] for n in
                 [description.get('recommendedName', {})] + description.get('alternativeNames', [])
                 if 'fullName' in n]
        names += biology.get('extra_product_aliases', {}).get(family, [])
        if tax_id == 9606:
            genes += split_names(impi[family]['Symbol'])
            names += split_names(impi[family]['Name'])
        embl = []
        for x in crossrefs:
            if x['database'] != 'EMBL':
                continue
            properties = {p['key']: p['value'] for p in x.get('properties', [])}
            protein_id = properties.get('ProteinId', '')
            if protein_id not in ('', '-'):
                embl.append({'protein_id': protein_id, 'parent_accession': x['id'],
                             'cross_reference_status': properties.get('Status', '')})
        selected.append({'family': family, 'tax_id': tax_id, 'uniprot': accession,
                         'reference_protein': sequence, 'genes': sorted(set(genes)),
                         'products': sorted(set(names)), 'embl': sorted(embl, key=lambda x: (
                             x['protein_id'], x['parent_accession'], x['cross_reference_status'])),
                         'translation_table': biology['families'][family]['genetic_code']})
    selected.sort(key=lambda r: (r['tax_id'], r['family'], r['uniprot']))
    if len({(r['tax_id'], r['family']) for r in selected}) != len(selected):
        raise ValueError('Multiple reviewed references per family/species need a declared resolution')
    if {r['family'] for r in selected if r['tax_id'] == 9606} != set(biology['families']):
        raise ValueError('Missing declared human references')
    return selected, excluded


def index_query(tax_id, references):
    # Quoted exact query terms provide discovery; normalized exact matching below
    # remains the eligibility gate if the portal performs broader text matching.
    terms = set()
    for ref in references:
        for gene in ref['genes']:
            terms.add(('gene', gene))
            terms.add(('gene_synonym', gene))
        terms.update(('product', product) for product in ref['products'])
    if not terms:
        raise ValueError('Cannot query an empty reference set')
    if any('"' in value or '\\' in value for _, value in terms):
        raise ValueError('Unsafe quoted ENA query term')
    expression = ' OR '.join(f'{field}="{value}"' for field, value in sorted(terms))
    return f'tax_eq({tax_id}) AND base_count<=303 AND ({expression})'


def metadata_reason(row, tax_id):
    if str(row['tax_id']) != str(tax_id):
        return 'wrong_taxon'
    if row.get('codon_start') not in ('', '1'):
        return 'partial_frame'
    for key in ['pseudo', 'pseudo_gene', 'exception', 'transl_except',
                'artificial_location', 'ribosomal_slippage']:
        if row.get(key, '').lower() not in ('', 'false', '0'):
            return key
    return None


def match_metadata(row, refs):
    genes = {normalize_name(n) for n in split_names(row['gene'] + ';' + row['gene_synonym'])}
    product = normalize_name(row['product'])
    matches = [r for r in refs if genes.intersection(map(normalize_name, r['genes']))
               or product in {normalize_name(n) for n in r['products'] + r['genes']}]
    if len(matches) > 1:
        raise ValueError(f'Ambiguous metadata family mapping: {row["protein_id"]}')
    return matches[0] if matches else None


def validate_secondary_record(record, reference, metadata=None, min_identity=0.95, max_aa=100):
    features = [f for f in record.features if f.type == 'CDS']
    if len(features) != 1:
        return None, 'not_one_cds'
    feature = features[0]
    q = feature.qualifiers
    if q.get('protein_id', [''])[0] != record.id or (
            metadata is not None and metadata['protein_id'] != record.id):
        return None, 'accession_mismatch'
    source_taxa = {v for f in record.features if f.type == 'source'
                   for v in f.qualifiers.get('db_xref', []) if v.startswith('taxon:')}
    if source_taxa != {f'taxon:{reference["tax_id"]}'}:
        return None, 'wrong_taxon'
    for key in ['pseudo', 'pseudogene', 'exception', 'transl_except',
                'ribosomal_slippage', 'artificial_location']:
        if key in q:
            return None, key
    if q.get('codon_start', ['1']) != ['1']:
        return None, 'partial_frame'
    if feature.location is None or any(not isinstance(p.start, ExactPosition) or
            not isinstance(p.end, ExactPosition) for p in feature.location.parts):
        return None, 'partial_location'
    dna = str(record.seq).upper()  # ENA coding endpoint has already oriented/spliced it.
    if not dna or set(dna) - set('ACGT'):
        return None, 'ambiguous_bases'
    if len(dna) % 3:
        return None, 'incomplete_codon'
    if metadata is not None:
        if len(dna) != int(metadata['base_count']):
            return None, 'length_metadata_mismatch'
        if metadata.get('sequence_md5') and metadata['sequence_md5'] != hashlib.md5(dna.encode()).hexdigest():
            return None, 'sequence_hash_mismatch'
    try:
        table = int(q.get('transl_table', ['1'])[0])
    except (ValueError, IndexError):
        return None, 'unexpected_genetic_code'
    if table != reference['translation_table']:
        return None, 'unexpected_genetic_code'
    try:
        protein = str(Seq(dna).translate(table=table, cds=True))
    except Exception:
        return None, 'invalid_complete_cds'
    if protein != ''.join(q.get('translation', [])):
        return None, 'translation_mismatch'
    if not 1 <= len(protein) <= max_aa:
        return None, 'not_microprotein'
    canonical = reference['reference_protein']
    if len(protein) != len(canonical):
        return None, 'reference_length_mismatch'
    identity = sum(a == b for a, b in zip(protein, canonical)) / len(canonical)
    if identity < min_identity:
        return None, 'reference_identity_below_threshold'
    rna = dna.replace('T', 'U')
    return {'rna': rna, 'sequence_sha256': hashlib.sha256(rna.encode()).hexdigest(),
            'protein': protein, 'family': reference['family'], 'tax_id': reference['tax_id'],
            'translation_table': table, 'reference_identity': identity}, None


def merge_records(records, panel):
    """Global nucleotide deduplication preserves cross-taxon accession evidence."""
    priority = {tax_id: i for i, tax_id in enumerate(panel)}
    grouped = {}
    for record in records:
        key = record['sequence_sha256']
        if key != hashlib.sha256(record['rna'].encode()).hexdigest():
            raise ValueError('Incorrect RNA hash before deduplication')
        if key not in grouped:
            grouped[key] = {k: v for k, v in record.items() if k != 'reference_identity'}
            grouped[key]['provenance'] = list(record['provenance'])
            grouped[key]['taxa'] = list(record.get('taxa', [record['tax_id']]))
        else:
            target = grouped[key]
            if any(target[k] != record[k] for k in ['family', 'translation_table', 'protein', 'rna']):
                raise ValueError('Ambiguous cross-family/code/protein RNA duplicate')
            target['provenance'].extend(record['provenance'])
            target['taxa'].extend(record.get('taxa', [record['tax_id']]))
    for record in grouped.values():
        record['taxa'] = sorted(set(record['taxa']), key=priority.__getitem__)
        record['tax_id'] = record['taxa'][0]
        unique = {json.dumps(p, sort_keys=True): p for p in record['provenance']}
        record['provenance'] = [unique[k] for k in sorted(unique)]
    return [grouped[k] for k in sorted(grouped)]


def fetch_records(download, accessions, missing):
    """Split failed/partially missing batches; never equate network failure to absence."""
    ids = sorted(set(accessions))
    if not ids:
        return
    digest = hashlib.sha256(','.join(ids).encode()).hexdigest()[:24]
    url = 'https://www.ebi.ac.uk/ena/browser/api/embl/' + ','.join(ids)
    try:
        path = download.get(f'cds/{digest}.embl', url)
    except requests.HTTPError as error:
        status = error.response.status_code
        if status not in (400, 404):
            raise
        if len(ids) == 1:
            missing.append({'accession': ids[0], 'url': url, 'status': status,
                            'retrieved_utc': now(), 'reason': 'browser_endpoint_unavailable'})
            return
        midpoint = len(ids) // 2
        yield from fetch_records(download, ids[:midpoint], missing)
        yield from fetch_records(download, ids[midpoint:], missing)
        return
    records = list(SeqIO.parse(path, 'embl'))
    observed = {r.id for r in records}
    if len(records) != len(observed) or observed - set(ids):
        raise ValueError(f'Duplicate/unrequested browser records: {path}')
    for record in records:
        yield record, path
    omitted = sorted(set(ids) - observed)
    if omitted and len(ids) == 1:
        missing.append({'accession': ids[0], 'url': url, 'status': 200,
                        'source_file': str(path), 'source_file_sha256': sha256(path),
                        'reason': 'browser_endpoint_omitted_record'})
    elif omitted:
        for accession in omitted:
            yield from fetch_records(download, [accession], missing)


def acquire_secondary(root):
    root = Path(root).resolve()
    biology_path = root / 'configs/secondary/biology.json'
    biology = read_json(biology_path)
    panel = biology['mammals'] + biology['nonmammalian_vertebrates']
    if panel[0] != 9606 or set(map(str, panel)) != set(biology['taxa']):
        raise ValueError('Declared taxonomy panel inconsistent')
    raw = root / 'data/raw/secondary'
    report_dir = root / 'reports/secondary'
    download = Download(raw)
    dependencies = {}
    for name in ['impi-2021-q4pre.xlsx', 'human-small-cds.tsv', 'human-small-cds-count.json']:
        path = root / 'data/raw' / name
        dependencies[name] = verify_receipt(path)
    impi = resolve_impi(root / 'data/raw/impi-2021-q4pre.xlsx', biology)
    query = '(' + ' OR '.join(f'organism_id:{tax_id}' for tax_id in panel) + ')'
    query += ' AND (length:[1 TO 100]) AND (reviewed:true) AND (protein_name:ATP synthase)'
    references_path = download.get('reviewed-panel-references.json', UNIPROT_URL,
                                   {'query': query, 'format': 'json'})
    refs, excluded_refs = reference_inventory(read_json(references_path)['results'], biology, impi)
    write_canonical_json(raw / 'references.json', refs)
    targets = defaultdict(list)
    metadata_stats, metadata_exclusions = [], Counter()
    for tax_id in panel:
        species_refs = [ref for ref in refs if ref['tax_id'] == tax_id]
        if not species_refs:
            continue
        if tax_id == 9606:
            path = root / 'data/raw/human-small-cds.tsv'
            count = int(read_json(root / 'data/raw/human-small-cds-count.json')['count'])
            query_text = 'tax_eq(9606) AND base_count<=303'
        else:
            query_text = index_query(tax_id, species_refs)
            params = {'result': 'coding', 'query': query_text, 'format': 'tsv', 'fields': FIELDS, 'limit': '0'}
            path = download.get(f'index/{tax_id}.tsv', ENA_SEARCH, params)
            count_path = download.get(f'index/{tax_id}-count.json',
                                     'https://www.ebi.ac.uk/ena/portal/api/count',
                                     {'result': 'coding', 'query': query_text, 'format': 'json'})
            count = int(read_json(count_path)['count'])
        index_sha256 = sha256(path)
        rows_seen, matched = 0, 0
        groups = defaultdict(list)
        with path.open(encoding='utf-8') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                rows_seen += 1
                ref = match_metadata(row, species_refs)
                if ref is None:
                    continue
                matched += 1
                reason = metadata_reason(row, tax_id)
                if reason:
                    metadata_exclusions[f'{tax_id}/{ref["family"]}/{reason}'] += 1
                    continue
                groups[(ref['uniprot'], row['sequence_md5'] or row['protein_id'])].append(row)
        if rows_seen != count:
            raise ValueError(f'Index count mismatch for {tax_id}: {rows_seen} != {count}')
        refs_by_id = {ref['uniprot']: ref for ref in species_refs}
        for (uniprot, _), rows in sorted(groups.items()):
            rows.sort(key=lambda r: r['protein_id'])
            representative = rows[0]
            targets[representative['protein_id']].append({
                'reference': refs_by_id[uniprot], 'metadata': representative, 'aliases': rows,
                'route': 'ena_index', 'index_file': path.relative_to(root).as_posix(),
                'index_file_sha256': index_sha256})
        metadata_stats.append({'tax_id': tax_id, 'query': query_text, 'rows': rows_seen,
                               'count_verified': True, 'matched_rows': matched,
                               'sequence_hash_groups': len(groups)})
        print(f'index {tax_id}: {matched} matched accessions / {len(groups)} CDS-hash groups', flush=True)
    reference_source_sha256 = sha256(references_path)
    for ref in refs:
        for xref in ref['embl']:
            targets[xref['protein_id']].append({'reference': ref, 'metadata': None,
                'aliases': [xref], 'route': 'uniprot_embl_cross_reference',
                'index_file': references_path.relative_to(root).as_posix(),
                'index_file_sha256': reference_source_sha256})
    target_list = sorted(targets)
    write_canonical_json(raw / 'targets.json', targets)
    rejected, missing, accepted = [], [], []
    for start in range(0, len(target_list), 50):
        for record, path in fetch_records(download, target_list[start:start + 50], missing):
            for context in targets[record.id]:
                ref = context['reference']
                result, reason = validate_secondary_record(record, ref, context['metadata'],
                    biology['min_reference_identity'], biology['max_amino_acids'])
                if reason:
                    rejected.append({'accession': record.id, 'family': ref['family'],
                        'tax_id': ref['tax_id'], 'route': context['route'], 'reason': reason,
                        'source_file': path.relative_to(root).as_posix(), 'source_file_sha256': sha256(path),
                        'metadata_aliases': len(context['aliases'])})
                    continue
                result['provenance'] = []
                parents = sorted({part.ref for feature in record.features if feature.type == 'CDS'
                                  for part in feature.location.parts if part.ref})
                for alias in context['aliases']:
                    result['provenance'].append({'accession': alias['protein_id'],
                        'retrieved_accession': record.id,
                        'parent_accession': alias['parent_accession'] if alias['protein_id'] != record.id
                            else ','.join(parents) or alias['parent_accession'],
                        'index_parent_accession': alias['parent_accession'],
                        'tax_id': ref['tax_id'], 'reference_uniprot': ref['uniprot'],
                        'reference_identity': result['reference_identity'], 'route': context['route'],
                        'source_file': path.relative_to(root).as_posix(), 'source_file_sha256': sha256(path),
                        'index_file': context['index_file'], 'index_file_sha256': context['index_file_sha256'],
                        'sequence_equivalent_metadata_alias': alias['protein_id'] != record.id,
                        'sequence_md5': alias.get('sequence_md5'),
                        'study_accession': alias.get('study_accession', ''),
                        'cross_reference_status': alias.get('cross_reference_status', '')})
                accepted.append(result)
        print(f'CDS {min(start + 50, len(target_list))}/{len(target_list)}; accepted routes {len(accepted)}', flush=True)
    records = merge_records(accepted, panel)
    output = raw / 'eligible-records.jsonl'
    payload = ''.join(json.dumps(r, sort_keys=True) + '\n' for r in records).encode('utf-8')
    write_canonical_json(raw / 'rejections.json', rejected)
    write_canonical_json(raw / 'unavailable-accessions.json', missing)
    manifest = {'created_utc': now(), 'stage': 'training_only', 'biology_sha256': sha256(biology_path),
        'references_file': (raw / 'references.json').relative_to(root).as_posix(),
        'references_sha256': sha256(raw / 'references.json'), 'references': len(refs),
        'impi_resolved_rows': impi, 'reference_exclusions': excluded_refs,
        'missing_reference_combinations': [{'tax_id': t, 'family': f} for t in panel for f in biology['families']
            if not any(r['tax_id'] == t and r['family'] == f for r in refs)],
        'metadata_inventory': metadata_stats, 'metadata_exclusions': dict(metadata_exclusions),
        'requested_distinct_accessions': len(targets), 'unavailable_accessions': missing,
        'record_rejections': dict(Counter(r['reason'] for r in rejected)),
        'eligible_records_file': output.relative_to(root).as_posix(), 'eligible_records_sha256': hashlib.sha256(payload).hexdigest(),
        'globally_unique_cds': len(records), 'unique_proteins': len({r['protein'] for r in records}),
        'by_primary_taxon_family': dict(Counter(f'{r["tax_id"]}/{r["family"]}' for r in records)),
        'by_observed_taxon_family': dict(Counter(f'{t}/{r["family"]}' for r in records for t in r['taxa'])),
        'source_dependencies': dependencies,
        'source_registry': [read_json(p) for p in sorted(raw.rglob('*.source.json'))],
        'deduplication': 'Exact RNA globally; human first then declared panel order for primary stratum; all accession/taxon provenance retained.',
        'representative_policy': 'Index metadata passes flags before one deterministic accession per nucleotide MD5; aliases are metadata-supported sequence equivalents, not independently CDS-validated records. All UniProt EMBL protein links are fetched separately.',
        'scope_limit': 'All eligible CDS recovered by the frozen panel/reference query and declared retrieval routes, not a globally exhaustive database census.'}
    report_path = report_dir / 'acquisition.json'
    identity_keys = ['biology_sha256', 'references_sha256', 'eligible_records_sha256', 'globally_unique_cds',
                     'by_primary_taxon_family', 'metadata_inventory', 'metadata_exclusions', 'record_rejections']
    preserve_identical_report(report_path, manifest, identity_keys, check_only=True)
    write_frozen_bytes(root, output, payload, report_path)
    preserve_identical_report(report_path, manifest, identity_keys)
    print(json.dumps({'globally_unique_cds': len(records), 'by_primary_taxon_family': manifest['by_primary_taxon_family'],
                      'unavailable_accessions': len(missing), 'record_rejections': manifest['record_rejections']}, indent=2), flush=True)
    return records


if __name__ == '__main__':
    acquire_secondary(Path.cwd())
