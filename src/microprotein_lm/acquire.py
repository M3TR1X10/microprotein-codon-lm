"""IMPI -> reviewed human microproteins -> ENA CDS metadata and records.

No reverse translation, inferred variants, cross-organism supplementation, or
synthetic sequences are used in acquisition. Cached requests are checksummed.
"""
import csv
import io
import json
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import openpyxl
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from Bio import SeqIO

from .io import now, read_json, sha256, write_json

IMPI_PAGE = 'https://www.mrc-mbu.cam.ac.uk/research-resources-and-facilities/impi'
IMPI_URL = 'https://www.mrc-mbu.cam.ac.uk/files/impi-2021-q4pre-20211001-dist_0.xlsx'
UNIPROT_URL = 'https://rest.uniprot.org/uniprotkb/stream'
ENA_SEARCH = 'https://www.ebi.ac.uk/ena/portal/api/search'
FIELDS = ('accession,protein_id,parent_accession,gene,gene_synonym,product,base_count,'
          'sequence_md5,tax_id,transl_table,codon_start,pseudo,pseudo_gene,exception,'
          'transl_except,artificial_location,ribosomal_slippage,study_accession')


class Download:
    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.session.headers['User-Agent'] = 'microprotein-codon-lm/0.1 scientific-data-acquisition'
        retry = Retry(total=4, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
        self.session.mount('https://', HTTPAdapter(max_retries=retry))

    def get(self, filename, url, params=None):
        path = self.root / filename
        manifest = path.with_name(path.name + '.source.json')
        request_url = requests.Request('GET', url, params=params).prepare().url
        if path.exists() and manifest.exists():
            source = read_json(manifest)
            if source['url'] != request_url or source['sha256'] != sha256(path):
                raise ValueError(f'Cache provenance mismatch: {path}; use a new raw directory')
            return path
        path.parent.mkdir(parents=True, exist_ok=True)
        with self.session.get(url, params=params, timeout=(20, 180), stream=True) as response:
            response.raise_for_status()
            temp = path.with_suffix(path.suffix + '.partial')
            with temp.open('wb') as f:
                for chunk in response.iter_content(1024 * 1024):
                    f.write(chunk)
            temp.replace(path)
            write_json(manifest, {'url': request_url, 'retrieved_utc': now(),
                'sha256': sha256(path), 'bytes': path.stat().st_size,
                'etag': response.headers.get('ETag'),
                'uniprot_release': response.headers.get('X-UniProt-Release')})
        time.sleep(0.12)
        return path


def fetch_insdseq_xml_with_failover(accession_ids: str, session: requests.Session = None) -> Tuple[bytes, str]:
    """Fetches INSDSeq XML from ENA XML API with failover to NCBI Entrez efetch.

    Returns:
        Tuple of (raw_xml_bytes, source_url)
    """
    if session is None:
        session = requests.Session()
        session.headers['User-Agent'] = 'microprotein-codon-lm/0.1 scientific-data-acquisition'

    # Primary: ENA XML API
    ena_url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession_ids}"
    try:
        resp = session.get(ena_url, timeout=(15, 60))
        if resp.status_code == 200 and (b'<INSD' in resp.content or b'<?xml' in resp.content):
            return resp.content, ena_url
    except Exception:
        pass

    # Failover 1: NCBI Entrez nuccore
    ncbi_nuccore_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession_ids}&retmode=xml"
    try:
        resp = session.get(ncbi_nuccore_url, timeout=(15, 60))
        if resp.status_code == 200 and (b'<INSD' in resp.content or b'<?xml' in resp.content):
            return resp.content, ncbi_nuccore_url
    except Exception:
        pass

    # Failover 2: NCBI Entrez protein
    ncbi_protein_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={accession_ids}&retmode=xml"
    try:
        resp = session.get(ncbi_protein_url, timeout=(15, 60))
        if resp.status_code == 200 and (b'<INSD' in resp.content or b'<?xml' in resp.content):
            return resp.content, ncbi_protein_url
    except Exception:
        pass

    raise RuntimeError(f"Failed to fetch INSDSeq XML for {accession_ids} from ENA and NCBI Entrez endpoints")



def discover(raw='data/raw'):
    download = Download(raw)
    workbook_path = download.get('impi-2021-q4pre.xlsx', IMPI_URL)
    uniprot_path = download.get('human-reviewed-small.json', UNIPROT_URL,
        {'query': '(organism_id:9606) AND (length:[1 TO 100]) AND (reviewed:true)', 'format': 'json'})
    workbook = openpyxl.load_workbook(workbook_path, read_only=True, data_only=True)
    rows = iter(workbook['IMPI-2021Q4pre'].values)
    header = next(rows)
    impi = {r[1]: dict(zip(header, r)) for r in rows}
    candidates = []
    for entry in read_json(uniprot_path)['results']:
        if entry.get('proteinDescription', {}).get('flag') == 'Fragment':
            continue
        for gene in entry.get('genes', []):
            name = gene.get('geneName', {}).get('value')
            row = impi.get(name)
            if not row or not str(row['IMPI Class']).startswith('Verified mitochondrial'):
                continue
            candidates.append({'gene': name, 'uniprot': entry['primaryAccession'],
                'protein': entry['proteinDescription']['recommendedName']['fullName']['value'],
                'product_names': sorted({row['Name'], entry['proteinDescription']['recommendedName']['fullName']['value'],
                    *[v['fullName']['value'] for v in entry['proteinDescription'].get('alternativeNames', [])]}),
                'reference_protein': entry['sequence']['value'],
                'amino_acids': entry['sequence']['length'],
                'synonyms': sorted({v['value'] for g in entry.get('genes', [])
                    for k in ['synonyms', 'orderedLocusNames', 'orfNames'] for v in g.get(k, [])}),
                'ensembl': row['Ensembl'], 'impi_class': row['IMPI Class']})
    candidates.sort(key=lambda c: (c['gene'], c['uniprot']))
    write_json(Path(raw) / 'candidates.json', candidates)
    # Retrieve all short human CDS metadata, so one query covers every candidate
    # fairly and is not truncated by a different per-family download budget.
    metadata_path = download.get('human-small-cds.tsv', ENA_SEARCH, {
        'result': 'coding', 'query': 'tax_eq(9606) AND base_count<=303',
        'fields': FIELDS, 'format': 'tsv', 'limit': '0'})
    aliases = defaultdict(set)
    products = defaultdict(set)
    product_aliases = read_json('configs/product_aliases.json')
    for i, candidate in enumerate(candidates):
        for alias in [candidate['gene'], *candidate['synonyms']]:
            aliases[normalize_name(alias)].add(i)
        for product in candidate['product_names'] + product_aliases.get(candidate['gene'], []):
            products[normalize_name(product)].add(i)
    groups = [dict() for _ in candidates]
    with metadata_path.open(encoding='utf-8') as f:
        total = 0
        matched = 0
        for row in csv.DictReader(f, delimiter='\t'):
            total += 1
            indices = set()
            for name in re.split('[;,]', row['gene'] + ';' + row['gene_synonym']):
                indices.update(aliases.get(normalize_name(name), []))
            indices.update(products.get(normalize_name(row['product']), []))
            indices.update(aliases.get(normalize_name(row['product']), []))
            if indices:
                matched += 1
            for idx in indices:
                groups[idx][row['accession']] = row
    count_path = download.get('human-small-cds-count.json',
        'https://www.ebi.ac.uk/ena/portal/api/count', {
            'result':'coding', 'query':'tax_eq(9606) AND base_count<=303', 'format':'json'})
    expected_count = int(read_json(count_path)['count'])
    if expected_count != total:
        raise ValueError(f'ENA inventory count mismatch: {total} downloaded, {expected_count} expected')
    write_json('reports/inventory-coverage.json', {'total_short_human_cds':total,
        'matched_to_candidates':matched, 'unmapped_records':total-matched,
        'candidate_count':len(candidates), 'count_verified':True,
        'scope':'Exact normalized gene/alias or curated product-name match; unmapped records are not eligible.'})
    counts = []
    for candidate, records in zip(candidates, groups):
        name = candidate['gene'] + '-' + candidate['uniprot']
        write_json(Path(raw) / 'metadata' / (name + '.json'), list(records.values()))
        counts.append({'gene': candidate['gene'], 'uniprot': candidate['uniprot'],
            'raw_records': len(records),
            'unique_cds_md5_upper_bound': len({r['sequence_md5'] or r['accession'] for r in records.values()})})
    write_json('reports/candidate-inventory.json', counts)
    print(json.dumps(sorted([c for c in counts if c['raw_records']],
        key=lambda x: -x['unique_cds_md5_upper_bound']), indent=2), flush=True)
    return candidates


def fetch_group(candidate, raw='data/raw'):
    """One representative per identical CDS hash; retain all aliases separately.

    Quality flags are screened on metadata before grouping. Raw counts never
    constitute a claim of usable training diversity.
    """
    name = candidate['gene'] + '-' + candidate['uniprot']
    records = read_json(Path(raw) / 'metadata' / (name + '.json'))
    groups = defaultdict(list)
    excluded = defaultdict(int)
    for r in records:
        reason = metadata_rejection(r)
        if reason:
            excluded[reason] += 1
            continue
        groups[r['sequence_md5'] or r['accession']].append(r)
    download = Download(raw)
    outputs = []
    # Per-hash representatives are deterministic; all source accessions remain
    # recorded, but are not interpreted as independent biological samples.
    representatives = [(sorted(rows, key=lambda r:r['protein_id'])[0], rows)
                       for _, rows in sorted(groups.items())]
    for start in range(0, len(representatives), 50):
        batch = representatives[start:start+50]
        ids = ','.join(row['protein_id'] for row, _ in batch)
        path = download.get(f'cds/{name}-{start:05d}.embl',
            'https://www.ebi.ac.uk/ena/browser/api/embl/' + ids)
        observed = {rec.id for rec in SeqIO.parse(path, 'embl')}
        if observed != {row['protein_id'] for row, _ in batch}:
            raise ValueError(f'ENA batch omitted or added records: {path}')
        for row, rows in batch:
            outputs.append({'path':str(path), 'metadata':row,
                'source_accessions':sorted(r['protein_id'] for r in rows)})
        print(f'{name}: retrieved {min(start+50,len(representatives))}/{len(representatives)} unique CDS records', flush=True)
    write_json(Path(raw) / 'metadata' / (name + '-fetched.json'),
        {'records': outputs, 'metadata_exclusions': dict(excluded)})
    return outputs, excluded


def metadata_rejection(row):
    if str(row['tax_id']) != '9606':
        return 'wrong_taxon'
    if row.get('codon_start') not in ('', '1'):
        return 'partial_frame'
    for key in ['pseudo', 'pseudo_gene', 'exception', 'transl_except',
                'artificial_location', 'ribosomal_slippage']:
        if row.get(key, '').lower() not in ('', 'false', '0'):
            return key
    return None


def normalize_name(value):
    return re.sub('[^a-z0-9]', '', value.lower())
