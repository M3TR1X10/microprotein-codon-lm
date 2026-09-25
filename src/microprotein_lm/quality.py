"""Fail-closed CDS checks on coding-oriented ENA records, not whole genomes."""
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SeqFeature import ExactPosition

from .insdseq_parser import INSDSeqParser
from .io import now, read_json, sha256, write_json
from .translation_engine import TranslationEngine
from .translation_evidence import calculate_tai, calculate_cai, determine_evidence_tier, annotate_record_evidence


class TaxonomyResolver:
    """Dynamic NCBI Taxonomy / UniProt API resolver for lineage trees and ranks."""

    FALLBACK_TAXONOMY = {
        '9606': {
            'tax_id': 9606,
            'scientific_name': 'Homo sapiens',
            'common_name': 'human',
            'rank': 'species',
            'lineage_tax_ids': [2759, 33154, 33208, 7711, 40674, 9443, 9604, 9605, 9606],
            'lineage': [
                {'tax_id': 2759, 'scientific_name': 'Eukaryota', 'rank': 'superkingdom'},
                {'tax_id': 33154, 'scientific_name': 'Opisthokonta', 'rank': 'no rank'},
                {'tax_id': 33208, 'scientific_name': 'Metazoa', 'rank': 'kingdom'},
                {'tax_id': 7711, 'scientific_name': 'Chordata', 'rank': 'phylum'},
                {'tax_id': 40674, 'scientific_name': 'Mammalia', 'rank': 'class'},
                {'tax_id': 9443, 'scientific_name': 'Primates', 'rank': 'order'},
                {'tax_id': 9604, 'scientific_name': 'Hominidae', 'rank': 'family'},
                {'tax_id': 9605, 'scientific_name': 'Homo', 'rank': 'genus'},
                {'tax_id': 9606, 'scientific_name': 'Homo sapiens', 'rank': 'species'}
            ],
            'is_eukaryota': True,
            'is_human': True,
            'source': 'fallback'
        },
        '10090': {
            'tax_id': 10090,
            'scientific_name': 'Mus musculus',
            'common_name': 'house mouse',
            'rank': 'species',
            'lineage_tax_ids': [2759, 33208, 7711, 40674, 9989, 10088, 10090],
            'lineage': [
                {'tax_id': 2759, 'scientific_name': 'Eukaryota', 'rank': 'superkingdom'},
                {'tax_id': 33208, 'scientific_name': 'Metazoa', 'rank': 'kingdom'},
                {'tax_id': 7711, 'scientific_name': 'Chordata', 'rank': 'phylum'},
                {'tax_id': 40674, 'scientific_name': 'Mammalia', 'rank': 'class'},
                {'tax_id': 9989, 'scientific_name': 'Rodentia', 'rank': 'order'},
                {'tax_id': 10088, 'scientific_name': 'Muridae', 'rank': 'family'},
                {'tax_id': 10090, 'scientific_name': 'Mus musculus', 'rank': 'species'}
            ],
            'is_eukaryota': True,
            'is_human': False,
            'source': 'fallback'
        },
        '10116': {
            'tax_id': 10116,
            'scientific_name': 'Rattus norvegicus',
            'common_name': 'Norway rat',
            'rank': 'species',
            'lineage_tax_ids': [2759, 33208, 7711, 40674, 9989, 10088, 10116],
            'lineage': [
                {'tax_id': 2759, 'scientific_name': 'Eukaryota', 'rank': 'superkingdom'},
                {'tax_id': 33208, 'scientific_name': 'Metazoa', 'rank': 'kingdom'},
                {'tax_id': 7711, 'scientific_name': 'Chordata', 'rank': 'phylum'},
                {'tax_id': 40674, 'scientific_name': 'Mammalia', 'rank': 'class'},
                {'tax_id': 9989, 'scientific_name': 'Rodentia', 'rank': 'order'},
                {'tax_id': 10088, 'scientific_name': 'Muridae', 'rank': 'family'},
                {'tax_id': 10116, 'scientific_name': 'Rattus norvegicus', 'rank': 'species'}
            ],
            'is_eukaryota': True,
            'is_human': False,
            'source': 'fallback'
        }
    }

    def __init__(self, cache_dir=None, session=None):
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.session = session or requests.Session()
        self.session.headers.setdefault('User-Agent', 'microprotein-codon-lm/0.1 taxonomy-resolver')
        retry = Retry(total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504])
        self.session.mount('https://', HTTPAdapter(max_retries=retry))
        self._memory_cache = {}

    def get_lineage(self, tax_id: Union[int, str]) -> Dict:
        tax_str = str(tax_id).strip()
        if tax_str in self._memory_cache:
            return self._memory_cache[tax_str]

        if self.cache_dir:
            cache_path = self.cache_dir / f"tax_{tax_str}.json"
            if cache_path.exists():
                try:
                    data = read_json(cache_path)
                    self._memory_cache[tax_str] = data
                    return data
                except Exception:
                    pass

        result = self._query_uniprot_taxonomy(tax_str)
        if not result:
            result = self._query_ncbi_taxonomy(tax_str)
        if not result:
            result = self._fallback_taxonomy(tax_str)

        if self.cache_dir and result:
            try:
                write_json(self.cache_dir / f"tax_{tax_str}.json", result)
            except Exception:
                pass

        self._memory_cache[tax_str] = result
        return result

    def _query_uniprot_taxonomy(self, tax_str: str) -> Optional[Dict]:
        url = f"https://rest.uniprot.org/taxonomy/{tax_str}"
        try:
            resp = self.session.get(url, timeout=(5, 10))
            if resp.status_code == 200:
                data = resp.json()
                lineage_items = []
                lineage_ids = []
                for node in data.get('lineage', []):
                    t_id = node.get('taxonId')
                    lineage_ids.append(t_id)
                    lineage_items.append({
                        'tax_id': t_id,
                        'scientific_name': node.get('scientificName', ''),
                        'rank': node.get('rank', 'no rank')
                    })
                curr_id = data.get('taxonId', int(tax_str) if tax_str.isdigit() else tax_str)
                lineage_ids.append(curr_id)
                lineage_items.append({
                    'tax_id': curr_id,
                    'scientific_name': data.get('scientificName', ''),
                    'rank': data.get('rank', 'species')
                })
                is_euk = any(node.get('scientific_name') == 'Eukaryota' or node.get('tax_id') == 2759 for node in lineage_items)
                return {
                    'tax_id': curr_id,
                    'scientific_name': data.get('scientificName', ''),
                    'common_name': data.get('commonName', ''),
                    'rank': data.get('rank', 'species'),
                    'lineage_tax_ids': lineage_ids,
                    'lineage': lineage_items,
                    'is_eukaryota': is_euk,
                    'is_human': str(curr_id) == '9606',
                    'source': 'uniprot_api'
                }
        except Exception:
            pass
        return None

    def _query_ncbi_taxonomy(self, tax_str: str) -> Optional[Dict]:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id={tax_str}&retmode=json"
        try:
            resp = self.session.get(url, timeout=(5, 10))
            if resp.status_code == 200:
                data = resp.json()
                res = data.get('result', {}).get(tax_str, {})
                if res and 'scientificname' in res:
                    s_name = res.get('scientificname', '')
                    c_name = res.get('commonname', '')
                    rank = res.get('rank', 'species')
                    return {
                        'tax_id': int(tax_str) if tax_str.isdigit() else tax_str,
                        'scientific_name': s_name,
                        'common_name': c_name,
                        'rank': rank,
                        'lineage_tax_ids': [int(tax_str)] if tax_str.isdigit() else [tax_str],
                        'lineage': [{'tax_id': tax_str, 'scientific_name': s_name, 'rank': rank}],
                        'is_eukaryota': True,
                        'is_human': str(tax_str) == '9606',
                        'source': 'ncbi_api'
                    }
        except Exception:
            pass
        return None

    def _fallback_taxonomy(self, tax_str: str) -> Dict:
        if tax_str in self.FALLBACK_TAXONOMY:
            return self.FALLBACK_TAXONOMY[tax_str]
        t_id = int(tax_str) if tax_str.isdigit() else tax_str
        return {
            'tax_id': t_id,
            'scientific_name': f'Organism_{tax_str}',
            'common_name': '',
            'rank': 'species',
            'lineage_tax_ids': [t_id],
            'lineage': [{'tax_id': t_id, 'scientific_name': f'Organism_{tax_str}', 'rank': 'species'}],
            'is_eukaryota': True,
            'is_human': str(tax_str) == '9606',
            'source': 'fallback'
        }


class CellularCompartmentResolver:
    """Dynamic EBI OLS & GO API resolver for cellular compartment annotations (mitochondrion, nucleus, cytosol)."""

    PRIMARY_GO_MAP = {
        'GO:0005739': 'mitochondrion',
        'GO:0005634': 'nucleus',
        'GO:0005829': 'cytosol',
        'GO:0005783': 'endoplasmic_reticulum',
        'GO:0005886': 'plasma_membrane',
        'GO:0005794': 'golgi_apparatus',
        'GO:0005576': 'extracellular',
    }

    KEYWORD_MAP = {
        'mitochondr': ('mitochondrion', 'GO:0005739'),
        'mitochondria': ('mitochondrion', 'GO:0005739'),
        'matrix': ('mitochondrion', 'GO:0005739'),
        'cristae': ('mitochondrion', 'GO:0005739'),
        'impi': ('mitochondrion', 'GO:0005739'),
        'nucle': ('nucleus', 'GO:0005634'),
        'chromatin': ('nucleus', 'GO:0005634'),
        'nucleolus': ('nucleus', 'GO:0005634'),
        'cytosol': ('cytosol', 'GO:0005829'),
        'cytoplas': ('cytosol', 'GO:0005829'),
        'endoplasmic': ('endoplasmic_reticulum', 'GO:0005783'),
        'sarcoplasmic': ('endoplasmic_reticulum', 'GO:0005783'),
        'plasma membrane': ('plasma_membrane', 'GO:0005886'),
        'cell membrane': ('plasma_membrane', 'GO:0005886'),
        'golgi': ('golgi_apparatus', 'GO:0005794'),
        'extracellular': ('extracellular', 'GO:0005576'),
        'secreted': ('extracellular', 'GO:0005576'),
    }

    def __init__(self, cache_dir=None, session=None):
        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.session = session or requests.Session()
        self.session.headers.setdefault('User-Agent', 'microprotein-codon-lm/0.1 compartment-resolver')
        retry = Retry(total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504])
        self.session.mount('https://', HTTPAdapter(max_retries=retry))
        self._memory_cache = {}

    def resolve_compartment(self, query: str) -> Dict:
        q_norm = str(query).strip().lower()
        if q_norm in self._memory_cache:
            return self._memory_cache[q_norm]

        if self.cache_dir:
            h = hashlib.md5(q_norm.encode()).hexdigest()
            cache_path = self.cache_dir / f"cmp_{h}.json"
            if cache_path.exists():
                try:
                    data = read_json(cache_path)
                    self._memory_cache[q_norm] = data
                    return data
                except Exception:
                    pass

        result = self._query_ebi_ols(q_norm)
        if not result:
            result = self._fallback_compartment(q_norm)

        if self.cache_dir and result:
            try:
                h = hashlib.md5(q_norm.encode()).hexdigest()
                write_json(self.cache_dir / f"cmp_{h}.json", result)
            except Exception:
                pass

        self._memory_cache[q_norm] = result
        return result

    def _query_ebi_ols(self, query: str) -> Optional[Dict]:
        url = f"https://www.ebi.ac.uk/ols4/api/search?q={requests.utils.quote(query)}&ontology=go"
        try:
            resp = self.session.get(url, timeout=(5, 10))
            if resp.status_code == 200:
                data = resp.json()
                docs = data.get('response', {}).get('docs', [])
                if docs:
                    top = docs[0]
                    obo_id = top.get('obo_id', '')
                    label = top.get('label', query)
                    primary, default_go = self._map_keyword(query)
                    go_id = obo_id if obo_id.startswith('GO:') else default_go
                    if obo_id in self.PRIMARY_GO_MAP:
                        primary = self.PRIMARY_GO_MAP[obo_id]
                    return {
                        'query': query,
                        'primary_compartment': primary,
                        'go_id': go_id,
                        'label': label,
                        'is_organelle': primary in ('mitochondrion', 'nucleus', 'cytosol', 'endoplasmic_reticulum', 'golgi_apparatus'),
                        'source': 'ebi_ols_api'
                    }
        except Exception:
            pass
        return None

    def _fallback_compartment(self, query: str) -> Dict:
        primary, go_id = self._map_keyword(query)
        return {
            'query': query,
            'primary_compartment': primary,
            'go_id': go_id,
            'label': primary.replace('_', ' ').title(),
            'is_organelle': primary in ('mitochondrion', 'nucleus', 'cytosol', 'endoplasmic_reticulum', 'golgi_apparatus'),
            'source': 'fallback'
        }

    def _map_keyword(self, query: str) -> Tuple[str, str]:
        q_lower = query.lower()
        for kw, (primary, go_id) in self.KEYWORD_MAP.items():
            if kw in q_lower:
                return primary, go_id
        return 'other', 'GO:0005575'


class ThreeTierQualityCategorizer:
    """Three-Tier Quality Categorization:
    - Tier 1 (Gold/Swiss-Prot): Curated Swiss-Prot / IMPI verified, high identity (>=0.95), canonical AUG.
    - Tier 2 (Genomic CDS): Standard GenBank/ENA annotated CDS feature, complete reading frame.
    - Tier 3 (sORF leads): Uncharacterized sORF prediction, non-canonical initiation (CTG, GTG, TTG), or novel lead.
    """

    TIER_1_GOLD = 1
    TIER_2_GENOMIC = 2
    TIER_3_SORF = 3

    TIER_NAMES = {
        1: "Tier 1 (Gold/Swiss-Prot)",
        2: "Tier 2 (Genomic CDS)",
        3: "Tier 3 (sORF leads)"
    }

    def categorize(self, record_data: Dict, metadata: Optional[Dict] = None, reference: Optional[Dict] = None) -> Dict:
        metadata = metadata or {}
        reference = reference or {}

        rationale = []

        is_swiss_prot = bool(
            reference.get('reviewed') or
            str(reference.get('impi_class', '')).startswith('Verified') or
            record_data.get('is_swiss_prot') or
            reference.get('uniprot', '').startswith(('P', 'Q', 'O'))
        )
        identity = record_data.get('reference_identity', 0.0)
        init_type = record_data.get('initiation_type', '')
        init_codon = record_data.get('initiation_codon', '')
        protein_len = len(record_data.get('protein', ''))

        # Check for Tier 1
        if is_swiss_prot and identity >= 0.95 and init_type == 'CANONICAL':
            rationale.append("Swiss-Prot/IMPI verified reference")
            rationale.append(f"High reference identity ({identity:.3f})")
            rationale.append(f"Canonical initiation codon ({init_codon})")
            rationale.append("Complete verified CDS")
            return {
                'tier': 1,
                'tier_name': self.TIER_NAMES[1],
                'rationale': rationale
            }

        # Check for Tier 3 (sORF lead / non-canonical TIS / uncharacterized sORF lead)
        if init_type == 'NEAR_COGNATE' or init_codon in ('CUG', 'GUG', 'UUG', 'ACG') or record_data.get('is_sorf_lead'):
            rationale.append(f"Non-canonical initiation codon ({init_codon})")
            rationale.append(f"sORF candidate microprotein ({protein_len} aa)")
            if not is_swiss_prot:
                rationale.append("Uncharacterized sORF prediction lead")
            return {
                'tier': 3,
                'tier_name': self.TIER_NAMES[3],
                'rationale': rationale
            }

        # Default to Tier 2 (Genomic CDS)
        rationale.append("Standard annotated genomic CDS feature")
        rationale.append(f"Valid reading frame ({protein_len} aa)")
        rationale.append(f"Reference identity ({identity:.3f})")
        if is_swiss_prot:
            rationale.append("Curated gene family candidate")

        tier_num = 1 if (is_swiss_prot and identity >= 0.95) else 2
        return {
            'tier': tier_num,
            'tier_name': self.TIER_NAMES[tier_num],
            'rationale': rationale
        }


def categorize_quality_tier(record_data: Dict, metadata: Optional[Dict] = None, reference: Optional[Dict] = None) -> Dict:
    categorizer = ThreeTierQualityCategorizer()
    return categorizer.categorize(record_data, metadata, reference)


class DynamicIdentityGate:
    def __init__(self, references: list, min_identity=0.95):
        self.references = references
        self.min_identity = min_identity
        self.aligner = PairwiseAligner()
        try:
            self.aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        except Exception:
            pass
        self.aligner.mode = 'global'
        try:
            self.aligner.end_gap_score = 0.0
        except AttributeError:
            pass

    def evaluate(self, protein):
        best_score = -1
        best_alignment = None
        best_identity = -1.0
        for ref in self.references:
            alignments = self.aligner.align(ref, protein)
            if alignments:
                aln = alignments[0]
                c = aln.counts()
                target_len = len(aln.target)
                identity = c.identities / target_len
                if aln.score > best_score or (aln.score == best_score and identity > best_identity):
                    best_score = aln.score
                    best_alignment = aln
                    best_identity = identity
        
        if not best_alignment:
            return None, 'no_alignment'
        
        c = best_alignment.counts()
        target_len = len(best_alignment.target)
        identity = c.identities / target_len
        positives = getattr(c, 'positives', None)
        if positives is None:
            positives = 0
        positive_rate = positives / target_len
        
        if identity < self.min_identity:
            return None, 'reference_identity_below_threshold'
            
        if c.internal_gaps > 0:
            variant = 'INDEL_VARIANT'
        elif c.gaps > 0:
            terminals = []
            if c.left_insertions > 0: terminals.append('N_TERMINAL_EXTENSION')
            if c.left_deletions > 0: terminals.append('N_TERMINAL_TRUNCATION')
            if c.right_insertions > 0: terminals.append('C_TERMINAL_EXTENSION')
            if c.right_deletions > 0: terminals.append('C_TERMINAL_TRUNCATION')
            variant = '+'.join(terminals) if terminals else 'INDEL_VARIANT'
        elif c.mismatches > 0:
            variant = 'SUBSTITUTION_ONLY'
        else:
            variant = 'EXACT_MATCH'
            
        return {
            'identity': identity,
            'positive_rate': positive_rate,
            'variant_class': variant
        }, None


def validate_record(record, metadata, reference, min_identity=0.95):
    engine = TranslationEngine()
    features = [f for f in record.features if f.type == 'CDS']
    if len(features) != 1:
        return None, 'not_one_cds'
    feature = features[0]
    q = feature.qualifiers
    if record.id != metadata['protein_id'] or q.get('protein_id', [''])[0] != record.id:
        return None, 'accession_mismatch'
    source_taxa = {v for f in record.features if f.type == 'source'
                   for v in f.qualifiers.get('db_xref', []) if v.startswith('taxon:')}
    if source_taxa != {'taxon:9606'}:
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
    # ENA's coding endpoint already splices and orients the nucleotide sequence.
    # CDS locations refer to the PARENT accession; do not extract them a second time.
    dna = str(record.seq).upper()
    if not dna or set(dna) - set('ACGT'):
        return None, 'ambiguous_bases'
    if len(dna) % 3:
        return None, 'incomplete_codon'
    if len(dna) != int(metadata['base_count']):
        return None, 'length_metadata_mismatch'
    md5 = hashlib.md5(dna.encode()).hexdigest()
    if metadata.get('sequence_md5') and metadata['sequence_md5'] != md5:
        return None, 'sequence_hash_mismatch'
    table = int(q.get('transl_table', ['1'])[0])
    expected_table = 2 if reference['gene'].startswith('MT-') else 1
    if table != expected_table:
        return None, 'unexpected_genetic_code'
    try:
        protein = engine.translate_cds(dna, table=table, cds=True)
    except Exception:
        return None, 'invalid_complete_cds'
    if protein != ''.join(q.get('translation', [''])):
        return None, 'translation_mismatch'
    if not 1 <= len(protein) <= 100:
        return None, 'not_microprotein'
    references = []
    if 'reference_protein' in reference:
        references.append(reference['reference_protein'])
    if 'isoforms' in reference:
        if isinstance(reference['isoforms'], dict):
            references.extend(reference['isoforms'].values())
        else:
            references.extend(reference['isoforms'])
    if not references:
        return None, 'no_reference_sequence'

    gate = DynamicIdentityGate(references, min_identity)
    eval_res, reason = gate.evaluate(protein)
    if reason:
        return None, reason
    
    identity = eval_res['identity']

    rna_seq = dna.replace('T', 'U')
    initiation_codon = rna_seq[:3]
    initiation_type = engine.get_initiation_type(initiation_codon, table)
    stop_codon = rna_seq[-3:]
    flank_5 = getattr(record, 'flanking_5prime_utr_30nt', metadata.get('flanking_5prime_utr_30nt', ''))
    flank_3 = getattr(record, 'flanking_3prime_utr_30nt', metadata.get('flanking_3prime_utr_30nt', ''))

    ribo_coverage = getattr(record, 'ribo_coverage', metadata.get('ribo_coverage', None))
    has_ms = metadata.get('has_ms_evidence', False)
    is_conserved = metadata.get('is_conserved', False)

    tai_score = calculate_tai(rna_seq)
    cai_score = calculate_cai(rna_seq)
    tier = determine_evidence_tier(
        has_ms_evidence=has_ms,
        ribo_coverage=ribo_coverage,
        is_conserved=is_conserved
    )

    tax_resolver = TaxonomyResolver()
    comp_resolver = CellularCompartmentResolver()
    categorizer = ThreeTierQualityCategorizer()

    tax_info = tax_resolver.get_lineage(metadata.get('tax_id', 9606))
    comp_query = reference.get('impi_class', metadata.get('product', reference.get('gene', 'mitochondrion')))
    comp_info = comp_resolver.resolve_compartment(comp_query)

    res = {
        'accession': record.id,
        'tax_id': 9606,
        'gene': reference['gene'],
        'uniprot': reference['uniprot'],
        'rna': rna_seq,
        'protein': protein,
        'translation_table': table,
        'reference_identity': identity,
        'positive_rate': eval_res['positive_rate'],
        'variant_class': eval_res['variant_class'],
        'sequence_sha256': hashlib.sha256(rna_seq.encode()).hexdigest(),
        'parent_accession': metadata['parent_accession'],
        'study_accession': metadata.get('study_accession', ''),
        'initiation_codon': initiation_codon,
        'initiation_type': initiation_type,
        'stop_codon': stop_codon,
        'flanking_5prime_utr_30nt': flank_5,
        'flanking_3prime_utr_30nt': flank_3,
        'tai_score': tai_score,
        'cai_score': cai_score,
        'evidence_tier': tier,
        'taxonomy_lineage': tax_info,
        'cellular_compartment': comp_info,
    }
    if ribo_coverage is not None:
        res['ribo_coverage'] = ribo_coverage

    tier_info = categorizer.categorize(res, metadata=metadata, reference=reference)
    res.update({
        'quality_tier': tier_info['tier'],
        'quality_tier_name': tier_info['tier_name'],
        'quality_tier_rationale': tier_info['rationale'],
    })

    return res, None


def select(raw='data/raw', processed='data/processed', min_identity=0.95):
    from .acquire import fetch_group
    if not 0 < min_identity <= 1:
        raise ValueError('min_identity must be in (0,1]')
    inventory = read_json('reports/candidate-inventory.json')
    candidates = {(c['gene'], c['uniprot']): c for c in read_json(Path(raw) / 'candidates.json')}
    rank = []
    cohort = {}
    for count in sorted(inventory, key=lambda c: (-c['unique_cds_md5_upper_bound'], c['gene'], c['uniprot'])):
        ref = candidates[count['gene'], count['uniprot']]
        fetched, excluded = fetch_group(ref, raw)
        cache = {}
        accepted = {}
        reasons = Counter()
        rejection_rows = []
        for item in fetched:
            if item['path'] not in cache:
                cache[item['path']] = {r.id: r for r in SeqIO.parse(item['path'], 'embl')}
            accession = item['metadata']['protein_id']
            result, reason = validate_record(cache[item['path']][accession], item['metadata'], ref, min_identity)
            if reason:
                reasons[reason] += 1
                rejection_rows.append({'accession': accession, 'reason': reason})
            else:
                result['source_accessions'] = item['source_accessions']
                result['source_file_sha256'] = sha256(item['path'])
                accepted[result['sequence_sha256']] = result
        tier_counts = Counter(r.get('quality_tier', 2) for r in accepted.values())
        row = {**count, 'usable_unique_cds': len(accepted),
            'unique_proteins': len({r['protein'] for r in accepted.values()}),
            'tier_distribution': {
                'tier_1_gold': tier_counts[1],
                'tier_2_genomic': tier_counts[2],
                'tier_3_sorf': tier_counts[3],
            },
            'metadata_rejections': dict(excluded), 'sequence_rejections': dict(reasons)}
        rank.append(row)
        cohort[ref['gene'], ref['uniprot']] = list(accepted.values())
        write_json(Path(raw) / 'metadata' / (ref['gene'] + '-' + ref['uniprot'] + '-rejections.json'), rejection_rows)
        if fetched:
            print(f"{ref['gene']}: {len(accepted)} usable distinct CDS; rejections={dict(reasons)}", flush=True)
    rank.sort(key=lambda c: (-c['usable_unique_cds'], c['gene'], c['uniprot']))
    write_json('reports/candidate-ranking.json', rank)
    if not rank or rank[0]['usable_unique_cds'] < 2:
        raise ValueError('Fewer than two distinct usable sequences: training cohort not created')
    winner = rank[0]
    records = sorted(cohort[winner['gene'], winner['uniprot']], key=lambda r: r['sequence_sha256'])
    root = Path(processed)
    root.mkdir(parents=True, exist_ok=True)
    out = root / 'cohort.jsonl'
    out.write_text(''.join(__import__('json').dumps(r, sort_keys=True) + '\n' for r in records), encoding='utf-8')
    tier_summary = Counter(r.get('quality_tier', 2) for r in records)
    manifest = {'created_utc': now(), 'stage': 'training_only', 'organism': 'Homo sapiens', 'tax_id': 9606,
        'candidate': candidates[winner['gene'], winner['uniprot']],
        'selection_rule': 'Maximum unique CDS passing fixed QC; alphabetical gene then UniProt accession breaks ties',
        'min_protein_identity': min_identity, 'function_evidence': 'Curated family annotation and reference sequence similarity; variant function is not experimentally confirmed',
        'cohort_sha256': sha256(out), 'sequences': len(records),
        'unique_proteins': len({r['protein'] for r in records}),
        'tier_distribution': {
            'tier_1_gold': tier_summary[1],
            'tier_2_genomic': tier_summary[2],
            'tier_3_sorf': tier_summary[3]
        },
        'taxonomy': records[0].get('taxonomy_lineage', {}) if records else {},
        'cellular_compartment': records[0].get('cellular_compartment', {}) if records else {},
        'training_sequences': len(records), 'validation_sequences': 0, 'test_sequences': 0,
        'deduplication': 'Exact RNA sequence; original accession multiplicity retained, never weighted',
        'source_sha256': {name: sha256(Path(raw) / name) for name in
            ['impi-2021-q4pre.xlsx', 'human-reviewed-small.json', 'human-small-cds.tsv', 'candidates.json']},
        'scope_limit': 'Winner among the explicitly mapped IMPI-verified, reviewed human <=100-aa candidates in this ENA snapshot, not all possible database sequences'}
    write_json(root / 'manifest.json', manifest)
    write_json('reports/data-manifest.json', manifest)
    return manifest

