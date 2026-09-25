"""Fail-closed CDS checks on coding-oriented ENA records, not whole genomes."""
import hashlib
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.SeqFeature import ExactPosition

from .acquire import fetch_group
from .insdseq_parser import INSDSeqParser
from .io import now, read_json, sha256, write_json
from .translation_engine import TranslationEngine


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
        for ref in self.references:
            alignments = self.aligner.align(ref, protein)
            if alignments:
                aln = alignments[0]
                if aln.score > best_score:
                    best_score = aln.score
                    best_alignment = aln
        
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

    return {
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
    }, None


def select(raw='data/raw', processed='data/processed', min_identity=0.95):
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
        row = {**count, 'usable_unique_cds': len(accepted),
            'unique_proteins': len({r['protein'] for r in accepted.values()}),
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
    manifest = {'created_utc': now(), 'stage': 'training_only', 'organism': 'Homo sapiens', 'tax_id': 9606,
        'candidate': candidates[winner['gene'], winner['uniprot']],
        'selection_rule': 'Maximum unique CDS passing fixed QC; alphabetical gene then UniProt accession breaks ties',
        'min_protein_identity': min_identity, 'function_evidence': 'Curated family annotation and reference sequence similarity; variant function is not experimentally confirmed',
        'cohort_sha256': sha256(out), 'sequences': len(records),
        'unique_proteins': len({r['protein'] for r in records}),
        'training_sequences': len(records), 'validation_sequences': 0, 'test_sequences': 0,
        'deduplication': 'Exact RNA sequence; original accession multiplicity retained, never weighted',
        'source_sha256': {name: sha256(Path(raw) / name) for name in
            ['impi-2021-q4pre.xlsx', 'human-reviewed-small.json', 'human-small-cds.tsv', 'candidates.json']},
        'scope_limit': 'Winner among the explicitly mapped IMPI-verified, reviewed human <=100-aa candidates in this ENA snapshot, not all possible database sequences'}
    write_json(root / 'manifest.json', manifest)
    write_json('reports/data-manifest.json', manifest)
    return manifest
