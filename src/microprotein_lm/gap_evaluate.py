"""Case-level evaluation, isolated from training and prospective data selection."""
import hashlib
import json
import math
from pathlib import Path

import torch

from .gap_baselines import GapBaselines, METHODS
from .gap_inference import GapDecoder
from .gap_metrics import reconstruction, TrainingPositionSupport
from .io import read_json, sha256, now
from .model import Decoder, ModelConfig
from .secondary_io import write_json


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def atomic_json(path, value):
    import time
    path = Path(path)
    temporary = path.with_suffix('.partial.json')
    write_json(temporary, value)
    for attempt in range(10):
        try:
            temporary.replace(path)
            return
        except (PermissionError, OSError):
            if attempt == 9:
                raise
            time.sleep(0.1)


def state_hash(model):
    h = hashlib.sha256()
    for key, value in sorted(model.state_dict().items()):
        tensor = value.detach().contiguous()
        h.update(f'{key}|{tensor.dtype}|{tuple(tensor.shape)}\n'.encode())
        h.update(tensor.cpu().numpy().tobytes())
    return h.hexdigest()


def load_model(root, model_spec):
    root = Path(root)
    info = model_spec['checkpoint']
    path = root / info['file']
    if sha256(path) != info['sha256'] or path.stat().st_size != info['bytes']:
        raise ValueError('Audited final checkpoint bytes changed')
    metadata = read_json(root / model_spec['run_metadata_file'])
    checkpoint = torch.load(path, weights_only=True, map_location='cpu')
    if checkpoint['identity'] != model_spec['checkpoint_identity'] or checkpoint['update'] != 2000:
        raise ValueError('Checkpoint does not match frozen evaluation identity')
    cfg = metadata.get('config', metadata)
    from .tokenization import Tokenizer
    tok = Tokenizer(model_spec['mode'])
    maximum = cfg['max_nt_context'] // tok.width - 1
    model_cfg = ModelConfig(vocab_size=len(tok.vocabulary), max_tokens=maximum, token_width=tok.width, **cfg['model'])
    model = Decoder(model_cfg)
    model.load_state_dict(checkpoint['model'], strict=True)
    if state_hash(model) != info['model_state_sha256']:
        raise ValueError('Model tensor state disagrees with audited identity')
    del checkpoint
    return GapDecoder(model, model_spec['mode'])


def contexts(record, case):
    """Split known flanks and hidden truth; only the caller handles all three."""
    rna = record['rna']
    if hashlib.sha256(rna.encode()).hexdigest() != case['sequence_sha256']:
        raise ValueError('Case RNA identity mismatch')
    s, g = case['start_codon'], case['gap_codons']
    left = case['left_context_codons']
    right = case['right_context_codons']
    m = len(rna) // 3 - 1
    if type(s) is not int or type(g) is not int or s < 1 or g < 1 or s + g > m:
        raise ValueError('Gap must be internal to sense codons after initiation')
    if left != 'full' and (type(left) is not int or not 1 <= left <= s):
        raise ValueError('Requested left context is unavailable')
    if type(right) is not int or right < 0 or s + g + right > m:
        raise ValueError('Requested right sense-codon context is unavailable')
    offset = 0 if left == 'full' else (s - left) * 3
    return {'prefix': rna[offset:s * 3], 'truth': rna[s * 3:(s + g) * 3],
            'right': rna[(s + g) * 3:(s + g + right) * 3], 'offset_nt': offset}


def conservation_outcomes(support, model_support, record, case, truth, prediction):
    common = support.describe(record, case['start_codon'], truth, prediction)
    # Allele masks use the identical all-training union for every comparator.
    # Whether an organism/family appeared in this model is a different property.
    visible = model_support.describe(record, case['start_codon'], truth, prediction)
    return {**common, 'conservation_reference': 'all_excluded_training_union',
            **{key: visible[key] for key in ('family_seen', 'taxon_seen', 'family_taxon_seen')}}


def bank_identity(case, model_spec, plan_hash):
    return {'plan_sha256': plan_hash, 'model_id': model_spec['model_id'],
            'sequence_sha256': case['sequence_sha256'], 'anchor_id': case['anchor_id'],
            'start_codon': case['start_codon'], 'gap_codons': case['gap_codons'],
            'left_context_codons': case['left_context_codons']}


def evaluate_case(decoder, case, record, model_spec, union_support, model_support,
                  baselines, plan_hash, bank_directory, beam_width=8):
    parts = contexts(record, case)
    prefix, truth, right, offset = [parts[k] for k in ('prefix', 'truth', 'right', 'offset_nt')]
    result = {'schema_version': 1, 'plan_sha256': plan_hash,
              'model_id': model_spec['model_id'], 'arm': model_spec['arm'],
              'seed': model_spec['seed'], 'mode': model_spec['mode'],
              'dataset': model_spec['dataset'], 'case': case,
              'truth': truth, 'initial_left_codons': len(prefix) // 3,
              'absolute_offset_nt': offset, 'novelty': record.get('novelty', {})}
    if case['stage'] == 'right_flank':
        identity = bank_identity(case, model_spec, plan_hash)
        bank_path = Path(bank_directory) / (digest(identity) + '.json')
        if bank_path.exists():
            bank = read_json(bank_path)
            if bank.get('identity') != identity or digest(bank['candidates']) != bank['candidates_sha256']:
                raise ValueError('Frozen candidate bank identity changed')
        else:
            # This function receives no hidden gap or suffix. Commit the bank
            # before suffix scoring or measuring candidate correctness.
            candidates = decoder.proposals(prefix, case['gap_codons'], width=beam_width, offset_nt=offset)
            bank = {'identity': identity, 'created_utc': now(), 'beam_width': beam_width,
                    'candidates': candidates, 'candidates_sha256': digest(candidates)}
            atomic_json(bank_path, bank)
        scored = decoder.rerank(prefix, bank['candidates'], right, offset)
        prediction = scored[0]['rna']
        candidate_metrics = [reconstruction(truth, c['rna'], record['translation_table'])
                             for c in bank['candidates']]
        result.update(method='beam_rerank', candidate_bank_sha256=sha256(bank_path),
                      candidates_sha256=bank['candidates_sha256'], candidate_count=len(candidate_metrics),
                      candidate_exact_recall=int(any(c['exact_span'] for c in candidate_metrics)),
                      candidate_max_correct_codons=max(c['correct_codons'] for c in candidate_metrics),
                      candidate_max_correct_amino_acids=max(c['correct_amino_acids'] for c in candidate_metrics),
                      selected_gap_log_probability=scored[0]['log_probability'],
                      selected_right_log_probability=scored[0]['right_log_probability'],
                      selected_joint_log_score=scored[0]['joint_log_score'])
    else:
        # Score truth and generate in different calls; generation accepts only
        # the observed prefix and length, never the scored true continuation.
        logp = decoder.continuation_log_probability(prefix, truth, offset)
        prediction = decoder.proposals(prefix, case['gap_codons'], width=1, offset_nt=offset)[0]['rna']
        result.update(method='codon_synchronous_greedy', gap_log_probability=logp,
                      gap_bits_per_base=-logp / (len(truth) * math.log(2)))
        baseline_results = {}
        for method in METHODS:
            b = baselines.score_and_generate(method, prefix, truth, case['start_codon'], record)
            b.update(reconstruction(truth, b['prediction'], record['translation_table']))
            baseline_results[method] = b
        result['baselines'] = baseline_results
    result.update(prediction=prediction, **reconstruction(truth, prediction, record['translation_table']))
    result.update(conservation_outcomes(union_support, model_support, record, case, truth, prediction))
    return result


def validate_completed_result(result, plan_hash, model_spec, case):
    if result.get('plan_sha256') != plan_hash or result.get('model_id') != model_spec['model_id'] or result.get('case') != case:
        raise ValueError('Saved evaluation result does not match the frozen plan')
    for key in ('prediction', 'truth'):
        if len(result[key]) != 3 * case['gap_codons'] or set(result[key]) - set('AUCG'):
            raise ValueError('Saved gap result has an invalid sequence')
    measured = reconstruction(result['truth'], result['prediction'], case['translation_table'])
    if any(result.get(k) != v for k, v in measured.items()):
        raise ValueError('Saved reconstruction metrics disagree with saved sequences')
    if 'gap_bits_per_base' in result:
        if not math.isfinite(result['gap_bits_per_base']) or result['gap_bits_per_base'] < 0:
            raise ValueError('Invalid saved probability score')
