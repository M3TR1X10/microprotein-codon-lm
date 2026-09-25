"""Frozen-cohort secondary experiments; every reported outcome is training fit.

This module deliberately leaves the first pilot's implementation untouched.
The supplied CDS have already passed biological QC; this reader independently
checks identity, alphabet, reading frame, labels, and training-only boundaries.
"""
import copy
import ctypes
import hashlib
import json
import math
import os
import platform
import time
from collections import Counter, defaultdict
from contextlib import nullcontext
from dataclasses import asdict
from pathlib import Path

import numpy as np
import torch
from torch.nn import functional as F

from .io import now, read_json, sha256
from .model import Decoder, ModelConfig
from .secondary_io import write_json
from .tokenization import Tokenizer
from .train import git_state, lr_at


SCOPE = 'TRAINING ONLY; no validation, test, gap completion, or structure scoring'


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def engine_source_hashes():
    root = Path(__file__).parent
    return {name: sha256(root / name) for name in
            ('secondary_train.py', 'secondary_io.py', 'model.py', 'tokenization.py', 'train.py', 'io.py')}


class SecondaryCorpus:
    """Complete independent CDS, with no concatenation, truncation or special token."""

    def __init__(self, cohort_dir, mode):
        self.root = Path(cohort_dir)
        self.mode = mode
        self.tokenizer = Tokenizer(mode)
        self.manifest = read_json(self.root / 'manifest.json')
        cohort_hash = sha256(self.root / 'cohort.jsonl')
        if cohort_hash != self.manifest['cohort_sha256']:
            raise ValueError('Cohort checksum does not match its frozen manifest')
        if self.manifest.get('stage', 'training_only') != 'training_only':
            raise ValueError('Only training-only cohorts are supported')
        if any(self.manifest.get(k, 0) for k in ('validation_sequences', 'test_sequences')):
            raise ValueError('Validation and test cohorts are not supported')
        self.records = [json.loads(line) for line in
                        (self.root / 'cohort.jsonl').read_text(encoding='utf-8').splitlines()
                        if line.strip()]
        if not self.records:
            raise ValueError('Cohort must contain at least one CDS')
        self.sequences = []
        for r in self.records:
            rna = r['rna']
            if len(rna) < 6 or len(rna) % 3:
                raise ValueError('Every CDS requires complete codons and at least one target codon')
            if hashlib.sha256(rna.encode()).hexdigest() != r['sequence_sha256']:
                raise ValueError('Sequence checksum mismatch')
            if not isinstance(r['family'], str) or not r['family'].strip():
                raise ValueError('Each CDS requires a nonempty family label')
            if any(type(r[k]) is not int or r[k] <= 0 for k in ('tax_id', 'translation_table')):
                raise ValueError('Taxon and translation table must be positive integers')
            self.sequences.append(np.asarray(self.tokenizer.encode(rna), dtype=np.int64))
        ids = [r['sequence_sha256'] for r in self.records]
        if len(set(ids)) != len(ids):
            raise ValueError('Cohort must contain distinct RNA sequences; duplicates are not new data')
        self.families = [r['family'] for r in self.records]
        self.taxa = [str(r['tax_id']) for r in self.records]
        self.codes = [str(r['translation_table']) for r in self.records]
        self.first_target = 3 // self.tokenizer.width
        self.target_counts = np.asarray([len(s) - self.first_target for s in self.sequences])
        self.family_target_tokens = dict(Counter())
        for family, count in zip(self.families, self.target_counts):
            self.family_target_tokens[family] = self.family_target_tokens.get(family, 0) + int(count)
        self.family_counts = dict(Counter(self.families))
        self.variable_masks = self._variable_masks()
        self.identity = {'cohort_sha256': cohort_hash,
                         'manifest_sha256': sha256(self.root / 'manifest.json'),
                         'ordered_sequence_metadata_sha256': _digest([
                             {k: r[k] for k in ('sequence_sha256', 'family', 'tax_id', 'translation_table')}
                             for r in self.records])}

    def __len__(self):
        return len(self.sequences)

    def _variable_masks(self):
        # A variable codon has >1 observed triplet at an aligned CDS position
        # within the same family AND taxon. Absent tail positions are not alleles.
        # The base arm scores all three bases of those same variable codons.
        alleles = defaultdict(lambda: defaultdict(set))
        for r in self.records:
            for pos in range(1, len(r['rna']) // 3):
                alleles[(r['family'], r['tax_id'])][pos].add(r['rna'][3*pos:3*pos+3])
        masks = []
        width = self.tokenizer.width
        for r, seq in zip(self.records, self.sequences):
            mask = np.zeros(len(seq) - 1, dtype=bool)
            for pos, values in alleles[(r['family'], r['tax_id'])].items():
                if len(values) > 1 and 3*pos < len(r['rna']):
                    # Target seq[j] sits at index j-1 in the shifted loss array.
                    mask[3*pos//width-1:3*(pos+1)//width-1] = True
            masks.append(mask)
        return masks

    def batch(self, indices, device='cpu'):
        selected = [self.sequences[int(i)] for i in indices]
        length = max(len(s) for s in selected) - 1
        x = torch.zeros((len(selected), length), dtype=torch.long)
        y = torch.full_like(x, -100)
        for row, seq in enumerate(selected):
            x[row, :len(seq)-1] = torch.from_numpy(seq[:-1])
            y[row, :len(seq)-1] = torch.from_numpy(seq[1:])
            y[row, :self.first_target-1] = -100
        return x.to(device), y.to(device)


def build_secondary_model(config, corpus):
    maximum = config['max_nt_context'] // corpus.tokenizer.width - 1
    if max(map(len, corpus.sequences)) - 1 > maximum:
        raise ValueError('Context is too short; truncation is not permitted')
    c = ModelConfig(vocab_size=len(corpus.tokenizer.vocabulary), max_tokens=maximum,
                    token_width=corpus.tokenizer.width, **config['model'])
    model = Decoder(c)
    if not config.get('position_encoding', True):
        model.position.zero_()
    return model


def objective_loss(logits, targets, indices, corpus, objective, effective_batch_size, target_count):
    """One microbatch's contribution; sum these before the optimizer update.

    For family_macro, uniform-CDS sampling gives the unbiased estimator
    N/(F*B) sum_i L_i/T_family(i). There is no random minibatch-weight denominator.
    """
    losses = F.cross_entropy(logits.reshape(-1, logits.shape[-1]), targets.reshape(-1),
                             ignore_index=-100, reduction='none').reshape_as(targets)
    if objective == 'token':
        return losses.sum() / target_count
    if objective != 'family_macro':
        raise ValueError('loss must be token or family_macro')
    weights = torch.tensor([1.0/corpus.family_target_tokens[corpus.families[int(i)]]
                            for i in indices], device=logits.device, dtype=losses.dtype)
    return (losses.sum(dim=1)*weights).sum()*len(corpus)/(len(corpus.family_counts)*effective_batch_size)


def _bucket():
    return {'nll': 0.0, 'correct': 0, 'tokens': 0}


def _add(bucket, nll, correct, tokens):
    bucket['nll'] += float(nll)
    bucket['correct'] += int(correct)
    bucket['tokens'] += int(tokens)


def _metrics(bucket, width):
    n = bucket['tokens']
    return {'training_nll_per_token': bucket['nll']/n if n else None,
            'training_bits_per_base': bucket['nll']/(n*width*math.log(2)) if n else None,
            'training_token_accuracy': bucket['correct']/n if n else None,
            'scored_training_tokens': n, 'scored_training_bases': n*width}


@torch.no_grad()
def training_diagnostics(model, corpus, batch_size, device):
    was_training = model.training
    model.eval()
    totals, variable = _bucket(), _bucket()
    families, codes, taxa = defaultdict(_bucket), defaultdict(_bucket), defaultdict(_bucket)
    family_variables = defaultdict(_bucket)
    try:
        for start in range(0, len(corpus), batch_size):
            indices = list(range(start, min(start+batch_size, len(corpus))))
            x, y = corpus.batch(indices, device)
            logits = model(x)
            loss = F.cross_entropy(logits.reshape(-1, logits.shape[-1]), y.reshape(-1),
                                   ignore_index=-100, reduction='none').reshape_as(y)
            values = loss.float().cpu().numpy()
            correct = ((logits.argmax(-1) == y) & (y != -100)).cpu().numpy()
            for row, i in enumerate(indices):
                nll, right, count = float(values[row].sum(dtype=np.float64)), int(correct[row].sum()), int(corpus.target_counts[i])
                for bucket in (totals, families[corpus.families[i]], codes[corpus.codes[i]], taxa[corpus.taxa[i]]):
                    _add(bucket, nll, right, count)
                mask = corpus.variable_masks[i]
                vnll = float(values[row, :len(mask)][mask].sum(dtype=np.float64))
                vright, vcount = int(correct[row, :len(mask)][mask].sum()), int(mask.sum())
                _add(variable, vnll, vright, vcount)
                _add(family_variables[corpus.families[i]], vnll, vright, vcount)
    finally:
        model.train(was_training)
    width = corpus.tokenizer.width
    per_family = {k: _metrics(v, width) for k, v in sorted(families.items())}
    return {**_metrics(totals, width),
            'family_macro_bits_per_base': sum(v['training_bits_per_base'] for v in per_family.values())/len(per_family),
            'per_family': per_family,
            'per_code': {k: _metrics(v, width) for k, v in sorted(codes.items())},
            'per_taxon': {k: _metrics(v, width) for k, v in sorted(taxa.items())},
            'variable_positions': {**_metrics(variable, width),
                'definition': 'All bases of variable codons within family and taxon; aligned from CDS start; first codon excluded',
                'per_family': {k: _metrics(v, width) for k, v in sorted(family_variables.items())}}}


def fit_training_baselines(corpus, cache_path=None):
    """Fit and score training-only categorical distributions; cache by full identity.

    family_position has privileged family metadata. No model receives this label.
    Empirical probabilities are scored only on the events used to fit them.
    """
    path = Path(cache_path) if cache_path is not None else corpus.root/f'baselines-{corpus.mode}.json'
    identity = {'data': corpus.identity, 'mode': corpus.mode, 'source_hashes': engine_source_hashes(),
                'prior_total_mass': 2.0, 'version': 1}
    if path.exists():
        old = read_json(path)
        if old.get('cache_identity') == identity:
            return old
    vocab, width = len(corpus.tokenizer.vocabulary), corpus.tokenizer.width
    maximum = max(map(len, corpus.sequences))
    unigram, bigram = np.zeros(vocab), np.zeros((vocab, vocab))
    position = np.zeros((maximum, vocab))
    family_position = {f: np.zeros((maximum, vocab)) for f in corpus.family_counts}
    events = []
    for i, seq in enumerate(corpus.sequences):
        positions = np.arange(corpus.first_target, len(seq))
        targets, previous = seq[positions], seq[positions-1]
        np.add.at(unigram, targets, 1)
        np.add.at(bigram, (previous, targets), 1)
        np.add.at(position, (positions, targets), 1)
        np.add.at(family_position[corpus.families[i]], (positions, targets), 1)
        events.append((positions, previous, targets))
    report = {'scope': SCOPE, 'cache_identity': identity, 'mode': corpus.mode,
              'prior_total_mass': 2.0, 'prior_per_category': 2.0/vocab,
              'oracle_warning': 'family_position knows the family label; models receive sequence tokens only',
              'metrics': {'uniform': {'training_bits_per_base': math.log2(vocab)/width}},
              'target_bases': int(corpus.target_counts.sum())*width}
    for label, alpha in (('empirical', 0.0), ('equal_total_prior_2', 2.0/vocab)):
        result = {}
        for name in ('unigram', 'bigram', 'position', 'family_position'):
            totals, families = _bucket(), defaultdict(_bucket)
            for i, (positions, previous, targets) in enumerate(events):
                if name == 'unigram':
                    counts = np.broadcast_to(unigram, (len(targets), vocab))
                elif name == 'bigram':
                    counts = bigram[previous]
                elif name == 'position':
                    counts = position[positions]
                else:
                    counts = family_position[corpus.families[i]][positions]
                denominator = counts.sum(axis=1) + alpha*vocab
                probability = (counts[np.arange(len(targets)), targets]+alpha)/denominator
                nll = float(-np.log(probability).sum())
                correct = int((counts.argmax(axis=1) == targets).sum())
                _add(totals, nll, correct, len(targets))
                _add(families[corpus.families[i]], nll, correct, len(targets))
            per_family = {f: _metrics(b, width) for f, b in sorted(families.items())}
            result[name] = {**_metrics(totals, width), 'per_family': per_family,
                           'family_macro_bits_per_base': sum(b['training_bits_per_base'] for b in per_family.values())/len(per_family)}
        report['metrics'][label] = result
    write_json(path, report)
    return report


def process_peak_rss_bytes():
    """Process lifetime peak, including earlier runs; explicitly not incremental RSS."""
    try:
        if os.name == 'nt':
            from ctypes import wintypes
            class Counters(ctypes.Structure):
                _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD),
                            ('PeakWorkingSetSize', ctypes.c_size_t), ('WorkingSetSize', ctypes.c_size_t),
                            ('QuotaPeakPagedPoolUsage', ctypes.c_size_t), ('QuotaPagedPoolUsage', ctypes.c_size_t),
                            ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t), ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                            ('PagefileUsage', ctypes.c_size_t), ('PeakPagefileUsage', ctypes.c_size_t)]
            counters = Counters()
            counters.cb = ctypes.sizeof(counters)
            kernel = ctypes.WinDLL('kernel32', use_last_error=True)
            kernel.GetCurrentProcess.restype = wintypes.HANDLE
            psapi = ctypes.WinDLL('psapi', use_last_error=True)
            psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
            if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
                return None
            return int(counters.PeakWorkingSetSize)
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return int(peak if platform.system() == 'Darwin' else peak*1024)
    except (ImportError, AttributeError, OSError):
        return None


def _atomic_checkpoint(path, payload):
    temporary = path.with_suffix('.partial')
    torch.save(payload, temporary)
    temporary.replace(path)


def _configuration(config):
    config = copy.deepcopy(config)
    for key, value in {'loss': 'token', 'position_encoding': True, 'precision': 'float32',
                       'device': 'auto', 'deterministic': True}.items():
        config.setdefault(key, value)
    if config.get('stage') != 'training_only':
        raise ValueError('Secondary experiments support training_only only')
    for key in ('batch_size', 'gradient_accumulation', 'updates', 'diagnostic_interval',
                'checkpoint_interval', 'cpu_threads', 'max_nt_context'):
        if type(config[key]) is not int or config[key] < 1:
            raise ValueError(f'{key} must be a positive integer')
    if config['max_nt_context'] < 6 or config['max_nt_context'] % 3:
        raise ValueError('max_nt_context must be a multiple of three >= 6')
    if config['loss'] not in ('token', 'family_macro') or type(config['position_encoding']) is not bool:
        raise ValueError('Invalid loss or position_encoding setting')
    if config['precision'] not in ('float32', 'bfloat16'):
        raise ValueError('precision must be float32 or bfloat16')
    if not math.isfinite(config['gradient_clip']) or config['gradient_clip'] <= 0:
        raise ValueError('gradient_clip must be positive and finite')
    lr_at(0, config)
    return config


def train_secondary(config, mode, seed, cohort_dir, output, resume=False, stop_after=None, max_seconds=None):
    """Train one frozen secondary arm. stop_after is an absolute update index.

    Only checkpoint.pt (rolling latest) and final.pt (completed run) are retained.
    A partial stop records diagnostics but never represents a completed experiment.
    max_seconds is an invocation wall-time guard, not an experimental setting.
    It stops at an optimizer boundary and then writes diagnostics/checkpoints;
    callers must reserve time for that orderly shutdown. Setup is included.
    """
    setup_started = time.monotonic()
    config = _configuration(config)
    if type(seed) is not int or type(resume) is not bool:
        raise ValueError('seed must be an integer and resume must be boolean')
    if stop_after is not None and (type(stop_after) is not int or stop_after < 1):
        raise ValueError('stop_after must be a positive absolute optimizer-update index')
    if max_seconds is not None and (not isinstance(max_seconds, (int, float))
                                    or not math.isfinite(max_seconds) or max_seconds <= 0):
        raise ValueError('max_seconds must be positive and finite')
    path = Path(output)
    if path.exists() and any(path.iterdir()) and not resume:
        raise FileExistsError(f'Run exists: {path}; choose a new output or resume')
    corpus = SecondaryCorpus(cohort_dir, mode)
    source = engine_source_hashes()
    identity = {'config': config, 'mode': mode, 'seed': seed, 'data_identity': corpus.identity,
                'source_hashes': source}
    checkpoint = None
    if resume:
        checkpoint = torch.load(path/'checkpoint.pt', map_location='cpu', weights_only=True)
        if checkpoint['identity'] != identity:
            raise ValueError('Resume requires identical config, seed, mode, cohort, metadata and source code')
    torch.set_num_threads(config['cpu_threads'])
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(config['deterministic'])
    requested = config['device']
    device = torch.device(('cuda' if torch.cuda.is_available() else 'cpu') if requested == 'auto' else requested)
    mixed = config['precision'] == 'bfloat16'
    if mixed and (device.type != 'cuda' or not torch.cuda.is_bf16_supported()):
        raise ValueError('bfloat16 requires a supported CUDA device')
    if device.type == 'cuda':
        torch.cuda.reset_peak_memory_stats(device)
    model = build_secondary_model(config, corpus).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=config['learning_rate'],
                                  betas=(0.9, 0.95), weight_decay=config['weight_decay'])
    sample_rng = torch.Generator().manual_seed(seed+100000)
    history, update_start, elapsed, bases_seen, sequences_seen = [], 0, 0.0, 0, 0
    family_exposures = {f: {'sequences': 0, 'target_bases': 0} for f in corpus.family_counts}
    sampling_trace = '0'*64
    peak_cuda_allocated, peak_cuda_reserved = 0, 0
    if checkpoint is not None:
        model.load_state_dict(checkpoint['model'])
        optimizer.load_state_dict(checkpoint['optimizer'])
        torch.set_rng_state(checkpoint['torch_rng'])
        sample_rng.set_state(checkpoint['sample_rng'])
        if device.type == 'cuda':
            torch.cuda.set_rng_state_all(checkpoint['cuda_rng'])
        history = checkpoint['history']
        update_start, elapsed = checkpoint['update'], checkpoint['elapsed_seconds']
        bases_seen, sequences_seen = checkpoint['bases_seen'], checkpoint['sequences_seen']
        family_exposures, sampling_trace = checkpoint['per_family_exposures'], checkpoint['sampling_trace_sha256']
        peak_cuda_allocated = checkpoint.get('peak_cuda_allocated_bytes', 0)
        peak_cuda_reserved = checkpoint.get('peak_cuda_reserved_bytes', 0)
    end_update = min(config['updates'], stop_after) if stop_after is not None else config['updates']
    if end_update <= update_start:
        raise ValueError('No optimizer updates remain to run')
    baselines = fit_training_baselines(corpus)
    path.mkdir(parents=True, exist_ok=True)
    write_json(path/'baselines.json', baselines)
    if checkpoint is None:
        write_json(path/'run.json', {'created_utc': now(), **identity, 'scope': SCOPE,
            'model': asdict(model.config), 'parameters': sum(p.numel() for p in model.parameters()),
            'python': platform.python_version(), 'torch': str(torch.__version__),
            'platform': platform.platform(), 'device': str(device), 'git': git_state(),
            'sequences': len(corpus), 'family_counts': corpus.family_counts,
            'target_bases_per_corpus_pass': int(corpus.target_counts.sum())*corpus.tokenizer.width,
            'vocabulary': corpus.tokenizer.vocabulary,
            'checkpoint_retention': 'rolling checkpoint.pt and final.pt upon completion',
            'diagnostic_precision': 'float32',
            'sampling': 'Uniform distinct CDS with replacement; dedicated seeded index generator',
            'variable_position_definition': 'Within family and taxon, all bases of observed variable codons'} )
        history.append({'update': 0, 'sequences_seen': 0, 'bases_seen': 0,
                        **training_diagnostics(model, corpus, config['batch_size'], device)})
    setup_seconds = time.monotonic()-setup_started
    started = time.monotonic()
    last_progress_at = started
    effective_batch = config['batch_size']*config['gradient_accumulation']
    stop_reason = 'completed' if end_update == config['updates'] else 'stop_after'
    completed_update = update_start
    for update in range(update_start, end_update):
        model.train()
        learning_rate = lr_at(update, config)
        for group in optimizer.param_groups:
            group['lr'] = learning_rate
        optimizer.zero_grad(set_to_none=True)
        draws = torch.randint(len(corpus), (config['gradient_accumulation'], config['batch_size']), generator=sample_rng)
        flat = draws.flatten().tolist()
        target_count = int(corpus.target_counts[flat].sum())
        # Materialize one microbatch on the device at a time.
        for indices in draws.tolist():
            x, y = corpus.batch(indices, device)
            context = torch.autocast('cuda', dtype=torch.bfloat16) if mixed else nullcontext()
            with context:
                loss = objective_loss(model(x), y, indices, corpus, config['loss'], effective_batch, target_count)
            if not torch.isfinite(loss):
                raise FloatingPointError('Non-finite secondary training loss')
            loss.backward()
        norm = torch.nn.utils.clip_grad_norm_(model.parameters(), config['gradient_clip'], error_if_nonfinite=True)
        optimizer.step()
        bases_seen += target_count*corpus.tokenizer.width
        sequences_seen += effective_batch
        for i in flat:
            family_exposures[corpus.families[i]]['sequences'] += 1
            family_exposures[corpus.families[i]]['target_bases'] += int(corpus.target_counts[i])*corpus.tokenizer.width
        sampling_trace = hashlib.sha256(bytes.fromhex(sampling_trace)+draws.numpy().tobytes()).hexdigest()
        current = update+1
        completed_update = current
        step_finished = time.monotonic()
        budget_hit = (max_seconds is not None and step_finished-setup_started >= max_seconds
                      and current < config['updates'])
        if budget_hit:
            stop_reason = 'max_seconds'
        diagnostic_due = current % config['diagnostic_interval'] == 0 or current == end_update or budget_hit
        if not diagnostic_due and step_finished-last_progress_at >= 60:
            print(f"PROGRESS secondary {config.get('dataset', corpus.root.name)} {mode} seed={seed} "
                  f"update={current}/{config['updates']} elapsed={elapsed+step_finished-started:.0f}s",
                  flush=True)
            last_progress_at = step_finished
        if diagnostic_due:
            history.append({'update': current, 'learning_rate': learning_rate, 'gradient_norm': float(norm),
                            'sequences_seen': sequences_seen, 'bases_seen': bases_seen,
                            **training_diagnostics(model, corpus, config['batch_size'], device)})
            print(f"secondary {config.get('dataset', corpus.root.name)} {mode} seed={seed} update={current}: "
                  f"TRAIN bits/base={history[-1]['training_bits_per_base']:.5f}", flush=True)
            last_progress_at = time.monotonic()
        if current % config['checkpoint_interval'] == 0 or current == end_update or budget_hit:
            if device.type == 'cuda':
                peak_cuda_allocated = max(peak_cuda_allocated, int(torch.cuda.max_memory_allocated(device)))
                peak_cuda_reserved = max(peak_cuda_reserved, int(torch.cuda.max_memory_reserved(device)))
            payload = {'identity': identity, 'update': current, 'model': model.state_dict(),
                'optimizer': optimizer.state_dict(), 'torch_rng': torch.get_rng_state(),
                'sample_rng': sample_rng.get_state(),
                'cuda_rng': torch.cuda.get_rng_state_all() if device.type == 'cuda' else [],
                'history': history, 'elapsed_seconds': elapsed+time.monotonic()-started,
                'bases_seen': bases_seen, 'sequences_seen': sequences_seen,
                'per_family_exposures': family_exposures, 'sampling_trace_sha256': sampling_trace,
                'peak_cuda_allocated_bytes': peak_cuda_allocated, 'peak_cuda_reserved_bytes': peak_cuda_reserved}
            _atomic_checkpoint(path/'checkpoint.pt', payload)
            if current == config['updates']:
                _atomic_checkpoint(path/'final.pt', payload)
            write_json(path/'metrics.json', history)
        if budget_hit:
            break
    elapsed += time.monotonic()-started
    summary = {'scope': SCOPE, 'dataset': config.get('dataset', corpus.root.name),
        'mode': mode, 'seed': seed, 'config': config, 'updates': completed_update,
        'complete': completed_update == config['updates'], 'stop_reason': stop_reason,
        'runtime_max_seconds_this_invocation': max_seconds, 'initial': history[0], 'final': history[-1],
        'parameters': sum(p.numel() for p in model.parameters()), 'data_identity': corpus.identity,
        'cohort_sha256': corpus.identity['cohort_sha256'], 'source_hashes': source,
        'sequences': len(corpus), 'bases_seen': bases_seen, 'sequences_seen': sequences_seen,
        'sequence_exposures_per_unique_cds': sequences_seen/len(corpus),
        'per_family_exposures': family_exposures, 'sampling_trace_sha256': sampling_trace,
        'elapsed_seconds': elapsed, 'setup_seconds_this_invocation': setup_seconds,
        'timing_scope': 'Training loop, interval diagnostics and checkpoint writes; setup and initial diagnostics excluded',
        'process_peak_rss_bytes': process_peak_rss_bytes(),
        'rss_scope': 'Process-lifetime high-water mark, including earlier runs in this process',
        'peak_cuda_allocated_bytes': peak_cuda_allocated if device.type == 'cuda' else None,
        'peak_cuda_reserved_bytes': peak_cuda_reserved if device.type == 'cuda' else None}
    write_json(path/'summary.json', summary)
    return summary
