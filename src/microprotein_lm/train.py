"""Training and training-set diagnostics only. No validation/test data path."""
import math
import platform
import subprocess
import time
from contextlib import nullcontext
from dataclasses import asdict
from pathlib import Path

import torch

from .data import Corpus
from .io import read_json, write_json, now, sha256
from .model import ModelConfig, Decoder, summed_loss


def lr_at(update, config):
    warmup, total = config['warmup_updates'], config['updates']
    peak, floor = config['learning_rate'], config['min_learning_rate']
    if not 0 <= warmup < total or not 0 < floor <= peak:
        raise ValueError('Require 0 <= warmup < updates and 0 < minimum LR <= peak LR')
    if update < warmup:
        return peak*(update+1)/warmup
    progress = (update-warmup)/max(1,total-warmup-1)
    return floor+(peak-floor)*0.5*(1+math.cos(math.pi*progress))


@torch.no_grad()
def training_diagnostics(model, corpus, batch_size, device):
    model.eval()
    nll, correct, tokens = 0.,0,0
    for start in range(0,len(corpus),batch_size):
        x,y = corpus.batch(list(range(start,min(start+batch_size,len(corpus)))),device)
        logits = model(x)
        nll += float(summed_loss(logits,y))
        mask = y != -100
        tokens += int(mask.sum())
        correct += int(((logits.argmax(-1) == y)&mask).sum())
    model.train()
    return {'training_nll_per_token':nll/tokens,
        'training_bits_per_base':nll/(tokens*corpus.tokenizer.width*math.log(2)),
        'training_token_accuracy':correct/tokens,
        'scored_training_tokens':tokens,'scored_training_bases':tokens*corpus.tokenizer.width}


def git_state():
    def read(*args):
        try:
            return subprocess.check_output(['git',*args],stderr=subprocess.DEVNULL,text=True).strip()
        except (OSError,subprocess.CalledProcessError):
            return None
    return {'commit':read('rev-parse','HEAD'),'dirty':read('status','--porcelain')}


def source_hashes():
    return {p.name:sha256(p) for p in Path(__file__).parent.glob('*.py')}


def train(config, mode, seed, processed='data/processed', output=None, resume=False, stop_after=None):
    config = dict(config)
    if mode not in ('base','codon') or config.get('stage') != 'training_only':
        raise ValueError('This implementation supports training-only base/codon experiments')
    for key in ['batch_size','gradient_accumulation','updates','diagnostic_interval','checkpoint_interval','cpu_threads']:
        if not isinstance(config[key],int) or config[key] < 1:
            raise ValueError(f'{key} must be a positive integer')
    if config['max_nt_context'] < 6 or config['max_nt_context'] % 3:
        raise ValueError('max_nt_context must be a multiple of three >= 6')
    lr_at(0,config)
    torch.set_num_threads(config['cpu_threads'])
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(config.get('deterministic',True))
    requested = config.get('device','auto')
    device = torch.device(('cuda' if torch.cuda.is_available() else 'cpu') if requested == 'auto' else requested)
    if config.get('precision','float32') not in ('float32','bfloat16'):
        raise ValueError('Supported precision: float32, bfloat16')
    mixed = config.get('precision') == 'bfloat16'
    if mixed and (device.type != 'cuda' or not torch.cuda.is_bf16_supported()):
        raise ValueError('bfloat16 requires a supported CUDA device; use float32 here')
    corpus = Corpus(processed,mode)
    c = ModelConfig(vocab_size=len(corpus.tokenizer.vocabulary),
        max_tokens=config['max_nt_context']//corpus.tokenizer.width-1,
        token_width=corpus.tokenizer.width, **config['model'])
    if corpus.metadata['max_sequence_tokens']-1 > c.max_tokens:
        raise ValueError('Context too short; truncation is not permitted')
    model = Decoder(c).to(device)
    optimizer = torch.optim.AdamW(model.parameters(),lr=config['learning_rate'],
        betas=(0.9,0.95),weight_decay=config['weight_decay'])
    sample_rng = torch.Generator().manual_seed(seed+100000)
    path = Path(output or f'runs/{mode}-seed{seed}')
    if path.exists() and any(path.iterdir()) and not resume:
        raise FileExistsError(f'Run exists: {path}; choose a new output or --resume')
    path.mkdir(parents=True,exist_ok=True)
    source = source_hashes()
    start_update, history, elapsed, bases_seen, sequences_seen = 0,[],0.,0,0
    if resume:
        checkpoint = torch.load(path/'checkpoint.pt',map_location='cpu',weights_only=True)
        if (checkpoint['config'] != config or checkpoint['mode'] != mode or checkpoint['seed'] != seed
                or checkpoint['cohort_sha256'] != corpus.metadata['cohort_sha256']
                or checkpoint['source_hashes'] != source):
            raise ValueError('Resume requires identical config, seed, mode, cohort and source code')
        model.load_state_dict(checkpoint['model'])
        optimizer.load_state_dict(checkpoint['optimizer'])
        torch.set_rng_state(checkpoint['torch_rng'])
        if device.type == 'cuda':
            torch.cuda.set_rng_state_all(checkpoint['cuda_rng'])
        sample_rng.set_state(checkpoint['sample_rng'])
        start_update,history = checkpoint['update'],checkpoint['history']
        elapsed,bases_seen,sequences_seen = checkpoint['elapsed_seconds'],checkpoint['bases_seen'],checkpoint['sequences_seen']
    else:
        write_json(path/'run.json',{'created_utc':now(),'config':config,'mode':mode,'seed':seed,
            'model':asdict(c),'parameters':sum(p.numel() for p in model.parameters()),
            'cohort_sha256':corpus.metadata['cohort_sha256'],'source_hashes':source,
            'python':platform.python_version(),'torch':str(torch.__version__),
            'device':str(device),'platform':platform.platform(),'git':git_state(),
            'scope':'Training-only pilot; no hyperparameter winner or predictive accuracy is inferred'})
        history.append({'update':0,**training_diagnostics(model,corpus,config['batch_size'],device)})
    end_update = min(config['updates'],stop_after) if stop_after is not None else config['updates']
    if end_update <= start_update:
        raise ValueError('No updates left to run')
    started = time.monotonic()
    for update in range(start_update,end_update):
        model.train()
        learning_rate = lr_at(update,config)
        for group in optimizer.param_groups:
            group['lr'] = learning_rate
        optimizer.zero_grad(set_to_none=True)
        indices = torch.randint(len(corpus),(config['gradient_accumulation'],config['batch_size']),generator=sample_rng)
        batches = [corpus.batch(idx.tolist(),device) for idx in indices]
        target_count = sum(int((y != -100).sum()) for _,y in batches)
        for x,y in batches:
            context = torch.autocast('cuda',dtype=torch.bfloat16) if mixed else nullcontext()
            with context:
                loss = summed_loss(model(x),y)/target_count
            if not torch.isfinite(loss):
                raise FloatingPointError('Non-finite training loss')
            loss.backward()
        norm = torch.nn.utils.clip_grad_norm_(model.parameters(),config['gradient_clip'],error_if_nonfinite=True)
        optimizer.step()
        bases_seen += target_count*corpus.tokenizer.width
        sequences_seen += indices.numel()
        current = update+1
        if current % config['diagnostic_interval'] == 0 or current == end_update:
            entry = {'update':current,'learning_rate':learning_rate,'gradient_norm':float(norm),
                **training_diagnostics(model,corpus,config['batch_size'],device)}
            history.append(entry)
            print(f"{mode} seed={seed} update={current}: TRAIN bits/base={entry['training_bits_per_base']:.4f}",flush=True)
        if current % config['checkpoint_interval'] == 0 or current == end_update:
            checkpoint = {'config':config,'mode':mode,'seed':seed,'update':current,
                'cohort_sha256':corpus.metadata['cohort_sha256'],'source_hashes':source,
                'model':model.state_dict(),'optimizer':optimizer.state_dict(),
                'torch_rng':torch.get_rng_state(), 'sample_rng':sample_rng.get_state(),
                'cuda_rng':torch.cuda.get_rng_state_all() if device.type == 'cuda' else [],
                'history':history,'elapsed_seconds':elapsed+time.monotonic()-started,
                'bases_seen':bases_seen,'sequences_seen':sequences_seen}
            torch.save(checkpoint,path/'checkpoint.partial')
            (path/'checkpoint.partial').replace(path/'checkpoint.pt')
            write_json(path/'metrics.json',history)
    elapsed += time.monotonic()-started
    summary = {'mode':mode,'seed':seed,'updates':end_update,'complete':end_update == config['updates'],
        'scope':'TRAINING ONLY; no validation, test, gap completion, or structure scoring',
        'initial':history[0],'final':history[-1],'elapsed_seconds':elapsed,
        'parameters':sum(p.numel() for p in model.parameters()),'bases_seen':bases_seen,
        'sequences_seen':sequences_seen,'sequence_exposures_per_unique_cds':sequences_seen/len(corpus),
        'cohort_sha256':corpus.metadata['cohort_sha256']}
    write_json(path/'summary.json',summary)
    return summary
