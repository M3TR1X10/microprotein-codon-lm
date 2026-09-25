"""Compile completed pilot outputs; never select a best model."""
import csv
import json
import statistics
from collections import Counter
from pathlib import Path

import numpy as np

from microprotein_lm.io import read_json,write_json,sha256
from microprotein_lm.tokenization import Tokenizer

raw = Path('data/raw')
sources = []
for path in sorted(raw.rglob('*.source.json')):
    receipt = read_json(path)
    data_path = Path(str(path).removesuffix('.source.json'))
    if sha256(data_path) != receipt['sha256']:
        raise ValueError(f'Source file has changed: {data_path}')
    sources.append({'file':data_path.relative_to(raw).as_posix(),**receipt})
write_json('reports/source-registry.json',{'sources':sources,
    'note':'Raw archives are required if the provider changes content at these URLs.',
    'protocol_inputs':{'configs/product_aliases.json':sha256('configs/product_aliases.json'),
        'configs/pilot.json':sha256('configs/pilot.json')}})

records = [json.loads(l) for l in Path('data/processed/cohort.jsonl').read_text().splitlines()]
matrix = np.array([list(r['rna']) for r in records])
distances = [(matrix[i] != matrix[j]).sum()/matrix.shape[1]
    for i in range(len(matrix)) for j in range(i)]
diversity = {'sequences':len(records),'nucleotides_per_sequence':matrix.shape[1],
    'variable_positions':int(sum(len(set(matrix[:,i])) > 1 for i in range(matrix.shape[1]))),
    'pairwise_nucleotide_identity_mean':float(1-np.mean(distances)),
    'pairwise_nucleotide_identity_min':float(1-np.max(distances)),
    'exact_reference_protein_sequences':sum(r['reference_identity'] == 1 for r in records),
    'min_reference_protein_identity':min(r['reference_identity'] for r in records)}
write_json('reports/cohort-diversity.json',diversity)
with open('reports/codon-coverage.csv','w',newline='',encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['token_id','codon','total_occurrences','target_occurrences'])
    t = Tokenizer('codon')
    counts = Counter(c for r in records for c in t.encode(r['rna']))
    targets = Counter(c for r in records for c in t.encode(r['rna'])[1:])
    for i,token in enumerate(t.vocabulary):
        writer.writerow([i,token,counts[i],targets[i]])

results = read_json('reports/pilot-results.json')
history = []
for r in results:
    folder = Path('runs/pilot')/f"{r['mode']}-seed{r['seed']}"
    history.append({'mode':r['mode'],'seed':r['seed'],'metrics':read_json(folder/'metrics.json'),
        'run':read_json(folder/'run.json')})
write_json('reports/pilot-runs.json',history)
lines = ['# Pilot results — training only','',
    'All six prespecified runs completed 200 optimizer updates on the same 129 human MT-ATP8 coding sequences. '
    'Every run reduced training loss from initialization. These are training-set diagnostics, not unseen-sequence accuracy.','',
    '| Arm | Seed | Parameters | Initial bits/base | Final bits/base | Seconds |',
    '|---|---:|---:|---:|---:|---:|']
for r in results:
    lines.append(f"| {r['mode']} | {r['seed']} | {r['parameters']:,} | {r['initial']['training_bits_per_base']:.4f} | {r['final']['training_bits_per_base']:.4f} | {r['elapsed_seconds']:.1f} |")
lines += ['', 'Each run saw 3,200 sequence presentations (24.81 presentations per unique CDS) and 652,800 target nucleotides. '
    'The two arms used the same biological examples, seed-specific sampling order and update budget. Attention cost and output-class count differ.','',
    '![Training curves](training-curves.svg)','',
    'The curves show each seed separately. They are not uncertainty intervals over biological populations. '
    'Timing includes training-corpus diagnostics after initialization and checkpoint writing.','',
    '## Interpretation','']
for mode in ('base','codon'):
    values = [r['final']['training_bits_per_base'] for r in results if r['mode']==mode]
    lines.append(f"- {mode}: final training bits/base mean **{statistics.mean(values):.4f}**, sample SD **{statistics.stdev(values):.4f}** across three initialization seeds.")
baselines = read_json('reports/training-baselines.json')
lines += ['', 'The codon arm fit the training corpus more quickly under this fixed update/exposure budget. '
    'This is an optimization observation, not evidence that codon tokenization predicts new biological sequences better.','',
    '| Training-fit baseline | Base bits/base | Codon bits/base |','|---|---:|---:|']
for name in ['uniform','unigram','bigram','position']:
    b = baselines['arms']['base']['metrics'][name]['bits_per_base']
    c = baselines['arms']['codon']['metrics'][name]['bits_per_base']
    lines.append(f'| {name} | {b:.4f} | {c:.4f} |')
lines += ['', '**Neither transformer surpassed its corresponding smoothed positional baseline.** '
    'A close single-gene cohort is largely predictable by nucleotide/codon position; transformer loss reduction alone is weak evidence of useful biological modeling. '
    'No hyperparameter winner, functional claim, gap-completion accuracy or structural accuracy is established.','',
    f"The cohort has **{diversity['variable_positions']} variable nucleotide positions** and **{diversity['pairwise_nucleotide_identity_mean']:.2%} mean pairwise nucleotide identity**. "
    f"{diversity['exact_reference_protein_sequences']} sequences encode exactly the reviewed reference protein. "
    'These observations explain the strong position baseline and limit the effective diversity of 129 distinct CDS.','',
    'All 64 codon classes are represented in the encoder and model head; `codon-coverage.csv` records which classes actually occurred. '
    'Unobserved classes have no positive training examples. Vocabulary coverage must not be confused with observed-data coverage.','',
    '## Artifacts and checks','',
    '- `pilot-results.json`: per-run initial/final metrics, exposures and timings.',
    '- `pilot-runs.json`: full loss histories, executable settings, source hashes and environment.',
    '- `source-registry.json` and `data-manifest.json`: input retrieval and cohort identity.',
    '- Local `runs/pilot/<arm>-seed<seed>/checkpoint.pt`: model, optimizer and RNG states.',
    '- Offline software tests verify vocabulary, biological QC, boundaries, causality, corruption rejection and exact CPU resume.',
    '', 'The pilot ran on CPU with Python 3.12.14 and PyTorch 2.14.0. Validation, final testing, gap filling, AlphaFold and PyMOL were not run.']
Path('reports/pilot-report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

# Static research figure: each seed visible; no implied biological confidence interval.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
fig,axes = plt.subplots(1,2,figsize=(10,3.8),sharey=True,layout='constrained')
colors = ['#215a8a','#ce6b29','#48856a']
for ax,mode in zip(axes,['base','codon']):
    for color,item in zip(colors,[h for h in history if h['mode']==mode]):
        ax.plot([m['update'] for m in item['metrics']],
            [m['training_bits_per_base'] for m in item['metrics']],
            marker='o',markersize=3,color=color,label=f"Seed {item['seed']}")
    pos = baselines['arms'][mode]['metrics']['position']['bits_per_base']
    ax.axhline(pos,color='#666666',linestyle='--',linewidth=1,label='Position baseline (train fit)')
    ax.set_title('4-base tokens' if mode=='base' else '64-codon tokens')
    ax.set_xlabel('Optimizer updates')
    ax.set_ylim(0,2.25)
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
axes[0].set_ylabel('Training negative log likelihood (bits/base)')
fig.suptitle('Human MT-ATP8: training feasibility only',fontsize=13)
fig.savefig('reports/training-curves.svg')
fig.savefig('reports/training-curves.png',dpi=180)
print('Wrote source registry, coverage, diversity, pilot report and training curves')
