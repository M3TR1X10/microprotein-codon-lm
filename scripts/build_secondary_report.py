"""Report frozen secondary training runs; no training, selection or held-out scoring.

Run after scripts/run_secondary.py. Incomplete runs stay explicitly incomplete.
Counts are read from final manifests, never hardcoded from an acquisition preview.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from microprotein_lm.io import read_json, sha256
from microprotein_lm.secondary_io import write_json
from microprotein_lm.tokenization import Tokenizer

SCOPE = 'TRAINING ONLY; all diagnostics are resubstitution, with no held-out observations.'
PRIMARY = ['human_atp8_base', 'human_atp8_codon', 'human_complex_codon',
           'mammal_complex_codon', 'vertebrate_complex_codon']
CONTROLS = ['vertebrate_small', 'vertebrate_family_macro', 'vertebrate_no_position']
LABELS = {
    'human_atp8_base': 'Human ATP8 / bases',
    'human_atp8_codon': 'Human ATP8 / codons',
    'human_complex_codon': 'Human Complex V / codons',
    'mammal_complex_codon': 'Mammal Complex V / codons',
    'vertebrate_complex_codon': 'Vertebrate Complex V / codons',
    'vertebrate_small': 'Vertebrate / small capacity',
    'vertebrate_family_macro': 'Vertebrate / family-macro objective',
    'vertebrate_no_position': 'Vertebrate / no explicit positions',
    'full_vertebrate_capacity': 'Full pool / larger capacity',
}


def number(value, places=4):
    if value is None:
        return 'n/a'
    if not math.isfinite(value):
        raise ValueError('Nonfinite number in scientific report')
    return f'{value:,.{places}f}'


def mean_sd(values, places=4):
    values = [float(x) for x in values if x is not None]
    if not values:
        return 'n/a'
    if len(values) == 1:
        return f'{number(values[0], places)} (n=1; SD n/a)'
    return f'{number(statistics.mean(values), places)} ± {number(statistics.stdev(values), places)}'


def support(values):
    values = sorted({int(x) for x in values if x is not None})
    return ', '.join(f'{x:,}' for x in values) if values else 'n/a'


def linear_diversity(records):
    """O(total sequence length), avoiding all-pairs biological comparisons."""
    vocabulary = Tokenizer('codon').vocabulary
    total, targets = Counter(), Counter()
    alleles = defaultdict(lambda: defaultdict(Counter))
    groups, peptides = Counter(), defaultdict(set)
    sequence_ids, all_peptides = set(), set()
    for row in records:
        rna = row['rna']
        if hashlib.sha256(rna.encode()).hexdigest() != row['sequence_sha256']:
            raise ValueError('A reported sequence hash does not match its RNA')
        if row['sequence_sha256'] in sequence_ids:
            raise ValueError('Duplicate RNA in a reported cohort')
        sequence_ids.add(row['sequence_sha256'])
        if len(rna) % 3 or set(rna) - set('AUCG'):
            raise ValueError('Invalid RNA reading frame/alphabet')
        codons = [rna[i:i+3] for i in range(0, len(rna), 3)]
        total.update(codons)
        targets.update(codons[1:])
        group = f"{row['family']}:{row['tax_id']}"
        groups[group] += 1
        peptides[group].add(row['protein'])
        all_peptides.add(row['protein'])
        for pos, codon in enumerate(codons[1:], 1):
            alleles[group][pos][codon] += 1
    group_metrics = {}
    for group in sorted(groups):
        positions = alleles[group]
        variable = [counts for counts in positions.values() if len(counts) > 1]
        group_metrics[group] = {
            'n_sequences': groups[group], 'n_unique_peptides': len(peptides[group]),
            'target_codon_sites': len(positions), 'variable_codon_sites': len(variable),
            'variable_codon_occurrences': sum(sum(counts.values()) for counts in variable),
            'variable_target_bases': 3 * sum(sum(counts.values()) for counts in variable),
            'nucleotide_sites_variable_within_target_codons': sum(
                len({codon[i] for codon in counts}) > 1
                for counts in positions.values() for i in range(3)),
        }
    coverage = [{'token_id': i, 'codon': codon, 'total_occurrences': total[codon],
                 'target_occurrences': targets[codon]} for i, codon in enumerate(vocabulary)]
    return {
        'n_sequences': len(records), 'n_unique_peptides': len(all_peptides),
        'observed_total_codon_classes': sum(total[c] > 0 for c in vocabulary),
        'observed_target_codon_classes': sum(targets[c] > 0 for c in vocabulary),
        'total_nt': sum(len(r['rna']) for r in records),
        'target_nt': sum(len(r['rna'])-3 for r in records),
        'variable_codon_sites_summed_over_family_taxon': sum(g['variable_codon_sites'] for g in group_metrics.values()),
        'variable_target_bases': sum(g['variable_target_bases'] for g in group_metrics.values()),
        'family_taxon': group_metrics, 'codons': coverage,
    }


def source_registry(root):
    raw, sources = root/'data/raw', []
    for path in sorted(raw.rglob('*.source.json')):
        receipt = read_json(path)
        data_path = Path(str(path).removesuffix('.source.json'))
        actual = sha256(data_path)
        if actual != receipt['sha256']:
            raise ValueError(f'Source checksum mismatch: {data_path}')
        sources.append({**receipt, 'file': data_path.relative_to(raw).as_posix(),
                        'receipt_file': path.relative_to(raw).as_posix(),
                        'receipt_sha256': sha256(path), 'actual_sha256': actual})
    return {'scope': SCOPE, 'path_base': 'data/raw', 'sources': sources,
            'coverage': 'All receipt-bearing raw files present at report generation, including prior-pilot inputs; individual sequence provenance identifies used files.',
            'note': 'URLs may change. Preserve the raw archives; hashes alone cannot reconstruct a changed provider response.'}


def load_runs(root, plan, results, views):
    jobs = {(j['arm'], j['seed']): j for j in plan['jobs']}
    seen, runs, baselines = set(), [], {}
    for summary in results:
        key = summary['arm'], summary['seed']
        if key in seen or key not in jobs:
            raise ValueError(f'Duplicate or unplanned result: {key}')
        seen.add(key)
        job = jobs[key]
        if summary['config'] != job['config'] or summary['source_hashes'] != plan['engine_source_hashes']:
            raise ValueError(f'Result/frozen engine identity mismatch: {key}')
        if summary['data_identity']['manifest_sha256'] != job['manifest_sha256']:
            raise ValueError(f'Result/frozen manifest identity mismatch: {key}')
        if summary['role'] != job['role'] or summary['dataset'] != job['dataset'] or summary['mode'] != job['mode']:
            raise ValueError(f'Result/plan identity mismatch: {key}')
        if summary['cohort_sha256'] != job['cohort_sha256'] or summary['cohort_sha256'] != views[job['dataset']]['cohort_sha256']:
            raise ValueError(f'Result/cohort identity mismatch: {key}')
        directory = root/'runs/secondary'/f"{key[0]}-seed{key[1]}"
        metadata, metrics, baseline = [read_json(directory/f) for f in ('run.json', 'metrics.json', 'baselines.json')]
        if any(metadata[k] != summary[k] for k in ('config', 'mode', 'seed')):
            raise ValueError(f'Run configuration mismatch: {key}')
        if metadata['data_identity'] != summary['data_identity'] or metadata['source_hashes'] != summary['source_hashes']:
            raise ValueError(f'Run metadata mismatch: {key}')
        if baseline['cache_identity']['data'] != summary['data_identity'] or baseline['mode'] != summary['mode']:
            raise ValueError(f'Baseline/cohort identity mismatch: {key}')
        if baseline['cache_identity']['source_hashes'] != plan['engine_source_hashes']:
            raise ValueError(f'Baseline/frozen engine identity mismatch: {key}')
        if metrics[0] != summary['initial'] or metrics[-1] != summary['final']:
            raise ValueError(f'Metrics/summary mismatch: {key}')
        if summary['complete'] and summary['updates'] != job['config']['updates']:
            raise ValueError(f'Completed run has wrong update budget: {key}')
        baseline_key = f"{job['dataset']}:{job['mode']}"
        if baseline_key in baselines and baselines[baseline_key] != baseline:
            raise ValueError(f'Different baseline for same cohort/mode: {baseline_key}')
        baselines[baseline_key] = baseline
        runs.append({'arm': job['arm'], 'role': job['role'],
                     'run_directory': directory.relative_to(root).as_posix(),
                     'metadata': metadata, 'summary': summary, 'metrics': metrics,
                     'baseline_key': baseline_key})
    return runs, baselines


def create_curves(root, runs, baselines):
    # Standard plotting for the exported scientific figure; no image generator.
    os.environ.setdefault('MPLCONFIGDIR', str(root/'work/matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'svg.fonttype': 'none'})
    colors = {17: '#17618f', 29: '#bc632b', 43: '#347d56'}
    arms = PRIMARY + CONTROLS + ['full_vertebrate_capacity']
    fig, axes = plt.subplots(3, 3, figsize=(15, 12), layout='constrained')
    for ax, arm in zip(axes.flat, arms):
        chosen = [r for r in runs if r['arm'] == arm]
        for run in sorted(chosen, key=lambda r: r['summary']['seed']):
            s, history = run['summary'], run['metrics']
            ax.plot([m['update'] for m in history], [m['training_bits_per_base'] for m in history],
                    color=colors.get(s['seed'], '#584996'), linewidth=1.7, marker='o', markersize=3,
                    label=f"Seed {s['seed']}" + ('' if s['complete'] else ' (partial)'))
        if chosen:
            baseline = baselines[chosen[0]['baseline_key']]['metrics']
            for style, label, value in [(':', 'Empirical pooled position', baseline['empirical']['position']['training_bits_per_base']),
                ('--', 'Family-position, prior mass 2', baseline['equal_total_prior_2']['family_position']['training_bits_per_base'])]:
                ax.axhline(value, color='#686868', linestyle=style, linewidth=1, label=label)
            ax.legend(fontsize=7, loc='best')
        else:
            ax.text(.5, .5, 'No recorded run', transform=ax.transAxes, ha='center', va='center')
        ax.set_title(LABELS[arm], fontsize=10)
        ax.set_xlabel('Optimizer updates')
        ax.set_ylabel('Training NLL (bits/base)')
        ax.set_ylim(bottom=0)
        ax.grid(alpha=.18)
    fig.suptitle('Secondary pilot: training fit only; each curve is one initialization\nFull-pool panel changes both dataset and capacity', fontsize=14)
    fig.savefig(root/'reports/secondary/curves.svg')
    fig.savefig(root/'reports/secondary/curves.png', dpi=180)
    plt.close(fig)


def data_card(root, cohorts, diversity):
    biology = read_json(root/'configs/secondary/biology.json')
    lines = ['# Secondary dataset card', '',
        'This snapshot contains observed ATP-synthase coding sequences for training only. '
        'Counts below are generated from frozen cohort files and checked against their manifests; '
        'they are not an estimate of independent individuals or all available vertebrate sequences.', '',
        '## Scope and sources', '',
        'The declared panel and exact family identifiers are in [biology.json](../configs/secondary/biology.json). '
        'IMPI establishes human mitochondrial eligibility; UniProt supplies reviewed family/species references; '
        'ENA supplies observed nucleotide records. The [acquisition report](../reports/secondary/acquisition.json) '
        'records the CDS-index/cross-reference routes, exclusions and coverage. '
        'When present, the [genome-acquisition report](../reports/secondary/genome-acquisition.json) separately records '
        'retrieval of mitochondrial parent records in the declared 14,000–20,000-nt range and extraction of explicitly annotated ATP8 CDS. '
        'Parent length does not establish whole-genome completeness; the extracted CDS must pass the complete-CDS checks. '
        'The original acquisition snapshot is preserved; the final cohort manifests identify which enriched pool was selected.', '',
        'The families share ATP synthase membership; they are not all ATP8 homologs. '
        'Cross-species reference coverage is incomplete and availability-based. '
        'Nonmammalian coverage and actual family-by-taxon counts must be read from this snapshot, not inferred from the species panel.', '',
        '## Quality and sequence identity', '',
        f"The operational maximum is {biology['max_amino_acids']} amino acids excluding stop, with "
        f"at least {biology['min_reference_identity']:.0%} amino-acid identity to the sequence's own-species reviewed reference and the reference length. "
        'Only complete, unambiguous CDS with allowed start, full stop, matching annotated translation and no internal stop or translation exceptions are eligible. '
        'Mitochondrial ATP8 uses genetic code 2; nuclear families use code 1. T is transcribed to U after QC. '
        'No reverse translation, invented variants or repaired bases are added.', '',
        'Identical RNA is globally deduplicated. Original accessions, versions, source files, hashes and taxonomic observations remain in per-sequence provenance in the enriched raw pool, joined from compact training views by sequence SHA-256. '
        'A primary selection stratum represents each shared sequence; taxon-stratum counts do not erase its other observed taxa. '
        'Annotation and high sequence identity do not experimentally establish equivalent molecular function.', '',
        '## Frozen cohort counts', '',
        '| View | Distinct CDS | Distinct peptides | Families | Taxon strata | Target bases | Observed target codons / 64 |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for name, manifest in cohorts['views'].items():
        d = diversity[name]
        lines.append(f"| {name} | {manifest['n_sequences']:,} | {d['n_unique_peptides']:,} | "
            f"{len(manifest['family_counts'])} | {len(manifest['taxon_counts'])} | {manifest['target_nt']:,} | {d['observed_target_codon_classes']} |")
    matched_n = cohorts['views']['human_atp8']['n_sequences']
    family_replacements = matched_n - cohorts['overlap_distinct_rna']['human_atp8__human_complex']
    vertebrate_replacements = matched_n - cohorts['overlap_distinct_rna']['mammal_complex__vertebrate_complex']
    lines += ['', f'The matched family expansion replaces **{family_replacements}/{matched_n} CDS ({family_replacements/matched_n:.2%})**. '
        f'The nonmammalian expansion replaces **{vertebrate_replacements}/{matched_n} CDS ({vertebrate_replacements/matched_n:.2%})** in the mammal view. '
        'These small interventions limit what aggregate training differences can reveal about broader family or taxonomic effects.', '',
        'The base and codon human-ATP8 arms share exactly one cohort. Matched Complex V views share family quotas and total CDS count. '
        'Equal sequence count does not mean equal nucleotide exposure, independent biological sampling, or identical effective diversity. '
        'The full pools retain all eligible distinct CDS recovered within the declared acquisition scope.', '',
        '## Family and taxon composition', '',
        '| View | Family | Primary taxon stratum | Distinct CDS |', '|---|---|---|---:|']
    for name, manifest in cohorts['views'].items():
        for group, count in manifest['family_taxon_counts'].items():
            family, tax = group.rsplit(':', 1)
            lines.append(f"| {name} | {family} | {biology['taxa'].get(tax, tax)} ({tax}) | {count:,} |")
    lines += ['', '## Intended interpretation and limitations', '',
        'The entire snapshot is available to training. It must not later be relabeled as an untouched validation/test set. '
        'The matched views replace records and overlap; they are not independent cohorts. '
        'Quotas may retain sparse families only in the human stratum, so adding organisms does not expand every family equally.', '',
        'Variable-codon positions are identified separately within each family and primary taxon stratum. '
        'The base arm scores all three bases of those same variable codons. Singleton or invariant strata contribute no variable-position targets. '
        'Different cohorts have different variable-position supports; those metrics are descriptive across cohorts.', '',
        'All 64 RNA codons are vocabulary classes even when observed counts are zero. '
        'The [coverage/diversity artifact](../reports/secondary/codon-coverage.json) reports every class, source of positional variation and its support. '
        'No pairwise-identity estimate or effective-sample-size estimate is computed.', '',
        'See the [protocol](secondary-protocol.md), [cohort manifests](../reports/secondary/cohorts.json), '
        '[retrieval registry](../reports/secondary/source-registry.json) and [training report](../reports/secondary/report.md). '
        'Raw archives and checkpoints remain outside Git tracking.']
    (root/'docs/secondary-data-card.md').write_bytes(('\n'.join(lines)+'\n').encode('utf-8'))


def build_markdown(root, plan, cohorts, runs, baselines, diversity):
    results = [r['summary'] for r in runs]
    completed = [s for s in results if s['complete']]
    groups = {arm: [s for s in completed if s['arm'] == arm] for arm in PRIMARY + CONTROLS}
    codon_pairs = {s['seed']: s for s in groups['human_atp8_codon']}
    for base in groups['human_atp8_base']:
        if base['seed'] in codon_pairs:
            codon = codon_pairs[base['seed']]
            for key in ('cohort_sha256', 'sampling_trace_sha256', 'bases_seen', 'sequences_seen', 'updates'):
                if base[key] != codon[key]:
                    raise ValueError(f'Base/codon exposure mismatch in {key}, seed {base["seed"]}')
    expected = Counter(j['arm'] for j in plan['jobs'])
    execution_path = root/'reports/secondary/execution.json'
    execution = read_json(execution_path) if execution_path.exists() else {}
    lines = ['# Secondary pilot — training only', '', SCOPE, '',
        f"**{len(completed)} of {len(plan['jobs'])} planned runs completed their fixed update budgets.** "
        f"{len(results)-len(completed)} recorded runs are partial; {len(plan['jobs'])-len(results)} have no reported result. "
        'Partial runs are excluded from completed-run averages. A pending seed is not a failed scientific hypothesis.', '',
        f"The aggregate execution budget is {plan['aggregate_budget_seconds']/3600:g} hours. "
        + (f"Recorded scheduler wall time is {execution['elapsed_seconds']/60:.1f} minutes. " if 'elapsed_seconds' in execution else '')
        + 'Run-loop timings below include interval diagnostics and checkpoints but exclude setup/initial diagnostics; scheduler time covers the broader execution.', '',
        '## Data and controlled comparisons', '',
        '| Frozen view | Distinct CDS | Distinct peptides | Families | Taxon strata | Target bases/pass | Variable target bases | Target codons observed / 64 |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for name, manifest in cohorts['views'].items():
        d = diversity[name]
        lines.append(f"| {name} | {manifest['n_sequences']:,} | {d['n_unique_peptides']:,} | {len(manifest['family_counts'])} | "
                     f"{len(manifest['taxon_counts'])} | {manifest['target_nt']:,} | {d['variable_target_bases']:,} | {d['observed_target_codon_classes']} |")
    lines += ['', 'The five primary arms vary tokenization or cohort composition. The three additional controls change capacity, '
        'family-macro loss weighting, or explicit positional encoding on exactly the same vertebrate cohort. '
        'Uniform sequence draws and matched update budgets control presentations; natural lengths make nucleotide exposure differ. '
        'The larger full-pool run changes both capacity and data and is an engineering demonstration.', '',
        '## Completed matched runs', '',
        'Values are mean ± sample standard deviation across completed initialization seeds, not confidence intervals or biological replicates. '
        'The completed/planned column exposes missing seeds. No model is selected by minimum training loss.', '',
        '| Arm | Completed/planned | Parameters | Updates | Initial bits/base | Final bits/base | Final family-macro bits/base | Final variable-codon bits/base | Variable target bases |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for arm in PRIMARY + CONTROLS:
        g = groups[arm]
        lines.append(f"| {LABELS[arm]} | {len(g)}/{expected[arm]} | {support(s['parameters'] for s in g)} | "
            f"{support(s['updates'] for s in g)} | {mean_sd(s['initial']['training_bits_per_base'] for s in g)} | "
            f"{mean_sd(s['final']['training_bits_per_base'] for s in g)} | {mean_sd(s['final']['family_macro_bits_per_base'] for s in g)} | "
            f"{mean_sd(s['final']['variable_positions']['training_bits_per_base'] for s in g)} | "
            f"{support(s['final']['variable_positions']['scored_training_bases'] for s in g)} |")
    lines += ['', '![Training curves](curves.svg)', '',
        'Every plotted curve represents one initialization. Dashed family-position baselines receive privileged family labels. '
        'There is no held-out curve. Base/codon token accuracy and raw token perplexity are not compared.', '',
        'Variable-codon masks are cohort-specific. The base/codon pair shares the same biological support; '
        'the three vertebrate controls share their reference arm\'s support. Across different cohorts, changing support prevents a paired accuracy interpretation.', '',
        '## One-factor contrasts on the vertebrate cohort', '',
        'Differences are control minus the standard vertebrate codon arm, paired by initialization seed. '
        'Positive values mean higher training negative log likelihood for that metric. These are descriptive optimization contrasts, not significance tests.', '',
        '| Control | Complete seed pairs | Overall delta bits/base | Family-macro delta bits/base | Variable-codon delta bits/base |',
        '|---|---:|---:|---:|---:|']
    reference = {s['seed']: s for s in groups['vertebrate_complex_codon']}
    for arm in CONTROLS:
        pairs = [(s, reference[s['seed']]) for s in groups[arm] if s['seed'] in reference]
        overall, macro, variable = [], [], []
        for control, ref in pairs:
            if control['sampling_trace_sha256'] != ref['sampling_trace_sha256'] or control['cohort_sha256'] != ref['cohort_sha256']:
                raise ValueError(f'Unmatched sequence exposure in control {arm}')
            overall.append(control['final']['training_bits_per_base']-ref['final']['training_bits_per_base'])
            macro.append(control['final']['family_macro_bits_per_base']-ref['final']['family_macro_bits_per_base'])
            a, b = control['final']['variable_positions'], ref['final']['variable_positions']
            if a['training_bits_per_base'] is not None and b['training_bits_per_base'] is not None:
                variable.append(a['training_bits_per_base']-b['training_bits_per_base'])
        lines.append(f"| {LABELS[arm]} | {len(pairs)} | {mean_sd(overall)} | {mean_sd(macro)} | {mean_sd(variable)} |")
    lines += ['', '## Training-fit baselines', '',
        'Each distribution is fit and scored on the same training observations. Empirical rows use no smoothing. '
        'Smoothed rows use total prior mass 2 per categorical distribution (2/4 per base category, 2/64 per codon category). '
        'Equal total prior removes the first pilot\'s prior-mass imbalance; different factorizations and conditional contexts still differ. '
        'The family-position baseline is privileged: the transformer receives no family or organism token. Uniform entropy is 2 bits/base.', '',
        '| Cohort / tokens | Baseline | Empirical bits/base | Prior-mass-2 bits/base | Prior-mass-2 family-macro bits/base |',
        '|---|---|---:|---:|---:|']
    for key, baseline in sorted(baselines.items()):
        for name in ('unigram', 'bigram', 'position', 'family_position'):
            empirical = baseline['metrics']['empirical'][name]
            smooth = baseline['metrics']['equal_total_prior_2'][name]
            lines.append(f"| {key} | {name}{' (privileged)' if name == 'family_position' else ''} | "
                         f"{number(empirical['training_bits_per_base'])} | {number(smooth['training_bits_per_base'])} | {number(smooth['family_macro_bits_per_base'])} |")
    lines += ['', '## Per-family training fit', '',
        'Completed matched runs only; support is the number of target bases per full training-corpus pass, not multiplied by initialization seeds.', '',
        '| Arm | Family | Initializations | Target bases | Final bits/base | Variable target bases | Final variable-codon bits/base |',
        '|---|---|---:|---:|---:|---:|---:|']
    for arm in PRIMARY + CONTROLS:
        g = groups[arm]
        for family in sorted({f for s in g for f in s['final']['per_family']}):
            metrics = [s['final']['per_family'][family] for s in g]
            variables = [s['final']['variable_positions']['per_family'][family] for s in g]
            lines.append(f"| {LABELS[arm]} | {family} | {len(metrics)} | {support(m['scored_training_bases'] for m in metrics)} | "
                f"{mean_sd(m['training_bits_per_base'] for m in metrics)} | {support(m['scored_training_bases'] for m in variables)} | "
                f"{mean_sd(m['training_bits_per_base'] for m in variables)} |")
    lines += ['', '## Full-pool capacity demonstration', '',
        'This separate run uses the full recovered vertebrate pool and a larger codon model. '
        'Both factors change, so its fit is excluded from matched contrasts; one initialization has no across-seed SD.', '',
        '| Arm | State | CDS | Parameters | Updates | Initial bits/base | Last bits/base | Last family-macro bits/base |',
        '|---|---|---:|---:|---:|---:|---:|---:|']
    engineering = [s for s in results if s['role'] != 'matched']
    for s in engineering:
        lines.append(f"| {s['arm']} | {'Complete' if s['complete'] else 'Partial'} | {s['sequences']:,} | {s['parameters']:,} | "
            f"{s['updates']:,} | {number(s['initial']['training_bits_per_base'])} | {number(s['final']['training_bits_per_base'])} | "
            f"{number(s['final']['family_macro_bits_per_base'])} |")
    if not engineering:
        lines.append('| Full-pool run | No reported result | — | — | — | — | — | — |')
    lines += ['', '## Exposure, runtime and memory', '',
        '| Arm | Seed | State | Updates | Sequence presentations | Presentations / unique CDS | Target bases presented | Loop minutes | Peak process MiB |',
        '|---|---:|---|---:|---:|---:|---:|---:|---:|']
    for s in results:
        rss = s.get('process_peak_rss_bytes')
        lines.append(f"| {s['arm']} | {s['seed']} | {'Complete' if s['complete'] else 'Partial: '+s['stop_reason']} | {s['updates']:,} | "
            f"{s['sequences_seen']:,} | {number(s['sequence_exposures_per_unique_cds'], 2)} | {s['bases_seen']:,} | "
            f"{number(s['elapsed_seconds']/60, 2)} | {number(rss/1024**2 if rss is not None else None, 1)} |")
    improvements = sum(s['final']['training_bits_per_base'] < s['initial']['training_bits_per_base'] for s in completed)
    lines += ['', f'{improvements} of {len(completed)} completed runs reduced overall training bits/base from initialization. '
        'This is an optimization check; it does not establish unseen-sequence prediction.', '',
        '## Interpretation limits and next decisions', '',
        '- The first pilot used a different acquisition snapshot, model capacity, training budget and smoothing rule. '
        'Changes from its reported loss cannot be attributed to one factor; the new matched controls provide the interpretable within-round contrasts.',
        '- More varied cohorts can have higher inherent entropy. Higher training loss across cohorts does not by itself mean poorer modeling. '
        'Compare baseline gaps, per-family fit and support before interpreting aggregate changes.',
        '- Sparse families may have one or a few distinct CDS; a family-macro objective gives those sequences substantial weight and can encourage memorization. '
        'Zero variable-position support yields n/a, not zero error.',
        '- The species panel and reviewed references are availability-based. Missing family/species combinations and global sequence sharing constrain taxonomic interpretation. '
        'Distinct CDS and initialization seeds are not independent biological replicates.',
        '- No explicit-position control removes the positional encoding only; causal order and prefix length remain available.',
        '- Capacity changes width, depth and head count together. It is a bundled architecture comparison, not an isolated estimate of any one dimension.',
        '- Maximum data means the eligible recovered pool under declared routes and QC; maximum model means the selected profiled candidate under the desktop reserve and time budget. '
        'Neither claim describes an absolute global maximum.',
        '- A future predictive decision requires a separately designed untouched evaluation cohort. '
        'This round performs no validation, testing, gap completion, structural prediction, or molecular-function assessment.', '',
        '## Reproducibility artifacts', '',
        '- [Frozen protocol](../../docs/secondary-protocol.md), [data card](../../docs/secondary-data-card.md), [plan](plan.json), '
        '[cohort manifests and overlaps](cohorts.json), [acquisition](acquisition.json).',
        '- [Per-run results](results.json), [histories, metadata and baselines](runs.json), [all 64 codon counts and diversity supports](codon-coverage.json), '
        '[source receipts and hashes](source-registry.json), [hardware profile](hardware-benchmark.json).',
        '- [PNG figure](curves.png) and [SVG figure](curves.svg). Raw archives and local rolling/final checkpoints remain outside Git tracking.',
        '- Software check outcomes are recorded separately by the project; generating this report does not certify an unexecuted check.']
    (root/'reports/secondary/report.md').write_bytes(('\n'.join(lines)+'\n').encode('utf-8'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default=str(ROOT))
    parser.add_argument('--data-only', action='store_true', help='Export audited data artifacts before any model fitting')
    args = parser.parse_args()
    root = Path(args.root).resolve()
    output = root/'reports/secondary'
    plan, cohorts = read_json(output/'plan.json'), read_json(output/'cohorts.json')
    results = read_json(output/'results.json') if (output/'results.json').exists() else []
    if plan['scope'] != 'training_only' or cohorts['stage'] != 'training_only':
        raise ValueError('Report requires training-only plan/cohorts')
    for path, checksum in plan['input_hashes'].items():
        if sha256(root/path) != checksum:
            raise ValueError(f'Frozen plan input changed before reporting: {path}')
    diversity = {}
    for name, manifest in cohorts['views'].items():
        path = root/'data/processed/secondary'/name/'cohort.jsonl'
        if sha256(path) != manifest['cohort_sha256']:
            raise ValueError(f'Frozen cohort changed: {name}')
        with path.open(encoding='utf-8') as handle:
            d = linear_diversity([json.loads(line) for line in handle if line.strip()])
        if d['n_sequences'] != manifest['n_sequences'] or d['target_nt'] != manifest['target_nt']:
            raise ValueError(f'Manifest counts disagree: {name}')
        diversity[name] = {'cohort_sha256': manifest['cohort_sha256'], **d}
    stamp = datetime.now(timezone.utc).isoformat()
    acquisition_reports = {}
    for name in ('acquisition.json', 'genome-acquisition.json'):
        path = output/name
        if path.exists():
            acquisition = read_json(path)
            acquisition_reports[name] = {'file': path.relative_to(root).as_posix(),
                'sha256': sha256(path),
                'top_level_scalar_fields': {k: v for k, v in acquisition.items()
                    if v is None or isinstance(v, (str, int, float, bool))},
                'note': 'The linked report retains full inventories, exclusions and nested counts.'}
    write_json(output/'source-registry.json', {**source_registry(root), 'acquisition_reports': acquisition_reports})
    write_json(output/'codon-coverage.json', {'scope': SCOPE, 'generated_utc': stamp,
        'vocabulary': Tokenizer('codon').vocabulary,
        'definition': 'All64 codon classes including zeros. Target counts omit first codon, include full stop. Variation is within family and primary taxon stratum; all3 bases of variable codons count. No pairwise identity computation.',
        'cohorts': diversity})
    data_card(root, cohorts, diversity)
    if args.data_only:
        print('Wrote audited data card, source registry and codon coverage before model fitting.')
        return
    runs, baselines = load_runs(root, plan, results, cohorts['views'])
    write_json(output/'runs.json', {'scope': SCOPE, 'generated_utc': stamp,
        'plan_sha256': sha256(output/'plan.json'), 'report_builder_sha256': sha256(Path(__file__)),
        'runs': runs, 'baselines': baselines})
    create_curves(root, runs, baselines)
    build_markdown(root, plan, cohorts, runs, baselines, diversity)
    print(f'Wrote secondary report, plots, run histories, sources, codon coverage and data card for {len(runs)} recorded runs.')


if __name__ == '__main__':
    main()
