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


def verify_paired_exposure(control, reference, label):
    for key in ('cohort_sha256', 'sampling_trace_sha256', 'updates', 'sequences_seen', 'bases_seen'):
        if control[key] != reference[key]:
            raise ValueError(f'Paired exposure mismatch in {key}: {label}')
    for stage in ('initial', 'final'):
        a, b = control[stage], reference[stage]
        if a['scored_training_bases'] != b['scored_training_bases']:
            raise ValueError(f'Paired scored-region mismatch: {label}')
        if a['variable_positions']['scored_training_bases'] != b['variable_positions']['scored_training_bases']:
            raise ValueError(f'Paired variable-position support mismatch: {label}')
        for category in ('per_family', 'per_taxon', 'per_code'):
            if ({k: v['scored_training_bases'] for k, v in a[category].items()} !=
                    {k: v['scored_training_bases'] for k, v in b[category].items()}):
                raise ValueError(f'Paired {category} support mismatch: {label}')
    if control['per_family_exposures'] != reference['per_family_exposures']:
        raise ValueError(f'Paired family exposure mismatch: {label}')


def paired_contrasts(groups, plan):
    """Only completed, seed-paired outcomes; directions come from the frozen protocol."""
    designs = [
        ('Tokenization: codon minus base', 'human_atp8_codon', 'human_atp8_base',
         'overall_delta', -1, 'Codon lower overall bits/base'),
        ('Capacity: small minus standard', 'vertebrate_small', 'vertebrate_complex_codon',
         'overall_delta', 1, 'Small architecture higher overall bits/base'),
        ('Family weighting: macro minus token', 'vertebrate_family_macro', 'vertebrate_complex_codon',
         'macro_delta', -1, 'Family-macro objective lower family-macro bits/base'),
        ('Explicit position: removed minus supplied', 'vertebrate_no_position', 'vertebrate_complex_codon',
         'overall_delta', 1, 'Removed positions higher overall bits/base'),
    ]
    contrasts = []
    for label, arm, ref_arm, metric, direction, expectation in designs:
        refs = {s['seed']: s for s in groups[ref_arm]}
        pairs = []
        for control in sorted(groups[arm], key=lambda s: s['seed']):
            if control['seed'] not in refs:
                continue
            ref = refs[control['seed']]
            verify_paired_exposure(control, ref, f'{label}, seed {control["seed"]}')
            a, b = control['final'], ref['final']
            av, bv = a['variable_positions'], b['variable_positions']
            variable = (av['training_bits_per_base']-bv['training_bits_per_base']
                        if av['training_bits_per_base'] is not None and bv['training_bits_per_base'] is not None else None)
            pairs.append({'seed': control['seed'],
                'overall_delta': a['training_bits_per_base']-b['training_bits_per_base'],
                'macro_delta': a['family_macro_bits_per_base']-b['family_macro_bits_per_base'],
                'variable_delta': variable, 'variable_target_bases': av['scored_training_bases'],
                'sequences_seen': control['sequences_seen'], 'bases_seen': control['bases_seen']})
        planned = len({j['seed'] for j in plan['jobs'] if j['arm'] == arm} &
                      {j['seed'] for j in plan['jobs'] if j['arm'] == ref_arm})
        contrasts.append({'label': label, 'pairs': pairs, 'planned': planned,
                          'expectation': expectation, 'metric': metric, 'direction': direction})
    return contrasts


def expectation_outcome(value, direction, complete=True):
    if value is None:
        return 'Pending'
    outcome = 'Met' if value*direction > 0 else 'Reversed' if value*direction < 0 else 'Unchanged'
    return outcome if complete else f'Provisional: {outcome.lower()} (missing seed pairs)'


def validate_execution(execution, plan_hash):
    if not execution:
        return
    if execution.get('plan_sha256') != plan_hash:
        raise ValueError('Execution accounting belongs to a different frozen plan')
    seconds = execution.get('elapsed_seconds')
    if not isinstance(seconds, (int, float)) or not math.isfinite(seconds) or seconds < 0:
        raise ValueError('Execution accounting contains invalid elapsed time')
    if execution.get('status') not in ('running', 'complete', 'budget_stopped', 'timed_out',
                                        'worker_failed', 'invalid_output'):
        raise ValueError('Unknown scheduler execution status')


def applied_timing_amendment(root, plan_hash):
    directory = root/'reports/secondary/interruption'
    receipt_path = directory/'application.json'
    if not receipt_path.exists():
        return None
    receipt = read_json(receipt_path)
    amendment = read_json(directory/'amendment.json')
    if (receipt['amendment_sha256'] != sha256(directory/'amendment.json')
            or amendment['plan_sha256'] != plan_hash
            or receipt['original_execution_sha256'] != sha256(directory/'original-execution.json')
            or amendment['evidence_sha256'] != sha256(directory/'power-evidence.json')):
        raise ValueError('Timing amendment evidence/plan identity mismatch')
    raw = amendment['original_elapsed_wall_seconds']
    excluded = amendment['excluded_hardware_deep_idle_seconds']
    charged = amendment['charged_elapsed_before_resume_seconds']
    if (not all(math.isfinite(v) for v in (raw, excluded, charged))
            or not 0 < excluded < raw or abs(raw-excluded-charged) > 1e-6):
        raise ValueError('Invalid timing amendment arithmetic')
    return amendment


def charged_loop_seconds(summary, amendment):
    raw = summary['elapsed_seconds']
    if amendment and (summary['arm'], summary['seed']) == (amendment['affected_arm'], amendment['affected_seed']):
        charged = raw-amendment['excluded_hardware_deep_idle_seconds']
        if charged < 0:
            raise ValueError('Interrupted run timing is smaller than its documented exclusion')
        return charged
    return raw


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
        if (type(summary['complete']) is not bool or type(summary['updates']) is not int
                or not 0 < summary['updates'] <= job['config']['updates']
                or summary['complete'] != (summary['updates'] == job['config']['updates'])
                or summary['final']['update'] != summary['updates']):
            raise ValueError(f'Run has inconsistent completion/update accounting: {key}')
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
    panels, plotted_values = {}, []
    for arm in arms:
        chosen = sorted((r for r in runs if r['arm'] == arm), key=lambda r: r['summary']['seed'])
        baseline_lines = []
        if chosen:
            baseline = baselines[chosen[0]['baseline_key']]['metrics']
            baseline_lines = [(':', 'Empirical pooled position', baseline['empirical']['position']['training_bits_per_base']),
                ('--', 'Family-position, prior mass 2', baseline['equal_total_prior_2']['family_position']['training_bits_per_base'])]
        panels[arm] = chosen, baseline_lines
        plotted_values.extend(m['training_bits_per_base'] for r in chosen for m in r['metrics'])
        plotted_values.extend(value for _, _, value in baseline_lines)
    if any(not math.isfinite(value) or value < 0 for value in plotted_values):
        raise ValueError('Training curve NLL and baselines must be finite and nonnegative')
    # Use one scale and range across panels. Zero entropy is a valid baseline:
    # never suppress it, clip it to epsilon, or silently switch individual panels.
    use_log = bool(plotted_values) and min(plotted_values) > 0
    maximum = max(plotted_values, default=2.0)
    upper = max(maximum * 1.15, 1e-12)
    log_lower = min(plotted_values) / 1.25 if use_log else None
    modes = [(False, 'curves-linear'), (use_log, 'curves')]
    for logarithmic, stem in modes:
        fig, axes = plt.subplots(3, 3, figsize=(15, 12), layout='constrained')
        for ax, arm in zip(axes.flat, arms):
            chosen, baseline_lines = panels[arm]
            for run in chosen:
                s, history = run['summary'], run['metrics']
                ax.plot([m['update'] for m in history], [m['training_bits_per_base'] for m in history],
                        color=colors.get(s['seed'], '#584996'), linewidth=1.7, marker='o', markersize=3,
                        label=f"Seed {s['seed']}" + ('' if s['complete'] else ' (partial)'))
            for style, label, value in baseline_lines:
                ax.axhline(value, color='#686868', linestyle=style, linewidth=1, label=label)
            if chosen:
                ax.legend(fontsize=7, loc='best')
            else:
                ax.text(.5, .5, 'No recorded run', transform=ax.transAxes, ha='center', va='center')
            ax.set_title(LABELS[arm], fontsize=10)
            ax.set_xlabel('Optimizer updates')
            ax.set_ylabel('Training NLL (bits/base)\n' + ('Logarithmic y-axis' if logarithmic else 'Linear y-axis'))
            if logarithmic:
                ax.set_yscale('log')
                ax.set_ylim(log_lower, upper)
                ax.grid(which='minor', axis='y', alpha=.08)
            else:
                ax.set_ylim(0, upper)
            ax.grid(which='major', alpha=.18)
        if logarithmic:
            scale_note = 'Logarithmic y-axis; original bits/base values; identical linear view: curves-linear.svg'
        elif plotted_values and min(plotted_values) == 0:
            scale_note = 'Linear y-axis: a plotted value is zero; all observations and baselines retained without offsets'
        else:
            scale_note = 'Linear y-axis; identical observations and baselines; common range across panels'
        fig.suptitle('Secondary pilot: training fit only; each curve is one initialization\n'
                     'Full-pool panel changes both dataset and capacity\n' + scale_note, fontsize=12)
        fig.savefig(root/f'reports/secondary/{stem}.svg')
        fig.savefig(root/f'reports/secondary/{stem}.png', dpi=180)
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
        'Different cohorts have different variable-position supports; those metrics are descriptive across cohorts. '
        'A position with even one alternative codon is variable, so this mask does not isolate minority-codon observations. '
        'When all target codon positions vary somewhere in a cohort, variable-position fit equals overall fit.', '',
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
    full_support_views = [name for name, manifest in cohorts['views'].items()
                          if manifest['target_nt'] > 0 and diversity[name]['variable_target_bases'] == manifest['target_nt']]
    support_note = (' Every target codon position varies at least once in '+', '.join(f'`{name}`' for name in full_support_views)+
                    '; variable-position and overall fit are identical there, so this diagnostic does not resolve rare-variant performance.'
                    if full_support_views else '')
    contrasts = paired_contrasts(groups, plan)
    expected = Counter(j['arm'] for j in plan['jobs'])
    execution_path = root/'reports/secondary/execution.json'
    execution = read_json(execution_path) if execution_path.exists() else {}
    validate_execution(execution, sha256(root/'reports/secondary/plan.json'))
    timing_amendment = applied_timing_amendment(root, sha256(root/'reports/secondary/plan.json'))
    if execution.get('complete') and (len(completed) != len(plan['jobs']) or execution['status'] != 'complete'):
        raise ValueError('Scheduler completion disagrees with reported completed runs')
    lines = ['# Secondary pilot — training only', '', SCOPE, '',
        f"**{len(completed)} of {len(plan['jobs'])} planned runs completed their fixed update budgets.** "
        f"{len(results)-len(completed)} recorded runs are partial; {len(plan['jobs'])-len(results)} have no reported result. "
        'Partial runs are excluded from completed-run averages. A pending seed is not a failed scientific hypothesis.', '',
        f"The aggregate execution budget is {plan['aggregate_budget_seconds']/3600:g} hours. "
        + ((f"Budget-chargeable execution time is {execution['elapsed_seconds']/60:.1f} minutes. " if timing_amendment else
            f"Recorded scheduler wall time is {execution['elapsed_seconds']/60:.1f} minutes. ") if 'elapsed_seconds' in execution else '')
        + 'Run-loop timings below include interval diagnostics and checkpoints but exclude setup/initial diagnostics; scheduler time covers the broader execution.', '',
        '## Data and controlled comparisons', '',
        '| Frozen view | Distinct CDS | Distinct peptides | Families | Taxon strata | Target bases/pass | Variable target bases | Target codons observed / 64 |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    if execution:
        descriptions = {'running': 'Work is still running or its last recorded controller state is active.',
            'complete': 'Every planned run completed its update budget.',
            'budget_stopped': 'The controller stopped under its remaining-time guard.',
            'timed_out': 'The outer aggregate deadline terminated a worker; its latest checkpoint may precede termination.',
            'worker_failed': 'A worker failed; unreported work must not be described as merely pending.',
            'invalid_output': 'A worker produced missing or inconsistent output that failed identity checks.'}
        # Put status before the cohort table, preserving valid Markdown table adjacency.
        insert_at = lines.index('## Data and controlled comparisons')
        note = [f"Scheduler status: **{execution['status']}**. {descriptions[execution['status']]}", '']
        recovery = execution.get('recovery_accounting')
        if recovery and recovery.get('kind') != 'standby_budget_amendment':
            charged = recovery.get('extra_wall_seconds_charged', 0)
            note += [f"Recovery accounting conservatively charged {charged/60:.2f} additional minutes since an interrupted worker launch, "
                     'including downtime. Scheduler wall time therefore is not necessarily active training time.', '']
        if timing_amendment:
            excluded = timing_amendment['excluded_hardware_deep_idle_seconds']
            note += [f"A documented lid-triggered standby interruption stopped the original execution after eight complete runs. "
                f"The [post-launch runtime amendment](../../docs/secondary-runtime-amendment.md) excludes only {excluded/60:.2f} minutes "
                f"of Windows-recorded hardware deep idle, once. All other time remains charged. Including that preserved exclusion, "
                f"recorded execution elapsed is {(execution['elapsed_seconds']+excluded)/60:.1f} minutes. "
                'Budget-chargeable time is conservative accounting, not exact CPU-compute time. Model settings, data and update budgets did not change.', '',
                'The seed-29 base run retained its checkpoint and extra interruption diagnostic. Its model outcome remains in the matched contrasts; '
                'its interrupted timing is excluded from hardware-speed comparisons. Raw loop times below retain standby; the separate charged column applies the single documented exclusion.', '']
        lines[insert_at:insert_at] = note
    for name, manifest in cohorts['views'].items():
        d = diversity[name]
        lines.append(f"| {name} | {manifest['n_sequences']:,} | {d['n_unique_peptides']:,} | {len(manifest['family_counts'])} | "
                     f"{len(manifest['taxon_counts'])} | {manifest['target_nt']:,} | {d['variable_target_bases']:,} | {d['observed_target_codon_classes']} |")
    matched_n = cohorts['views']['human_atp8']['n_sequences']
    overlap = cohorts['overlap_distinct_rna']
    replacement_counts = [matched_n-overlap[key] for key in
        ('human_atp8__human_complex', 'human_complex__mammal_complex', 'mammal_complex__vertebrate_complex')]
    family_n, mammal_n, vertebrate_n = replacement_counts
    lines += ['', f'The family expansion replaces **{family_n}/{matched_n} CDS ({family_n/matched_n:.2%})**; '
        f'the mammalian expansion replaces **{mammal_n}/{matched_n} ({mammal_n/matched_n:.2%})**; '
        f'the nonmammalian expansion replaces **{vertebrate_n}/{matched_n} ({vertebrate_n/matched_n:.2%})**. '
        'Nonmammalian additions are ATP8 only. Small replacement fractions limit what aggregate losses reveal about broader diversity.', '',
        'The five primary arms vary tokenization or cohort composition. The three additional controls change capacity, '
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
        'the three vertebrate controls share their reference arm\'s support. Across different cohorts, changing support prevents a paired accuracy interpretation. '
        'This mask counts every observation at a variable position, not only minority codons.'+support_note, '',
        '## Controlled contrasts and prespecified expectations', '',
        'Differences follow the subtraction order in each row and are paired by initialization seed. '
        'Positive values mean higher training negative log likelihood for that metric. These are descriptive optimization contrasts, not significance tests.', '',
        'Met/reversed describes only the frozen directional training expectation evaluated on the mean paired difference. '
        'Missing seed pairs make that interpretation provisional. No direction was asserted for raw loss across different taxonomic cohorts. '
        'The capacity contrast changes width, depth and heads together; it is not a single-dimension intervention.', '',
        '| Contrast | Complete/planned pairs | Overall delta bits/base | Family-macro delta bits/base | Variable-codon delta bits/base | Prespecified expectation | Mean-pair outcome |',
        '|---|---:|---:|---:|---:|---|---|']
    for contrast in contrasts:
        pairs = contrast['pairs']
        expected_delta = statistics.mean(p[contrast['metric']] for p in pairs) if pairs else None
        outcome = expectation_outcome(expected_delta, contrast['direction'], len(pairs) == contrast['planned'])
        lines.append(f"| {contrast['label']} | {len(pairs)}/{contrast['planned']} | "
            f"{mean_sd((p['overall_delta'] for p in pairs), 6)} | {mean_sd((p['macro_delta'] for p in pairs), 6)} | "
            f"{mean_sd((p['variable_delta'] for p in pairs), 6)} | {contrast['expectation']} | {outcome} |")
    lines += ['', 'Each matched pair below has verified identical cohort, sampling trace, updates, sequence/base presentations, '
        'family exposures and scored supports. Variable support is comparable within each pair only.', '',
        '| Contrast | Seed | Overall delta bits/base | Family-macro delta bits/base | Variable-codon delta bits/base | Variable target bases | Sequence presentations | Target bases presented | Seed expectation |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---|']
    for contrast in contrasts:
        for pair in contrast['pairs']:
            lines.append(f"| {contrast['label']} | {pair['seed']} | {number(pair['overall_delta'], 6)} | "
                f"{number(pair['macro_delta'], 6)} | {number(pair['variable_delta'], 6)} | {pair['variable_target_bases']:,} | "
                f"{pair['sequences_seen']:,} | {pair['bases_seen']:,} | "
                f"{expectation_outcome(pair[contrast['metric']], contrast['direction'])} |")
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
    biology = read_json(root/'configs/secondary/biology.json')
    for key, heading in [('per_taxon', 'Per-taxon training fit'), ('per_code', 'Per-genetic-code training fit')]:
        lines += ['', f'## {heading}', '',
            'Completed matched runs only. Support counts target bases per corpus pass, not independent individuals. '
            'Taxa are the frozen primary selection strata; genetic codes describe translation conventions, not added model tokens.', '',
            '| Arm | Stratum | Initializations | Target bases | Final bits/base |', '|---|---|---:|---:|---:|']
        for arm in PRIMARY + CONTROLS:
            g = groups[arm]
            for stratum in sorted({f for s in g for f in s['final'][key]}):
                metrics = [s['final'][key][stratum] for s in g]
                label = (f"{biology['taxa'].get(stratum, stratum)} ({stratum})" if key == 'per_taxon' else
                         {'1': '1 — standard nuclear code', '2': '2 — vertebrate mitochondrial code'}.get(stratum, stratum))
                lines.append(f"| {LABELS[arm]} | {label} | {len(metrics)} | "
                    f"{support(m['scored_training_bases'] for m in metrics)} | {mean_sd(m['training_bits_per_base'] for m in metrics)} |")
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
        '| Arm | Seed | State | Updates | Sequence presentations | Presentations / unique CDS | Target bases presented | Raw loop minutes | Budget-charged loop minutes | Peak process MiB |',
        '|---|---:|---|---:|---:|---:|---:|---:|---:|---:|']
    for s in results:
        rss = s.get('process_peak_rss_bytes')
        lines.append(f"| {s['arm']} | {s['seed']} | {'Complete' if s['complete'] else 'Partial: '+s['stop_reason']} | {s['updates']:,} | "
            f"{s['sequences_seen']:,} | {number(s['sequence_exposures_per_unique_cds'], 2)} | {s['bases_seen']:,} | "
            f"{number(s['elapsed_seconds']/60, 2)} | {number(charged_loop_seconds(s, timing_amendment)/60, 2)} | "
            f"{number(rss/1024**2 if rss is not None else None, 1)} |")
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
        '- [Per-run results](results.json), [histories, metadata and baselines](runs.json), [saved-model integrity inventory](checkpoints.json), [all 64 codon counts and diversity supports](codon-coverage.json), '
        '[source receipts and hashes](source-registry.json), [hardware profile](hardware-benchmark.json).',
        '- [PNG figure](curves.png) and [SVG figure](curves.svg), with the identical measurements on [linear axes](curves-linear.svg). Raw archives and local rolling/final checkpoints remain outside Git tracking.',
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
