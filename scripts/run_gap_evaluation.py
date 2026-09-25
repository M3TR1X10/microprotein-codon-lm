"""Prospective gap evaluation execution runner.

Evaluates audited completed pretraining model checkpoints on frozen holdout CDS
across nested gaps, left contexts, and right-flank candidate reranking.
Enforces torch.inference_mode(), zero state mutation, atomic checkpointing,
and a 90-minute local execution cap.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from microprotein_lm.gap_baselines import GapBaselines
from microprotein_lm.gap_design import load_verified_holdout, make_cases
from microprotein_lm.gap_evaluate import atomic_json, evaluate_case, load_model, validate_completed_result
from microprotein_lm.gap_metrics import TrainingPositionSupport
from microprotein_lm.io import now, read_json, sha256


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def find_model_specs(root, eval_config):
    root = Path(root)
    checkpoints_path = root / 'reports/secondary/checkpoints.json'
    if not checkpoints_path.exists():
        raise FileNotFoundError('Secondary checkpoint audit record missing; run audit_secondary_checkpoints.py')
    checkpoints_report = read_json(checkpoints_path)
    if checkpoints_report.get('status') != 'passed':
        raise ValueError('Secondary checkpoint audit did not pass')
    
    verified_map = {}
    for run in checkpoints_report.get('verified_runs', []):
        key = (run['arm'], run['seed'])
        verified_map[key] = run
        
    plan = read_json(root / 'reports/secondary/plan.json')
    runs_by_arm_seed = {}
    for r in plan.get('jobs', []):
        runs_by_arm_seed[(r['arm'], r['seed'])] = r

    matched_arms = [
        'human_atp8_base', 'human_atp8_codon', 'human_complex_codon',
        'mammal_complex_codon', 'vertebrate_complex_codon', 'vertebrate_small',
        'vertebrate_family_macro', 'vertebrate_no_position'
    ]
    primary_seeds = eval_config['primary']['model_seeds'] # [17, 29, 43]
    
    specs = []
    # Primary matched arms & seeds
    for arm in matched_arms:
        for seed in primary_seeds:
            key = (arm, seed)
            if key not in verified_map:
                print(f"[Notice] Checkpoint for {arm}-seed{seed} not yet complete/verified. Excluding.", flush=True)
                continue
            v_info = verified_map[key]
            run_info = runs_by_arm_seed[key]
            model_id = f"{arm}-seed{seed}"
            meta_path = f"runs/secondary/{model_id}/summary.json"
            ckpt_file = v_info['files']['final.pt']['file']
            ckpt_path = root / ckpt_file
            ckpt_obj = torch.load(ckpt_path, weights_only=True, map_location='cpu')
            specs.append({
                'model_id': model_id,
                'arm': arm,
                'seed': seed,
                'mode': run_info['mode'],
                'dataset': run_info['dataset'],
                'role': 'primary_matched',
                'checkpoint': v_info['files']['final.pt'],
                'checkpoint_identity': ckpt_obj['identity'],
                'run_metadata_file': meta_path,
                'cohort_path': f"data/processed/secondary/{run_info['dataset']}/cohort.jsonl"
            })

    # Exploratory full_vertebrate_capacity seed 17
    exp_arm = eval_config['exploratory']['full_pool_arm']
    exp_seed = eval_config['exploratory']['full_pool_seed']
    if (exp_arm, exp_seed) in verified_map:
        v_info = verified_map[(exp_arm, exp_seed)]
        run_info = runs_by_arm_seed[(exp_arm, exp_seed)]
        model_id = f"{exp_arm}-seed{exp_seed}"
        ckpt_file = v_info['files']['final.pt']['file']
        ckpt_path = root / ckpt_file
        ckpt_obj = torch.load(ckpt_path, weights_only=True, map_location='cpu')
        specs.append({
            'model_id': model_id,
            'arm': exp_arm,
            'seed': exp_seed,
            'mode': run_info['mode'],
            'dataset': run_info['dataset'],
            'role': 'exploratory_full_pool',
            'checkpoint': v_info['files']['final.pt'],
            'checkpoint_identity': ckpt_obj['identity'],
            'run_metadata_file': f"runs/secondary/{model_id}/summary.json",
            'cohort_path': f"data/processed/secondary/{run_info['dataset']}/cohort.jsonl"
        })

    return specs


def read_jsonl(path):
    with path.open(encoding='utf-8') as f:
        return [json.loads(line) for line in f if line.strip()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='configs/gap-evaluation.json')
    parser.add_argument('--output', default='reports/gap/results.json')
    parser.add_argument('--bank-dir', default='reports/gap/banks')
    parser.add_argument('--max-seconds', type=int, default=5400)
    parser.add_argument('--dry-run', action='store_true', help='Evaluate a 2-case dry run profile')
    args = parser.parse_args()

    start_wall_time = time.time()
    config_path = ROOT / args.config
    eval_config = read_json(config_path)
    plan_hash = sha256(config_path)

    torch.set_num_threads(eval_config['runtime']['cpu_threads'])
    torch.set_num_interop_threads(eval_config['runtime']['interop_threads'])

    bank_dir = ROOT / args.bank_dir
    bank_dir.mkdir(parents=True, exist_ok=True)
    out_path = ROOT / args.output
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("Loading verified holdout...", flush=True)
    holdout = load_verified_holdout(ROOT, evaluation_config=args.config)
    records = holdout['records']
    records_by_sha256 = {r['sequence_sha256']: r for r in records}
    training_records = holdout['training_records']
    union_support = TrainingPositionSupport(training_records)

    print("Building deterministic test cases...", flush=True)
    cases_plan = make_cases(records, eval_config)
    
    stages = cases_plan['support']['stages']
    print(f"Case breakdown: Primary={len(cases_plan['primary'])}, Focused={len(cases_plan['focused_context'])}, RightFlank={len(cases_plan['right_flank'])}, Exploratory={len(cases_plan['exploratory'])}", flush=True)

    model_specs = find_model_specs(ROOT, eval_config)
    print(f"Found {len(model_specs)} verified completed model checkpoints.", flush=True)
    if not model_specs:
        raise ValueError("No verified completed checkpoints found to evaluate!")

    # Cache per-dataset cohorts & supports to avoid redundant I/O
    model_cohort_supports = {}
    model_baselines = {}
    for spec in model_specs:
        ds = spec['dataset']
        if ds not in model_cohort_supports:
            cohort_records = read_jsonl(ROOT / spec['cohort_path'])
            # Ensure family label is present for all records
            for r in cohort_records:
                if not r.get('family') and r.get('gene') in ('MT-ATP8', 'ATP8'):
                    r['family'] = 'ATP8'
            model_cohort_supports[ds] = TrainingPositionSupport(cohort_records)
            model_baselines[ds] = GapBaselines(cohort_records, prior_total_mass=eval_config['baselines']['prior_total_mass'])

    # Load existing results if resuming
    existing_results = {}
    if out_path.exists():
        try:
            prev = read_json(out_path)
            if prev.get('plan_sha256') == plan_hash:
                for res in prev.get('results', []):
                    key = (res['model_id'], res['case']['case_id'])
                    existing_results[key] = res
                print(f"Resuming existing execution with {len(existing_results)} completed results.", flush=True)
        except Exception as e:
            print(f"Warning: Failed to load previous results ({e}). Starting fresh.", flush=True)

    new_results = dict(existing_results)
    budget_exceeded = False

    beam_width = eval_config['right_flank']['beam_width']

    total_evaluations = 0
    skipped_evaluations = 0

    for model_spec in model_specs:
        if budget_exceeded:
            break

        model_id = model_spec['model_id']
        ds = model_spec['dataset']
        model_support = model_cohort_supports[ds]
        baselines = model_baselines[ds]

        # Determine cases applicable to this model
        applicable_cases = []
        # Primary cases apply to all matched arms
        if model_spec['role'] == 'primary_matched':
            applicable_cases.extend(cases_plan['primary'])

        # Focused context and Right flank apply to specified arms
        if model_spec['arm'] in eval_config['focused_context']['arms']:
            applicable_cases.extend(cases_plan['focused_context'])
            applicable_cases.extend(cases_plan['right_flank'])

        # Exploratory cases apply to exploratory arm or primary arms
        if model_spec['role'] == 'exploratory_full_pool':
            applicable_cases.extend(cases_plan['exploratory'])
        elif model_spec['role'] == 'primary_matched':
            applicable_cases.extend(cases_plan['exploratory'])

        if args.dry_run:
            applicable_cases = applicable_cases[:2]

        # Filter out already computed cases
        uncomputed = [c for c in applicable_cases if (model_id, c['case_id']) not in new_results]
        if not uncomputed:
            skipped_evaluations += len(applicable_cases)
            continue

        print(f"Loading checkpoint {model_id} ({len(uncomputed)} cases to run)...", flush=True)
        decoder = load_model(ROOT, model_spec)

        for case in uncomputed:
            elapsed = time.time() - start_wall_time
            if elapsed >= args.max_seconds:
                print(f"[Cap Reached] Reached 90-minute limit ({elapsed:.1f}s / {args.max_seconds}s). Stopping cleanly.", flush=True)
                budget_exceeded = True
                break

            rec = records_by_sha256[case['sequence_sha256']]
            res = evaluate_case(
                decoder=decoder,
                case=case,
                record=rec,
                model_spec=model_spec,
                union_support=union_support,
                model_support=model_support,
                baselines=baselines,
                plan_hash=plan_hash,
                bank_directory=bank_dir,
                beam_width=beam_width
            )
            validate_completed_result(res, plan_hash, model_spec, case)
            new_results[(model_id, case['case_id'])] = res
            total_evaluations += 1

            # Save periodically
            if total_evaluations % 20 == 0:
                snapshot = {
                    'schema_version': 'gap-results-v1',
                    'created_utc': now(),
                    'plan_sha256': plan_hash,
                    'eval_config_sha256': plan_hash,
                    'manifest_identity': holdout['identity'],
                    'n_evaluations': len(new_results),
                    'budget_exceeded': budget_exceeded,
                    'elapsed_seconds': time.time() - start_wall_time,
                    'results': list(new_results.values())
                }
                atomic_json(out_path, snapshot)

        del decoder
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    total_elapsed = time.time() - start_wall_time
    final_snapshot = {
        'schema_version': 'gap-results-v1',
        'created_utc': now(),
        'plan_sha256': plan_hash,
        'eval_config_sha256': plan_hash,
        'manifest_identity': holdout['identity'],
        'n_evaluations': len(new_results),
        'budget_exceeded': budget_exceeded,
        'elapsed_seconds': total_elapsed,
        'models_evaluated': len(model_specs),
        'results': list(new_results.values())
    }
    atomic_json(out_path, final_snapshot)

    print(f"\nGap evaluation run finished!")
    print(f"Total evaluated cases: {total_evaluations} (Skipped cached: {skipped_evaluations})")
    print(f"Total stored results: {len(new_results)}")
    print(f"Elapsed wall time: {total_elapsed:.2f}s / {args.max_seconds}s limit")
    print(f"Output saved to: {out_path}")


if __name__ == '__main__':
    main()
