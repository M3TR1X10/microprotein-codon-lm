# Resume checklist — 25 September 2026

The owner requested this handoff because account usage is nearly exhausted. **Leave the existing training process running. Do not start a second scheduler.** The owner explicitly said not to stop training.

## Current state

- Workspace: `C:\Users\jddub\OneDrive\Desktop\Small microprotein next codon prediction modeling using LLM trainging techniques`.
- Private repository: <https://github.com/M3TR1X10/microprotein-codon-lm>.
- Published/local `main`: `284e9eb26617a0533e0cbd586454b3380793b79b`. Its CI passed: <https://github.com/M3TR1X10/microprotein-codon-lm/actions/runs/36151472709>.
- Immutable pretraining commit: `70941dc89abd40803a153b9b407a8d3a8f779b1d`. All 18 frozen input/engine hashes and the plan were verified unchanged after resumption.
- **Training is active in native execution session `43882`.** Last observed at approximately 15:19 UTC: **9/25 runs complete**; job 10 is `human_atp8_codon`, seed 29. Human ATP8/base seed 29 completed its 2,000 updates, final training bits/base approximately 0.06309. There were about 142 minutes left after the shutdown reserve when job 10 began; this is a timestamped observation, not a live estimate.
- Original four-hour training allowance remains enforced by the existing scheduler. One documented, published, one-time standby correction has already been applied. **Do not apply it again.** Read `docs/secondary-runtime-amendment.md` and the immutable evidence/archive plus application receipt in `reports/secondary/interruption/`.
- The original training plan, code, cohorts, model settings, seed order and endpoints must remain unchanged. Current reporting files may still describe only the first eight runs and are not final.
- New user authorization: develop and execute a separate progressive gap test, **original nine species only**, comparing prefix-only completion with right-flank candidate reranking. No validation-based tuning or AlphaFold/PyMOL work is authorized yet.
- New test acquisition found **20 candidate untouched CDS**: 17 cattle ATP8, one human ATP8, one mouse ATP5F1E and one mouse ATP5ME. Six encode a peptide already present in training. Other panel species currently contribute no novel eligible sequence. These are a small, cattle-dominated holdout, not a balanced nine-species benchmark.
- Candidate test data have **not been frozen**, and **no biological gap inference has run**. Acquisition processes have finished. A final cached QC rebuild is required after the last acquisition-source changes.
- New test-stage files and ongoing training outputs are local, mostly untracked/uncommitted. Preserve them. The latest published commit does not yet contain this test work or final training outcomes.

## 1. Let training finish and verify it

Owner: coordinator; use an independent audit agent after completion.

1. Read `reports/secondary/execution.json` and poll existing session 43882 if available. Do not launch `scripts/run_secondary.py` while the original scheduler is alive. A process inventory attempt through CIM was denied by the sandbox; the native session itself confirmed active training.
2. If training finishes, record completed/partial/missing jobs and the final charged elapsed time. Do not extend the update or time budget to hide an incomplete matrix. Preserve partial checkpoints.
3. If the app/session is interrupted, inspect execution state and running processes before any resume. The original runner conservatively charges interruption downtime. No additional sleep credit is automatic.
4. Run `.venv\Scripts\python.exe work/audit_secondary_commit.py` to verify frozen sources, then `.venv\Scripts\python.exe scripts/audit_secondary_checkpoints.py`. The existing inventory contains eight verified completed runs. The audit must extend it while preserving their established hashes; never delete that inventory to bypass a mismatch.
5. Run `.venv\Scripts\python.exe scripts/build_secondary_report.py`, inspect both logarithmic and linear scientific figures, and update the README/notebook with actual final results. Keep training fit distinct from held-out accuracy. Explain the 403/412 shared CDS and only nine ATP8 replacements between matched mammal/vertebrate views.
6. Report the interrupted base/seed29 model in scientific contrasts, but exclude its raw duration from hardware-speed comparisons. Preserve raw wall time and separately label budget-chargeable time.

## 2. Finish and independently audit test-data preparation

Owner: biology/acquisition agent; reviewer: experimental-design agent.

1. Read `docs/gap-acquisition-handoff.md`, `docs/gap-data-card.md`, `configs/gap-biology.json` and `src/microprotein_lm/gap_acquire.py`.
2. Complete the final cached rebuild: `.venv\Scripts\python.exe scripts/build_gap_data.py --include-outer-parents`. This must not alter any training file. Confirm current-source hashes, raw receipt/data hashes, coordinate/translation agreement and the acquisition report's counts.
3. Verify exclusion against the original pilot and the entire 533-CDS training union, including the full-pool capacity arm. Verify global RNA hashes and stable source protein/parent accessions. A hidden gap in a training sequence is not an unseen test sequence.
4. Audit taxonomy, family, own-species reviewed reference, complete CDS, <=100 aa, >=95% reference identity, genetic code and all rejection reasons. Unreviewed UniProt entries are retrieval leads only; they do not replace reviewed QC references.
5. Retain both retrieval routes and their ascertainment limits: 1,859 unreviewed-entry EMBL links plus the bounded, metadata-selected outside-length mitochondrial-parent route. Preserve the exact final inventory counts from the report rather than this preliminary handoff.
6. Verify nearest-training nucleotide/peptide descriptors and source qualifiers. Identical peptides and close homologs can remain after exact RNA exclusion; do not call records independent people or lineages.
7. After review, freeze with `.venv\Scripts\python.exe scripts/build_gap_data.py --include-outer-parents --freeze`. The manifest must refuse future silent changes. Do not relax eligibility to increase test size.

## 3. Complete the prospective test plan and implementation

Owner: evaluation-engine agent; reviewer: experimental-design agent.

Read `docs/gap-evaluation-protocol.md`, `configs/gap-evaluation.json` and `docs/gap-design-handoff.md`. These are drafts until reviewed, frozen and published before biological inference.

- Primary matrix: eight completed matched arms, seeds 17/29/43, same held-out sequences and nested gaps of 1, 2, 4, 8 and 16 codons. Score true-gap joint bits/base and free-running codon-synchronous greedy reconstruction. All 64 codons remain allowed; stops are outcomes, not repaired errors.
- Focused context matrix: human ATP8 base/codon and standard vertebrate codon; initial left contexts 1/2/4/8/16/32/full. Preserve absolute nucleotide positions. **Only one human test CDS was found so far**, so this panel cannot support population claims or meaningful biological uncertainty bands.
- Right-flank matrix: same focused models, left 8/full, right 0/1/2/4/8, deterministic beam width eight. Freeze each prefix-only candidate bank before looking at the suffix or correctness; reuse the bank across right lengths. Compare against right=0 from the same bank. These are restricted-candidate joint ranking scores, not a normalized bidirectional posterior.
- Optional long gaps 32/48/64 use unchanged anchors only where supported; compare shorter lengths on the same support intersection. The full-pool larger model is a separate, confounded engineering extension.
- Proposed separate inference cap is 90 minutes, after training, with four CPU threads and one loaded checkpoint at a time. Synthetic profiling must precede final plan freeze. Do not run substantial inference concurrently with training.

Existing new implementation:

| File | State / remaining work |
|---|---|
| `src/microprotein_lm/gap_inference.py` | External KV cache, absolute positions, exact base-to-64-codon probabilities, greedy/beam generation and suffix reranking. Eighteen targeted tests passed, including causality, checkpoint nonmutation and full-path tie handling. |
| `src/microprotein_lm/gap_metrics.py` | Genetic-code-aware reconstruction, minority/unseen/rare/conserved/variable support. Four initial tests passed; the latest rare/conservation extension still needs a rerun and targeted coverage. |
| `src/microprotein_lm/gap_baselines.py` | Uniform, unigram, bigram, position and privileged family/taxon/position baselines fit only to each model's training data. Four targeted tests passed. |
| `src/microprotein_lm/gap_evaluate.py` | New case-level checkpoint loading, flank splitting, probability/generation, bank persistence, baselines and result validation. **Not yet tested or integrated.** Review identities, resume behavior, capacity boundaries, nonmutation and truth separation before use. |
| `src/microprotein_lm/gap_design.py` | Assigned to design agent; consult its handoff for final existence/completion and tests. Must provide deterministic panel/anchor/case selection, leakage verification and 99% within-stratum sequence clusters. |
| `scripts/run_gap_evaluation.py` | **Still to implement.** Frozen plan, audited final checkpoint allowlist, source/data hashes, bounded serialized execution, durable complete-case results, candidate-bank reuse, conservative resume accounting and explicit missing-cell reporting. |
| `scripts/build_gap_report.py` | **Still to implement.** Paired support-aware summaries and scientific figures described below. |

The coordinator and reviewer must verify that all common conservation masks use the **same complete training union**, while per-model family/taxon visibility and baseline fitting use that model's own training cohort.

## 4. Build and verify the graph/report collection

Owner: reporting/statistics agent; reviewer: biology plus experimental-design agents.

1. Produce leakage/acquisition flow and actual support counts by family/taxon/code; display empty strata and every missing or unexecuted cell.
2. Plot teacher-forced bits/base, exact-gap match, base/codon accuracy and amino-acid accuracy versus gap length, with training-only baselines and separate individual-seed traces.
3. Add left-context-by-gap heatmaps, within-gap error trajectories, right-context change from same-bank right=0, candidate recall/attainable-accuracy ceilings, premature-stop rates and conservation/allele-frequency results.
4. Pair models/conditions on identical sequence/anchor support. Distinguish primary versus optional long-gap/full-pool results. Never pool unsupported taxa into a claim about all nine species.
5. Treat distinct CDS as repeated-measure units; gaps, codons and model seeds do not multiply biological sample size. Cluster near-identical CDS before any optional bootstrap. The protocol requires at least five clusters per contributing stratum for bands; the current human/mouse strata will fail that threshold. Show counts and seed ranges instead of invented confidence intervals.
6. Save raw per-case predictions, truth, coordinates, candidate-bank identities and machine-readable summaries locally; publish compact reproducibility records and figures. Inspect every figure for correct labels, support and legibility.

## 5. Freeze, publish, execute and close

Owner: coordinator, after explicit independent review reports.

1. Run the completed offline test suite. Before this stage, 118 prior tests passed and the amendment's CI passed. Subsequent targeted checks passed as described above; this is **not yet a final full-suite result**. Fixtures are not research data.
2. Finish synthetic inference profiling and freeze the executable test plan with exact source/config/protocol/test/case/checkpoint hashes. Allow only completed integrity-audited final checkpoints, never minimum-loss checkpoint selection.
3. **Publish the test protocol, implementation and frozen plan before any biological model inference.** Preserve the already published training preregistration and operational amendment.
4. Run the bounded evaluator once training has stopped. Resume only unchanged identities and complete-case outputs. Preserve failures and partials; do not tune methods on test outcomes.
5. Generate final test reports/graphs, obtain independent statistical/biological interpretation reviews, update README navigation, publish results to the private repository and confirm CI plus local/remote commit identity.

GitHub publication uses the authenticated connector because shell Git HTTPS credentials are unavailable. Existing ignored helpers are `work/github_payload.py` and `work/import_remote_commit.py`; inspect them before reuse. Stage exact reviewed files, verify the remote parent and resulting tree SHA, update main non-forcibly, then import the exact connector-created commit locally. `.git` writes require sandbox escalation. Raw archives/checkpoints remain ignored. Do not commit credentials or large model files.

## Multi-agent coordination rules for the next session

- Coordinator owns the live training process, frozen-plan checks, scope/budget decisions, integration and publication. **No agent may kill/restart training or change a frozen file.**
- Biology agent owns acquisition/QC/test manifests and reviews biological claims.
- Experimental-design agent owns cases, leakage/replication checks and prospective-method review; it independently reviews the engine's information boundaries.
- Evaluation-engine agent owns inference runner/resume logic and meaningful correctness fixtures; it does not select data by performance.
- Reporting agent owns paired aggregation, uncertainty and figures. With limited concurrency, reuse an agent only after its preceding deliverable is complete.
- Assign disjoint file ownership, state exact APIs before coding, and send concise cross-review findings. Independent review must name evidence, unresolved limitations and concrete fixes; it is not a claim of certification.
- Do not spawn more work merely to spend time. Advance bounded questions: data eligibility, leakage, decoding correctness, equal support, timing, then interpretation. Keep training and testing artifacts separate.
- User requested a wrap-up now. No automatic follow-up or monitoring automation was created. Continue when the user resumes work; the already launched training process should keep running on its own while the computer stays awake.
