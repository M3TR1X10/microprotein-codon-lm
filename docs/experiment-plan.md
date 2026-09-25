# Small experiments and decision gates

## Completed foundation

| Subproblem | Deliverable | Acceptance criterion |
|---|---|---|
| Resolve organism/family scope | Recorded owner decisions | One human family, close variants, no validation/test |
| Distinguish gene catalog from nucleotide data | Source registry and data card | Every nucleotide comes from a real CDS record |
| Choose family by quantity | Candidate ranking and QC ledger | Largest count after fixed gates and deduplication |
| Encode two views | Base and codon binary datasets | Exact round trips, same records, 4/64 fixed classes |
| Protect biological boundaries | Offset-indexed loader and masks | No training transition from one CDS to another |
| Verify causality and reproducibility | Offline software checks | Future tokens cannot affect earlier logits; resume matches |
| Establish training feasibility | Prespecified six-run pilot | Finite optimization, saved state, honest training diagnostics |

## Initial hyperparameters

`configs/pilot.json` defines the executable experiment; this table explains the choices.

| Parameter | Value | Reason |
|---|---:|---|
| Vocabulary | 4 bases / 64 triplets | Exhaustive biological symbols, no learned merges |
| Embedding width | 48 | Small CPU-feasible initial representation, independent of vocabulary size |
| Layers / heads | 2 / 4 | Compact causal decoder; 12 dimensions per attention head |
| Positional encoding | Fixed sinusoidal | Removes unequal learned position-table parameter counts |
| Context capacity | 303 nt | Supports a 100-aa CDS plus full stop in either representation |
| Dropout | 0.1 | Conservative starting regularization, not tuned on held-out data |
| Batch / accumulation | 16 / 1 | Same sequence exposure per update in both arms |
| AdamW LR / minimum | 0.001 / 0.0001 | Explicit positive warmup/cosine schedule |
| Warmup / updates | 20 / 200 | Bounded optimization pilot rather than open-ended training |
| Weight decay / clipping | 0.01 / 1.0 | Initial optimizer stabilization settings |
| Seeds | 17, 29, 43 | Fixed seed list, report every run |
| Precision | float32 | Portable, inspectable CPU starting point |

These are starting values, not scientifically established optimal hyperparameters. Lower training loss does not establish that a configuration is better at prediction.

## Training-only follow-on stages

Do not launch these automatically as part of the first pilot. Record a protocol amendment and resource budget before each stage.

1. **Optimizer sensitivity:** Keep cohort, architecture and exposures fixed. Compare peak learning rates `0.0003` and `0.001` using the same seeds and both arms. Inspect numerical stability, loss curves and resource cost. Do not report the lower-loss configuration as the best generalizer.
2. **Capacity sensitivity:** Fix the stable optimizer setting; compare widths 24 and 48, then one and two layers. Change one factor at a time. Measure parameters, time and training fit. A positional categorical baseline remains a required reference.
3. **Compute scaling:** Increase accumulation before increasing GPU memory demands; keep effective batch size and exposure budget explicit. Benchmark float32 versus supported CUDA bfloat16 separately. Record the device and library versions. Do not compare loss at unmatched exposures and call it a tokenization effect.
4. **Data scaling:** Quantify new distinct sequences and provenance first. A future multi-organism expansion requires a new versioned cohort and explicit authorization; it is not an augmentation of the existing human cohort.

Distributed training, experiment-tracking services and cloud schedulers are intentionally unnecessary for this proof of concept. The module boundaries and immutable run directories support later additions without rewriting the biological filters or tokenizers.

## Deferred decision

Determining whether the model shows predictive promise requires validation or unseen-sequence evaluation, which the owner has deferred. No amount of training-only optimization resolves that question. Keep the final test/gap/structure work as a separate future stage; do not consume unseen examples while tuning this phase.
