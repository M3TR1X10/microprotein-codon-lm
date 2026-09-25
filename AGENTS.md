# Working rules for this scientific project

Read `docs/protocol.md` and the current dataset/run manifests before changing the experiment.

- Preserve the archived first pilot and frozen ATP-synthase-only secondary training in `docs/secondary-protocol.md`, with equal-count comparisons and a four-hour local training budget. The owner subsequently authorized a separate prospective test stage in `docs/gap-evaluation-protocol.md`: original nine species only, progressive codon gaps, prefix-only generation and right-flank reranking. Keep its untouched data and inference outputs separate from training. Validation-based tuning and AlphaFold/PyMOL remain deferred.
- Preserve the documented one-time standby correction in `docs/secondary-runtime-amendment.md` and its original archives. Budget-chargeable time differs from raw elapsed time; never silently reapply the credit or change the frozen scientific settings.
- Separate hypotheses, operational assumptions, measurements and interpretations. Record protocol changes before running the affected experiment. Never turn training fit into a claim of generalization or molecular function.
- Break new work into a specific biological or engineering question, a bounded change, an appropriate check and a recorded result.
- Retain real sequence provenance, accession versions, retrieval receipts and hashes. Do not reverse-translate proteins, invent variants, repair sequences or silently supplement organisms outside the declared panel.
- Select the cohort by declared quality gates and usable distinct counts. Never use model performance to choose biological data after the fact.
- Keep base/codon comparisons aligned on data, biological target region and exposure. Preserve exact 4/64 vocabularies and reading frames.
- Preserve source/model configurations and run outputs. New experiments use new run directories. Repeated reads/checks are not independent experiments.
- Run appropriate offline software checks for changes to biological QC, tokenization, data boundaries, objectives or resume behavior. Software fixtures are never training data.
- Store raw data and checkpoints outside Git tracking. Commit compact methods, manifests and honest reports; retain upstream attribution.
- Do not introduce cloud spending, external tracking accounts, structural-model downloads or broad hyperparameter sweeps without an agreed need and budget.
