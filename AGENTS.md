# Working rules for this scientific project

Read `docs/protocol.md` and the current dataset/run manifests before changing the experiment.

- Preserve the owner's current scope: one human protein family with close annotated variants; training only. Validation, test-set work, gap completion and AlphaFold/PyMOL are deferred until requested.
- Separate hypotheses, operational assumptions, measurements and interpretations. Record protocol changes before running the affected experiment. Never turn training fit into a claim of generalization or molecular function.
- Break new work into a specific biological or engineering question, a bounded change, an appropriate check and a recorded result.
- Retain real sequence provenance, accession versions, retrieval receipts and hashes. Do not reverse-translate proteins, invent variants, repair sequences or silently supplement other organisms.
- Select the cohort by declared quality gates and usable distinct counts. Never use model performance to choose biological data after the fact.
- Keep base/codon comparisons aligned on data, biological target region and exposure. Preserve exact 4/64 vocabularies and reading frames.
- Preserve source/model configurations and run outputs. New experiments use new run directories. Repeated reads/checks are not independent experiments.
- Run appropriate offline software checks for changes to biological QC, tokenization, data boundaries, objectives or resume behavior. Software fixtures are never training data.
- Store raw data and checkpoints outside Git tracking. Commit compact methods, manifests and honest reports; retain upstream attribution.
- Do not introduce cloud spending, external tracking accounts, structural-model downloads or broad hyperparameter sweeps without an agreed need and budget.
