# Laboratory notebook

## 2026-09-24 America/Toronto / 2026-09-25 UTC

**Question:** Build the first reproducible stage of a microprotein next-token modeling experiment, using the requested small-language-model scaffold and IMPI as the initial biological source.

**Owner decisions:** Start with one family from one organism. Include close annotated protein variants rather than requiring identical amino-acid sequences. Train only; defer validation as well as testing. Preserve AlphaFold/PyMOL as a future objective.

**Operational decisions before training:** Use the human scope of the IMPI workbook, reviewed references up to 100 aa, complete CDS, same reference length and at least 95% amino-acid identity. Select by distinct usable CDS count only. Use all eligible sequences for training. Retain all 64 codon classes. Report a common predicted region after the first codon and a common loss unit of bits per base.

**Source investigation:** The IMPI page was unavailable through the text browser but downloaded directly from its official server. Its current linked workbook is IMPI-2021-Q4pre, with gene annotations rather than nucleotide sequences. UniProt supplies reviewed references and ENA supplies original CDS. This source limitation was communicated during development.

**Ascertainment correction before selection:** Initial gene-only matching missed records without gene symbols. Added exact normalized curated product names and two explicit aliases supported by ENA annotations. Retrieved 300,582 short human CDS metadata records and verified ENA's total. No training metric was used during ascertainment.

**Selection:** 61 reference candidates were screened. MT-ATP8 had 129 usable distinct CDS; MT-ND4L had 109; other mapped families had at most one. Chose MT-ATP8 by count. Frozen cohort hash and all source hashes appear in the reports.

**Implementation:** Two fixed encoders, boundary-preserving binary datasets, a small causal decoder, training-only diagnostics, empirical baselines and resumable state. The current environment has CPU PyTorch; no CUDA device is available to this runtime. Synthetic software-test fixtures are isolated from biological training data.

**Pilot specification:** The six-run, 200-update design in `configs/pilot.json` was written before the pilot ran. All three seeds are reported for both arms. No post-hoc early stopping or winner selection is performed. Source-file hashes in each run capture the implementation used, including the initial uncommitted development state.

**Interpretation:** See `reports/pilot-report.md` for measured outcomes. Training fit is evidence of working optimization, not unseen-sequence accuracy or functional preservation. This small, closely related cohort has a strong positional baseline.
