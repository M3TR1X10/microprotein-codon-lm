# Microprotein codon language models

A reproducible, **training-only** proof of concept comparing a four-base RNA language model with a 64-codon language model. The first cohort is **human MT-ATP8**, selected by usable distinct coding-sequence count within an explicitly defined IMPI candidate screen.

This experiment asks whether a small causal transformer can learn the training distribution and whether both representations can be trained reproducibly. It does **not** yet measure prediction on unseen sequences. Validation, test sets, gap completion, AlphaFold and PyMOL are deferred at the project owner's request.

## Current experiment

| Property | Value |
|---|---|
| Organism | Homo sapiens; NCBI taxonomy 9606 |
| Family | MT-ATP8 / ATP synthase F(0) complex subunit 8 |
| Reference | UniProt P03928; 68 amino acids |
| Eligible coding sequences | 129 distinct sequences, each 207 nucleotides including stop |
| Distinct translated proteins | 53 close variants |
| Quality threshold | Same reference length, at least 95% amino-acid identity, valid complete CDS |
| Genetic code | Vertebrate mitochondrial, NCBI translation table 2 |
| Allocation | 129 training; 0 validation; 0 test |
| Models | Two-layer causal transformers; 48-dimensional embeddings; 4 attention heads |
| Pilot | 200 optimizer updates per arm for each of seeds 17, 29, 43 |

The quality threshold is an operational definition of a close annotated family, **not experimental proof that every variant retains equivalent function**. All sequences are genuine database CDS; no reverse translation or synthetic augmentation is used.

IMPI supplies gene identity and mitochondrial evidence. It does not supply nucleotide sequences: **all training nucleotides come from ENA**, with UniProt used for reviewed reference sequences and names. This differs from the requested bulk source because the IMPI download contains annotations rather than sequence data. See [the data card](docs/data-card.md) for the complete source distinction and selection limits.

## Run the pipeline

Python 3.11+ is supported; the recorded pilot environment uses Python 3.12.14 and PyTorch 2.14.0 on CPU. Commands below assume the repository root as the current directory.

```powershell
python -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements-lock.txt
.venv\Scripts\python.exe -m pip install --no-deps -e .

.venv\Scripts\python.exe -m microprotein_lm.cli discover
.venv\Scripts\python.exe -m microprotein_lm.cli select
.venv\Scripts\python.exe -m microprotein_lm.cli prepare
.venv\Scripts\python.exe -m microprotein_lm.cli baselines
.venv\Scripts\python.exe scripts/run_pilot.py
```

On Linux/macOS, use `.venv/bin/python`. For a compatible CUDA installation, install the appropriate PyTorch build using [PyTorch's installation instructions](https://pytorch.org/get-started/locally/); record that environment separately. Default precision is float32; CUDA bfloat16 is optional when supported. GPU availability is checked rather than assumed.

To reproduce the **original frozen source snapshot** rather than query current databases, first run `python scripts/restore_snapshot.py`, then the commands above. Snapshot restoration verifies every SHA-256 and fails if an upstream record has changed. Keep an archival copy of `data/raw` for long-term reproducibility; an accession URL is not a permanent content guarantee.

Run one arm or resume an interrupted run:

```powershell
.venv\Scripts\python.exe -m microprotein_lm.cli train --mode codon --seed 17 --output runs/my-codon-run
.venv\Scripts\python.exe -m microprotein_lm.cli train --mode codon --seed 17 --output runs/my-codon-run --resume
.venv\Scripts\python.exe -m pytest -q
```

Resume requires unchanged configuration, source code and cohort. It resumes the original update budget; it does not extend a completed run. Existing outputs are protected against accidental replacement. The software tests use synthetic fixtures to check correctness; they are not biological model test sets.

## Where to look

| Path | Purpose |
|---|---|
| [docs/protocol.md](docs/protocol.md) | Questions, hypotheses, controls, fixed decisions and interpretation rules |
| [docs/data-card.md](docs/data-card.md) | Eligibility, provenance, filtering, count-based selection and limitations |
| [docs/experiment-plan.md](docs/experiment-plan.md) | Small experimental stages and hyperparameter decisions |
| [docs/scaffold.md](docs/scaffold.md) | How the requested SLM scaffold was adapted |
| [docs/lab-notebook.md](docs/lab-notebook.md) | Decision and execution record |
| [configs/pilot.json](configs/pilot.json) | Executable initial hyperparameters |
| [reports/pilot-report.md](reports/pilot-report.md) | Measured training results and their interpretation |
| `src/microprotein_lm/` | Acquisition, QC, tokenization, data loading, model, training and baselines |
| `data/raw/` | Cached original downloads and source receipts; ignored by Git |
| `data/processed/` | One common cohort plus the two encoded representations; ignored by Git |
| `runs/` | Checkpoints, exact run settings and training diagnostics; ignored by Git |
| `reports/` | Compact reproducibility records and results committed to the private repository |

The setup follows [ChaitanyaK77's Building-a-Small-Language-Model-SLM-](https://github.com/ChaitanyaK77/Building-a-Small-Language-Model-SLM-): preparation, fixed token IDs in binary storage, a compact causal decoder, and configurable training. The upstream MIT notice is retained in [references/SLM-LICENSE](references/SLM-LICENSE). No external experiment-tracking account or paid training service is required.
