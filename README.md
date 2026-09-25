# Microprotein codon language models

A reproducible, **training-only** project comparing four-base RNA and 64-codon language models. The archived first pilot uses human MT-ATP8. The secondary pilot expands to seven small ATP synthase subunit families across human, five other mammals, and three nonmammalian vertebrates.

## Architectural Upgrades

The system includes key architectural upgrades for comprehensive microprotein modeling:

1. **Extended Sequence Window (5' UTR + Kozak + CDS + 3' UTR)**:
   - Supports 5' UTR (-50 nt), Kozak context window (-6 to +4 nt around TIS), CDS region, and 3' UTR (+30 nt).
   - `WindowConfig` and `ExtendedWindowExtractor` extract, align, and pad extended genomic contexts for models.

2. **Non-Canonical Initiation Codon (TIS) Dynamics**:
   - Expands initiation handling beyond canonical `AUG` to near-cognate start codons (`CUG`, `GUG`, `UUG`, `ACG`, `AUU`).
   - Introduces empirical initiation efficiency priors tensor (`E_init`) mapping start codons to initiation weights (e.g. `AUG`: 1.0, `CUG`: 0.6, `GUG`: 0.5, `UUG`: 0.4, `ACG`: 0.3, `AUU`: 0.2).
   - Generates reading-phase coordinate tracks relative to active TIS.

3. **Multi-Track Dataset Representation**:
   - `MultiTrackDataset` and `collate_multitrack_batch` encode samples into 4 parallel tracks:
     * **Track 1**: Sequence Tokens (base or codon token IDs)
     * **Track 2**: Segment Type IDs (0 = 5' UTR, 1 = Kozak, 2 = CDS, 3 = 3' UTR)
     * **Track 3**: Frame Offsets relative to active TIS (-1, 0, 1, 2)
     * **Track 4**: Ribo-seq P-site footprint coverage vector (float values)

4. **Multi-Tier Experimental Validation & Evidence-Weighted Loss**:
   - 4-Tier evidence hierarchy:
     * **Tier 1**: MS/MS + Ribo-seq ($w_1 = 1.0$)
     * **Tier 2**: Ribo-seq TIS ($w_2 = 0.75$)
     * **Tier 3**: Conservation ($w_3 = 0.50$)
     * **Tier 4**: Computational / sORF prediction ($w_4 = 0.10$)
   - `EvidenceWeightedLoss` computes sample-weighted cross-entropy loss $L_{\text{weighted}}$, scaling token loss by evidence tier weight $w_i$ and optional $E_{\text{init}}$ initiation priors.

## Secondary experiment

This pretraining commit fixes the [protocol](docs/secondary-protocol.md), [data card](docs/secondary-data-card.md), and [run plan](reports/secondary/plan.json) before model fitting. Training outcomes will be added in a subsequent commit. The [pipeline guide](docs/secondary-pipeline.md) explains the commands and artifacts; the [laboratory notebook](docs/secondary-lab-notebook.md) records the questions, acquisition corrections and decisions.

The five primary comparisons are human ATP8/base, human ATP8/codon, human ATP-synthase/codon, mammal ATP-synthase/codon, and vertebrate ATP-synthase/codon. They use equal distinct-CDS counts; the three expanded-family views also share family quotas. Three additional controls change capacity, family weighting or explicit positional encoding on the same vertebrate cohort. Each matched arm uses three initialization seeds. A separate larger codon model trains the full recovered pool as a capacity demonstration.

Biological eligibility limits useful data size. The pipeline combines ENA's coding index, all reviewed-reference EMBL protein links, and annotated ATP8 extraction from the declared mitochondrial parent-sequence panel. It preserves source accessions, original annotation coordinates, download receipts, hashes and rejection reasons, then deduplicates exact RNA. The original pilot's source omissions are documented; its data and results are retained.

Hardware profiling selected a common 192-wide/four-layer transformer and a separate 384-wide/six-layer codon model. Training is serialized on four CPU threads under a four-hour aggregate limit, with identity-checked checkpoints and an outer deadline. These are measured practical candidates for the installed environment, not a claim to have found the computer's absolute maximum.

```powershell
.venv\Scripts\python.exe scripts/build_secondary_data.py
.venv\Scripts\python.exe -m pytest -q
.venv\Scripts\python.exe scripts/run_secondary.py --freeze-only
.venv\Scripts\python.exe scripts/build_secondary_report.py --data-only
.venv\Scripts\python.exe scripts/run_secondary.py
.venv\Scripts\python.exe scripts/build_secondary_report.py
```

See the [pipeline guide](docs/secondary-pipeline.md) for artifact locations and snapshot constraints. Raw archives and checkpoints remain local and are excluded from Git; compact manifests, results and figures are versioned. Software tests are correctness fixtures, not biological validation or test sets.

## Initial experiment — archived

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
