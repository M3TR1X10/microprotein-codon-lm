# Observed gap-test data

This benchmark is acquired after the secondary training matrix was frozen and before any gap predictions. It contains real observed CDS from the original nine species and seven declared ATP-synthase families only. It does not remove or relabel training observations as test data. The exact test counts and source receipts are recorded in [the acquisition report](../reports/gap/acquisition.json); a final manifest is written only by the explicit freeze command.

## Acquisition and biological gates

The [gap acquisition configuration](../configs/gap-biology.json) binds the full **533-CDS** training union by SHA-256. This includes the original 129-CDS pilot and all sequences available to the full-pool capacity run, even if that run is still pending when test acquisition begins. Reviewed own-species references, InterPro family mappings, genetic codes and CDS quality gates are inherited without alteration from the frozen secondary experiment.

The first route queries unreviewed UniProt records from the same families and species and follows all real EMBL protein-accession links. Unreviewed entries are discovery leads only: they do not replace reviewed biological references or establish experimentally verified function. Each source CDS must independently satisfy the original complete-CDS, correct-taxonomy, unambiguous-base, allowed-start/full-stop, translation-agreement, canonical-length and >=95% own-species protein-identity criteria. There are no synthetic, reverse-translated or repaired sequences.

The second route queries ENA mitochondrial parent sequences outside the training route's 14,000–20,000-nt interval, bounded to 150–100,000 nt. A metadata-only inventory found 186,099 records. Broad `complete sequence` text predominantly identified control-region/rRNA fragments, so the declared payload screen uses case-insensitive description terms `ATP`, `A6L` or `genome`. This selects 1,379 parent records containing 3,120,778 nucleotide bases. That screen was chosen from metadata and retrieval cost before any test prediction. It is an ascertainment rule: eligible ATP8 in records without those description terms can be missed. Complete CDS are extracted only from explicit source annotations; genomic completeness is not inferred from record length or title.

Source files, download receipts, source and index hashes, accession versions, parent coordinates and available source-feature qualifiers are retained in `data/raw/gap` and each accepted record's provenance. Dataset construction neither runs nor ranks models.

## Exclusion, dependence and grouping

Exact RNA SHA-256 duplicates of any of the 533 training-pool sequences are excluded, regardless of organism or accession. Records sharing a stable protein or parent accession with training provenance are also excluded; a changed version or mirrored record is not a new independent source. The exclusion ledger records matching accessions and rejection reasons. New records are globally deduplicated by RNA; contradictory cross-family/code/protein assignments fail rather than being resolved silently.

Exact-sequence exclusion does **not** establish distant homology or independent individuals. The own-species reference gate deliberately permits close annotated variants. Identical amino-acid sequences encoded by different RNA can remain. Per-record novelty annotations include the nearest edit distance among all training CDS in the same family, and separate same-family/species/length nucleotide and peptide Hamming distances. These annotations describe relatedness; no model score selects or rejects a sequence. Unknown sample provenance stays unknown.

Many gap locations, context widths or gap sizes from one CDS are repeated measurements on the same sequence. They are not new biological replicates. RNA hashes identify repeated measurements; peptide hashes identify exact peptide groups. Additional near-duplicate clustering and case-selection rules belong to the separately frozen gap-evaluation protocol. Report sparse strata and zero-support cells explicitly. A single human sequence cannot establish human population accuracy, and support from one family cannot establish performance for all ATP-synthase families.

## Freeze and reproducibility

Run `.venv/Scripts/python.exe scripts/build_gap_data.py --include-outer-parents` to acquire and audit candidates. After code and data review, append `--freeze` to write `data/processed/gap/manifest.json`. Once that manifest exists, changed cohort/configuration/reference/source identity is refused. The compact final cohort remains at `data/processed/gap/cohort.jsonl`, with its full source provenance and novelty descriptors. Rebuilding does not change frozen training files or checkpoints.

This is a sequence-disjoint, same-species test assembled through newly declared retrieval routes. It is not a prospective clinical or population study, an independent-individual cohort, a random vertebrate sample, or proof that every usable public record has been recovered. Gap completion, nucleotide matching and translation matching do not themselves establish molecular function or a correct folded structure.
