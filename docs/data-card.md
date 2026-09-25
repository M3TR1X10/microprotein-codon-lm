# Dataset card — human MT-ATP8 pilot

## Source roles

- **Eligibility:** [IMPI 2021-Q4pre workbook](https://www.mrc-mbu.cam.ac.uk/files/impi-2021-q4pre-20211001-dist_0.xlsx), linked from the [requested IMPI page](https://www.mrc-mbu.cam.ac.uk/research-resources-and-facilities/impi). It supplies human Ensembl identifiers, gene names and mitochondrial evidence. Its 27,366 body rows include entries outside the mitochondrial class and blank rows; membership in the workbook alone is not eligibility.
- **Reference:** Reviewed human UniProt entries with sequence length 1–100 aa. Gene and product annotations are mapped to IMPI records classified `Verified mitochondrial`, including that class's explanatory suffixes.
- **Original nucleotides:** [ENA coding records](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/browser-api.html). Query `tax_eq(9606) AND base_count<=303`; retrieve every metadata record with `limit=0`, and verify the returned count against ENA's count endpoint.

The experiment is **IMPI-guided, ENA-sourced**. It cannot honestly be described as containing mostly sequences downloaded from IMPI, which provides no coding-sequence column. No amino-acid sequence is reverse-translated. The selected reference is [UniProt P03928](https://www.uniprot.org/uniprotkb/P03928/entry).

## Candidate ascertainment

The frozen screen contains 61 reviewed small human reference proteins linked to verified IMPI genes. Match ENA gene symbols, gene synonyms and exact normalized reference product names; normalization removes punctuation and letter case, not meaningful words. Two explicit product aliases are recorded in `configs/product_aliases.json`: the ENA annotation `ATP synthase F0 subunit 8` is linked to ATP8 by record WAB69273.1; the ND4L product name appears with the ND4L gene in the retrieved metadata. This closes an observed naming gap without substring matching unrelated ATP synthase subunits.

For alternative short products from a larger gene (for example alternative MIEF1 or ND4 reading frames), a gene-name match alone is insufficient: complete translation and reference length/identity checks must also pass.

Coverage is limited to records mapped by this declared procedure. Unmapped records, unreviewed reference proteins, unknown microproteins, other organisms, and valid alleles not deposited or indexed in this ENA snapshot are outside the selection universe. A zero count means no mapped usable record under this procedure, not biological absence. This is a count-based winner within a reproducible screen, not proof of the globally largest possible family in every database.

## Quality rules

Apply the same gates to all candidates before selecting the winner:

1. Human taxonomy 9606 and one CDS; record and protein identifiers must agree.
2. No pseudogene, translation exception, programmed frameshift, ribosomal slippage, or artificial-location flag. Require codon start 1 and exact, complete feature endpoints.
3. Use the coding-oriented, spliced sequence returned by ENA's **coding accession** endpoint. Its feature coordinates refer to the parent nucleotide record, so applying those coordinates to the returned CDS again would be incorrect.
4. Require unambiguous A/C/G/T, a length divisible by three, and agreement with length and nucleotide MD5 metadata.
5. Translate a complete CDS with the appropriate NCBI code: table 2 for the human mitochondrial genes, table 1 for the nuclear genes. Require an allowed initiator, a terminal stop, no internal stop, and agreement with the deposited protein translation. Do not append missing stops or repair incomplete mitochondrial stops.
6. Require at most 100 amino acids, the same reference length, and at least 95% amino-acid identity to the reviewed reference. The terminal stop is retained in RNA but omitted from the protein sequence.
7. Deduplicate the resulting RNA exactly; preserve accession aliases and source multiplicity without using them as training weights.

Sequence-identity and CDS integrity establish membership/quality proxies, not functional assays. Disease-associated or function-altering close variants can still be included because their functional status is not established by these filters. That limitation is consistent with the owner's choice to include close annotated variants.

The genetic code matters: UGA codes for tryptophan in vertebrate mitochondria, and AGA/AGG are assigned as stops in that table ([NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)). All codons remain vocabulary classes regardless of their translation semantics.

## Selection outcome

| Candidate | Mapped records | Distinct CDS before QC | Usable distinct CDS |
|---|---:|---:|---:|
| MT-ATP8 | 5,318 | 137 | 129 |
| MT-ND4L | 5,318 | 112 | 109 |
| Other mapped eligible families | See ranking | See ranking | At most 1 per family |

The MT-ATP8 metadata gate excludes one record with a partial reading frame before downloads; seven of the remaining distinct sequences contain ambiguous nucleotides. The final 129 RNA sequences encode 53 distinct 68-aa proteins. Each RNA is 207 nt: 204 protein-coding bases plus one complete stop codon. Both models condition on the first 3 nt and learn targets covering 204 nt per example.

The maximum usable count is the only selection score. Ties use gene symbol and then UniProt accession alphabetically. No structure quality, model loss, expected ease of prediction or preferred gene function enters the selection. Full counts, rejection summaries, raw-query coverage and source hashes are in `reports/candidate-ranking.json`, `reports/inventory-coverage.json`, `reports/data-manifest.json` and `reports/source-registry.json`.

## Reproducibility and limitations

Each download has a URL, UTC timestamp, SHA-256, size and available server release/ETag. Cache mismatches fail closed; failed requests are retried and then fail visibly. Selection never turns a network error into a zero count. Query totals are checked to prevent silent truncation. Protein records sharing an identical nucleotide hash are represented once after metadata quality screening; source-accession lists are retained. Where identical sequences have conflicting annotations, this representative policy may be conservative; it does not infer independent samples from duplicates.

This is a highly related single-gene human cohort, not a broad sample of genomic diversity. Uniform weighting targets the distribution of observed distinct variants rather than human population frequencies. No sample independence, healthy status, haplogroup balance, geographic representativeness, or equivalent variant function is claimed. IMPI's 2021 annotations and current UniProt/ENA content also have different release dates.

Raw source files and checkpoints remain local and ignored by Git. Source inventories and hashes are versioned so an authorized researcher can reconstruct and audit the input. Upstream data terms remain applicable; retaining the code scaffold's MIT notice does not relicense third-party biological data.
