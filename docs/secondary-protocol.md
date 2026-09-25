# Secondary pilot protocol — fixed before model fitting

The owner authorizes an ATP-synthase-only expansion, five comparative dataset/model views, additional controlled variables, independent agent review, and up to four hours of local training. Validation, unseen-sequence testing, gap filling and structure work remain deferred. Here "pilot test" means a training experiment, not a new held-out test set.

## Biological expansion

Correct the first acquisition's incomplete handling of semicolon-separated IMPI symbols and supplement ENA search results with the complete set of UniProt EMBL protein cross-references for each eligible reference. These are ascertainment corrections, not biological quality relaxations. Preserve every original pilot file and report the corrected coverage explicitly.

**Coverage amendment before any secondary model fit:** the two routes recovered only 203 distinct CDS, whereas ENA's parent-sequence inventory identified approximately 71,000 mitochondrial genome records across the same nine species. Add a third, symmetric retrieval route for all nine taxa: `tax_eq(TAXON) AND organelle="mitochondrion" AND base_count>=14000 AND base_count<=20000`. Download all returned parent records in bounded batches, verify inventory completeness and source identity, and extract only explicitly annotated ATP8 CDS. This length-bounded parent screen is an ascertainment rule, not proof that every parent is a complete genome. Apply the same CDS quality and own-species reference gates, retain original parent coordinates/accession versions, and globally deduplicate after combining routes. Do not infer an ORF from an unannotated genome. Preserve the two-route snapshot and quantify the third route's actual additions. This amendment uses source coverage evidence, with no secondary training outcomes observed.

Use seven functionally related Complex V subunit families: ATP8, ATP5ME/e, ATP5MF/f, ATP5MJ/j, ATP5MK/k, ATP5F1E/epsilon and ATP5MGL/g-like. They share ATP synthase membership; they are not all ATP8 homologs. Human reference membership must trace to the verified IMPI row, resolving compound symbols and synonyms. Cross-organism references must be reviewed UniProt entries mapped through the declared InterPro family identifier. Exclude fragments and ATP5F1EP2, a pseudogene-related entry.

The declared species panel is human, mouse, rat, cattle, pig, Sumatran orangutan, chicken, zebrafish and African clawed frog. Five other mammals provide relatively broad reviewed-reference coverage; the three nonmammals add bird, fish and amphibian representatives. This availability-based panel is fixed before viewing model outcomes; it is not a random or comprehensive sample of vertebrate diversity. In the reviewed <=100-aa reference screen, nonmammalian coverage is ATP8 only. Record missing family/species combinations without inventing sequences.

Retain complete, unambiguous observed CDS with an allowed start, full stop, no internal stop or translation exceptions, <=100 translated amino acids, own-species reference length and >=95% amino-acid identity to that reference. Do not require 95% identity to the human protein. Translation must agree with the source annotation. Apply table 2 to vertebrate mitochondrial ATP8 and table 1 to the nuclear families. Verify source taxonomy. Deduplicate identical RNA globally, preserving all accession/species provenance; a sequence shared with human is assigned to the human selection stratum. Do not create synthetic variants, repair stops, or count duplicated accessions as new data.

## Maximum pools and matched views

Keep all eligible unique CDS in separate human, mammal and vertebrate pools. Dataset storage is streamed/disk-backed where appropriate; filling RAM is not the definition of useful dataset size. The maximal dataset means all eligible records recovered by the declared retrieval routes and species panel, not all sequences that could ever be discovered.

The five primary model views are:

1. Human ATP8, base vocabulary.
2. The identical human ATP8 cohort, codon vocabulary.
3. Human ATP-synthase families, codon vocabulary.
4. Human plus the selected other mammals, ATP-synthase families, codon vocabulary.
5. Those mammals plus the selected nonmammalian vertebrates, codon vocabulary.

Set common distinct-CDS count N to the available human ATP8 count. Fail if any declared matched view cannot fill its quotas. For the three Complex V views, freeze identical family quotas: include the available sparse human nuclear-family sequences, and use ATP8 for the remaining slots. If nuclear sequences alone exceed half of N, distribute at most floor(N/2) slots among those families by deterministic round robin, with the remainder ATP8. This ensures both original-family and added-family representation without changing inclusion after model results.

Select using a separate fixed cohort seed, shared across all optimizer seeds. Rank sequences and the species round-robin order by seeded SHA-256. In mammal and vertebrate views retain a common human subset within each family, approximately one third of its quota. Fill remaining mammal slots by round robin among eligible nonhuman species. For the vertebrate view replace approximately half of those nonhuman mammal slots with nonmammalian records from the same family when available; preserve the shared human subset and remaining mammal selections. Redistribute shortages within the same family deterministically. Require at least two nonhuman mammal species and at least one nonmammalian species in the final corresponding views. Report all actual counts and overlap hashes. Quota-one families remain human-only; nonmammalian expansion affects ATP8 only under the available reviewed references.

Equal-N views are subsets that replace records; they are not strict supersets. They have equal distinct-sequence counts and, under matched training budgets, equal sequence presentations. Natural CDS lengths differ, so target nucleotide counts and computation are measured and reported rather than declared equal. Within-cohort architectural/objective controls use exactly the same sequence draws.

## Experimental variables

Record the following directional expectations before fitting. They are falsifiable optimization expectations at this fixed budget, not statistical hypotheses about a biological population. Report the mean paired difference and each seed; no significance test or model selection uses these outcomes.

- Engineering: every completed run has finite updates and lower final overall training bits/base than initialization. Failure triggers investigation and a disclosed new version, not silent replacement.
- Tokenization: the codon arm is expected to reach lower final bits/base than the base arm on the identical human ATP8 cohort at equal sequence presentations. Reverse ordering contradicts this expectation under this budget.
- Capacity: the 192-wide body is expected to fit the frozen vertebrate cohort better than the original 48-wide body at equal updates. This measures an architecture-and-optimization response, not optimal model size.
- Family weighting: family-macro training is expected to lower the unweighted mean of per-family bits/base relative to token-weighted training on the same cohort. Overall token-weighted fit may move in the opposite direction; report both.
- Explicit position: removing sinusoidal position features is expected to make fitting this position-conserved corpus harder. The no-position arm still has causal order and prefix-length information, so the result cannot isolate every possible positional cue.
- Taxonomic expansion: no raw-loss direction is asserted across different cohorts. Compare descriptive changes in diversity, categorical-baseline fit and model fit with support counts. A harder expanded corpus may have higher loss even when optimization is successful.

The five primary arms use one common hardware-profiled transformer body and seeds 17, 29 and 43. Add three controls on the frozen vertebrate view: the original 48-wide/two-layer/four-head architecture, family-macro loss weighting, and removal of explicit positional encoding. The architecture comparison changes width, depth and head count together; it is a bundled capacity comparison and cannot isolate any one of those dimensions. The weighting and positional controls each change one setting. The latter still exposes causal order/prefix length and is not a claim that all position information has been removed. All arms retain exactly 4 or 64 output symbols; no family/species special tokens are added.

The family-macro loss changes the objective while retaining identical uniform sequence draws. Let T_f be total valid training target tokens in family f, F family count, N sequence count and B total sampled sequences across accumulation. The unbiased estimator is N/(F*B) times the sum of sequence cross-entropy sums divided by the corresponding T_f. Do not normalize by the minibatch's random weight sum.

Choose the common optimizer-update budget from measured base/codon timings and the four-hour aggregate limit before fitting. Use four CPU threads unless profiling supports a change; serialize substantial training processes to avoid memory contention and distorted timings. The arms are parallel scientific comparisons, not competing simultaneous memory allocations. Save the plan, hardware criteria and hashes before launch. Hardware profiling supports a common 192-wide, four-layer body, with 2,000 updates per run and 10% learning-rate warmup. The base arm of a 384-wide/six-layer body was skipped by the memory reserve gate. A separate one-seed run on the full vertebrate pool uses the largest successfully profiled codon body (384-wide/six-layer) and 2,000 updates; it is excluded from equal-size causal comparisons because both capacity and dataset change. These are the largest selected profiled candidates, not a measured absolute memory limit.

## Learning from the first pilot

The first codon runs improved loss but all ended with the same token accuracy; conserved positions dominated. The position baseline outperformed every transformer. Therefore report overall and family-macro bits/base, per-family/per-taxon/per-genetic-code fit, and within-family/species variable-position fit. Variable positions are defined by codon variation; the base arm scores all three bases of exactly those same codons. Masks are cohort-specific, so variable-fit comparisons across different cohorts have different supports and are descriptive; report the support counts. Baselines include pooled position and privileged family-position distributions. Use empirical training entropy plus total prior mass 2, avoiding the original cross-vocabulary smoothing-mass difference. Privileged family labels are given to those baselines, not the transformer, and are labeled explicitly.

Record exact exposure, wall time, peak working memory, finite gradients, initial/final diagnostics and complete run state. Checkpoint retention is one rolling checkpoint plus the completed final checkpoint, not every historical state. Use offline software tests for invariants; do not confuse those with biological validation data.

Training-only differences describe optimization and fit to these observed cohorts. Higher loss on a more diverse cohort is not necessarily worse modeling. Neither lower training loss nor three initialization seeds establish generalization, optimal hyperparameters, preserved function or organism-level statistical significance.
