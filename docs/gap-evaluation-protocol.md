# Gap evaluation — prospective protocol

## Scope and registration boundary

The owner now authorizes unseen-sequence testing, progressive gaps, prefix-only completion and right-flank reranking, restricted to the **same nine species and seven ATP synthase families** declared in `configs/secondary/biology.json`. No additional organism, model fitting, validation-based tuning, structure prediction or functional assay is part of this stage. The historical training protocols remain unchanged.

This design is informed by already observed training outcomes. It is prospective with respect to the new held-out model outcomes; it is not a claim that those training outcomes were unseen. Publish this protocol, `configs/gap-evaluation.json`, acquisition rules and evaluation implementation before any model scores the biological test set. Freeze the test sequences, exact gaps, checkpoint identities and executable evaluation plan before inference. Software fixtures may be used beforehand and must never become research examples.

Separate four questions:

1. **Data feasibility:** Are any eligible, genuinely untrained CDS recoverable within the approved panel after exclusion of every training sequence?
2. **Predictive probability:** Which frozen models assign greater conditional probability to the same missing nucleotides?
3. **Reconstruction:** How does exact and amino-acid reconstruction change as a known-length gap grows?
4. **Information and decoding:** How do left-context length and use of observed right context change results under fixed decoding rules?

The reference is the observed annotated CDS, not an experimentally validated structure or function. Correct amino acids, an uninterrupted reading frame and agreement with an annotation are separate measurements; none establishes molecular function.

## Test acquisition, exclusions and freeze

Retain the existing biological gates: real accession-versioned complete CDS, permitted taxonomy/family, unambiguous RNA after DNA-to-RNA transcription, complete start/stop, no internal stop or translation exception, at most 100 amino acids, reviewed own-species reference length, at least 95% amino-acid identity and agreement with annotated translation. Preserve parent accessions, coordinates, strand, taxonomy, translation table, raw retrieval receipts and hashes. An added accession, later retrieval date or new database route does not by itself create a new sequence.

Construct the exclusion set from the **union of the original pilot cohort and every frozen secondary cohort, including all 533 CDS in `full_vertebrate`**. Exclude a matching full-RNA SHA-256 globally, even if that sequence was absent from one particular model's 412-CDS view or shared across species. Compare linked source/parent accessions as a second leakage flag. A newly revised record linked to an already used source individual/parent is excluded from the primary panel when that relationship can be established; preserve the reason. Do not claim all unlinked records represent independent individuals. No training record may be reclassified as test because its codon is rare, its gap was hidden at inference, or a particular optimizer may not have sampled it.

Deduplicate eligible test RNA globally. Retain all taxonomic provenance, with one deterministic primary family/taxon stratum per distinct CDS. Report nearest training nucleotide identity within compatible family/taxon/length and exact peptide overlap. These describe similarity, not an assurance against shared ancestry or source overlap.

Freeze all eligible strict holdout in a separate pool. For the bounded primary panel, select at most **96 distinct CDS and 16 per family/taxon stratum**, using SHA-256 ranks with selection seed 20260926 and deterministic round robin across ranked strata. No score influences selection. Keep an exclusion ledger and counts at each gate. If fewer than 96 qualify, use all that satisfy the fixed quota. Missing taxa/families remain missing. If none qualify, report infeasibility and stop biological testing; do not relax the gates or substitute training resubstitution results.

The test manifest records the test/source hashes, exclusion-union hash and counts, QC configuration, selection algorithm, actual family/taxon/code support, protocol/config/implementation hashes and the frozen checkpoint inventory. Once scored, this snapshot is a used test set. Any later method development informed by it must be labeled exploratory and needs a fresh untouched set for confirmation.

## Models and comparative support

Use only completed, integrity-audited **final 2,000-update** checkpoints from the frozen secondary plan. The primary matrix comprises its eight matched arms and seeds 17, 29 and 43. An incomplete or unavailable checkpoint is excluded explicitly; do not choose an earlier or lower-loss checkpoint. The one-seed full-pool/larger-capacity model is a separately labeled engineering extension because both capacity and training data differ.

All primary models receive the same selected test sequences and gap locations. Report species/family results separately, including whether each stratum existed in that model's training view. Applying a human ATP8 model to another species or family is a deliberate extrapolation measurement. Do not present a pooled result dominated by human ATP8 as evidence for every family or species.

The tokenizer comparison is primary on the held-out **human ATP8** subset. Capacity, weighting and position comparisons pair the vertebrate controls with the standard vertebrate model on the identical test support. The three dataset expansions compare the same test targets, while disclosing their limited training-cohort replacement sizes. Training seeds describe model/optimizer variability, not biological replication.

## Gap and context construction

Use zero-based codon coordinates. Let `s` be the first missing codon, `G` the gap length and `M` the number of non-stop codons in the complete CDS. Retain the initiation codon and terminal stop as observed sequence; primary gaps contain only sense positions. Their correct length is supplied to the algorithm. This is a known-length reconstruction task, not discovery of gap length, start, stop or reading frame.

Primary lengths are **1, 2, 4, 8 and 16 codons**. For each eligible CDS, rank starts satisfying `1 <= s` and `s + 16 <= M` using the gap seed, sequence hash and start. Select one start and reuse it for all five lengths. The gaps are nested and their left prefix is identical, making the length comparison directly paired. Never choose a location because a model performs well or badly there. Report omitted short sequences. One selected position per sequence avoids making many correlated gaps look like independent observations.

For the focused context experiments select at most **16 human ATP8 CDS** from the frozen primary panel, without model outcomes. Require a common anchor with `s >= 32` and `s + 16 + 8 <= M`; rank valid starts deterministically and reuse the selected anchor across all conditions. The focused anchor can differ from the primary anchor and is identified separately. If fewer than 16 qualify, use their actual support; if none qualify, omit this panel. Do not silently shorten the requested context.

The left-context grid is **1, 2, 4, 8, 16, 32 and full-prefix codons**. Cropping removes earlier observed codons, while preserving the retained tokens' original nucleotide coordinates in the sinusoidal position buffer. Implement an external inference wrapper; do not modify frozen `model.py`, its state or any training file. The no-position checkpoint keeps its zero position buffer. Resetting cropped context to position zero would combine context ablation with positional misalignment and is not part of this protocol. Absolute position remains information given to position-enabled models; the context experiment is not an ablation of all position information.

Retained/generated prefixes grow naturally during completion; the stated left-context length refers to the observed prefix initially supplied, not a rolling window that discards generated codons. No hidden gap token, separator or new vocabulary class is introduced.

Optional secondary 32-, 48- and 64-codon gaps use the same primary anchor only where `s + G <= M`. Report each smaller support and recompute shorter-gap results on that intersection before making a length curve. These lengths are excluded from the primary context/reranking grid, which may have no support under its flank requirements. A missing long-gap panel is an eligibility result, not a reason to move its anchors after scoring.

## Probability scoring and prefix-only generation

For the true missing RNA `g` and observed left context `l`, the probability endpoint is

`gap_bits_per_base = -log P(g | l) / (3 * G * log(2))`.

Evaluate this with teacher forcing: later true gap tokens may condition on earlier true gap tokens only for this probability calculation. Never feed true gap tokens into free-running reconstruction. The base model's joint probability multiplies the three conditional base probabilities for each codon; the codon model supplies one probability for the triplet. Both score exactly the same biological target and condition on the same observed bases. Raw token perplexity and token accuracy are not cross-vocabulary endpoints.

Primary free-running generation uses **codon-synchronous greedy decoding**. At each step obtain all 64 next-codon joint log probabilities. For a base model, enumerate the four possible first bases and their possible second bases, then compute third-base probabilities; the 64 joint probabilities must sum to one within numerical tolerance. Select the highest-probability triplet, append that complete triplet, and repeat for exactly `G` steps. This is identical search granularity across vocabularies; it does not imply identical compute or parameter count. Resolve ties by the frozen codon vocabulary order. Temperature is one, with no top-k/top-p sampling, stop suppression or genetic-code mask.

All 64 codons remain available. An internal stop is recorded as an outcome; do not repair, discard or resample it. The generated sequence is precisely `3*G` nucleotides even after a predicted stop. Models run in evaluation mode under inference/no-gradient execution. No model or optimizer state changes. A KV cache is allowed only after cached and uncached distributions, absolute offsets and causal behavior pass offline equivalence tests; record the inference implementation hash.

## Right-flank reranking

The causal models cannot natively attend to a suffix across an empty gap. The authorized right-flank method is **candidate generation followed by reranking**, not a newly trained bidirectional model.

Use the focused 16-sequence panel, its fixed anchors, the five primary gap lengths, observed left lengths **8 and full**, and right lengths **0, 1, 2, 4 and 8 codons**. Evaluate only human ATP8/base, human ATP8/codon and standard vertebrate/codon, with all three training seeds. The full-prefix left-only context grid uses those same three arms.

For each model/seed/sequence/anchor/gap/left-context combination, construct a deterministic **beam of width eight** using the exact 64-way next-codon distributions. Retain the eight prefixes with highest summed gap log probability at each depth; break ties lexicographically in the frozen vocabulary. The candidate generator receives the left context and known gap length only. It receives neither the true gap, suffix, translation-derived answer nor candidate correctness. Freeze/hash its final candidate bank before scoring any suffix. Candidate banks differ between models, which limits direct interpretation of cross-model reranking results.

Reuse that identical candidate bank for every right-context length. Select the candidate maximizing

`log P(candidate | left) + log P(observed_right | left, candidate)`.

For a fixed cell every candidate and suffix has the same length, so no length normalization, tuned mixing coefficient or amino-acid filter is applied. The `right=0` result selects the highest-scoring candidate from the same beam and is the matched control. It can differ from primary greedy completion because the search differs. Report these separately rather than attributing that difference to right context.

For positive right length, only the observed suffix following the gap is supplied. Scoring retains original coordinates and uses generated candidate tokens, never the true hidden gap. This score is an unnormalized joint score useful for ranking a restricted candidate set; it is **not** a normalized probability of the gap given both flanks. If displaying a softmax over eight candidates, label it candidate-set-normalized confidence, not full posterior confidence.

After the bank is frozen, measure whether it contains the exact true gap and the best attainable nucleotide/amino-acid reconstruction. This gives a candidate-recall ceiling. Never inject the truth into the beam for scoring or reranking. Report changes from `right=0` on the same cases/bank, including harms as well as improvements.

## Baselines and metrics

Fit baseline counts only from each model's own frozen training cohort; never from test targets. The baseline hierarchy is uniform 64-codon, pooled codon unigram/bigram, absolute codon-position categorical, and privileged family-plus-taxon-plus-position categorical. Use fixed total Dirichlet prior mass two over 64 codons. A missing privileged stratum/position falls back to the pooled position distribution, and missing pooled position to pooled unigram. An absent bigram predecessor falls back to pooled unigram. Record fallbacks. The privileged baseline receives family and taxon labels the transformer does not receive. Use codon probabilities for both tokenizer arms so its joint gap score is identical when their training data are identical. Uniform scoring is exactly two bits/base.

Position baselines are essential because the training pilot showed strong position memorization. Score and greedily complete the same gaps with the same coordinate and information rules. Right-independent categorical baselines cannot gain suffix information; show them as fixed controls, not a reranking algorithm that somehow learns from the right flank. An optional nearest-training-template baseline must be registered separately before scoring, including its permitted family labels and flank-distance rule; it is not required by this initial matrix.

Report at minimum:

- Teacher-forced gap bits/base, averaged per sequence and by gap length; a token-weighted pooled value may accompany it.
- Exact nucleotide gap match, nucleotide accuracy and whole-codon accuracy from free-running generation.
- Codon error rate by offset within the generated gap, with support counts at every offset, to expose error accumulation.
- Exact amino-acid gap match and residue accuracy using each record's annotated genetic code; internal stop count/rate and stop-free reconstruction rate.
- For each gap/context, distinct CDS, families, taxa, genetic codes, gaps, target bases and available model seeds; no pooled percentage without its denominator.
- Paired changes versus the relevant model/baseline or `right=0`, with each training seed shown.

Translate internal gap codons with the appropriate [NCBI genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi): table 2 for vertebrate mitochondrial ATP8 and table 1 for these nuclear families. In particular, UGA is not an internal stop in table 2; AGA/AGG have different assignments from table 1. Do not apply initiation-codon remapping to an internal gap. Genetic-code-valid output is not evidence of correct function.

For conservation diagnostics, define masks **from the complete frozen training union**, without test scores: within family/taxon/position, conserved versus variable codons and target-allele frequency. Use these common masks for paired model comparisons. Separate majority alleles (including ties for the largest count), observed nonmajority alleles, train-unseen alleles, and observed alleles at frequency at most 1%; the rare-frequency flag may overlap the majority flag in an unusually diverse stratum. Do not interpret database frequency as population frequency. A stratum absent from the reference training union is `unsupported`, not automatically an unseen allele. If every position varies, say that variable-position accuracy duplicates overall support and does not isolate rare alleles. Also show distance-to-training strata and exact peptide overlap.

## Expectations and uncertainty

The prespecified directional expectations are descriptive: exact gap recovery will decrease with longer nested gaps; more observed left context will improve mean gap bits/base; reranking with right context will improve exact recovery relative to its same-bank `right=0` control. Report reversals and non-monotonic changes. No training-fit ranking guarantees held-out ranking, and no best setting is selected from this test matrix. The tokenizer and architectural contrasts are reported as paired held-out differences without asserting that the training winner must win again.

The replication unit is a distinct CDS, with repeated conditions and model seeds paired within it. Summarize per-CDS results before cohort aggregation; report an unweighted mean over supported family/taxon strata alongside the human ATP8 endpoint. Different gaps, contexts and three optimizer seeds are not additional biological samples.

Display individual-seed estimates and their range. Optional 95% descriptive uncertainty bands use 2,000 paired bootstrap draws, with all conditions for a sequence kept together. Build connected components at at least 99% full-CDS nucleotide identity within compatible family/taxon/length strata before scoring, and resample those sequence clusters rather than individual near-duplicate records. Resample within strata, aggregate with the original fixed stratum weights, and average the three paired model seeds inside each draw. Only display such bands when each contributing stratum has at least five clusters; otherwise show estimates/counts/seed ranges and mark bands unavailable. Clusters are a similarity control, not verified independent individuals. No p-values or biological-population confidence claims are planned.

## Bounded execution and reporting

Use a separate **90-minute local inference cap**, independent of the original four-hour training cap. Serialize inference; CPU, four intra-op threads, one inter-op thread, float32, one checkpoint loaded at a time. Perform synthetic-fixture profiling before freezing the executable plan. Do not run substantial evaluation concurrently with training. No cloud compute or external tracking is introduced.

Priority is (1) the primary full-prefix probability and greedy matrix, (2) the focused prefix-context grid, (3) focused beam/reranking, then (4) the optional long-gap and full-pool engineering extensions. Use frozen model/seed and case order, preserve complete case records and identities on interruption, and report unexecuted cells. A runtime limit does not permit selecting favorable cases, dropping errors, changing algorithms or claiming a partial matrix complete. Model outputs, timings and candidate banks remain separate from training outputs. The inference time is neither a new training budget nor part of model-quality scoring.

The compact report should provide a leakage/eligibility flow, a model/test-support table, paired model contrasts, accuracy and bits/base versus gap length, left-context-by-gap heatmaps, right-context gain curves with candidate recall, and code-aware validity/conservation summaries. Use the exact same support for each paired curve; display a separate curve or missing cell when support changes. Show primary and exploratory panels separately. Save per-case machine-readable results, plan/config/source/checkpoint/data hashes, candidate-bank hashes and environment. Retain negative outcomes and failures.

Before scoring, software checks must cover training/test exclusion, deterministic gap selection, nested/matched support, true coordinate offsets, cache equivalence, causal independence from future tokens, exact base-to-codon normalization, independent probability versus generation paths, truth/suffix-blind beam construction, fixed candidate reuse across suffix lengths, known-length output, genetic-code differences, unchanged checkpoint tensors and accurate partial/resume accounting. These are offline implementation checks, never biological test results.
