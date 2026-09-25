# Protocol v0.1 — training feasibility

## Scope and questions

The experimental unit is a **distinct, observed, complete coding sequence** of one annotated human microprotein family. Multiple accession records of the same nucleotide sequence are technical database duplicates for this experiment. They are not additional independent training examples. A distinct sequence is also not guaranteed to represent an independent person or evolutionary lineage.

The owner specified one organism, one functional family, and close protein variants; broader organism coverage is deferred. We operationalize a microprotein as a translated sequence of at most 100 amino acids, excluding the terminal stop. Close variants must match the reviewed reference length and have at least 95% positionwise amino-acid identity. This conservative, explicit threshold was fixed before cohort selection and training; it is not a validated function predictor.

1. **Data feasibility:** How many distinct usable sequences remain after defining the eligible family pool and applying fixed quality checks?
2. **Implementation feasibility:** Can both vocabularies represent the same original sequences exactly, preserve reading frame and sequence boundaries, and support causal next-token training?
3. **Optimization feasibility:** Does a prespecified small model train with finite gradients and decreasing training negative log likelihood across three seeds?
4. **Future predictive hypothesis, deferred:** Does either representation improve performance on unseen biological sequences under a separately approved evaluation design?

## Hypotheses and falsifiability

For the present pilot, the descriptive expectation is that each arm's final training loss is lower than its initialization, without numerical failures. A numerical failure, incorrect reading-frame treatment, irreproducible resume, data hash mismatch, or inability to outperform its own initialization falsifies the corresponding engineering expectation. These outcomes justify debugging, not changing biological inclusion rules until results look favorable.

No generalization null hypothesis is tested in this phase. No significance test, confidence interval across biological samples, best model, best hyperparameter or functional-accuracy claim will be reported. Three initialization seeds characterize optimizer variability; they are not independent biological replicates.

## Comparison design

| Controlled property | Rule |
|---|---|
| Biological data | Identical ordered cohort, RNA sequences and source hashes |
| Unit weighting | Uniform sampling over distinct CDS; duplicates carry no extra weight |
| Predicted region | All sequence after the first codon, including the complete stop codon |
| Conditioning | The first complete codon supplied in both arms |
| Context | Complete CDS, within a configured 303-nucleotide capacity; no truncation |
| Architecture | Same transformer width, depth, heads, dropout and fixed positional encoding |
| Optimization | Same optimizer, learning-rate schedule, batches, accumulation and update count |
| Sampling | Dedicated seeded sequence-index generator, shared sampling rule across arms |
| Reporting | Full training-corpus diagnostics at initialization and fixed intervals |

Base targets corresponding to nucleotides 2 and 3 are ignored so both arms score nucleotides 4 through 207. The base model factorizes a codon's probability into three conditional base probabilities; the codon model predicts the entire next triplet. No token represents a gap, unknown base, sequence separator, beginning, end or padding. Fill values used to form a rectangular batch are masked from the loss. Causal attention prevents right-side fill values from influencing earlier sequence positions.

Vocabulary size is 4 or 64. Embedding width is an independent architectural choice, not a synonym for vocabulary size. All 64 codon classes exist even when some have zero observed training counts. DNA T is transcribed to RNA U only after biological checks; reverse complementation and random frame shifts are not augmentations.

The two arms have the same transformer body but slightly different embedding/output parameter counts (weight-tied, 60 additional rows in the codon arm). Positional encodings are fixed and indexed in nucleotide coordinates. Equal sequence exposure is the primary compute control; attention cost and elapsed time differ. This pilot does not claim a perfectly parameter- and FLOP-matched comparison.

## Metrics and controls

Loss uses summed cross entropy divided by the actual count of valid target tokens, including across accumulation microbatches. Schedule steps are optimizer updates. Report bits per nucleotide as `total negative log likelihood / (ln(2) * number of predicted nucleotides)`. This permits a common unit; raw token perplexities and base/codon token accuracies are not directly comparable.

Uniform, training-fit unigram, bigram and per-position categorical baselines reveal how much of the training distribution can be explained without a transformer. Smoothing is fixed at 0.5 per category. Positional baselines have different support sizes across arms; their smoothed losses reflect that. All these are **resubstitution** metrics: the same training data fit and score the baselines.

Save every scheduled checkpoint and the final state; never choose a checkpoint by minimum training loss and call it the best predictor. Persist config, seed, data hash, source hashes, environment, optimizer, RNG states, exposures and training diagnostics. CPU resume is checked against uninterrupted training with exact tensor equality. Cross-device/library bitwise equality is not promised ([PyTorch reproducibility notes](https://docs.pytorch.org/docs/2.14/notes/randomness.html)).

## Explicitly deferred

No validation or test partition is created, no unseen sequence is scored, and no gap-completion experiment is implemented. The current left-to-right objective conditions on a prefix; using both sides of a missing gap is a future objective-design decision. This entire snapshot has been available to training and must not later be relabeled an untouched test set.

No AlphaFold job, generated protein folding, PyMOL comparison, structure metric or structure benchmark is implemented. Predicted structure similarity alone would not establish molecular function, nor would an AlphaFold prediction of the reference be experimental ground truth. Those questions remain for a later approved phase.
