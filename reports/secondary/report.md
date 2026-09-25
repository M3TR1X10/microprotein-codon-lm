# Secondary pilot — training only

TRAINING ONLY; all diagnostics are resubstitution, with no held-out observations.

**25 of 25 planned runs completed their fixed update budgets.** 0 recorded runs are partial; 0 have no reported result. Partial runs are excluded from completed-run averages. A pending seed is not a failed scientific hypothesis.

The aggregate execution budget is 4 hours. Budget-chargeable execution time is 222.1 minutes. Run-loop timings below include interval diagnostics and checkpoints but exclude setup/initial diagnostics; scheduler time covers the broader execution.

Scheduler status: **complete**. Every planned run completed its update budget.

A documented lid-triggered standby interruption stopped the original execution after eight complete runs. The [post-launch runtime amendment](../../docs/secondary-runtime-amendment.md) excludes only 262.00 minutes of Windows-recorded hardware deep idle, once. All other time remains charged. Including that preserved exclusion, recorded execution elapsed is 484.2 minutes. Budget-chargeable time is conservative accounting, not exact CPU-compute time. Model settings, data and update budgets did not change.

The seed-29 base run retained its checkpoint and extra interruption diagnostic. Its model outcome remains in the matched contrasts; its interrupted timing is excluded from hardware-speed comparisons. Raw loop times below retain standby; the separate charged column applies the single documented exclusion.

## Data and controlled comparisons

| Frozen view | Distinct CDS | Distinct peptides | Families | Taxon strata | Target bases/pass | Variable target bases | Target codons observed / 64 |
|---|---:|---:|---:|---:|---:|---:|---:|
| full_human | 422 | 185 | 7 | 1 | 86,262 | 84,096 | 64 |
| full_mammal | 524 | 244 | 7 | 6 | 106,884 | 88,803 | 64 |
| full_vertebrate | 533 | 250 | 7 | 9 | 108,345 | 88,929 | 64 |
| human_atp8 | 412 | 178 | 1 | 1 | 84,048 | 84,048 | 57 |
| human_complex | 412 | 182 | 7 | 1 | 84,222 | 82,056 | 64 |
| mammal_complex | 412 | 201 | 7 | 6 | 83,859 | 70,125 | 64 |
| vertebrate_complex | 412 | 206 | 7 | 9 | 83,532 | 69,432 | 64 |

The family expansion replaces **10/412 CDS (2.43%)**; the mammalian expansion replaces **85/412 (20.63%)**; the nonmammalian expansion replaces **9/412 (2.18%)**. Nonmammalian additions are ATP8 only. Small replacement fractions limit what aggregate losses reveal about broader diversity.

The five primary arms vary tokenization or cohort composition. The three additional controls change capacity, family-macro loss weighting, or explicit positional encoding on exactly the same vertebrate cohort. Uniform sequence draws and matched update budgets control presentations; natural lengths make nucleotide exposure differ. The larger full-pool run changes both capacity and data and is an engineering demonstration.

## Completed matched runs

Values are mean ± sample standard deviation across completed initialization seeds, not confidence intervals or biological replicates. The completed/planned column exposes missing seeds. No model is selected by minimum training loss.

| Arm | Completed/planned | Parameters | Updates | Initial bits/base | Final bits/base | Final family-macro bits/base | Final variable-codon bits/base | Variable target bases |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Human ATP8 / bases | 3/3 | 1,780,608 | 2,000 | 2.0323 ± 0.0760 | 0.0632 ± 0.0001 | 0.0632 ± 0.0001 | 0.0632 ± 0.0001 | 84,048 |
| Human ATP8 / codons | 3/3 | 1,792,128, 1,793,668 | 2,000 | 1.9978 ± 0.0136 | 0.0630 ± 0.0000 | 0.0630 ± 0.0000 | 0.0630 ± 0.0000 | 84,048 |
| Human Complex V / codons | 3/3 | 1,792,128, 1,793,668 | 2,000 | 1.9984 ± 0.0131 | 0.0615 ± 0.0015 | 0.0574 ± 0.0036 | 0.0619 ± 0.0016 | 82,056 |
| Mammal Complex V / codons | 3/3 | 1,792,128, 1,793,668 | 2,000 | 2.0019 ± 0.0132 | 0.0598 ± 0.0003 | 0.0624 ± 0.0049 | 0.0647 ± 0.0003 | 70,125 |
| Vertebrate Complex V / codons | 3/3 | 1,792,128, 1,793,668 | 2,000 | 2.0020 ± 0.0127 | 0.0599 ± 0.0004 | 0.0622 ± 0.0028 | 0.0648 ± 0.0004 | 69,432 |
| Vertebrate / small capacity | 3/3 | 59,712, 60,100 | 2,000 | 2.0071 ± 0.0048 | 0.1216 ± 0.0014 | 1.3955 ± 0.0380 | 0.0786 ± 0.0005 | 69,432 |
| Vertebrate / family-macro objective | 3/3 | 1,792,128, 1,793,668 | 2,000 | 2.0020 ± 0.0127 | 0.1205 ± 0.0007 | 0.0469 ± 0.0037 | 0.1012 ± 0.0001 | 69,432 |
| Vertebrate / no explicit positions | 3/3 | 1,792,128, 1,793,668 | 2,000 | 2.0228 ± 0.0299 | 0.0579 ± 0.0001 | 0.0716 ± 0.0019 | 0.0619 ± 0.0001 | 69,432 |

![Training curves](curves.svg)

Every plotted curve represents one initialization. Dashed family-position baselines receive privileged family labels. There is no held-out curve. Base/codon token accuracy and raw token perplexity are not compared.

Variable-codon masks are cohort-specific. The base/codon pair shares the same biological support; the three vertebrate controls share their reference arm's support. Across different cohorts, changing support prevents a paired accuracy interpretation. This mask counts every observation at a variable position, not only minority codons. Every target codon position varies at least once in `human_atp8`; variable-position and overall fit are identical there, so this diagnostic does not resolve rare-variant performance.

## Controlled contrasts and prespecified expectations

Differences follow the subtraction order in each row and are paired by initialization seed. Positive values mean higher training negative log likelihood for that metric. These are descriptive optimization contrasts, not significance tests.

Met/reversed describes only the frozen directional training expectation evaluated on the mean paired difference. Missing seed pairs make that interpretation provisional. No direction was asserted for raw loss across different taxonomic cohorts. The capacity contrast changes width, depth and heads together; it is not a single-dimension intervention.

| Contrast | Complete/planned pairs | Overall delta bits/base | Family-macro delta bits/base | Variable-codon delta bits/base | Prespecified expectation | Mean-pair outcome |
|---|---:|---:|---:|---:|---|---|
| Tokenization: codon minus base | 3/3 | -0.000161 ± 0.000081 | -0.000161 ± 0.000081 | -0.000161 ± 0.000081 | Codon lower overall bits/base | Met |
| Capacity: small minus standard | 3/3 | 0.061708 ± 0.001720 | 1.333334 ± 0.035806 | 0.013796 ± 0.000906 | Small architecture higher overall bits/base | Met |
| Family weighting: macro minus token | 3/3 | 0.060582 ± 0.000594 | -0.015300 ± 0.001695 | 0.036391 ± 0.000584 | Family-macro objective lower family-macro bits/base | Met |
| Explicit position: removed minus supplied | 3/3 | -0.001985 ± 0.000317 | 0.009394 ± 0.002499 | -0.002925 ± 0.000534 | Removed positions higher overall bits/base | Reversed |

Each matched pair below has verified identical cohort, sampling trace, updates, sequence/base presentations, family exposures and scored supports. Variable support is comparable within each pair only.

| Contrast | Seed | Overall delta bits/base | Family-macro delta bits/base | Variable-codon delta bits/base | Variable target bases | Sequence presentations | Target bases presented | Seed expectation |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Tokenization: codon minus base | 17 | -0.000218 | -0.000218 | -0.000218 | 84,048 | 32,000 | 6,528,000 | Met |
| Tokenization: codon minus base | 29 | -0.000068 | -0.000068 | -0.000068 | 84,048 | 32,000 | 6,528,000 | Met |
| Tokenization: codon minus base | 43 | -0.000197 | -0.000197 | -0.000197 | 84,048 | 32,000 | 6,528,000 | Met |
| Capacity: small minus standard | 17 | 0.061231 | 1.337671 | 0.013271 | 69,432 | 32,000 | 6,488,652 | Met |
| Capacity: small minus standard | 29 | 0.063616 | 1.366773 | 0.014842 | 69,432 | 32,000 | 6,490,143 | Met |
| Capacity: small minus standard | 43 | 0.060278 | 1.295557 | 0.013274 | 69,432 | 32,000 | 6,485,778 | Met |
| Family weighting: macro minus token | 17 | 0.061259 | -0.014070 | 0.036261 | 69,432 | 32,000 | 6,488,652 | Met |
| Family weighting: macro minus token | 29 | 0.060346 | -0.014597 | 0.037029 | 69,432 | 32,000 | 6,490,143 | Met |
| Family weighting: macro minus token | 43 | 0.060142 | -0.017233 | 0.035883 | 69,432 | 32,000 | 6,485,778 | Met |
| Explicit position: removed minus supplied | 17 | -0.002071 | 0.008943 | -0.002908 | 69,432 | 32,000 | 6,488,652 | Reversed |
| Explicit position: removed minus supplied | 29 | -0.001635 | 0.007151 | -0.002399 | 69,432 | 32,000 | 6,490,143 | Reversed |
| Explicit position: removed minus supplied | 43 | -0.002251 | 0.012088 | -0.003467 | 69,432 | 32,000 | 6,485,778 | Reversed |

## Training-fit baselines

Each distribution is fit and scored on the same training observations. Empirical rows use no smoothing. Smoothed rows use total prior mass 2 per categorical distribution (2/4 per base category, 2/64 per codon category). Equal total prior removes the first pilot's prior-mass imbalance; different factorizations and conditional contexts still differ. The family-position baseline is privileged: the transformer receives no family or organism token. Uniform entropy is 2 bits/base.

| Cohort / tokens | Baseline | Empirical bits/base | Prior-mass-2 bits/base | Prior-mass-2 family-macro bits/base |
|---|---|---:|---:|---:|
| full_vertebrate:codon | unigram | 1.7012 | 1.7012 | 2.2739 |
| full_vertebrate:codon | bigram | 0.7902 | 0.7913 | 1.8084 |
| full_vertebrate:codon | position | 0.4034 | 0.4054 | 1.9155 |
| full_vertebrate:codon | family_position (privileged) | 0.2705 | 0.2809 | 0.3890 |
| human_atp8:base | unigram | 1.7810 | 1.7810 | 1.7810 |
| human_atp8:base | bigram | 1.7452 | 1.7452 | 1.7452 |
| human_atp8:base | position | 0.0624 | 0.0660 | 0.0660 |
| human_atp8:base | family_position (privileged) | 0.0624 | 0.0660 | 0.0660 |
| human_atp8:codon | unigram | 1.6039 | 1.6039 | 1.6039 |
| human_atp8:codon | bigram | 0.5021 | 0.5036 | 0.5036 |
| human_atp8:codon | position | 0.0622 | 0.0643 | 0.0643 |
| human_atp8:codon | family_position (privileged) | 0.0622 | 0.0643 | 0.0643 |
| human_complex:codon | unigram | 1.6331 | 1.6331 | 2.4324 |
| human_complex:codon | bigram | 0.5454 | 0.5470 | 1.7953 |
| human_complex:codon | position | 0.1269 | 0.1297 | 2.0417 |
| human_complex:codon | family_position (privileged) | 0.0603 | 0.0716 | 0.3601 |
| mammal_complex:codon | unigram | 1.6767 | 1.6767 | 2.3640 |
| mammal_complex:codon | bigram | 0.7498 | 0.7513 | 1.8994 |
| mammal_complex:codon | position | 0.3338 | 0.3365 | 2.1255 |
| mammal_complex:codon | family_position (privileged) | 0.2681 | 0.2793 | 0.4692 |
| vertebrate_complex:codon | unigram | 1.6752 | 1.6752 | 2.3630 |
| vertebrate_complex:codon | bigram | 0.7637 | 0.7652 | 1.8960 |
| vertebrate_complex:codon | position | 0.3524 | 0.3550 | 2.1081 |
| vertebrate_complex:codon | family_position (privileged) | 0.2872 | 0.2984 | 0.4720 |

## Per-family training fit

Completed matched runs only; support is the number of target bases per full training-corpus pass, not multiplied by initialization seeds.

| Arm | Family | Initializations | Target bases | Final bits/base | Variable target bases | Final variable-codon bits/base |
|---|---|---:|---:|---:|---:|---:|
| Human ATP8 / bases | ATP8 | 3 | 84,048 | 0.0632 ± 0.0001 | 84,048 | 0.0632 ± 0.0001 |
| Human ATP8 / codons | ATP8 | 3 | 84,048 | 0.0630 ± 0.0000 | 84,048 | 0.0630 ± 0.0000 |
| Human Complex V / codons | ATP5F1E | 3 | 306 | 0.0625 ± 0.0038 | 6 | 0.3474 ± 0.0092 |
| Human Complex V / codons | ATP5ME | 3 | 414 | 0.0500 ± 0.0039 | 6 | 0.3651 ± 0.0024 |
| Human Complex V / codons | ATP5MF | 3 | 846 | 0.0413 ± 0.0033 | 36 | 0.2316 ± 0.0135 |
| Human Complex V / codons | ATP5MGL | 3 | 300 | 0.0483 ± 0.0046 | 0 | n/a |
| Human Complex V / codons | ATP5MJ | 3 | 174 | 0.0700 ± 0.0027 | 0 | n/a |
| Human Complex V / codons | ATP5MK | 3 | 174 | 0.0680 ± 0.0115 | 0 | n/a |
| Human Complex V / codons | ATP8 | 3 | 82,008 | 0.0618 ± 0.0016 | 82,008 | 0.0618 ± 0.0016 |
| Mammal Complex V / codons | ATP5F1E | 3 | 306 | 0.0826 ± 0.0136 | 0 | n/a |
| Mammal Complex V / codons | ATP5ME | 3 | 420 | 0.0714 ± 0.0155 | 0 | n/a |
| Mammal Complex V / codons | ATP5MF | 3 | 810 | 0.0507 ± 0.0055 | 0 | n/a |
| Mammal Complex V / codons | ATP5MGL | 3 | 300 | 0.0440 ± 0.0067 | 0 | n/a |
| Mammal Complex V / codons | ATP5MJ | 3 | 174 | 0.0693 ± 0.0118 | 0 | n/a |
| Mammal Complex V / codons | ATP5MK | 3 | 174 | 0.0594 ± 0.0039 | 0 | n/a |
| Mammal Complex V / codons | ATP8 | 3 | 81,675 | 0.0598 ± 0.0002 | 70,125 | 0.0647 ± 0.0003 |
| Vertebrate Complex V / codons | ATP5F1E | 3 | 306 | 0.0921 ± 0.0080 | 0 | n/a |
| Vertebrate Complex V / codons | ATP5ME | 3 | 420 | 0.0819 ± 0.0205 | 0 | n/a |
| Vertebrate Complex V / codons | ATP5MF | 3 | 810 | 0.0473 ± 0.0019 | 0 | n/a |
| Vertebrate Complex V / codons | ATP5MGL | 3 | 300 | 0.0382 ± 0.0001 | 0 | n/a |
| Vertebrate Complex V / codons | ATP5MJ | 3 | 174 | 0.0590 ± 0.0006 | 0 | n/a |
| Vertebrate Complex V / codons | ATP5MK | 3 | 174 | 0.0570 ± 0.0016 | 0 | n/a |
| Vertebrate Complex V / codons | ATP8 | 3 | 81,348 | 0.0599 ± 0.0005 | 69,432 | 0.0648 ± 0.0004 |
| Vertebrate / small capacity | ATP5F1E | 3 | 306 | 1.5282 ± 0.0737 | 0 | n/a |
| Vertebrate / small capacity | ATP5ME | 3 | 420 | 1.4928 ± 0.0457 | 0 | n/a |
| Vertebrate / small capacity | ATP5MF | 3 | 810 | 1.5098 ± 0.0838 | 0 | n/a |
| Vertebrate / small capacity | ATP5MGL | 3 | 300 | 1.7934 ± 0.1108 | 0 | n/a |
| Vertebrate / small capacity | ATP5MJ | 3 | 174 | 1.6115 ± 0.0840 | 0 | n/a |
| Vertebrate / small capacity | ATP5MK | 3 | 174 | 1.7505 ± 0.0969 | 0 | n/a |
| Vertebrate / small capacity | ATP8 | 3 | 81,348 | 0.0826 ± 0.0010 | 69,432 | 0.0786 ± 0.0005 |
| Vertebrate / family-macro objective | ATP5F1E | 3 | 306 | 0.0405 ± 0.0032 | 0 | n/a |
| Vertebrate / family-macro objective | ATP5ME | 3 | 420 | 0.0609 ± 0.0294 | 0 | n/a |
| Vertebrate / family-macro objective | ATP5MF | 3 | 810 | 0.0301 ± 0.0012 | 0 | n/a |
| Vertebrate / family-macro objective | ATP5MGL | 3 | 300 | 0.0240 ± 0.0081 | 0 | n/a |
| Vertebrate / family-macro objective | ATP5MJ | 3 | 174 | 0.0222 ± 0.0028 | 0 | n/a |
| Vertebrate / family-macro objective | ATP5MK | 3 | 174 | 0.0278 ± 0.0071 | 0 | n/a |
| Vertebrate / family-macro objective | ATP8 | 3 | 81,348 | 0.1228 ± 0.0008 | 69,432 | 0.1012 ± 0.0001 |
| Vertebrate / no explicit positions | ATP5F1E | 3 | 306 | 0.0767 ± 0.0060 | 0 | n/a |
| Vertebrate / no explicit positions | ATP5ME | 3 | 420 | 0.0805 ± 0.0033 | 0 | n/a |
| Vertebrate / no explicit positions | ATP5MF | 3 | 810 | 0.0632 ± 0.0028 | 0 | n/a |
| Vertebrate / no explicit positions | ATP5MGL | 3 | 300 | 0.0685 ± 0.0036 | 0 | n/a |
| Vertebrate / no explicit positions | ATP5MJ | 3 | 174 | 0.0827 ± 0.0064 | 0 | n/a |
| Vertebrate / no explicit positions | ATP5MK | 3 | 174 | 0.0720 ± 0.0127 | 0 | n/a |
| Vertebrate / no explicit positions | ATP8 | 3 | 81,348 | 0.0576 ± 0.0000 | 69,432 | 0.0619 ± 0.0001 |

## Per-taxon training fit

Completed matched runs only. Support counts target bases per corpus pass, not independent individuals. Taxa are the frozen primary selection strata; genetic codes describe translation conventions, not added model tokens.

| Arm | Stratum | Initializations | Target bases | Final bits/base |
|---|---|---:|---:|---:|
| Human ATP8 / bases | Homo sapiens (9606) | 3 | 84,048 | 0.0632 ± 0.0001 |
| Human ATP8 / codons | Homo sapiens (9606) | 3 | 84,048 | 0.0630 ± 0.0000 |
| Human Complex V / codons | Homo sapiens (9606) | 3 | 84,222 | 0.0615 ± 0.0015 |
| Mammal Complex V / codons | Mus musculus (10090) | 3 | 2,412 | 0.0609 ± 0.0011 |
| Mammal Complex V / codons | Rattus norvegicus (10116) | 3 | 1,608 | 0.0543 ± 0.0026 |
| Mammal Complex V / codons | Pongo abelii (9601) | 3 | 408 | 0.0850 ± 0.0048 |
| Mammal Complex V / codons | Homo sapiens (9606) | 3 | 66,774 | 0.0596 ± 0.0004 |
| Mammal Complex V / codons | Sus scrofa (9823) | 3 | 5,691 | 0.0587 ± 0.0036 |
| Mammal Complex V / codons | Bos taurus (9913) | 3 | 6,966 | 0.0616 ± 0.0013 |
| Vertebrate Complex V / codons | Mus musculus (10090) | 3 | 2,412 | 0.0558 ± 0.0015 |
| Vertebrate Complex V / codons | Rattus norvegicus (10116) | 3 | 1,608 | 0.0549 ± 0.0026 |
| Vertebrate Complex V / codons | Danio rerio (7955) | 3 | 162 | 0.0604 ± 0.0028 |
| Vertebrate Complex V / codons | Xenopus laevis (8355) | 3 | 165 | 0.0637 ± 0.0054 |
| Vertebrate Complex V / codons | Gallus gallus (9031) | 3 | 1,134 | 0.0615 ± 0.0013 |
| Vertebrate Complex V / codons | Pongo abelii (9601) | 3 | 408 | 0.0822 ± 0.0078 |
| Vertebrate Complex V / codons | Homo sapiens (9606) | 3 | 66,774 | 0.0604 ± 0.0002 |
| Vertebrate Complex V / codons | Sus scrofa (9823) | 3 | 5,289 | 0.0573 ± 0.0018 |
| Vertebrate Complex V / codons | Bos taurus (9913) | 3 | 5,580 | 0.0577 ± 0.0032 |
| Vertebrate / small capacity | Mus musculus (10090) | 3 | 2,412 | 0.1064 ± 0.0068 |
| Vertebrate / small capacity | Rattus norvegicus (10116) | 3 | 1,608 | 0.1478 ± 0.0128 |
| Vertebrate / small capacity | Danio rerio (7955) | 3 | 162 | 1.3073 ± 0.1063 |
| Vertebrate / small capacity | Xenopus laevis (8355) | 3 | 165 | 1.2920 ± 0.0555 |
| Vertebrate / small capacity | Gallus gallus (9031) | 3 | 1,134 | 0.1024 ± 0.0029 |
| Vertebrate / small capacity | Pongo abelii (9601) | 3 | 408 | 0.7460 ± 0.1759 |
| Vertebrate / small capacity | Homo sapiens (9606) | 3 | 66,774 | 0.1003 ± 0.0010 |
| Vertebrate / small capacity | Sus scrofa (9823) | 3 | 5,289 | 0.1492 ± 0.0079 |
| Vertebrate / small capacity | Bos taurus (9913) | 3 | 5,580 | 0.2386 ± 0.0053 |
| Vertebrate / family-macro objective | Mus musculus (10090) | 3 | 2,412 | 0.2540 ± 0.0335 |
| Vertebrate / family-macro objective | Rattus norvegicus (10116) | 3 | 1,608 | 0.3344 ± 0.0039 |
| Vertebrate / family-macro objective | Danio rerio (7955) | 3 | 162 | 2.3393 ± 0.2014 |
| Vertebrate / family-macro objective | Xenopus laevis (8355) | 3 | 165 | 2.3207 ± 0.0337 |
| Vertebrate / family-macro objective | Gallus gallus (9031) | 3 | 1,134 | 0.4259 ± 0.0514 |
| Vertebrate / family-macro objective | Pongo abelii (9601) | 3 | 408 | 1.3736 ± 0.0438 |
| Vertebrate / family-macro objective | Homo sapiens (9606) | 3 | 66,774 | 0.0857 ± 0.0009 |
| Vertebrate / family-macro objective | Sus scrofa (9823) | 3 | 5,289 | 0.1328 ± 0.0052 |
| Vertebrate / family-macro objective | Bos taurus (9913) | 3 | 5,580 | 0.1227 ± 0.0089 |
| Vertebrate / no explicit positions | Mus musculus (10090) | 3 | 2,412 | 0.0508 ± 0.0025 |
| Vertebrate / no explicit positions | Rattus norvegicus (10116) | 3 | 1,608 | 0.0520 ± 0.0016 |
| Vertebrate / no explicit positions | Danio rerio (7955) | 3 | 162 | 0.0707 ± 0.0025 |
| Vertebrate / no explicit positions | Xenopus laevis (8355) | 3 | 165 | 0.0675 ± 0.0065 |
| Vertebrate / no explicit positions | Gallus gallus (9031) | 3 | 1,134 | 0.0666 ± 0.0048 |
| Vertebrate / no explicit positions | Pongo abelii (9601) | 3 | 408 | 0.0623 ± 0.0047 |
| Vertebrate / no explicit positions | Homo sapiens (9606) | 3 | 66,774 | 0.0585 ± 0.0001 |
| Vertebrate / no explicit positions | Sus scrofa (9823) | 3 | 5,289 | 0.0536 ± 0.0011 |
| Vertebrate / no explicit positions | Bos taurus (9913) | 3 | 5,580 | 0.0567 ± 0.0020 |

## Per-genetic-code training fit

Completed matched runs only. Support counts target bases per corpus pass, not independent individuals. Taxa are the frozen primary selection strata; genetic codes describe translation conventions, not added model tokens.

| Arm | Stratum | Initializations | Target bases | Final bits/base |
|---|---|---:|---:|---:|
| Human ATP8 / bases | 2 — vertebrate mitochondrial code | 3 | 84,048 | 0.0632 ± 0.0001 |
| Human ATP8 / codons | 2 — vertebrate mitochondrial code | 3 | 84,048 | 0.0630 ± 0.0000 |
| Human Complex V / codons | 1 — standard nuclear code | 3 | 2,214 | 0.0512 ± 0.0037 |
| Human Complex V / codons | 2 — vertebrate mitochondrial code | 3 | 82,008 | 0.0618 ± 0.0016 |
| Mammal Complex V / codons | 1 — standard nuclear code | 3 | 2,184 | 0.0604 ± 0.0057 |
| Mammal Complex V / codons | 2 — vertebrate mitochondrial code | 3 | 81,675 | 0.0598 ± 0.0002 |
| Vertebrate Complex V / codons | 1 — standard nuclear code | 3 | 2,184 | 0.0607 ± 0.0032 |
| Vertebrate Complex V / codons | 2 — vertebrate mitochondrial code | 3 | 81,348 | 0.0599 ± 0.0005 |
| Vertebrate / small capacity | 1 — standard nuclear code | 3 | 2,184 | 1.5753 ± 0.0509 |
| Vertebrate / small capacity | 2 — vertebrate mitochondrial code | 3 | 81,348 | 0.0826 ± 0.0010 |
| Vertebrate / family-macro objective | 1 — standard nuclear code | 3 | 2,184 | 0.0358 ± 0.0054 |
| Vertebrate / family-macro objective | 2 — vertebrate mitochondrial code | 3 | 81,348 | 0.1228 ± 0.0008 |
| Vertebrate / no explicit positions | 1 — standard nuclear code | 3 | 2,184 | 0.0714 ± 0.0012 |
| Vertebrate / no explicit positions | 2 — vertebrate mitochondrial code | 3 | 81,348 | 0.0576 ± 0.0000 |

## Full-pool capacity demonstration

This separate run uses the full recovered vertebrate pool and a larger codon model. Both factors change, so its fit is excluded from matched contrasts; one initialization has no across-seed SD.

| Arm | State | CDS | Parameters | Updates | Initial bits/base | Last bits/base | Last family-macro bits/base |
|---|---|---:|---:|---:|---:|---:|---:|
| full_vertebrate_capacity | Complete | 533 | 10,675,204 | 2,000 | 2.0417 | 0.0610 | 0.0643 |

## Exposure, runtime and memory

| Arm | Seed | State | Updates | Sequence presentations | Presentations / unique CDS | Target bases presented | Raw loop minutes | Budget-charged loop minutes | Peak process MiB |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| human_atp8_base | 17 | Complete | 2,000 | 32,000 | 77.67 | 6,528,000 | 21.71 | 21.71 | 706.3 |
| human_atp8_codon | 17 | Complete | 2,000 | 32,000 | 77.67 | 6,528,000 | 5.34 | 5.34 | 407.5 |
| human_complex_codon | 17 | Complete | 2,000 | 32,000 | 77.67 | 6,541,776 | 5.79 | 5.79 | 473.9 |
| mammal_complex_codon | 17 | Complete | 2,000 | 32,000 | 77.67 | 6,515,085 | 5.67 | 5.67 | 476.6 |
| vertebrate_complex_codon | 17 | Complete | 2,000 | 32,000 | 77.67 | 6,488,652 | 5.64 | 5.64 | 476.2 |
| vertebrate_small | 17 | Complete | 2,000 | 32,000 | 77.67 | 6,488,652 | 0.89 | 0.89 | 315.9 |
| vertebrate_family_macro | 17 | Complete | 2,000 | 32,000 | 77.67 | 6,488,652 | 5.76 | 5.76 | 478.6 |
| vertebrate_no_position | 17 | Complete | 2,000 | 32,000 | 77.67 | 6,488,652 | 5.75 | 5.75 | 475.0 |
| human_atp8_base | 29 | Complete | 2,000 | 32,000 | 77.67 | 6,528,000 | 300.19 | 38.19 | 715.4 |
| human_atp8_codon | 29 | Complete | 2,000 | 32,000 | 77.67 | 6,528,000 | 7.20 | 7.20 | 408.1 |
| human_complex_codon | 29 | Complete | 2,000 | 32,000 | 77.67 | 6,538,464 | 5.54 | 5.54 | 479.3 |
| mammal_complex_codon | 29 | Complete | 2,000 | 32,000 | 77.67 | 6,512,283 | 6.72 | 6.72 | 475.4 |
| vertebrate_complex_codon | 29 | Complete | 2,000 | 32,000 | 77.67 | 6,490,143 | 6.05 | 6.05 | 477.6 |
| vertebrate_small | 29 | Complete | 2,000 | 32,000 | 77.67 | 6,490,143 | 1.13 | 1.13 | 314.0 |
| vertebrate_family_macro | 29 | Complete | 2,000 | 32,000 | 77.67 | 6,490,143 | 6.93 | 6.93 | 477.0 |
| vertebrate_no_position | 29 | Complete | 2,000 | 32,000 | 77.67 | 6,490,143 | 6.98 | 6.98 | 476.0 |
| human_atp8_base | 43 | Complete | 2,000 | 32,000 | 77.67 | 6,528,000 | 27.55 | 27.55 | 703.9 |
| human_atp8_codon | 43 | Complete | 2,000 | 32,000 | 77.67 | 6,528,000 | 4.17 | 4.17 | 453.9 |
| human_complex_codon | 43 | Complete | 2,000 | 32,000 | 77.67 | 6,541,563 | 5.12 | 5.12 | 523.4 |
| mammal_complex_codon | 43 | Complete | 2,000 | 32,000 | 77.67 | 6,512,544 | 6.14 | 6.14 | 521.5 |
| vertebrate_complex_codon | 43 | Complete | 2,000 | 32,000 | 77.67 | 6,485,778 | 5.95 | 5.95 | 521.4 |
| vertebrate_small | 43 | Complete | 2,000 | 32,000 | 77.67 | 6,485,778 | 0.96 | 0.96 | 359.3 |
| vertebrate_family_macro | 43 | Complete | 2,000 | 32,000 | 77.67 | 6,485,778 | 4.93 | 4.93 | 524.1 |
| vertebrate_no_position | 43 | Complete | 2,000 | 32,000 | 77.67 | 6,485,778 | 4.58 | 4.58 | 519.0 |
| full_vertebrate_capacity | 17 | Complete | 2,000 | 32,000 | 60.04 | 6,506,472 | 24.58 | 24.58 | 838.7 |

25 of 25 completed runs reduced overall training bits/base from initialization. This is an optimization check; it does not establish unseen-sequence prediction.

## Interpretation limits and next decisions

- The first pilot used a different acquisition snapshot, model capacity, training budget and smoothing rule. Changes from its reported loss cannot be attributed to one factor; the new matched controls provide the interpretable within-round contrasts.
- More varied cohorts can have higher inherent entropy. Higher training loss across cohorts does not by itself mean poorer modeling. Compare baseline gaps, per-family fit and support before interpreting aggregate changes.
- Sparse families may have one or a few distinct CDS; a family-macro objective gives those sequences substantial weight and can encourage memorization. Zero variable-position support yields n/a, not zero error.
- The species panel and reviewed references are availability-based. Missing family/species combinations and global sequence sharing constrain taxonomic interpretation. Distinct CDS and initialization seeds are not independent biological replicates.
- No explicit-position control removes the positional encoding only; causal order and prefix length remain available.
- Capacity changes width, depth and head count together. It is a bundled architecture comparison, not an isolated estimate of any one dimension.
- Maximum data means the eligible recovered pool under declared routes and QC; maximum model means the selected profiled candidate under the desktop reserve and time budget. Neither claim describes an absolute global maximum.
- A future predictive decision requires a separately designed untouched evaluation cohort. This round performs no validation, testing, gap completion, structural prediction, or molecular-function assessment.

## Reproducibility artifacts

- [Frozen protocol](../../docs/secondary-protocol.md), [data card](../../docs/secondary-data-card.md), [plan](plan.json), [cohort manifests and overlaps](cohorts.json), [acquisition](acquisition.json).
- [Per-run results](results.json), [histories, metadata and baselines](runs.json), [saved-model integrity inventory](checkpoints.json), [all 64 codon counts and diversity supports](codon-coverage.json), [source receipts and hashes](source-registry.json), [hardware profile](hardware-benchmark.json).
- [PNG figure](curves.png) and [SVG figure](curves.svg), with the identical measurements on [linear axes](curves-linear.svg). Raw archives and local rolling/final checkpoints remain outside Git tracking.
- Software check outcomes are recorded separately by the project; generating this report does not certify an unexecuted check.
