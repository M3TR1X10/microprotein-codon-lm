# Secondary pilot — training only

TRAINING ONLY; all diagnostics are resubstitution, with no held-out observations.

**16 of 25 planned runs completed their fixed update budgets.** 0 recorded runs are partial; 9 have no reported result. Partial runs are excluded from completed-run averages. A pending seed is not a failed scientific hypothesis.

The aggregate execution budget is 4 hours. Budget-chargeable execution time is 137.1 minutes. Run-loop timings below include interval diagnostics and checkpoints but exclude setup/initial diagnostics; scheduler time covers the broader execution.

Scheduler status: **running**. Work is still running or its last recorded controller state is active.

A documented lid-triggered standby interruption stopped the original execution after eight complete runs. The [post-launch runtime amendment](../../docs/secondary-runtime-amendment.md) excludes only 262.00 minutes of Windows-recorded hardware deep idle, once. All other time remains charged. Including that preserved exclusion, recorded execution elapsed is 399.1 minutes. Budget-chargeable time is conservative accounting, not exact CPU-compute time. Model settings, data and update budgets did not change.

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
| Human ATP8 / bases | 2/3 | 1,780,608 | 2,000 | 2.0476 ± 0.1007 | 0.0632 ± 0.0001 | 0.0632 ± 0.0001 | 0.0632 ± 0.0001 | 84,048 |
| Human ATP8 / codons | 2/3 | 1,792,128 | 2,000 | 1.9944 ± 0.0173 | 0.0630 ± 0.0000 | 0.0630 ± 0.0000 | 0.0630 ± 0.0000 | 84,048 |
| Human Complex V / codons | 2/3 | 1,792,128 | 2,000 | 1.9950 ± 0.0166 | 0.0612 ± 0.0019 | 0.0592 ± 0.0024 | 0.0615 ± 0.0020 | 82,056 |
| Mammal Complex V / codons | 2/3 | 1,792,128 | 2,000 | 2.0000 ± 0.0180 | 0.0598 ± 0.0004 | 0.0647 ± 0.0042 | 0.0648 ± 0.0004 | 70,125 |
| Vertebrate Complex V / codons | 2/3 | 1,792,128 | 2,000 | 2.0001 ± 0.0173 | 0.0597 ± 0.0003 | 0.0630 ± 0.0035 | 0.0646 ± 0.0004 | 69,432 |
| Vertebrate / small capacity | 2/3 | 59,712 | 2,000 | 2.0047 ± 0.0035 | 0.1222 ± 0.0014 | 1.4152 ± 0.0240 | 0.0787 ± 0.0008 | 69,432 |
| Vertebrate / family-macro objective | 2/3 | 1,792,128 | 2,000 | 2.0001 ± 0.0173 | 0.1205 ± 0.0010 | 0.0486 ± 0.0031 | 0.1012 ± 0.0002 | 69,432 |
| Vertebrate / no explicit positions | 2/3 | 1,792,128 | 2,000 | 2.0400 ± 0.0027 | 0.0579 ± 0.0000 | 0.0710 ± 0.0022 | 0.0619 ± 0.0000 | 69,432 |

![Training curves](curves.svg)

Every plotted curve represents one initialization. Dashed family-position baselines receive privileged family labels. There is no held-out curve. Base/codon token accuracy and raw token perplexity are not compared.

Variable-codon masks are cohort-specific. The base/codon pair shares the same biological support; the three vertebrate controls share their reference arm's support. Across different cohorts, changing support prevents a paired accuracy interpretation. This mask counts every observation at a variable position, not only minority codons. Every target codon position varies at least once in `human_atp8`; variable-position and overall fit are identical there, so this diagnostic does not resolve rare-variant performance.

## Controlled contrasts and prespecified expectations

Differences follow the subtraction order in each row and are paired by initialization seed. Positive values mean higher training negative log likelihood for that metric. These are descriptive optimization contrasts, not significance tests.

Met/reversed describes only the frozen directional training expectation evaluated on the mean paired difference. Missing seed pairs make that interpretation provisional. No direction was asserted for raw loss across different taxonomic cohorts. The capacity contrast changes width, depth and heads together; it is not a single-dimension intervention.

| Contrast | Complete/planned pairs | Overall delta bits/base | Family-macro delta bits/base | Variable-codon delta bits/base | Prespecified expectation | Mean-pair outcome |
|---|---:|---:|---:|---:|---|---|
| Tokenization: codon minus base | 2/3 | -0.000143 ± 0.000106 | -0.000143 ± 0.000106 | -0.000143 ± 0.000106 | Codon lower overall bits/base | Provisional: met (missing seed pairs) |
| Capacity: small minus standard | 2/3 | 0.062424 ± 0.001687 | 1.352222 ± 0.020579 | 0.014057 ± 0.001111 | Small architecture higher overall bits/base | Provisional: met (missing seed pairs) |
| Family weighting: macro minus token | 2/3 | 0.060802 ± 0.000645 | -0.014333 ± 0.000372 | 0.036645 ± 0.000543 | Family-macro objective lower family-macro bits/base | Provisional: met (missing seed pairs) |
| Explicit position: removed minus supplied | 2/3 | -0.001853 ± 0.000308 | 0.008047 ± 0.001267 | -0.002653 ± 0.000360 | Removed positions higher overall bits/base | Provisional: reversed (missing seed pairs) |

Each matched pair below has verified identical cohort, sampling trace, updates, sequence/base presentations, family exposures and scored supports. Variable support is comparable within each pair only.

| Contrast | Seed | Overall delta bits/base | Family-macro delta bits/base | Variable-codon delta bits/base | Variable target bases | Sequence presentations | Target bases presented | Seed expectation |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Tokenization: codon minus base | 17 | -0.000218 | -0.000218 | -0.000218 | 84,048 | 32,000 | 6,528,000 | Met |
| Tokenization: codon minus base | 29 | -0.000068 | -0.000068 | -0.000068 | 84,048 | 32,000 | 6,528,000 | Met |
| Capacity: small minus standard | 17 | 0.061231 | 1.337671 | 0.013271 | 69,432 | 32,000 | 6,488,652 | Met |
| Capacity: small minus standard | 29 | 0.063616 | 1.366773 | 0.014842 | 69,432 | 32,000 | 6,490,143 | Met |
| Family weighting: macro minus token | 17 | 0.061259 | -0.014070 | 0.036261 | 69,432 | 32,000 | 6,488,652 | Met |
| Family weighting: macro minus token | 29 | 0.060346 | -0.014597 | 0.037029 | 69,432 | 32,000 | 6,490,143 | Met |
| Explicit position: removed minus supplied | 17 | -0.002071 | 0.008943 | -0.002908 | 69,432 | 32,000 | 6,488,652 | Reversed |
| Explicit position: removed minus supplied | 29 | -0.001635 | 0.007151 | -0.002399 | 69,432 | 32,000 | 6,490,143 | Reversed |

## Training-fit baselines

Each distribution is fit and scored on the same training observations. Empirical rows use no smoothing. Smoothed rows use total prior mass 2 per categorical distribution (2/4 per base category, 2/64 per codon category). Equal total prior removes the first pilot's prior-mass imbalance; different factorizations and conditional contexts still differ. The family-position baseline is privileged: the transformer receives no family or organism token. Uniform entropy is 2 bits/base.

| Cohort / tokens | Baseline | Empirical bits/base | Prior-mass-2 bits/base | Prior-mass-2 family-macro bits/base |
|---|---|---:|---:|---:|
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
| Human ATP8 / bases | ATP8 | 2 | 84,048 | 0.0632 ± 0.0001 | 84,048 | 0.0632 ± 0.0001 |
| Human ATP8 / codons | ATP8 | 2 | 84,048 | 0.0630 ± 0.0000 | 84,048 | 0.0630 ± 0.0000 |
| Human Complex V / codons | ATP5F1E | 2 | 306 | 0.0633 ± 0.0049 | 6 | 0.3506 ± 0.0105 |
| Human Complex V / codons | ATP5ME | 2 | 414 | 0.0521 ± 0.0019 | 6 | 0.3638 ± 0.0013 |
| Human Complex V / codons | ATP5MF | 2 | 846 | 0.0432 ± 0.0002 | 36 | 0.2256 ± 0.0122 |
| Human Complex V / codons | ATP5MGL | 2 | 300 | 0.0498 ± 0.0054 | 0 | n/a |
| Human Complex V / codons | ATP5MJ | 2 | 174 | 0.0716 ± 0.0007 | 0 | n/a |
| Human Complex V / codons | ATP5MK | 2 | 174 | 0.0731 ± 0.0101 | 0 | n/a |
| Human Complex V / codons | ATP8 | 2 | 82,008 | 0.0614 ± 0.0020 | 82,008 | 0.0614 ± 0.0020 |
| Mammal Complex V / codons | ATP5F1E | 2 | 306 | 0.0864 ± 0.0168 | 0 | n/a |
| Mammal Complex V / codons | ATP5ME | 2 | 420 | 0.0759 ± 0.0188 | 0 | n/a |
| Mammal Complex V / codons | ATP5MF | 2 | 810 | 0.0507 ± 0.0078 | 0 | n/a |
| Mammal Complex V / codons | ATP5MGL | 2 | 300 | 0.0467 ± 0.0067 | 0 | n/a |
| Mammal Complex V / codons | ATP5MJ | 2 | 174 | 0.0719 ± 0.0154 | 0 | n/a |
| Mammal Complex V / codons | ATP5MK | 2 | 174 | 0.0614 ± 0.0026 | 0 | n/a |
| Mammal Complex V / codons | ATP8 | 2 | 81,675 | 0.0598 ± 0.0003 | 70,125 | 0.0648 ± 0.0004 |
| Vertebrate Complex V / codons | ATP5F1E | 2 | 306 | 0.0886 ± 0.0076 | 0 | n/a |
| Vertebrate Complex V / codons | ATP5ME | 2 | 420 | 0.0917 ± 0.0163 | 0 | n/a |
| Vertebrate Complex V / codons | ATP5MF | 2 | 810 | 0.0467 ± 0.0022 | 0 | n/a |
| Vertebrate Complex V / codons | ATP5MGL | 2 | 300 | 0.0381 ± 0.0001 | 0 | n/a |
| Vertebrate Complex V / codons | ATP5MJ | 2 | 174 | 0.0587 ± 0.0006 | 0 | n/a |
| Vertebrate Complex V / codons | ATP5MK | 2 | 174 | 0.0570 ± 0.0023 | 0 | n/a |
| Vertebrate Complex V / codons | ATP8 | 2 | 81,348 | 0.0597 ± 0.0004 | 69,432 | 0.0646 ± 0.0004 |
| Vertebrate / small capacity | ATP5F1E | 2 | 306 | 1.5431 ± 0.0977 | 0 | n/a |
| Vertebrate / small capacity | ATP5ME | 2 | 420 | 1.5186 ± 0.0136 | 0 | n/a |
| Vertebrate / small capacity | ATP5MF | 2 | 810 | 1.5561 ± 0.0343 | 0 | n/a |
| Vertebrate / small capacity | ATP5MGL | 2 | 300 | 1.8129 ± 0.1492 | 0 | n/a |
| Vertebrate / small capacity | ATP5MJ | 2 | 174 | 1.6163 ± 0.1182 | 0 | n/a |
| Vertebrate / small capacity | ATP5MK | 2 | 174 | 1.7768 ± 0.1210 | 0 | n/a |
| Vertebrate / small capacity | ATP8 | 2 | 81,348 | 0.0824 ± 0.0013 | 69,432 | 0.0787 ± 0.0008 |
| Vertebrate / family-macro objective | ATP5F1E | 2 | 306 | 0.0414 ± 0.0039 | 0 | n/a |
| Vertebrate / family-macro objective | ATP5ME | 2 | 420 | 0.0747 ± 0.0244 | 0 | n/a |
| Vertebrate / family-macro objective | ATP5MF | 2 | 810 | 0.0308 ± 0.0005 | 0 | n/a |
| Vertebrate / family-macro objective | ATP5MGL | 2 | 300 | 0.0227 ± 0.0110 | 0 | n/a |
| Vertebrate / family-macro objective | ATP5MJ | 2 | 174 | 0.0223 ± 0.0040 | 0 | n/a |
| Vertebrate / family-macro objective | ATP5MK | 2 | 174 | 0.0257 ± 0.0088 | 0 | n/a |
| Vertebrate / family-macro objective | ATP8 | 2 | 81,348 | 0.1227 ± 0.0011 | 69,432 | 0.1012 ± 0.0002 |
| Vertebrate / no explicit positions | ATP5F1E | 2 | 306 | 0.0771 ± 0.0084 | 0 | n/a |
| Vertebrate / no explicit positions | ATP5ME | 2 | 420 | 0.0789 ± 0.0025 | 0 | n/a |
| Vertebrate / no explicit positions | ATP5MF | 2 | 810 | 0.0631 ± 0.0040 | 0 | n/a |
| Vertebrate / no explicit positions | ATP5MGL | 2 | 300 | 0.0667 ± 0.0023 | 0 | n/a |
| Vertebrate / no explicit positions | ATP5MJ | 2 | 174 | 0.0860 ± 0.0038 | 0 | n/a |
| Vertebrate / no explicit positions | ATP5MK | 2 | 174 | 0.0677 ± 0.0146 | 0 | n/a |
| Vertebrate / no explicit positions | ATP8 | 2 | 81,348 | 0.0575 ± 0.0000 | 69,432 | 0.0619 ± 0.0000 |

## Per-taxon training fit

Completed matched runs only. Support counts target bases per corpus pass, not independent individuals. Taxa are the frozen primary selection strata; genetic codes describe translation conventions, not added model tokens.

| Arm | Stratum | Initializations | Target bases | Final bits/base |
|---|---|---:|---:|---:|
| Human ATP8 / bases | Homo sapiens (9606) | 2 | 84,048 | 0.0632 ± 0.0001 |
| Human ATP8 / codons | Homo sapiens (9606) | 2 | 84,048 | 0.0630 ± 0.0000 |
| Human Complex V / codons | Homo sapiens (9606) | 2 | 84,222 | 0.0612 ± 0.0019 |
| Mammal Complex V / codons | Mus musculus (10090) | 2 | 2,412 | 0.0605 ± 0.0012 |
| Mammal Complex V / codons | Rattus norvegicus (10116) | 2 | 1,608 | 0.0536 ± 0.0032 |
| Mammal Complex V / codons | Pongo abelii (9601) | 2 | 408 | 0.0832 ± 0.0052 |
| Mammal Complex V / codons | Homo sapiens (9606) | 2 | 66,774 | 0.0596 ± 0.0005 |
| Mammal Complex V / codons | Sus scrofa (9823) | 2 | 5,691 | 0.0606 ± 0.0023 |
| Mammal Complex V / codons | Bos taurus (9913) | 2 | 6,966 | 0.0609 ± 0.0010 |
| Vertebrate Complex V / codons | Mus musculus (10090) | 2 | 2,412 | 0.0556 ± 0.0020 |
| Vertebrate Complex V / codons | Rattus norvegicus (10116) | 2 | 1,608 | 0.0536 ± 0.0016 |
| Vertebrate Complex V / codons | Danio rerio (7955) | 2 | 162 | 0.0588 ± 0.0009 |
| Vertebrate Complex V / codons | Xenopus laevis (8355) | 2 | 165 | 0.0618 ± 0.0060 |
| Vertebrate Complex V / codons | Gallus gallus (9031) | 2 | 1,134 | 0.0616 ± 0.0018 |
| Vertebrate Complex V / codons | Pongo abelii (9601) | 2 | 408 | 0.0851 ± 0.0084 |
| Vertebrate Complex V / codons | Homo sapiens (9606) | 2 | 66,774 | 0.0603 ± 0.0001 |
| Vertebrate Complex V / codons | Sus scrofa (9823) | 2 | 5,289 | 0.0564 ± 0.0015 |
| Vertebrate Complex V / codons | Bos taurus (9913) | 2 | 5,580 | 0.0576 ± 0.0046 |
| Vertebrate / small capacity | Mus musculus (10090) | 2 | 2,412 | 0.1096 ± 0.0056 |
| Vertebrate / small capacity | Rattus norvegicus (10116) | 2 | 1,608 | 0.1406 ± 0.0045 |
| Vertebrate / small capacity | Danio rerio (7955) | 2 | 162 | 1.2484 ± 0.0420 |
| Vertebrate / small capacity | Xenopus laevis (8355) | 2 | 165 | 1.2825 ± 0.0750 |
| Vertebrate / small capacity | Gallus gallus (9031) | 2 | 1,134 | 0.1040 ± 0.0014 |
| Vertebrate / small capacity | Pongo abelii (9601) | 2 | 408 | 0.7270 ± 0.2444 |
| Vertebrate / small capacity | Homo sapiens (9606) | 2 | 66,774 | 0.1008 ± 0.0007 |
| Vertebrate / small capacity | Sus scrofa (9823) | 2 | 5,289 | 0.1532 ± 0.0052 |
| Vertebrate / small capacity | Bos taurus (9913) | 2 | 5,580 | 0.2412 ± 0.0040 |
| Vertebrate / family-macro objective | Mus musculus (10090) | 2 | 2,412 | 0.2472 ± 0.0443 |
| Vertebrate / family-macro objective | Rattus norvegicus (10116) | 2 | 1,608 | 0.3340 ± 0.0054 |
| Vertebrate / family-macro objective | Danio rerio (7955) | 2 | 162 | 2.2804 ± 0.2456 |
| Vertebrate / family-macro objective | Xenopus laevis (8355) | 2 | 165 | 2.3060 ± 0.0312 |
| Vertebrate / family-macro objective | Gallus gallus (9031) | 2 | 1,134 | 0.4486 ± 0.0467 |
| Vertebrate / family-macro objective | Pongo abelii (9601) | 2 | 408 | 1.3545 ± 0.0406 |
| Vertebrate / family-macro objective | Homo sapiens (9606) | 2 | 66,774 | 0.0859 ± 0.0013 |
| Vertebrate / family-macro objective | Sus scrofa (9823) | 2 | 5,289 | 0.1351 ± 0.0047 |
| Vertebrate / family-macro objective | Bos taurus (9913) | 2 | 5,580 | 0.1211 ± 0.0119 |
| Vertebrate / no explicit positions | Mus musculus (10090) | 2 | 2,412 | 0.0496 ± 0.0021 |
| Vertebrate / no explicit positions | Rattus norvegicus (10116) | 2 | 1,608 | 0.0510 ± 0.0002 |
| Vertebrate / no explicit positions | Danio rerio (7955) | 2 | 162 | 0.0721 ± 0.0012 |
| Vertebrate / no explicit positions | Xenopus laevis (8355) | 2 | 165 | 0.0649 ± 0.0066 |
| Vertebrate / no explicit positions | Gallus gallus (9031) | 2 | 1,134 | 0.0643 ± 0.0039 |
| Vertebrate / no explicit positions | Pongo abelii (9601) | 2 | 408 | 0.0645 ± 0.0038 |
| Vertebrate / no explicit positions | Homo sapiens (9606) | 2 | 66,774 | 0.0585 ± 0.0001 |
| Vertebrate / no explicit positions | Sus scrofa (9823) | 2 | 5,289 | 0.0541 ± 0.0011 |
| Vertebrate / no explicit positions | Bos taurus (9913) | 2 | 5,580 | 0.0568 ± 0.0028 |

## Per-genetic-code training fit

Completed matched runs only. Support counts target bases per corpus pass, not independent individuals. Taxa are the frozen primary selection strata; genetic codes describe translation conventions, not added model tokens.

| Arm | Stratum | Initializations | Target bases | Final bits/base |
|---|---|---:|---:|---:|
| Human ATP8 / bases | 2 — vertebrate mitochondrial code | 2 | 84,048 | 0.0632 ± 0.0001 |
| Human ATP8 / codons | 2 — vertebrate mitochondrial code | 2 | 84,048 | 0.0630 ± 0.0000 |
| Human Complex V / codons | 1 — standard nuclear code | 2 | 2,214 | 0.0531 ± 0.0018 |
| Human Complex V / codons | 2 — vertebrate mitochondrial code | 2 | 82,008 | 0.0614 ± 0.0020 |
| Mammal Complex V / codons | 1 — standard nuclear code | 2 | 2,184 | 0.0626 ± 0.0061 |
| Mammal Complex V / codons | 2 — vertebrate mitochondrial code | 2 | 81,675 | 0.0598 ± 0.0003 |
| Vertebrate Complex V / codons | 1 — standard nuclear code | 2 | 2,184 | 0.0618 ± 0.0036 |
| Vertebrate Complex V / codons | 2 — vertebrate mitochondrial code | 2 | 81,348 | 0.0597 ± 0.0004 |
| Vertebrate / small capacity | 1 — standard nuclear code | 2 | 2,184 | 1.6047 ± 0.0021 |
| Vertebrate / small capacity | 2 — vertebrate mitochondrial code | 2 | 81,348 | 0.0824 ± 0.0013 |
| Vertebrate / family-macro objective | 1 — standard nuclear code | 2 | 2,184 | 0.0385 ± 0.0038 |
| Vertebrate / family-macro objective | 2 — vertebrate mitochondrial code | 2 | 81,348 | 0.1227 ± 0.0011 |
| Vertebrate / no explicit positions | 1 — standard nuclear code | 2 | 2,184 | 0.0708 ± 0.0007 |
| Vertebrate / no explicit positions | 2 — vertebrate mitochondrial code | 2 | 81,348 | 0.0575 ± 0.0000 |

## Full-pool capacity demonstration

This separate run uses the full recovered vertebrate pool and a larger codon model. Both factors change, so its fit is excluded from matched contrasts; one initialization has no across-seed SD.

| Arm | State | CDS | Parameters | Updates | Initial bits/base | Last bits/base | Last family-macro bits/base |
|---|---|---:|---:|---:|---:|---:|---:|
| Full-pool run | No reported result | — | — | — | — | — | — |

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

16 of 16 completed runs reduced overall training bits/base from initialization. This is an optimization check; it does not establish unseen-sequence prediction.

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
