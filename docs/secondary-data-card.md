# Secondary dataset card

This snapshot contains observed ATP-synthase coding sequences for training only. Counts below are generated from frozen cohort files and checked against their manifests; they are not an estimate of independent individuals or all available vertebrate sequences.

## Scope and sources

The declared panel and exact family identifiers are in [biology.json](../configs/secondary/biology.json). IMPI establishes human mitochondrial eligibility; UniProt supplies reviewed family/species references; ENA supplies observed nucleotide records. The [acquisition report](../reports/secondary/acquisition.json) records the CDS-index/cross-reference routes, exclusions and coverage. When present, the [genome-acquisition report](../reports/secondary/genome-acquisition.json) separately records retrieval of mitochondrial parent records in the declared 14,000–20,000-nt range and extraction of explicitly annotated ATP8 CDS. Parent length does not establish whole-genome completeness; the extracted CDS must pass the complete-CDS checks. The original acquisition snapshot is preserved; the final cohort manifests identify which enriched pool was selected.

The families share ATP synthase membership; they are not all ATP8 homologs. Cross-species reference coverage is incomplete and availability-based. Nonmammalian coverage and actual family-by-taxon counts must be read from this snapshot, not inferred from the species panel.

## Quality and sequence identity

The operational maximum is 100 amino acids excluding stop, with at least 95% amino-acid identity to the sequence's own-species reviewed reference and the reference length. Only complete, unambiguous CDS with allowed start, full stop, matching annotated translation and no internal stop or translation exceptions are eligible. Mitochondrial ATP8 uses genetic code 2; nuclear families use code 1. T is transcribed to U after QC. No reverse translation, invented variants or repaired bases are added.

Identical RNA is globally deduplicated. Original accessions, versions, source files, hashes and taxonomic observations remain in per-sequence provenance in the enriched raw pool, joined from compact training views by sequence SHA-256. A primary selection stratum represents each shared sequence; taxon-stratum counts do not erase its other observed taxa. Annotation and high sequence identity do not experimentally establish equivalent molecular function.

## Frozen cohort counts

| View | Distinct CDS | Distinct peptides | Families | Taxon strata | Target bases | Observed target codons / 64 |
|---|---:|---:|---:|---:|---:|---:|
| full_human | 422 | 185 | 7 | 1 | 86,262 | 64 |
| full_mammal | 524 | 244 | 7 | 6 | 106,884 | 64 |
| full_vertebrate | 533 | 250 | 7 | 9 | 108,345 | 64 |
| human_atp8 | 412 | 178 | 1 | 1 | 84,048 | 57 |
| human_complex | 412 | 182 | 7 | 1 | 84,222 | 64 |
| mammal_complex | 412 | 201 | 7 | 6 | 83,859 | 64 |
| vertebrate_complex | 412 | 206 | 7 | 9 | 83,532 | 64 |

The matched family expansion replaces **10/412 CDS (2.43%)**. The nonmammalian expansion replaces **9/412 CDS (2.18%)** in the mammal view. These small interventions limit what aggregate training differences can reveal about broader family or taxonomic effects.

The base and codon human-ATP8 arms share exactly one cohort. Matched Complex V views share family quotas and total CDS count. Equal sequence count does not mean equal nucleotide exposure, independent biological sampling, or identical effective diversity. The full pools retain all eligible distinct CDS recovered within the declared acquisition scope.

## Family and taxon composition

| View | Family | Primary taxon stratum | Distinct CDS |
|---|---|---|---:|
| full_human | ATP5F1E | Homo sapiens (9606) | 2 |
| full_human | ATP5ME | Homo sapiens (9606) | 2 |
| full_human | ATP5MF | Homo sapiens (9606) | 3 |
| full_human | ATP5MGL | Homo sapiens (9606) | 1 |
| full_human | ATP5MJ | Homo sapiens (9606) | 1 |
| full_human | ATP5MK | Homo sapiens (9606) | 1 |
| full_human | ATP8 | Homo sapiens (9606) | 412 |
| full_mammal | ATP5F1E | Mus musculus (10090) | 1 |
| full_mammal | ATP5F1E | Rattus norvegicus (10116) | 1 |
| full_mammal | ATP5F1E | Homo sapiens (9606) | 2 |
| full_mammal | ATP5F1E | Bos taurus (9913) | 1 |
| full_mammal | ATP5ME | Mus musculus (10090) | 2 |
| full_mammal | ATP5ME | Rattus norvegicus (10116) | 1 |
| full_mammal | ATP5ME | Pongo abelii (9601) | 1 |
| full_mammal | ATP5ME | Homo sapiens (9606) | 2 |
| full_mammal | ATP5ME | Sus scrofa (9823) | 2 |
| full_mammal | ATP5ME | Bos taurus (9913) | 1 |
| full_mammal | ATP5MF | Mus musculus (10090) | 2 |
| full_mammal | ATP5MF | Pongo abelii (9601) | 1 |
| full_mammal | ATP5MF | Homo sapiens (9606) | 3 |
| full_mammal | ATP5MF | Sus scrofa (9823) | 2 |
| full_mammal | ATP5MF | Bos taurus (9913) | 1 |
| full_mammal | ATP5MGL | Homo sapiens (9606) | 1 |
| full_mammal | ATP5MJ | Mus musculus (10090) | 1 |
| full_mammal | ATP5MJ | Homo sapiens (9606) | 1 |
| full_mammal | ATP5MJ | Bos taurus (9913) | 1 |
| full_mammal | ATP5MK | Mus musculus (10090) | 1 |
| full_mammal | ATP5MK | Rattus norvegicus (10116) | 1 |
| full_mammal | ATP5MK | Homo sapiens (9606) | 1 |
| full_mammal | ATP5MK | Bos taurus (9913) | 1 |
| full_mammal | ATP8 | Mus musculus (10090) | 12 |
| full_mammal | ATP8 | Rattus norvegicus (10116) | 8 |
| full_mammal | ATP8 | Pongo abelii (9601) | 2 |
| full_mammal | ATP8 | Homo sapiens (9606) | 412 |
| full_mammal | ATP8 | Sus scrofa (9823) | 27 |
| full_mammal | ATP8 | Bos taurus (9913) | 32 |
| full_vertebrate | ATP5F1E | Mus musculus (10090) | 1 |
| full_vertebrate | ATP5F1E | Rattus norvegicus (10116) | 1 |
| full_vertebrate | ATP5F1E | Homo sapiens (9606) | 2 |
| full_vertebrate | ATP5F1E | Bos taurus (9913) | 1 |
| full_vertebrate | ATP5ME | Mus musculus (10090) | 2 |
| full_vertebrate | ATP5ME | Rattus norvegicus (10116) | 1 |
| full_vertebrate | ATP5ME | Pongo abelii (9601) | 1 |
| full_vertebrate | ATP5ME | Homo sapiens (9606) | 2 |
| full_vertebrate | ATP5ME | Sus scrofa (9823) | 2 |
| full_vertebrate | ATP5ME | Bos taurus (9913) | 1 |
| full_vertebrate | ATP5MF | Mus musculus (10090) | 2 |
| full_vertebrate | ATP5MF | Pongo abelii (9601) | 1 |
| full_vertebrate | ATP5MF | Homo sapiens (9606) | 3 |
| full_vertebrate | ATP5MF | Sus scrofa (9823) | 2 |
| full_vertebrate | ATP5MF | Bos taurus (9913) | 1 |
| full_vertebrate | ATP5MGL | Homo sapiens (9606) | 1 |
| full_vertebrate | ATP5MJ | Mus musculus (10090) | 1 |
| full_vertebrate | ATP5MJ | Homo sapiens (9606) | 1 |
| full_vertebrate | ATP5MJ | Bos taurus (9913) | 1 |
| full_vertebrate | ATP5MK | Mus musculus (10090) | 1 |
| full_vertebrate | ATP5MK | Rattus norvegicus (10116) | 1 |
| full_vertebrate | ATP5MK | Homo sapiens (9606) | 1 |
| full_vertebrate | ATP5MK | Bos taurus (9913) | 1 |
| full_vertebrate | ATP8 | Mus musculus (10090) | 12 |
| full_vertebrate | ATP8 | Rattus norvegicus (10116) | 8 |
| full_vertebrate | ATP8 | Danio rerio (7955) | 1 |
| full_vertebrate | ATP8 | Xenopus laevis (8355) | 1 |
| full_vertebrate | ATP8 | Gallus gallus (9031) | 7 |
| full_vertebrate | ATP8 | Pongo abelii (9601) | 2 |
| full_vertebrate | ATP8 | Homo sapiens (9606) | 412 |
| full_vertebrate | ATP8 | Sus scrofa (9823) | 27 |
| full_vertebrate | ATP8 | Bos taurus (9913) | 32 |
| human_atp8 | ATP8 | Homo sapiens (9606) | 412 |
| human_complex | ATP5F1E | Homo sapiens (9606) | 2 |
| human_complex | ATP5ME | Homo sapiens (9606) | 2 |
| human_complex | ATP5MF | Homo sapiens (9606) | 3 |
| human_complex | ATP5MGL | Homo sapiens (9606) | 1 |
| human_complex | ATP5MJ | Homo sapiens (9606) | 1 |
| human_complex | ATP5MK | Homo sapiens (9606) | 1 |
| human_complex | ATP8 | Homo sapiens (9606) | 402 |
| mammal_complex | ATP5F1E | Homo sapiens (9606) | 1 |
| mammal_complex | ATP5F1E | Bos taurus (9913) | 1 |
| mammal_complex | ATP5ME | Homo sapiens (9606) | 1 |
| mammal_complex | ATP5ME | Bos taurus (9913) | 1 |
| mammal_complex | ATP5MF | Homo sapiens (9606) | 1 |
| mammal_complex | ATP5MF | Sus scrofa (9823) | 1 |
| mammal_complex | ATP5MF | Bos taurus (9913) | 1 |
| mammal_complex | ATP5MGL | Homo sapiens (9606) | 1 |
| mammal_complex | ATP5MJ | Homo sapiens (9606) | 1 |
| mammal_complex | ATP5MK | Homo sapiens (9606) | 1 |
| mammal_complex | ATP8 | Mus musculus (10090) | 12 |
| mammal_complex | ATP8 | Rattus norvegicus (10116) | 8 |
| mammal_complex | ATP8 | Pongo abelii (9601) | 2 |
| mammal_complex | ATP8 | Homo sapiens (9606) | 321 |
| mammal_complex | ATP8 | Sus scrofa (9823) | 27 |
| mammal_complex | ATP8 | Bos taurus (9913) | 32 |
| vertebrate_complex | ATP5F1E | Homo sapiens (9606) | 1 |
| vertebrate_complex | ATP5F1E | Bos taurus (9913) | 1 |
| vertebrate_complex | ATP5ME | Homo sapiens (9606) | 1 |
| vertebrate_complex | ATP5ME | Bos taurus (9913) | 1 |
| vertebrate_complex | ATP5MF | Homo sapiens (9606) | 1 |
| vertebrate_complex | ATP5MF | Sus scrofa (9823) | 1 |
| vertebrate_complex | ATP5MF | Bos taurus (9913) | 1 |
| vertebrate_complex | ATP5MGL | Homo sapiens (9606) | 1 |
| vertebrate_complex | ATP5MJ | Homo sapiens (9606) | 1 |
| vertebrate_complex | ATP5MK | Homo sapiens (9606) | 1 |
| vertebrate_complex | ATP8 | Mus musculus (10090) | 12 |
| vertebrate_complex | ATP8 | Rattus norvegicus (10116) | 8 |
| vertebrate_complex | ATP8 | Danio rerio (7955) | 1 |
| vertebrate_complex | ATP8 | Xenopus laevis (8355) | 1 |
| vertebrate_complex | ATP8 | Gallus gallus (9031) | 7 |
| vertebrate_complex | ATP8 | Pongo abelii (9601) | 2 |
| vertebrate_complex | ATP8 | Homo sapiens (9606) | 321 |
| vertebrate_complex | ATP8 | Sus scrofa (9823) | 25 |
| vertebrate_complex | ATP8 | Bos taurus (9913) | 25 |

## Intended interpretation and limitations

The entire snapshot is available to training. It must not later be relabeled as an untouched validation/test set. The matched views replace records and overlap; they are not independent cohorts. Quotas may retain sparse families only in the human stratum, so adding organisms does not expand every family equally.

Variable-codon positions are identified separately within each family and primary taxon stratum. The base arm scores all three bases of those same variable codons. Singleton or invariant strata contribute no variable-position targets. Different cohorts have different variable-position supports; those metrics are descriptive across cohorts. A position with even one alternative codon is variable, so this mask does not isolate minority-codon observations. When all target codon positions vary somewhere in a cohort, variable-position fit equals overall fit.

All 64 RNA codons are vocabulary classes even when observed counts are zero. The [coverage/diversity artifact](../reports/secondary/codon-coverage.json) reports every class, source of positional variation and its support. No pairwise-identity estimate or effective-sample-size estimate is computed.

See the [protocol](secondary-protocol.md), [cohort manifests](../reports/secondary/cohorts.json), [retrieval registry](../reports/secondary/source-registry.json) and [training report](../reports/secondary/report.md). Raw archives and checkpoints remain outside Git tracking.
