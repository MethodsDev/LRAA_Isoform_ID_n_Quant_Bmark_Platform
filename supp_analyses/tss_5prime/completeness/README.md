# Read 5' completeness and the cap signature (PBMC vs ALS)

This is a measurement-only study, to review before any change to how LRAA calls TSS sites
or prunes them as degradation products. Nothing here changes LRAA or the site-usage
results.

- `make` builds the tables for each dataset, one dataset at a time, with `CPU` (12)
  workers. It takes about 4 minutes per dataset.
- `make report` also knits `tss_5prime_completeness.Rmd`.

The measurements come from LRAA's `util/sc/tss_5prime_completeness.py`. It uses LRAA's
read filters and transcribed-strand rule (as `alt_termini_read_check.py` does), and takes
each read's 5' clip, and whether LRAA strips it as untemplated G's, from `Pretty_alignment`.

## Inputs (per dataset)

- the aligned BAM;
- cell clusters;
- the site-usage site table (`<prefix>.sites.tsv`);
- the high-confidence site-level TSS events with their CAGE comparison
  (`<prefix>.TSS_events_vs_CAGE.events.tsv`);
- CTAT GENCODE v22 `ref_annot.gtf` and genome;
- the PBMC rDNA mask (same genome for both datasets);
- FANTOM5 CAGE peaks, from PBMC `__Site_Usage_Analysis/ext/`.

## Rules

**Single-TSS genes.** A gene qualifies if:

- it has exactly one TSS site in the dataset's site table;
- its symbol is a GENCODE `gene_name` on the same contig and strand;
- a GENCODE transcript of the gene has the site in its first exon, or within 50 nt
  upstream of it.

The **main transcript** is chosen among those transcripts: the one whose full intron chain
most reads reproduce. Ties go to more fitting reads, then to the longer model. Restricting
to transcripts that start at the site keeps a fragment model from making truncated reads
look full length.

**CAGE status.** Take the FANTOM5 peaks on the gene's strand over the main transcript plus
100 nt upstream, and merge peaks within 100 nt of each other into clusters. The status is:

- `single`: the top cluster holds at least 90% of the peak score and lies within 100 nt of
  the LRAA site;
- `multi`: there are peaks, but not like that;
- `none`: there are no peaks.

The analyses use genes with status `single` or `none`.

**Reads of a gene.** These are reads whose introns are all among the main transcript's
introns.

- A read is **polyA-anchored** if its 3' end lies within 50 nt of the observed 3' end: the
  most common read 3' end, in 10 nt bins, along the main transcript.
- Its **extent** is its distance along the transcript from its 5' end to that observed end.
- A read whose 5' end is within 50 nt of the TSS, or upstream of it, **reaches** it. Its
  extent is then the gene's length L.

**Completeness curves.** S(d) for a gene is the share of its anchored reads with extent of
at least d. K(d) is the mean of S(d) over genes with L >= d and at least 20 anchored reads,
so genes are weighted equally. It is computed for the whole dataset, per cell-type class
and per cluster.

**Cap signature.** Each read's transcript-5' clip is one of:

- `noclip`;
- `G`: LRAA strips the clip as untemplated G's (a whole clip of 1-3 G's, or a proximal run
  of at least 3 G's);
- `other`.

Because an untemplated G aligns when the genome has a G just upstream, the genomic base
upstream of the read's start is recorded too. Rates are compared at upstream bases other
than G.

**Event sites.** For each high-confidence TSS event's upstream (U) and downstream (D)
site, the tables record:

- the reads starting within 25 nt of the site (the site-counting tolerance);
- from the reads that start at U, whether D lies in their aligned exons;
- the median aligned distances from U and from D to those reads' 3' ends.

**Simulated threshold.** This is section D of the report. In each cluster, the reads
expected at D from truncated U molecules, per read reaching U, are

    r = (K(d_D - 25) - K(d_D + 25)) / K(d_U)

The event is flagged if either of these holds in the cluster where D's share is higher:

- D's reads number at most 2 × U's reads × r;
- after subtracting U × r from each cluster's D reads, the U-vs-D shift falls below 0.2, or
  below half the observed shift.

Events where D is not inside the U molecules (fewer than half the U-start reads cover D)
are exempt.

## Outputs

`<dataset>.` prefixed tables: `genes.tsv`, `extent.tsv.gz`, `capsig.tsv.gz`,
`event_pairs.tsv` and `event_sites.tsv.gz`. They are described in the script's header.

## Tandem TSS pairs (`make tandem`, `make tandem_report`)

Features that might tell a real downstream TSS from a truncation product where the two are
structurally identical: two TSSs of one splice pattern, the downstream one (D) inside the
upstream one's (U) first exon. No cap / untemplated-G feature. Code: LRAA
`util/sc/tandem_tss_features.py` (header documents every feature); report
`tandem_tss_features.Rmd`; completeness curves via `k_curves.R`.

- Likely real: tandem pairs of the dataset's LRAA models with distinct FANTOM5 CAGE peaks at
  both sites. False: read-start peaks in first exons of single-TSS genes that pass LRAA's own
  TSS peak calling (aggregation + peakiness). MCF7 (SG-NEx ONT, LRAA v0.40.0 de novo) adds
  PBMC's CAGE-supported pairs measured in its reads, and PBMC + ALS single-TSS genes with a
  single CAGE promoter (`pool_single_tss_genes.py`).
- Both sets were chosen by LRAA's own D/U degradation rule, so D/U separates them by
  construction; the report compares features at matched D/U.
- ViennaRNA (`pip install --user ViennaRNA`) gives the folding energies; refTSS v4.1 is
  downloaded into `ext/`.

## Sub-studies

- `cap_signature_other_datasets/`: the cap signature in the bulk benchmark datasets (PacBio,
  ONT, SIRVs, MORFs) as one G over a per-library baseline (2026-10-27).
  - `cap_signature_other_datasets/morf_untrimmed/`: MORF transcripts ARE capped (G>=2 84.7-89.9%
    at the vector TSSs vs 6.0% at internal starts) and the vector's EF-1alpha 5' UTR intron is
    spliced in 95.9% of decisive reads, while the cloned ORF is intronless; measured on the
    untrimmed MORF reads, which the benchmark BAMs lack (2026-11-11).
- `cap_calibration/`: per-peak capped fractions are bimodal in every library with the
  signature, and one mode at 0 without it; the basis of LRAA's library-level cap calibration
  (pre-pass, v0.47.0), with the FTD false TSS sites placed on it (2026-10-29).
