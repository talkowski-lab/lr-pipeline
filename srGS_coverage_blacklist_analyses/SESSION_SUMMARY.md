# srGS (gnomAD v3) genome-wide coverage

## Goal
Extract short-read (srGS) cohort coverage across GRCh38 in 100bp bins, in a
form directly comparable to the lrGS (HPRC+HGSVC mosdepth) bins in
`../lrGS_coverage_blacklist_analyses/`.

## Input
`gs://gcp-public-data--gnomad/release/3.0.1/coverage/genomes/gnomad.genomes.r3.0.1.coverage.ht`
- Hail **Table** (not MatrixTable), keyed by `locus` (GRCh38), one row per base;
  2,873,066,187 rows over chr1-22, X, Y (no chrM); 5,000 partitions, ~120 GB.
- Row fields: `mean`, `median_approx`, `total_DP`, `over_{1,5,10,15,20,25,30,50,100}`
  (fraction of samples with depth over X). Cohort aggregates only - no per-sample data.
- `total_DP / mean` ~ 71,700 samples (gnomAD v3).
- Identical content as text: `.../gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz` (75 GB).
- Loci absent from the table (N-gaps and other zero-coverage sequence, ~64 Mb
  beyond N-gaps) are treated as zero coverage in all samples.

## Scripts
- `extract_chrom_coverage_bins.py` - for one chromosome, downloads only the HT
  partitions overlapping it (via `rows/metadata.json.gz` range bounds), bins,
  and exports a bgzipped BED with every bin across the full chromosome length.
  `--hail-table-path --chrom --bin-size --low-thresholds 5,10 --out-bed X.bed.bgz`
- `ExtractHailCoverageBins.wdl` - scatters the script over chromosomes and
  concatenates. Script is staged at
  `gs://fc-a22a385b-ed3b-45ba-9d6c-87ca01c1e6b8/XZ/scripts/extract_chrom_coverage_bins.py`.
  Terra config `ExtractHailCoverageBins_agora` (Agora `xzhao_methods/ExtractHailCoverageBins`).
- `mask_low_dp_genotypes.py`, `bed_from_filtered_chr22.py`, `bed_from_masked_chr22.py` -
  earlier chr22 genotype-masking work (not part of this extraction).

## Output columns (no header)
`chrom, start, end, mean_cov, n_loci, low5_frac, low10_frac`
- `mean_cov` - sum of per-locus cohort mean depth / bin length (absent loci = 0)
- `n_loci` - loci in the bin present in the HT
- `low<X>_frac` - mean over the bin of (1 - over_X); absent loci count as 1.0.
  srGS analogue of the lrGS "fraction of samples low-coverage" (lrGS cutoff =
  20% of each sample's median; ~6x at 30x, so `low5_frac` is the closest match).

## Outputs
- `gnomAD_srGS.coverage_bins.low4.bed.gz` - earlier product (pre-dates this run).
- `gnomAD_srGS.r3.0.1.coverage_bins.100bp.bed.gz` - genome-wide 100bp bins
  (30,882,711 bins, contiguous, chr1-22/X/Y full length; 473 MB, local only -
  gitignored). Terra submission `31bc6f0b-9a5a-472c-a817-c633feebbc0d`; copy at
  `gs://fc-a22a385b-ed3b-45ba-9d6c-87ca01c1e6b8/submissions/31bc6f0b-9a5a-472c-a817-c633feebbc0d/ExtractHailCoverageBins/43978acd-8033-4ebb-9026-23c2796af0d8/call-CombineBeds/`.

## Findings
- chr22 test: 508,185 bins (full 50.82 Mb); 37.07 Mb of loci present; 137,442
  bins entirely absent; only 51 partially-absent bins (gaps are essentially
  bin-aligned); typical euchromatic bins ~33x with low5_frac ~0.

- Genome-wide: loci present 2,873.07 Mb (matches HT row count exactly); 215.17 Mb
  of bins entirely absent; bp-weighted mean coverage 27.98x (autosomes ~30x;
  chrX 22.6x; chrY 2.7x - XX/XY mix).
- Low-coverage regions (bins with low5_frac >= threshold, merged) vs lrGS
  blacklist (`../lrGS_coverage_blacklist_analyses/hgsvc_hprc.blacklist.p*.raw.bed.gz`):

  | threshold | srGS | lrGS | shared | srGS-only | lrGS-only |
  |---|---|---|---|---|---|
  | p90, all | 243.02 Mb | 172.41 Mb | 171.87 Mb | 71.16 Mb | 0.54 Mb |
  | p90, no chrY | 203.14 Mb | 138.77 Mb | 138.23 Mb | 64.91 Mb | 0.54 Mb |
  | p50, all | 270.53 Mb | 185.63 Mb | 184.39 Mb | 86.14 Mb | 1.24 Mb |
  | p50, no chrY | 213.99 Mb | 151.75 Mb | 150.51 Mb | 63.49 Mb | 1.24 Mb |

  The lrGS blacklist is almost entirely nested inside the srGS one; srGS has
  ~65 Mb of extra poorly-covered sequence outside chrY. Cutoffs differ (srGS:
  absolute <5x; lrGS: <20% of each sample's median), so treat as approximate.
  srGS chrY is low almost everywhere at p50 because ~half the cohort is XX.

## Caveats
- No per-sample data: per-sample relative cutoffs and sex-aware chrX/chrY rules
  used for lrGS cannot be reproduced. chrY `mean` is diluted by XX samples; chrX
  is an XX/XY mix.
- `over_X` boundary (> vs >=) follows gnomAD's definition; `low5_frac` is
  approximately "fraction of samples below 5x".
- The 2026-08-28 run (old script: mean over present loci only, no empty bins)
  outputs were deleted from the bucket; superseded by this run.

## Regenerate
Submit Terra config `ExtractHailCoverageBins_agora` in
LR_GNOMAD_1_CO-AoU_TALK/LR_GNOMAD-AoU_TALK_Annotation-Pipeline (default
`chroms` = chr1-22, X, Y; `bin_size` 100; `low_thresholds` "5,10").
