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
- `srGS_low_lrGS_ok_TR_overlap.sh <srGS_bins.bed.gz> <lrGS_coverage_counts.tsv> <out_dir> <SR_MIN> <LR_MAX>
  <MAX_STR_MOTIF> <exclude_chroms> <segdup.bed.gz> <simprep.bed.gz> <centromere.bed.gz> LABEL=catalog.bed.gz ...` -
  bins with srGS `low5_frac >= SR_MIN` and lrGS low-sample fraction `< LR_MAX`, merged into regions and
  annotated with TR catalogs (STR = shortest motif <= MAX_STR_MOTIF, VNTR otherwise), TRF simple
  repeats, SegDup and centromere. Writes `regions.bed.gz`, `tr_loci.bed.gz`, `summary.tsv`.
- `gene_TR_coverage_lrGS_vs_srGS.py (--gene G [G ...] | --gene-list F) --gtf GTF --sr-bins SR.bed.gz
  --lr-summary LR_RD.summary.bed.gz --lr-n-samples N --sr-norm X --lr-norm Y [--sr-norm-x X --lr-norm-x Y]
  --catalog LABEL=cat.bed.gz ... [--transcript GENE=TX ...] [--flank 1000] [--sr-max 0.5 --lr-min 0.7] --out-prefix P` -
  TR catalog loci in each gene (+flank) with per-locus and per-100bp-bin srGS vs lrGS depth (raw and
  normalized to the genome-wide mean of autosomes or chrX) and fraction of samples < 5x; one PDF page per
  gene, bins with srGS < sr-max and lrGS >= lr-min shaded. Compound loci are classed by the dominant motif
  (component spanning most bp). Exon track: Ensembl_canonical, else APPRIS principal transcript.
  Needs tabix-indexed bin files. Normalizers (bins fully present in the HT): autosomes srGS 30.657x /
  lrGS 54.963x; chrX srGS 23.312x / lrGS 40.029x. lrGS RD summary:
  `final_vcfs/low_cov_benchmark/gnomAD_lrGS.hgsvc_hprc.RD.summary.bed.gz` (292 samples; omits each
  chromosome's final partial bin).
- `coding_TR_coverage_lrGS_vs_srGS.py --gtf GTF --catalog LABEL=cat.bed.gz ... --sr-bins --lr-summary
  --lr-n-samples --sr-norm --lr-norm --sr-norm-x --lr-norm-x --segdup SegDup.bed.gz [--sr-max 0.5 --lr-min 0.7]
  --out-prefix P` - genome-wide coding TRs (catalog loci overlapping a protein_coding CDS, overlapping loci
  across catalogs merged; chrY excluded), flagged if any 100bp bin has srGS < sr-max and lrGS >= lr-min
  (normalized depth). Writes all regions, flagged list, gene list, size table and `size_counts.pdf` (total coding TRs by size, linear y from 0) and `size_proportion.pdf` (% flagged by size, 0-50%; bar labels = flagged counts).
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

- srGS-low / lrGS-well-covered regions (`srGS_low_lrGS_ok.sr90_lr10/`: srGS low5_frac >= 0.9,
  lrGS < 10% samples low, chrY excluded; catalogs TRExplorer v1.0.1 + Vamos v2.1 from
  `final_vcfs/STR_catalog/`, tracks from `final_vcfs/low_cov_benchmark/`):
  14,784 regions, 35.70 Mb. Most are small (10,155 < 1 kb) but bp is dominated by large ones.
  - Any overlap with a catalog TR: 6,032 regions (40.8%) - STR only 3,479, VNTR only 1,242,
    both 1,311; 8,752 (59.2%) overlap none. Catalog TR loci cover only 1.34 Mb (3.8%) of region bp.
  - Majority-bp class (priority centromere > catalog TR > TRF simple repeat > SegDup > other):
    centromere 1,364 regions / 19.29 Mb (54%); SegDup 8,539 / 12.82 Mb (36%);
    STR 1,490 / 0.61 Mb; VNTR 925 / 0.35 Mb; non-catalog TRF simple repeat 844 / 1.56 Mb;
    other 1,622 / 1.06 Mb.
  - By bp, the srGS-only deficit is centromeric satellite (absent from both TR catalogs) and SegDups.
    STR/VNTR-explained regions are ~16% of regions but < 3% of bp.

- CEL (`gene_TR_coverage/CEL.{TR_loci.tsv,bins.tsv,pdf}`; chr9:133,061,981-133,071,861, CEL-201):
  25 catalog TR loci in gene +/- 1 kb (20 STR, 5 VNTR); 3 exonic - CTG STR (exon 1), AAG STR (exon 8)
  and the exon 11 33-bp VNTR (MODY8 locus, chr9:133,071,240-133,071,608; TRExplorer variation
  cluster whose MOTIFS also list single-C components).
  - lrGS is flat across the whole gene, VNTR included: 47.6-48.1x (0.87x of genome mean), 0% of
    samples < 5x in every bin.
  - srGS: ~1.0x over most STRs, but the 33-bp VNTR averages 11.9x (0.39x), 31% of samples < 5x,
    with the worst bin (133,071,300-400) at 4.1x (0.13x), 66% of samples < 5x. Milder srGS dips at the
    exon 8 AAG STR (0.56x, 4.5% < 5x) and intronic GTG/GGA STRs (0.68-0.84x).

- Coding TRs (`coding_TR_coverage/`; TRExplorer + Vamos, GENCODE v39 CDS, srGS < 0.5x and lrGS >= 0.7x):
  - 48,679 coding catalog loci -> 40,948 merged coding TR regions in 14,707 genes
    (STR 37,285 regions / 14,073 genes; VNTR 3,663 / 2,976).
  - Flagged (srGS-poor, lrGS-ok): 576 regions in 341 genes (STR 410, VNTR 166). 453 of the 576 are
    >= 50% inside SegDup (NBPF, RIMBP3C, CT47A, TBC1D3, GOLGA8 families - srGS depth ~0 from
    multi-mapping); 123 regions in 111 genes are outside SegDup (e.g. EPPK1, ACAN, MUC5B, CFAP46, LYPD8,
    KLHDC4, NACAD, BOK, MLPH, INPP5E, UBC - mostly GC-rich VNTRs). CEL's exon 11 VNTR is flagged
    (worst bin 0.13x srGS vs 0.87x lrGS) but is labelled SegDup (CEL/CELP homology), so the SegDup
    split is context, not a clean repeat-vs-mapping separation.
  - Flagged fraction rises with size: ~1% for <= 100 bp, 2.6% at 101-200, 11.7% at 201-500,
    15.2% at 501-1000, 21.6% > 1 kb (`coding_TR.size_table.tsv`).
  - Per-gene plots of all 341 flagged genes: `flagged_genes.pdf` (order = `coding_TR.flagged_genes.ordered.txt`:
    non-SegDup worst-first, then SegDup).

## Caveats
- No per-sample data: per-sample relative cutoffs and sex-aware chrX/chrY rules
  used for lrGS cannot be reproduced. chrY `mean` is diluted by XX samples; chrX
  is an XX/XY mix.
- `over_X` boundary (> vs >=) follows gnomAD's definition; `low5_frac` is
  approximately "fraction of samples below 5x".
- The 2026-08-28 run (old script: mean over present loci only, no empty bins)
  outputs were deleted from the bucket; superseded by this run.

- Coverage is at 100 bp bin resolution: short STRs inherit their bin's depth, so a flagged 10 bp
  STR may reflect its surroundings (often a SegDup) rather than the STR itself.

## Regenerate
Submit Terra config `ExtractHailCoverageBins_agora` in
LR_GNOMAD_1_CO-AoU_TALK/LR_GNOMAD-AoU_TALK_Annotation-Pipeline (default
`chroms` = chr1-22, X, Y; `bin_size` 100; `low_thresholds` "5,10").
