# GD breakpoint profiles

Workflow `wdl_analyses/GDBreakpointProfiles.wdl` (Dockstore entry `GDBreakpointProfiles`, branches main and xz_analyses).
For each genomic-disorder (GD) locus, the segmental duplication (SD) at the left and right breakpoint, widened by half
its length on each side, is profiled with read depth, BAF and LD (1 kb copy-number bins by default, `cn_bin`). One figure per locus: left SD on the
left, right SD on the right; DEL carriers red to orange, DUP carriers blue to green; arrows at the BAF breakpoint
estimate of each DEL carrier.

Docker: `dockerfiles/Dockerfile.GDBreakpointProfiles` (scripts in this folder are copied to `/opt/gd_scripts`), image
`us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers/gd-breakpoint-profiles:<tag>` (current: 2026-11-09.9). The image also holds the two SD references in `/opt/gd_scripts/references/` (from `references/` here), used when `segdup_bed` / `bp_sd_bed` are not given.

## Steps
1. `gd_sd_windows.py` (run once, locally): GD call IDs -> one row per locus with the left/right breakpoint SDs and
   windows (`inputs/gd_loci.windows.tsv`, the workflow's `loci_tsv`).
2. Task `PrepareWindow` (two per locus) runs `prepare_region_data.sh`: depth slice of `rd_matrix`, VCF slice of the
   contig's VCF (gs:// streamed with pysam, requester pays), BAF (`vcf_to_baf_table.py`) and GT (`extract_gt.sh`)
   tables, samples renamed with `id_map` (`subset_rename_samples.py`).
3. Task `PlotLocus` runs `analyze_gd_locus.sh`: BAF loss scan (`find_baf_terminal_loh.py`; right side of the left SD,
   left side of the right SD), LD without the locus carriers (`compute_ld.py`) and the figure (`plot_two_regions.R`).

## Scripts
### gd_sd_windows.py
```bash
python3 gd_sd_windows.py --calls inputs/gd_calls.txt --bp-sd inputs/GD_NAHR.breakpoint_SDs.bed.gz \
  --segdup inputs/hg38.SegDup.sorted.merged.bed.gz --out inputs/gd_loci.windows.tsv [--flank-frac 0.5]
```
Left SD = breakpoint-SD row ending at the GD start, right SD = row starting at the GD end; loci missing from the
breakpoint-SD table use the merged SD block touching the GD edge, chained outward through gaps <= `--max-gap` (5 kb).
Output columns: locus, chrom, gd_start, gd_end, types, n_calls, highlight_pattern, left_sd_{name,start,end,source},
left_win_{start,end}, right_sd_{name,start,end,source}, right_win_{start,end} (0-based half-open).

### prepare_region_data.sh
```bash
bash prepare_region_data.sh --chr chr16 --start 30088531 --end 30343458 --name 16p11.2_BP5 \
  --rd aou_lrGS_GD_carriers.RD.txt.gz --vcf "gs://.../{chr}.g.vcf.bgz" --samples GD_carrier.samples \
  --id-map id_map.tsv [--requester-pays-project PROJECT] [--outdir .] [--keep-intermediate]
```
Outputs `<name>.renamed.{RD.txt,baf.tsv,GT.tsv}.gz`; warns when a table has rows off `--chr`.

### analyze_gd_locus.sh
```bash
bash analyze_gd_locus.sh --name 16p11.2-BP4-5 --highlight-pattern "GD_16p11\.2\-BP4\-5_(DEL|DUP)_chr16_29638676_30188531" \
  --left-prefix L --left-title 16p11:BP4 --left-window 29166925,29795926 \
  --right-prefix R --right-title 16p11:BP5 --right-window 30111067,30420922 --gd-start 29638676 --gd-end 30188531 \
  --medians GD_samples.RD_median --segdup hg38.SegDup.sorted.merged.bed.gz --bp-sd GD_NAHR.breakpoint_SDs.bed.gz --outdir out
```
Outputs `<name>.breakpoint_profile.pdf`, `<name>.{left,right}_loh.ranked.tsv` (per-sample BAF loss and breakpoint
estimate), clean-site BAF tables and LD tables.

### Helpers
`vcf_to_baf_table.py` (BAF from FORMAT/AD), `extract_gt.sh` (GT table), `subset_rename_samples.py` (keep + rename
samples), `find_baf_terminal_loh.py` (BAF loss at one end of a window, breakpoint estimate), `compute_ld.py` (r^2 of
dosages, binned), `plot_two_regions.R` (two-column figure; `--bin-size`/`--baf-bin` 0 = automatic per window).
Each prints its usage with `-h`.

## Built-in references (`references/`, baked into the image at /opt/gd_scripts/references)
- `gd_loci.windows.tsv`: the default `loci_tsv` (39 loci from `gd_calls.txt`, built by gd_sd_windows.py).
- `gd_calls.txt`: the 59 GD call IDs with carriers it was built from.
- `hg38.SegDup.sorted.merged.bed.gz`: merged hg38 segmental duplications (gatk-sv GenomicContext resource), SD track.
- `GD_NAHR.breakpoint_SDs.bed.gz`: the NAHR breakpoint SD blocks of each GD cluster (e.g. 16p11 BP4, BP5), with the GD
  coordinates that end or start at each block (left/right_GD_candidates); labels the SD breakpoint track and pairs each
  GD with its left and right SD in gd_sd_windows.py.

## Inputs
Data files are not tracked (repository .gitignore). Local copies: `inputs/` next to this README (`gd_calls.txt`, the 59
GD call IDs with carriers; `gd_loci.windows.tsv`, 39 loci; the two SD references) and the inputs template
`wdl_analyses/GDBreakpointProfiles.inputs.json` (fill in the workspace paths of loci_tsv, rd_matrix, rd_medians,
samples, id_map and the SD BEDs, and the requester-pays project). gvcfs: one VCF per contig, e.g.
`gs://prod-drc-broad/longreads/demo_group/work_in_progress/pacbio_unfiltered_small_variants_per_chromosome/chr16.g.vcf.bgz`.
Give all contig VCFs in one run (not one run per contig): SelectLoci reads the loci table, keeps only the contigs that
carry GD loci (`gd_contigs` output: contig, number of loci, matched VCF) and runs their loci; loci on a contig without a
VCF are listed in `skipped_loci`.

## Workflow outputs
- `breakpoint_profiles_tar`: all profile PDFs, `<prefix>_GD_breakpoint_profiles/<chrom>/<prefix>.<locus>.breakpoint_profile.pdf`.
- `breakpoint_tables_tar`: per locus, the ranked BAF-loss tables, the PlotLocus tables tar and the renamed window RD/BAF/GT
  and locus RD tables, `<prefix>_GD_breakpoint_tables/<chrom>/<locus>/`.
- `gd_contigs`, `skipped_loci`.

### fetch_vcf_window.py
One region of a VCF (local or gs://) into an indexed VCF, robust to long remote reads: for gs:// the file is reopened
with a fresh OAuth token every 40 min (tokens last about an hour) and, after a read error ("truncated file"), reopened
and resumed right after the last record written (at most 6 errors). Used by prepare_region_data.sh (step 2).
```
python3 fetch_vcf_window.py --vcf gs://.../chr1.g.vcf.bgz --region chr1:145541830-147568500 --out window.vcf.gz
```

### extract_gd_cluster_data.sh
Depth, BAF and GT for every breakpoint SD of one GD cluster (SD +/- half its length) and depth across the whole span
(first SD start minus its flank to last SD end plus its flank), samples renamed; e.g. 22q11.21 A-D:
```bash
bash extract_gd_cluster_data.sh --cluster 22q11 --breakpoints A,B,C,D --bp-sd references/GD_NAHR.breakpoint_SDs.bed.gz \
  --rd aou_lrGS_GD_carriers.RD.txt.gz --vcf "gs://.../{chr}.g.vcf.bgz" --id-map id_map.tsv --prefix aou_lrGS --outdir 22q11
```
Outputs `<prefix>.<cluster>_<BP>.renamed.{RD.txt,baf.tsv,GT.tsv}.gz`, `<prefix>.<cluster>_<first>-<last>.region.renamed.RD.txt.gz`
and `<prefix>.<cluster>.windows.bed`.
