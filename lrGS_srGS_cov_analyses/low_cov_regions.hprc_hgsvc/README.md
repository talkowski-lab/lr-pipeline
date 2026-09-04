# hgsvc_hprc long-read low-coverage blacklist

A cohort-level blacklist of genomic regions that are systematically poorly
covered by long-read genome sequencing (LRS), derived from 292 samples
across the HPRC and HGSVC long-read cohorts.

> **Data availability**: the large intermediate/output data files
> (`hgsvc_hprc.cohort.coverage_counts.tsv`, `hgsvc_hprc.cohort.sample_cutoffs.tsv`,
> and the `hgsvc_hprc.blacklist.p{90,70,50}.raw.bed` / `.annotated.bed` files)
> are excluded from this repo by the repo-wide `.gitignore` (`*.tsv`, `*.bed`
> patterns) and are only present in the local working copy — only the
> figures, scripts, and this README are version-controlled here.

## 1. Data

[mosdepth](https://github.com/brentp/mosdepth) was run on long-read
alignments for 292 samples from HPRC and HGSVC. Coverage was binned into
100bp windows genome-wide, giving mean coverage per 100bp bin per sample.

Per-sample coverage histograms are under [`cov_dist.per_sample/`](cov_dist.per_sample/)
(one PNG per sample, `{sample}.coverage_histogram.png`).

## 2. Per-sample low-coverage cutoff

For each sample, a "low coverage" cutoff was set at **20% of that sample's
own overall median coverage**. A 100bp bin is called low-coverage for a
sample if its mean coverage falls below that sample's cutoff.

Sex chromosomes are handled specially to account for ploidy:

- **Autosomes**: cutoff = 20% of overall median coverage (all samples).
- **chrX**: cutoff = 20% of overall median for female samples (diploid,
  same as autosomes); cutoff = **10%** of overall median for male samples
  (hemizygous — chrX coverage is expected to be roughly half the autosomal
  level in males, so the cutoff is halved to match).
- **chrY**: cutoff = 10% of overall median, **male samples only**. Female
  samples are excluded entirely from chrY (not just given a lenient cutoff)
  since they have no chrY to call coverage on.

Per-sample cutoffs and overall median coverage are in
[`hgsvc_hprc.cohort.sample_cutoffs.tsv`](hgsvc_hprc.cohort.sample_cutoffs.tsv)
(`sample_id`, `cutoff`, `median_coverage`).

Per-100bp-bin cohort summary (every bin genome-wide, low-coverage sample
count out of the number of samples eligible at that bin — 292 everywhere
except chrY, where only male samples are eligible) is in
[`hgsvc_hprc.cohort.coverage_counts.tsv`](hgsvc_hprc.cohort.coverage_counts.tsv)
(`chrom`, `start`, `end`, `low_coverage_sample_count`, `eligible_sample_count`).

## 3. Blacklist definition

A 100bp bin is added to the blacklist if **the fraction of eligible samples
showing low coverage at that bin meets or exceeds a threshold**. Adjacent
qualifying bins are merged (`bedtools merge`) into contiguous regions.

The blacklist is generated at three thresholds, since the "right" cutoff
depends on how conservative a user wants to be:

| threshold | file | regions | total size |
|---|---|---|---|
| ≥90% of eligible samples | `hgsvc_hprc.blacklist.p90.raw.bed` | 1,084 | 172.41 Mb |
| ≥70% of eligible samples | `hgsvc_hprc.blacklist.p70.raw.bed` | 1,472 | 178.71 Mb |
| ≥50% of eligible samples | `hgsvc_hprc.blacklist.p50.raw.bed` | 1,814 | 185.63 Mb |

The three thresholds are **perfectly nested** (p90 ⊆ p70 ⊆ p50 — verified
with `bedtools subtract`, zero bp difference either way) — see
[`hgsvc_hprc.blacklist.cutoff_venn.png`](hgsvc_hprc.blacklist.cutoff_venn.png).
Loosening the threshold from 90% to 50% only adds ~13 Mb (~7%) on top of the
90% blacklist, i.e. most 100bp bins genome-wide are either low-coverage in
nearly *no* samples or nearly *all* of them — there's relatively little
"in-between" territory.

### Annotation

Each blacklist region (at each threshold) is annotated against five
reference tracks, all in the parent directory
[`../`](..) (`low_cov_benchmark/`):

| track | file |
|---|---|
| telomere | `hg38.telomere.bed.gz` |
| centromere | `hg38.centromere.bed.gz` |
| N-masked (assembly gaps) | `hg38.n_masked_regions.bed.gz` |
| segmental duplication | `hg38.SegDup.sorted.merged.bed` |
| simple repeat | `hg38.SimpRep.sorted.merged.bed` |

`hgsvc_hprc.blacklist.p{90,70,50}.annotated.bed` columns:

- `chrom`, `start`, `end` — blacklist region coordinates
- `size` — region length (bp)
- `telomere_bp`, `centromere_bp`, `n_masked_bp`, `segdup_bp`, `simprep_bp` —
  overlap with each track, in bp, computed **independently per track** (a
  region can and often does overlap more than one track at once, so these
  are not mutually exclusive)
- `primary_annotation` — a single dominant-category label per region,
  assigned by priority order `centromere > telomere > N-mask >
  simple_repeat > seg_dup` (first category with nonzero overlap wins;
  `none` if the region matches none of the five tracks — true for only a
  handful of regions, <0.001% of total blacklist bp at every threshold)

`primary_annotation` breakdown (region count / total bp of regions whose
*dominant* category is X — note this is a per-merged-region majority vote,
not additive across thresholds, since loosening the threshold can merge
previously-separate regions of different dominant categories into one):

| category | p90 | p70 | p50 |
|---|---|---|---|
| centromere | 685 / 50.03 Mb | 982 / 54.58 Mb | 1,187 / 75.45 Mb |
| N-mask | 241 / 68.39 Mb | 233 / 68.76 Mb | 227 / 69.02 Mb |
| telomere | 47 / 53.09 Mb | 47 / 53.13 Mb | 46 / 37.21 Mb |
| simple_repeat | 82 / 0.85 Mb | 163 / 2.14 Mb | 265 / 3.82 Mb |
| seg_dup | 26 / 0.04 Mb | 37 / 0.09 Mb | 68 / 0.11 Mb |
| none | 3 / ~0 Mb | 10 / ~0 Mb | 21 / 0.02 Mb |

Essentially the entire blacklist at every threshold is explained by known
repeat/gap structure (N-masked assembly gaps, centromeric satellite, and
telomeric repeat dominate by total bp; centromere dominates by region
*count* since centromeric blacklist regions tend to be numerous but each
individually smaller).

## 4. Genome-wide visualization

Per-chromosome plots are under
[`blacklist_region.per_chr/`](blacklist_region.per_chr/)
(`{chrom}.blacklist.png`), one file per chromosome, each with three stacked
panels (top to bottom: ≥90%, ≥70%, ≥50%). Each panel shows the
low-coverage sample count across the chromosome (blue), the threshold line
for that panel, and the resulting merged blacklist regions shaded in red.

A Venn diagram of the overlap between the three thresholds is at
[`hgsvc_hprc.blacklist.cutoff_venn.png`](hgsvc_hprc.blacklist.cutoff_venn.png).

(The pre-existing [`low_cov_region.per_chr/`](low_cov_region.per_chr/) and
[`cov_dist.per_sample/`](cov_dist.per_sample/) plots show the raw
per-chromosome and per-sample coverage signal without any blacklist
overlay, for reference.)

## 5. Reproducing this pipeline

All generation scripts are under [`scripts/`](scripts/):

- `make_blacklist.sh` — orchestrates the full pipeline end to end
- `annotate_blacklist.py` — annotates a raw blacklist bed with the 5 tracks
  + `primary_annotation` (usage: `python3 annotate_blacklist.py 90 70 50`)
- `downsample_coverage.py` — downsamples `hgsvc_hprc.cohort.coverage_counts.tsv`
  (30.8M 100bp bins) to 5kb-max-aggregated windows for fast plotting
- `plot_blacklist_per_chr.py` — generates the 3-panel per-chromosome plots
- `plot_blacklist_venn.py` — generates the cutoff-overlap Venn diagram

```bash
cd low_cov_regions.hprc_hgsvc
bash scripts/make_blacklist.sh
python3 scripts/plot_blacklist_venn.py
```

Requires `bedtools`, `python3` (`matplotlib`, `matplotlib_venn`, `numpy`).
