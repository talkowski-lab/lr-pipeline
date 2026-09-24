# lrGS_meQTL

Cis-meQTL scanning with plink (LD-pruning, sparse-GRM marker set) + SAIGE
(sparse-GRM null model + SPA test), scattered per contig.

Two workflows share `MeQTLTasks.wdl`:

- **`GenotypeMeQTL_SAIGE.wdl`** - diploid genotype dosage vs. per-sample methylation.
- **`HaplotypeMeQTL_SAIGE.wdl`** - per-haplotype pseudo-genotype vs. per-haplotype
  methylation. Requires a fully phased input VCF.

## Inputs (both workflows)

`vcfs`, `methylation_files`, and `contigs` are parallel arrays - index `i` in
all three must describe the same contig (e.g. `vcfs[2]` and
`methylation_files[2]` are both chr22, `contigs[2] = "chr22"`).

| Input | Default | Notes |
|---|---|---|
| `prefix` | (required) | Run name, prefixed onto every output file: `~{prefix}.~{contig}...` |
| `min_call_rate` | 0.9 | Per-site minimum fraction of non-missing samples to test |
| `cis_window` | 2,000,000 | +/- bp window around each site's position for the association scan |
| `sites_per_shard` | 200 | How many qualifying sites each `RunCisMeQTLChunk` shard processes serially |
| `ld_prune_window_kb` / `ld_prune_step` / `ld_prune_r2` | 50 / 5 / 0.2 | `plink --indep-pairwise` params |
| `vcf_half_call` | `"missing"` | How plink treats GT half-calls (e.g. `0/.`); long-read phased VCFs can have these |
| `num_random_markers_for_grm` / `relatedness_cutoff` / `min_maf_for_grm` / `max_missing_rate_for_grm` | 2000 / 0.125 / 0.01 / 0.15 | `createSparseGRM.R` params |
| `covariates_file` | none | Optional TSV: `person_id` + covariate columns, tab-separated |
| `covar_col_list` / `qcovar_col_list` | `""` | Comma-separated covariate column names (qcovar = categorical); passed straight to SAIGE |
| `min_samples_per_site` | 20 | Sites with fewer non-missing, covariate-complete samples than this are skipped (logged, not fatal) |
| `inv_normalize` | true | Inverse-normalize the phenotype in step1 (`--invNormalize`) |

`plink_docker`, `bcftools_docker`, and `saige_docker` are required inputs (no
in-WDL default, per repo convention). Pre-verified images to pass in:
`quay.io/biocontainers/plink:1.90b6.21--h031d066_5`,
`quay.io/biocontainers/bcftools:1.19--h8b25389_1`, `wzhou88/saige:1.3.6`.

Every task also takes a `RuntimeAttr? runtime_attr_override`, exposed at the
workflow level as one `runtime_attr_<task>` input per task (e.g.
`runtime_attr_create_grm`), so resources can be tuned per task without
editing the WDL.

## Methylation file format

A gzipped, tab-separated, wide-format bed: `chrom`, `start`, `end`, then one
column per sample (or per haplotype, `sample_hap1`/`sample_hap2`, for
`HaplotypeMeQTL_SAIGE.wdl`), values missing as `.`/`NA`/empty. This matches
`hprc_methylated.chr22.combined.bed.gz` (genotype workflow) and
`hprc_methylated.chr22.haplotype.bed.gz` (haplotype workflow).

`FilterMethylationSites` computes each row's call rate as
`(non-missing sample values) / (total sample columns in the file)` and keeps
rows >= `min_call_rate`. **Note:** a site's call rate is defined against every
sample column present in that contig's methylation file - if the genotype VCF
covers a broader sample set than the methylation file (as in the exploratory
session this pipeline was built from: 292 genotyped samples vs. 231 with
methylation calls), that's fine and expected; per-site sample size is
whatever overlaps after joining phenotype, genotype, and (if supplied)
covariates.

## Design decisions worth knowing about

- **GRM scope is per contig, not genome-wide.** Since the only input
  parallelism is "one VCF per contig," each contig's sparse GRM is built from
  that contig's own LD-pruned markers. If you want a single genome-wide GRM
  shared across contigs instead, that requires restructuring (build the GRM
  once outside the contig scatter, e.g. from a genome-wide pruned marker
  set) - ask if you want that variant.
- **One null model (step1) per site, not per contig.** SAIGE's step1 null
  model is specific to one phenotype; the pruned marker set and sparse GRM
  (both contig-level) are reused unchanged across every site's step1 call.
- **Cis-window restriction is applied via SAIGE's own
  `--rangestoIncludeFile`** on the whole-contig VCF, rather than pre-slicing
  a per-site plink/VCF subset. Functionally equivalent, avoids an extra
  per-site file-prep step.
- **Sites are processed in a bash/python loop inside `RunCisMeQTLChunk`, not
  one WDL task per site.** A contig can have hundreds of thousands of
  qualifying sites; one Cromwell task per site would be impractical. Chunks
  of `sites_per_shard` sites scatter in parallel instead; tune
  `sites_per_shard` down for more parallelism (more, smaller shards) or up
  for less scheduling overhead.
- **Multiallelic sites are silently skipped by SAIGE step2** (its VCF reader
  only handles biallelic records) - this showed up as `Warning: skipping
  multiallelic variant` in interactive testing. If you need multiallelics
  included, they'd need to be split upstream (e.g. `bcftools norm -m-`)
  before this pipeline.
- **Failed or under-powered sites are skipped, not fatal.** Both
  insufficient-sample-count sites and SAIGE step1/step2 failures are logged
  to `skipped_sites.log` (per chunk, surfaced as
  `per_chunk_skipped_logs` in the workflow outputs) and excluded from
  `chunk.assoc.txt`, so one bad site doesn't fail the whole contig.
- **Haplotype encoding**: `SplitPhasedVcfToHaplotypes` rewrites each phased
  record's `sample` genotype (`a|b`) into two homozygous pseudo-diploid
  genotypes (`a/a` for `sample_hap1`, `b/b` for `sample_hap2`), so plink and
  SAIGE - both diploid-oriented - work unchanged on 2x the (pseudo-)samples.
  An unphased genotype (`a/b`) is still split, but the hap1/hap2 assignment
  in that case is arbitrary since there's no real phase information. A
  single-allele GT (e.g. hemizygous chrX/chrY in males) is duplicated onto
  both haplotypes rather than left truly haploid - a simplification worth
  revisiting if sex chromosomes matter for your analysis.
- **Output coordinate convention**: a site's "position" for both the cis
  window and the `pheno_pos` output column is its bed `start` (0-based),
  matching the ad hoc convention used in the interactive session this
  pipeline formalizes.

## Outputs

- `combined_assoc`: one TSV across all contigs, all tested sites. Columns:
  `pheno_site_id`, `pheno_chrom`, `pheno_pos`, `n_samples`, then SAIGE's own
  step2 columns (`CHR POS MarkerID Allele1 Allele2 AC_Allele2 AF_Allele2
  MissingRate BETA SE Tstat var p.value N`).
- `per_contig_assoc`: same, split per contig.
- `per_contig_pruned_bed` / `per_contig_sparse_grm`: intermediate plink/GRM
  files, in case you want to reuse them outside this workflow.
- `per_chunk_skipped_logs`: per-chunk skip reasons (nested array: contig x
  chunk).
