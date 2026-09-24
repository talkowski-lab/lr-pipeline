# lrGS_meQTL

Cis-meQTL scanning, scattered per contig, with two interchangeable backends
that share the same VCF + wide-format methylation table inputs:

- **SAIGE** (`MeQTLTasks.wdl`, sparse-GRM null model + SPA test) - handles
  related samples via a sparse GRM; tests one phenotype (site) at a time.
- **tensorQTL** (`MeQTLTensorQTLTasks.wdl`, permutation-based cis mapping) -
  wraps the core task from
  [AoU-Multiomics-Analysis/tensorQTL_cis_permutations](https://github.com/AoU-Multiomics-Analysis/tensorQTL_cis_permutations),
  which normally expects pre-made plink2 files + phenotype bed + covariates;
  here it builds those from a VCF + methylation table instead. Vectorized/GPU,
  no sparse-GRM relatedness correction, tests every qualifying site on a
  contig in one call.

Each backend has a genotype (diploid) and a haplotype workflow:

| Backend | Genotype | Haplotype |
|---|---|---|
| SAIGE | `GenotypeMeQTL_SAIGE.wdl` | `HaplotypeMeQTL_SAIGE.wdl` |
| tensorQTL | `GenotypeMeQTL_tensorQTL.wdl` | `HaplotypeMeQTL_tensorQTL.wdl` |

All four take the same `vcfs` / `methylation_files` / `contigs` / `prefix`
inputs (see below) and the haplotype variant of each backend requires a
fully phased input VCF.

## Common inputs

`vcfs`, `methylation_files`, and `contigs` are parallel arrays - index `i` in
all three must describe the same contig (e.g. `vcfs[2]` and
`methylation_files[2]` are both chr22, `contigs[2] = "chr22"`).

`prefix` (required, all four workflows) is the run name, prefixed onto every
output file: `~{prefix}.~{contig}...`.

## Methylation file format

A gzipped, tab-separated, wide-format bed: `chrom`, `start`, `end`, then one
column per sample (or per haplotype, `sample_hap1`/`sample_hap2`, for the
haplotype workflows), values missing as `.`/`NA`/empty. This matches
`hprc_methylated.chr22.combined.bed.gz` (genotype workflows) and
`hprc_methylated.chr22.haplotype.bed.gz` (haplotype workflows).

Both backends compute each row's call rate as `(non-missing sample values) /
(total sample columns in the file)` and keep rows >= `min_call_rate`.
**Note:** a site's call rate is defined against every sample column present
in that contig's methylation file - if the genotype VCF covers a broader
sample set than the methylation file (as in the exploratory session this
pipeline was built from: 292 genotyped samples vs. 231 with methylation
calls), that's fine and expected.

---

## SAIGE workflows

### Inputs

| Input | Default | Notes |
|---|---|---|
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

`plink_docker`, `bcftools_docker`, and `saige_docker` default to pre-verified
images (`quay.io/biocontainers/plink:1.90b6.21--h031d066_5`,
`quay.io/biocontainers/bcftools:1.19--h8b25389_1`, `wzhou88/saige:1.3.6`) but
can be overridden. Note this deviates from this repo's usual "docker is
never hardcoded" convention, at the user's request, for convenience filling
out the Terra UI.

Every task also takes a `RuntimeAttr? runtime_attr_override`, exposed at the
workflow level as one `runtime_attr_<task>` input per task (e.g.
`runtime_attr_create_grm`), so resources can be tuned per task without
editing the WDL.

### Design decisions worth knowing about

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
  multiallelic variant` in interactive testing. The tensorQTL workflows
  instead split multiallelics upstream (`bcftools norm -m -any`) rather than
  dropping them; the SAIGE workflows don't do this today.
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
  revisiting if sex chromosomes matter for your analysis. tensorQTL's
  haplotype workflow reuses this same task.
- **Output coordinate convention**: a site's "position" for both the cis
  window and the `pheno_pos` output column is its bed `start` (0-based),
  matching the ad hoc convention used in the interactive session this
  pipeline formalizes.

### Outputs

- `combined_assoc`: one TSV across all contigs, all tested sites. Columns:
  `pheno_site_id`, `pheno_chrom`, `pheno_pos`, `n_samples`, then SAIGE's own
  step2 columns (`CHR POS MarkerID Allele1 Allele2 AC_Allele2 AF_Allele2
  MissingRate BETA SE Tstat var p.value N`).
- `per_contig_assoc`: same, split per contig.
- `per_contig_pruned_bed` / `per_contig_sparse_grm`: intermediate plink/GRM
  files, in case you want to reuse them outside this workflow.
- `per_chunk_skipped_logs`: per-chunk skip reasons (nested array: contig x
  chunk).

---

## tensorQTL workflows

### Inputs

| Input | Default | Notes |
|---|---|---|
| `min_call_rate` | 0.9 | Per-site minimum fraction of non-missing samples to test |
| `cis_window` | 1,000,000 | tensorQTL's `--window`, matching the upstream repo's own default |
| `vcf_half_call` | `"missing"` | How plink2 treats GT half-calls |
| `covariates_file` | none | Optional TSV: `person_id` + covariate columns, tab-separated (same shape as the SAIGE workflows' `covariates_file`) - transposed internally into tensorQTL's `covariates x samples` layout. If omitted, an intercept-only (zero-covariate) file is generated automatically |
| `phenotype_groups` / `fdr` / `qvalue_lambda` / `pval_threshold` / `seed` / `flags` | none | Passed straight through to `python3 -m tensorqtl`, same as the upstream repo's workflow |

`plink2_docker`, `bcftools_docker`, `python_docker`, and `tensorqtl_docker`
default to pre-verified images but can be overridden (as with the SAIGE
workflows, this deviates from this repo's usual "docker is never hardcoded"
convention at the user's request): `quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0`,
`quay.io/biocontainers/bcftools:1.19--h8b25389_1`, `wzhou88/saige:1.3.6` for
`python_docker` (already has Python 3.8; any other dependency-free Python 3
image also works), and `gcr.io/broad-cga-francois-gtex/tensorqtl:latest`
(the upstream repo's own image - GPU-enabled, requires `nvidia-tesla-p100`
by default via `tensorqtl_gpu_type`/`tensorqtl_num_gpus`/`tensorqtl_gpu_zones`).

Same `RuntimeAttr? runtime_attr_<task>` pattern as the SAIGE workflows.

### Design decisions worth knowing about

- **No per-site chunking.** tensorQTL is vectorized/GPU-based and reads the
  whole phenotype matrix for a contig in one call, unlike SAIGE's one
  null-model-per-site loop - so there's no `sites_per_shard` equivalent here.
- **Multiallelic sites are split, not dropped**: `NormalizeVcf` runs
  `bcftools norm -m -any` on every contig VCF before plink2 ever sees it.
  This also sidesteps plink2's ~254-ALT-allele import limit on the rare
  complex multiallelic site. This is a deliberate improvement over the SAIGE
  workflows, which silently skip multiallelics instead (see above) - not
  something either pipeline's behavior depends on the other for.
  Chromosome naming (`chr22` vs `22`) is preserved through the pgen
  conversion via plink2's `--output-chr chrM`, since it defaults to
  stripping the `chr` prefix, which would otherwise silently break tensorQTL's
  cis-window matching against the methylation bed's `#chr` column.
- **Missing phenotype values are mean-imputed, not sample-subsetted.**
  tensorQTL needs one rectangular, complete phenotype x sample matrix for
  the whole contig (unlike SAIGE, which fits a fresh null model per site
  over whatever samples that site has); after the call-rate filter, any
  remaining missing values in a qualifying site are filled with that site's
  own mean across present samples.
- **Covariates file is auto-generated, including the no-covariates case.**
  tensorQTL's core task requires a `--covariates` file (unlike SAIGE, where
  `covarColList`/`qCovarColList` can simply be empty); `BuildCovariates`
  always produces one, transposed from `covariates_file` if supplied, or
  header-only (zero covariate rows, intercept-only model) if not.
- **`phenotype_groups` is exposed but not meQTL-specific.** It's tensorQTL's
  mechanism for grouping phenotypes that should share one permutation test
  (e.g. multiple splice junctions per gene, per the upstream repo's sQTL
  workflow) - left as a pass-through optional input in case you want to
  group CpG sites by region, but no methylation-specific grouping is applied
  by default.

### Outputs

- `combined_cis_qtl`: one gzipped TSV across all contigs, tensorQTL's own
  `cis_qtl.txt.gz` columns (one row per tested phenotype/site: top variant,
  permutation-derived p-value, q-value, etc.), headers de-duplicated across
  contigs.
- `per_contig_cis_qtl` / `per_contig_log`: same, split per contig, plus
  tensorQTL's own run log.
- `per_contig_phenotype_bed` / `per_contig_covariates`: the generated
  tensorQTL inputs, in case you want to reuse or inspect them.
- `per_contig_n_sites`: qualifying (call-rate-passing) site count per contig.
