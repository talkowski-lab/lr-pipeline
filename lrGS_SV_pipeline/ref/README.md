# Reference Files for Long-Read SV Discovery

This folder contains reference files required to run raw structural variant (SV) discovery methods from long-read whole-genome sequences (lrGS).

---

## Files

### `canonical.chromosomes.txt`
A plain-text list of the 25 canonical chromosomes (`chr1`–`chr22`, `chrX`, `chrY`, `chrM`) in GRCh38.  
Used to restrict variant discovery and downstream analyses to primary chromosomes, excluding alternate/unplaced contigs.

---

### `grch38-noalt.ref-bundle.json`
A JSON bundle pointing to all GRCh38 no-alt analysis set reference resources hosted on GCS.  
Includes:

| Key | Description |
|-----|-------------|
| `fasta` | GRCh38 no-alt reference FASTA (`GCA_000001405.15`) |
| `fai` | FASTA index |
| `dict` | Sequence dictionary |
| `tandem_repeat_bed` | Tandem repeat annotations (TRF) — used by tools such as PBSV and Sniffles to flag TR-overlapping SVs |
| `PAR_bed` | Pseudoautosomal region (PAR) coordinates for chrX/Y |
| `mt_chr_name` | Mitochondrial chromosome name (`chrM`) |
| `size_balanced_scatter_interval_ids` | IDs for size-balanced scatter intervals (for parallelization) |
| `size_balanced_scatter_intervallists_locators` | Paths to size-balanced interval list files |
| `intervallists_autosomes` | Canonical autosome interval list |
| `intervallists_allosomes` | Canonical allosome (chrX/Y) interval list |
| `chromosome_ploidy_priors` | Contig ploidy priors (for GATK-based CNV calling) |
| `haplotype_map` | Haplotype database for fingerprinting / sample QC |

All GCS paths point to `gs://fc-7f861a33-ddb4-4b2f-8d10-5679c9df6108/references/hg38/`.

---

### `pav.minimal.config.json`
Minimal configuration file for **PAV** (Phase Assembly Variant caller), which calls SVs and indels from phased haplotype assemblies.  
Specifies:

| Key | Description |
|-----|-------------|
| `reference` | Path to the reference FASTA used for alignment (`asm/ref.fa`) |
| `asm_pattern` | Glob pattern for per-sample, per-haplotype assembly FASTAs (`asm/{asm_name}/{hap}.fa.gz`) |

See the [PAV documentation](https://github.com/EichlerLab/pav) for full configuration options.
