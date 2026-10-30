# Long-Read CNV Pipeline
This workspace calls depth-based copy number variants from aligned long-read whole-genome sequencing cohorts, with three callers: GATK gCNV, in cohort or case-cohort mode, Sawfish and HiFiCNV.

The lr-pipeline [repository](https://github.com/talkowski-lab/lr-pipeline) documents every workflow. The [CNV pipeline page](https://github.com/talkowski-lab/lr-pipeline/blob/main/docs/pipeline_cnvs.md) describes each step and the choices to make before running it, such as the bin size and the gCNV mode, and its [Run Order](https://github.com/talkowski-lab/lr-pipeline/blob/main/docs/pipeline_cnvs.md#run-order) table matches the numbered method configs in this workspace.


## Using This Workspace
- **Data tables** - `sample` holds one row per sample, with its sex, aligned BAM and the unaligned BAM it was aligned from; per-sample outputs are written back to it. `sample_set` holds one row per cohort; cohort-level outputs are written to it. Samples are treated as unrelated.
- **Workspace data** - every Docker image, reference file and cohort input (PED) the method configs read. `truth_platinum_pedigree_sv_vcf` is a truth set for benchmarking the CNV calls.
- **Method configs** - prefixed with their Run Order step. Configs that share a number can run at the same time once the previous step has finished. Per-sample configs run on the `sample_set` with the `this.samples` expression, and cohort configs on the `sample_set` directly. The `00` configs are only needed for samples without an aligned BAM: `00-ResetAlignedReads` makes an unaligned BAM from reads aligned elsewhere, and `00-MinimapReadAlignment` aligns it. `04-gCNV-Cohort` and `04-gCNV-Case-Cohort` run gCNV in its two modes. `05-CreateCNVDepthPlots` plots `binned_coverage` across the cohort-mode gCNV calls within `cnv_plot_regions_bed`, for review.
- **Unnumbered configs** - the [Optional Analyses](https://github.com/talkowski-lab/lr-pipeline/blob/main/docs/pipeline_cnvs.md#optional-analyses), which no numbered step needs: `Automop` deletes intermediate files no data table references, which it only lists while `dry_run` is `true`.
- **Calling contigs** - `01-CreateDepthIntervals` sets the contigs gCNV calls on, which here are chr2, chr6 and chr22; set its `contigs` to `workspace.contigs` for a whole-genome callset. Sawfish and HiFiCNV always call the whole genome.

To run your own cohort, add your samples to `sample` with their sex and aligned BAMs, group them in a `sample_set`, point `ped` at your PED file, and run the configs in step order.
