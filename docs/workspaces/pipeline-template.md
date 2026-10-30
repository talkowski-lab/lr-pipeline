# Long-Read Pipeline
This workspace takes long-read whole-genome sequencing cohorts from unaligned reads to an annotated cohort VCF: alignment and assembly, SNV/indel, SV, tandem repeat and mobile element calling, cohort integration, phasing, annotation and release filtering. The same pipeline was used to build the HPRC/HGSVC and All of Us long-read callsets.

The lr-pipeline [repository](https://github.com/talkowski-lab/lr-pipeline) documents every workflow. The [pipeline page](https://github.com/talkowski-lab/lr-pipeline/blob/main/docs/pipeline.md) describes each step, and its [Run Order](https://github.com/talkowski-lab/lr-pipeline/blob/main/docs/pipeline.md#run-order) table matches the numbered method configs in this workspace.


## Using This Workspace
- **Data tables** - `sample` holds one row per sample, with its sex and unaligned BAM; per-sample outputs are written back to it. `sample_set` holds one row per cohort; cohort-level outputs are written to it. `contig` holds one row per contig, linked to its cohort; every workflow from step 13 runs on these rows; `contig_set` groups them for the configs that run over every contig at once. The trio-based workflows - `EvaluateNCRCutoffs` in step 27 and the inheritance summaries - read parent-child trios from the PED, so the cohort needs at least one.
- **Workspace data** - every Docker image, reference file and cohort input (PED, ancestry) the method configs read. The `truth_*` entries are truth sets for benchmarking the callset.
- **Method configs** - prefixed with their Run Order step. Configs that share a number can run at the same time once the previous step has finished. Per-sample configs run on the `sample_set` with the `this.samples` expression, cohort configs on the `sample_set` directly, and per-contig configs on the `contig` rows.
- **Unnumbered configs** - the [Optional Analyses](https://github.com/talkowski-lab/lr-pipeline/blob/main/docs/pipeline.md#optional-analyses), which no numbered step needs: `SummarizeAnnotations` and `SummarizeAnnotationsGATKSV` summarize the PASS calls and draw their plots, `CreateIGVPlots` takes IGV snapshots of the reads at `igv_variant_list`, and `Automop` deletes intermediate files no data table references, which it only lists while `dry_run` is `true`.

To run your own cohort, add your samples to `sample`, group them in a `sample_set`, point `ped` and `ancestry` at your files, and run the configs in step order.
