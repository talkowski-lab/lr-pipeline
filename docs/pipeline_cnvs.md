# CNV Pipeline
This document describes how to call depth-based CNVs on a new cohort, starting from aligned lrGS reads and finishing with a genotyped cohort CNV VCF. It is organized into two steps: computing each sample's depth, then calling CNVs across the cohort with GATK gCNV. `MosDepth` and `ConcatenateMosDepth` are shared with the [Depth Summary](pipeline.md#depth-summary) of the main pipeline, so their outputs can be reused if they already exist.

Numbered lists are sequential; bulleted lists run in parallel.


## Inputs
Everything below is cohort-specific and must be supplied per run. Reference files and Docker images are supplied separately via Terra workspace data - see [References](references.md) and [Docker images](dockers.md).
- **Aligned lrGS reads** - one aligned BAM and index per sample, plus the sample IDs, as produced by [MinimapReadAlignment](../wdl/tools/MinimapReadAlignment.wdl).
- **PED** - six-column pedigree covering every sample, providing sex (`1` male, `2` female). Required by depth genotyping.
- **Contigs** - `chr1` to `chr22`, `chrX` and `chrY`, in that order. `CreateCohortDepthFiles` estimates ploidy against this exact karyotype, so every step below uses the same list.
- **Bin size** - the width, in bp, of the read-count bins gCNV calls over, e.g. 2000. `CreateSampleReadCounts` and `CreateDepthIntervals` must use the same value.


## 1. Per-Sample Depth
Each sample's per-base depth is computed once, then reshaped into the two forms the cohort step needs: one combined BED for the cohort coverage matrix, and binned read counts for gCNV.
1. **[MosDepth](../wdl/tools/MosDepth.wdl)** - compute per-base coverage per sample, once per contig, with `stream_mode = true` so each contig is streamed straight from the BAM rather than split out first. Leave `bin_size` unset to keep per-base output.
2. Per-sample reshaping, in parallel on the per-contig BEDs:
   - **[ConcatenateMosDepth](../wdl/annotation_utils/ConcatenateMosDepth.wdl)** - concatenate the per-contig BEDs into one indexed per-base BED.
   - **[CreateSampleReadCounts](../wdl/annotation_utils/CreateSampleReadCounts.wdl)** - bin the per-contig BEDs at the bin size into a read-counts file, which is the sample's gCNV depth profile.

Step 1 ends with a combined per-base BED and a binned read-counts file per sample.


## 2. Cohort CNV Calling
The cohort files are built once, then `LongReadCNVs` calls, clusters and genotypes CNVs across every sample.
1. Cohort files, in parallel:
   - **[CreateDepthIntervals](../wdl/annotation_utils/CreateDepthIntervals.wdl)** - write the interval list matching the read-count bins, from the reference index, the contigs and the bin size.
   - **[CreateCohortDepthFiles](../wdl/annotation_utils/CreateCohortDepthFiles.wdl)** - build the cohort binned-coverage matrix, per-sample median coverage and ploidy estimates from the combined per-base BEDs.
2. **[LongReadCNVs](../wdl/tools/LongReadCNVs.wdl)** - call CNVs with GATK gCNV over the intervals from the read-count files, then merge, cluster and genotype them against the binned-coverage matrix and median coverage. Set `sort_depth_profiles = true`, pass the PED as `pedigree`, and pass the contigs as `contigs`. Clustering and genotyping run only on the contigs left in the filtered gCNV intervals, so a partial interval list, e.g. chr15 alone, runs no empty shards. Run in one of two modes:
   - **Cohort mode** - leave `num_training_samples` at `-1`, so every sample is used to fit the contig-ploidy and gCNV models.
   - **Case-cohort mode** - set `num_training_samples` below the cohort size, e.g. 100, so the models are fitted on that many randomly drawn samples and every remaining sample is called against them. Fewer than a few dozen training samples degrades the fitted models.

   GermlineCNVCaller memory grows with samples times intervals per shard, so raising `num_intervals_per_scatter` above its default also needs more memory in `runtime_attr_germline_cnv_caller`.

Step 2 ends with the merged cohort CNV VCF, the per-sample ploidy table, and the genotyped depth VCF with its read-depth table.
