# WDL standardization proposal

Status: proposal only; no WDL, Dockstore, or archive files changed.

Audit scope: all 83 active `wdl/**/*.wdl` files and all 53 `archive/wdl/**/*.wdl` files. Comparisons used workflow names, inputs/outputs, call graphs, task bodies/commands, repository references, documentation, and git history. Recommendations are confidence-ranked so mechanical cleanup stays separate from behavior-changing work.

## Naming rules

Use these verbs consistently:

- `Annotate`: add fields/tags to existing records.
- `Create`: derive a new artifact such as a matrix, summary, intervals file, metadata file, or plot set.
- `Convert`: change representation or file format.
- `Extract`: emit selected records as a separate artifact; `Subset`: retain the same representation with fewer records/columns.
- `Concatenate`: join non-overlapping shards/contigs in order.
- `Merge`: reconcile files of the same kind, usually samples or overlapping records.
- `Combine`: apply domain logic across caller/type-specific inputs.
- `Integrate`: combine different variant classes into one callset.
- `Normalize`, `Resolve`, `Filter`, `Summarize`, and `Evaluate`: use when these describe the actual operation better than generic `Process`, `Generate`, or `Plot`.
- Keep thin tool wrappers named for the underlying tool. Add a purpose suffix only when multiple workflows wrap the same tool.
- Keep file name and workflow name identical, PascalCase, with no underscores.

## Proposed renamings

### Active workflows: high confidence

| Current | Proposed | Reason |
|---|---|---|
| `AddEndTRs` | `AnnotateTREndTags` | Adds an `END` annotation; does not add TR records. |
| `CountAnnotations` | `SummarizeAnnotations` | Produces count/list/plotting summaries rather than annotations. |
| `CreateCoverageFile` | `CreateCohortCoverageSummary` | Output contains cohort mean/median/total/threshold statistics, not a generic coverage file. |
| `CreateDepthFiles` | `CreateCohortDepthFiles` | Makes cohort bincov, median coverage, estimated copy-number, and plot artifacts. |
| `CreateIntervalsFile` | `CreateDepthIntervals` | Intervals specifically match depth/read-count bins. |
| `CreateMetadataFile` | `CreateCohortMetadata` | Removes vague `File`; identifies cohort-level artifact. |
| `CreateReadCountsFile` | `CreateSampleReadCounts` | Makes one sample's binned read-count artifact. |
| `DiagnoseSingletons` | `SummarizeSingletonCalls` | Workflow deterministically tabulates calls; it does not run a diagnosis. |
| `DropGenotypes` | `StripGenotypes` | Matches existing reusable task `Helpers.StripGenotypes`. |
| `FlagLowCoverageRegions` | `IdentifyLowCoverageRegions` | Emits low-coverage BEDs/counts; does not flag an input VCF. |
| `GenerateTRGTJson` | `CreateTRGTHistograms` | Current output is a histogram TSV, not JSON. |
| `MergeBackbonePhased` | `FillBackbonePhasedGenotypes` | Pulls phased `GT`/`PS` into remaining unphased genotypes; not a general merge. |
| `PlotPhasingResults` | `EvaluateBackbonePhasing` | Produces accuracy/status tables and currently produces no plots. |
| `PostProcess` | `PostprocessCallset` | Removes globally vague name while preserving established “postprocess” terminology. |
| `TransformAlleleType` | `NormalizeAlleleTypes` | Reclassifies allele types into canonical classes and records subtypes. |
| `VcfToBed` | `ConvertVcfToBed` | Uses standard representation-change verb. |
| `PALMER` | `PALMERAssembly` | This variant specifically runs paired paternal/maternal assembly haplotypes. |
| `PALMERMerge` | `MergePALMERCallsets` | This is merge orchestration, not a PALMER caller run. |
| `HiPhaseMerge` | `MergeHiPhaseCallsets` | Same rationale; clearer distinction from `HiPhase`. |

Do not rename precise active names such as `ConcatenateMosDepth`, `ExtractSampleVcfs`, `FillPhasedGenotypes`, `ResolveHaplotypeOverlaps`, or tool-only wrappers such as `TRGT`, `MosDepth`, `VcfDist`, and `Whatshap` merely to force every output-producing workflow under `Create`.

### Archived workflows retained after deduplication

Apply these only to archived files that remain:

| Current | Proposed | Reason |
|---|---|---|
| `PhaseCallset_vBackbone` | `PhaseCallsetWithBackbone` | Removes forbidden underscore/version-style naming. |
| `PhaseCallset_vCommon` | `PhaseCallsetCommon` | Removes forbidden underscore/version-style naming. |
| `IndexCram` | `CreateCramIndex` | Describes produced artifact. |
| `AgglovarMerge` | `MergeWithAgglovar` | Standard merge verb first. |
| `TruvariMerge` | `MergeWithTruvari` | Standard merge verb first. |
| `StatisticalPhasingPreprocess` | `PreprocessStatisticalPhasing` | Standard action-first form. |
| `CompareBAMs` | `CompareBams` | PascalCase acronym treatment consistent with `Vcf`. |
| `CompareSamplesToVcf` / workflow `ExtractSamplesFromVcf` | `CompareVcfSamples` | Fixes file/workflow mismatch and reflects set comparison outputs. |
| `CreatePEDAncestry` | `CreatePedigreeAndAncestryFiles` | States both outputs; avoids unexplained compressed name. |
| `DownloadConvertBAM` | `CreateFastqFromS3Reads` | Workflow accepts BAM/FASTQ S3 inputs and emits merged FASTQ. |
| `ExtractRegionFromBAM` | `ExtractBamRegion` | Action-first, consistent acronym style. |
| `ParseAbsoluteOrigin` | `NormalizeDuplicationOrigins` | Resolves relative/absolute `ORIGIN` coordinates rather than merely parsing. |
| `RenameInfoFields` | `RenameVcfInfoFields` | Identifies affected representation. |
| `SVAddRawCallers` | `AnnotateSvCallerSupport` | Describes annotation semantics and uses action-first form. |

Four archived file/workflow mismatches must either be removed by deduplication or corrected atomically: `AnnotateExternalAF.wdl` defines `AnnotateExternalAFs`; `FillFormatFields_vBasic.wdl` and `FillFormatFields_vBcfTools.wdl` both define `FillFormatFields`; `CompareSamplesToVcf.wdl` defines `ExtractSamplesFromVcf`.

## Proposed restructurings and generalizations

### 1. Collapse `AnnotateVcf` and `AnnotateVcfCleared` (high confidence)

Keep one `AnnotateVcf` workflow. Add optional modes for:

- clearing specified destination INFO fields before annotation;
- restoring records/alleles from an optional untrimmed VCF;
- existing record sharding, input INFO stripping, TSV sorting/filtering, and sequential TSV annotation.

Both workflows already share contig subsetting, annotation metadata arrays, TSV subsetting, and final concatenation. `AnnotateVcfCleared` should disappear after its clear/swap behavior becomes an optional path in `AnnotateVcf`. Fold `FindUntrimmedAlleles` into this workflow or promote it as an active helper, because the current clear/swap workflow depends on a precomputed file whose creator is archived.

### 2. Extract PALMER shared tasks (high confidence)

Move `RunPALMERShard` and `MergePALMEROutputs` to one importable `wdl/utils/PALMERTasks.wdl` (or a better existing shared module). Their normalized task bodies are byte-identical in `PALMER.wdl` and `PALMERDiploid.wdl`. Generalize `ConvertPALMERToVcf` with a `haplotype` input; its two copies differ mainly because assembly mode passes `1|0`/`0|1`, while diploid-read mode hardcodes `1/1`.

Keep two top-level workflows because paired assembly haplotypes and one diploid read BAM have meaningfully different orchestration. Shared tasks remove duplication without forcing a risky all-purpose PALMER workflow.

### 3. Replace three CNV-local concat tasks with `Helpers.ConcatVcfs` (high confidence)

`GenotypeDepth.wdl`, `DepthPreprocessing.wdl`, and `DepthClustering.wdl` each define `ConcatVCFs`. Replace them with `Helpers.ConcatVcfs`, choosing existing `allow_overlaps` and `naive` inputs to preserve each command's current concat/sort behavior. This also moves these imported CNV modules toward current `RuntimeAttr` conventions.

### 4. Generalize `CombineTRs` for single-sample and cohort inputs (high confidence)

Active `CombineTRs` and archived `CombineCohortTRs` have the same four local tasks and orchestration. Differences are only workflow name, one `String sample_id` versus `Array[String] sample_ids`, and output names. Change active workflow to accept `Array[String] sample_ids` (optionally retain `String? sample_id` for a transition), use one canonical `combined_tr_vcf` output name, then delete archived `CombineCohortTRs`.

### 5. Share mosdepth median-binning implementation (medium confidence)

`CreateDepthFiles`, `CreateReadCountsFile`, and `FlagLowCoverageRegions` independently stream run-length-encoded mosdepth BEDs into fixed bins and calculate base-weighted medians. Extract one tested script or helper task with explicit options for:

- zero-based BED versus one-based interval output;
- dropping or retaining partial terminal bins;
- gzip/tabix output requirements;
- whole-file versus one-contig processing.

Do not combine `CreateCoverageFile`: it aggregates cohort mean/median/threshold fractions across samples and has different semantics. First add parity fixtures for gaps, partial bins, contig transitions, and even-sized medians; current implementations differ at these boundaries.

### 6. Split matching from reporting in old callset-overlap workflows (medium confidence)

Use active `AnnotateCallsetOverlap` as sole matching/annotation engine. If archived plot/summary outputs from `AnnotateCallsetOverlap_AF` remain useful, create a separate `EvaluateCallsetOverlap` workflow that consumes active match TSVs. Remove embedded copies of exact/Truvari/bedtools matching. Archived task-only `annotation/AnnotateCallsetOverlap.wdl` can then be deleted.

### 7. Consolidate archived FORMAT/GQ families before retaining any (medium confidence)

- Keep newest `FillFormatFields.wdl`; delete `_vBasic` and `_vBcfTools` prototypes after a small parity test covering GT, AD, GQ, PL, phasing, and gated field copying.
- Keep `GQCalculateCounts.wdl`; delete older `GQCalculation.wdl`. Git diff shows `GQCalculateCounts` adds caller-aware GQ emissions, fuzzy-match parameters, and revised bins.
- Replace duplicate trio discovery in retained GQ workflows with active `Helpers.FindTrios` (their commands are effectively identical).

### 8. Generalize common archive-only wrappers only if they still need direct execution (medium confidence)

If archived workflows must remain runnable, create one active, parameterized subset/concatenate workflow rather than retaining separate wrappers for contigs, sample lists, and per-sample extraction. If no direct execution requirement exists, delete wrappers and use active helper tasks from callers.

## Proposed archive deletions due to redundancy

### Delete after the active consolidations above (high confidence)

| Archived file | Replacement/evidence |
|---|---|
| `annotation_utils/CombineCohortTRs.wdl` | Generalized active `CombineTRs`; four local tasks are identical. |
| `annotation_utils/AnnotatePostProcess.wdl` | Behavior split into active `PostprocessCallset` and `AnnotateVRS`; its VRS task is identical to active task. |
| `annotation_utils/GQCalculation.wdl` | Superseded by newer, more general `GQCalculateCounts.wdl`. |
| `annotation_utils/FillFormatFields_vBasic.wdl` | Prototype superseded by `FillFormatFields.wdl`. |
| `annotation_utils/FillFormatFields_vBcfTools.wdl` | Prototype superseded by `FillFormatFields.wdl`. |
| `tools/SHAPEITPhase.wdl` | `PhaseCallset_vCommon.wdl` contains the same SHAPEIT task family plus newer shared phasing stages. |
| `tools/TRGTMerge.wdl` | `TRGTMergeContig` command is identical to active `HiPhaseMerge`; active orchestration includes its merge path. |
| `annotation/AnnotateCallsetOverlap.wdl` | Task-only legacy library, not imported; replace any wanted reports as described above. |
| `annotation_utils/PALMERToVcf.wdl` | Conversion is maintained inside active PALMER workflows and should move to shared PALMER tasks. |

### Delete as thin wrappers around active helpers (high confidence if no external Terra method depends on their entrypoint)

| Archived file | Active replacement |
|---|---|
| `annotation_utils/DownloadAWSFile.wdl` | `Helpers.TransferAWSToGCS`. |
| `annotation_utils/FilterTRGTCalls.wdl` | `Helpers.FilterTRGTVcf`. |
| `annotation_utils/SubsetVcfToContigs.wdl` | `Helpers.SubsetVcfToContig` + `Helpers.ConcatVcfs`. |
| `annotation_utils/CombineVcfsAcrossContigs.wdl` | `Helpers.StripGenotypes` when requested + `Helpers.ConcatVcfs`. |
| `annotation_utils/MergeTRs.wdl` | `Helpers.SubsetVcfToContig` + `Helpers.ConcatVcfs`; no TR-specific merge logic. |
| `annotation_utils/CombineVcfs.wdl` | Old helper composition and references obsolete helper naming; use `Helpers.MergeVcfs` or `Helpers.ConcatVcfs` according to overlap semantics. |
| `utils/MergeSplitVCF.wdl` | Old task library covered by active split/subset/concat helpers. |

### Delete only after explicit output-parity tests (medium confidence)

| Archived file | Likely replacement |
|---|---|
| `annotation/AnnotateSingletonReads.wdl` | Active `PostprocessCallset` singleton filtering plus `SummarizeSingletonCalls`. |
| `annotation_utils/AnnotateExternalVariants.wdl` | Active `BedtoolsClosestSV` and `AnnotateCallsetOverlap`. |
| `annotation_utils/AnnotateExternalAF.wdl` | Active overlap annotation TSV + unified `AnnotateVcf`; verify duplicate-DUP handling and final headers. |
| `annotation_utils/IntegrateHGSVCReference.wdl` | Active `IntegrateVcfs`; verify sample swapping and source/size annotations. |
| `annotation_utils/PreprocessGregorVcf.wdl` | Existing normalize/attribute/rename/strip helper chain; retain only if direct entrypoint matters. |
| `annotation_utils/SubsetVcfToPerSample.wdl` | Active `ExtractSampleVcfs` or a generalized subset workflow; verify field-dropping and multi-input behavior. |
| `annotation_utils/SubsetVcfToSamples.wdl` | Active `Helpers.SubsetVcfToSamples` + concat. |
| `annotation_utils/UpdateGenotypes.wdl` | Active `PostprocessCallset` genotype-transfer mode; task commands overlap heavily, but exact match rules and stats outputs need parity checks. |
| `annotation_utils/MakeDepthMetrics.wdl` | Current depth workflow family; git history moved it to archive during the recent depth reorganization. Verify all median-coverage outputs remain represented. |

Do not delete remaining archive workflows merely because they are old. `TransferMethylationTags`, statistical/backbone phasing, `GQCutoffs`, `MergeSites`, `MergeVcfs`, raw-caller FORMAT support, `BenchmarkSTRs`, `CreateDepthProfile`, `CompareBAMs`, `DownloadConvertBAM`, and similar files provide behavior not presently available from one active workflow.

## Implementation order

1. Record external entrypoints: Terra method configurations, saved input JSON keys, and scripts outside this repository. Workflow renames change fully qualified input keys and may break these consumers.
2. Add behavior/parity fixtures for proposed collapses.
3. Extract shared tasks and generalize active workflows without renaming them yet; validate all active WDLs.
4. Perform active renames atomically: file, workflow declaration, imports/calls, `.dockstore.yml`, `docs/workflows.md`, repository-structure docs, and any input templates.
5. Rename retained archive files/workflows and repair relative imports.
6. Delete high-confidence archive duplicates; update `archive/docs/workflows.md` so it does not contain dead links.
7. Run `womtool validate` over active WDLs, Dockstore sync checks, repository-wide old-name search, and focused fixture tests.

For live Terra configurations, choose one compatibility policy before implementation: either a clean breaking rename, or temporary deprecated wrapper workflows retaining old names for one release. Do not silently retain duplicate implementations.
