version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "MethylationQC.wdl" as MQC
import "CompareMethylationQC.wdl" as CMQC

# DownsampledMethylationQC — PacBio vs ONT methylation QC on one region, with every sample
# downsampled to the same read depth BEFORE the methylation pileup.
#
# MethylationQC's matched metrics downsample each CpG independently, because a bedMethyl
# has no read IDs. Adjacent CpGs then stop sharing molecules, so within-sample metrics such
# as neighbour discordance still depend on each sample's source depth. This workflow
# downsamples reads instead. Per sample it streams one region (default chr21) out of the
# BAM through its index, measures the region's median depth, keeps a random
# target_depth/depth fraction of reads, and reruns the platform's pileup on them (pb-CpG-tools for PacBio,
# modkit for ONT). The pileups go through MethylationQC once per platform and once per
# optional subgroup (e.g. methylation caller), and CompareMethylationQC compares the two
# platforms. Use the raw (not "matched") metrics from these runs, since the reads are
# already depth-matched.
#
# BAMs are read in place and never localized, so only the region is transferred. Samples
# whose region depth is below target_depth are dropped; downsample_summary lists every
# sample's source depth, fraction kept and resulting depth.
workflow DownsampledMethylationQC {
    input {
        Array[String] sample_ids
        Array[String] bams
        Array[String] platforms
        Array[String]? subgroups
        Array[String] subgroup_names = []
        String prefix

        String region = "chr21"
        Float target_depth = 20.0
        Int depth_min_mapq = 1
        Int downsample_seed = 42
        String? requester_pays_project

        File pacbio_ref_fa
        File pacbio_ref_fai
        Int pbcpg_min_mapq = 1
        File ont_ref_fa
        File ont_ref_fai
        String modkit_args = "--cpg --ignore h --combine-strands --mod-thresholds m:0.5"

        File qc_ref_fa
        File qc_ref_fai
        File? cpg_islands_bed
        String pacbio_beta_source = "counts"
        Int min_coverage = 10
        Int matched_coverage = 10
        Int min_samples_per_site = 10
        Int max_neighbor_distance = 50
        Int pca_site_stride = 10

        String samtools_docker = "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
        String pbcpg_docker = "quay.io/pacbio/pb-cpg-tools@sha256:afd5468a423fe089f1437d525fdc19c704296f723958739a6fe226caa01fba1c"
        String modkit_docker = "hangsuunc/modkit@sha256:3f8d7a245229b7207463b734fb5f6c95857a031481a9e98708f569d45d3780c0"
        String methylation_docker = "quay.io/ymostovoy/lr-methylation:latest"

        RuntimeAttr? runtime_attr_downsample
        RuntimeAttr? runtime_attr_pbcpg
        RuntimeAttr? runtime_attr_modkit
    }

    parameter_meta {
        sample_ids: "Sample IDs, unique, aligned by index with bams and platforms."
        bams: "gs:// paths of coordinate-sorted, indexed BAMs with MM/ML tags (index at <bam>.bai or <bam>.csi). Read by streaming the region through the index; never localized."
        platforms: "PacBio or ONT for each sample."
        subgroups: "Optional per-sample label (e.g. methylation caller). MethylationQC runs once more for each label in subgroup_names, on the samples carrying it. A subgroup must be single-platform. Use any other value (e.g. NA) for samples that belong to none."
        subgroup_names: "Subgroup labels to run MethylationQC for."
        prefix: "Basename for all outputs."
        region: "Region to analyse, as a contig (chr21) or contig:start-end. MethylationQC indexes the whole contig."
        target_depth: "Median read depth (over covered bases in region) every sample is downsampled to. Samples below it are dropped."
        depth_min_mapq: "Minimum MAPQ for reads counted in the region depth. The median over covered bases is used: collapsed repeats inflate a mean by up to 2x on chr21, and zero-depth bases (assembly gaps) are skipped."
        downsample_seed: "Seed for samtools --subsample. Reads are kept or dropped by read name, so every alignment of a read is kept together."
        requester_pays_project: "Google project billed when a BAM sits in a requester-pays bucket (e.g. this workspace's project)."
        pacbio_ref_fa: "Reference the PacBio pileup uses (the one the PacBio BAMs were aligned to)."
        pacbio_ref_fai: "Index for pacbio_ref_fa."
        pbcpg_min_mapq: "aligned_bam_to_cpg_scores --min-mapq. 1 matches the original cohort pileups."
        ont_ref_fa: "Reference the ONT pileup uses (the one the ONT BAMs were aligned to)."
        ont_ref_fai: "Index for ont_ref_fa."
        modkit_args: "modkit pileup options besides threads, ref, region and log. The default reproduces the original cohort pileups, including modkit's automatic filter threshold."
        qc_ref_fa: "Uncompressed reference for MethylationQC's CpG index. Use the same one for every platform so the comparison joins on the same sites."
        qc_ref_fai: "Index for qc_ref_fa."
        cpg_islands_bed: "Optional CpG-island BED for island/shore/shelf/open-sea stratification (see MethylationQC)."
        pacbio_beta_source: "MethylationQC pacbio_beta_source. counts (est_mod_count / cov) is computed the same way as ONT's modkit fraction."
        min_coverage: "MethylationQC min_coverage."
        matched_coverage: "MethylationQC matched_coverage. Its matched metrics redraw reads per CpG and are not needed here."
        min_samples_per_site: "MethylationQC and CompareMethylationQC min_samples_per_site."
        max_neighbor_distance: "MethylationQC max_neighbor_distance."
        pca_site_stride: "MethylationQC pca_site_stride. 10 keeps ~43k of chr21's ~430k CpGs."
        samtools_docker: "Image with samtools/htslib built with GCS support, and gcloud."
        pbcpg_docker: "pb-CpG-tools image. The default is 3.0.0_build1, the image used for the original PacBio pileups."
        modkit_docker: "modkit image. The default (modkit 0.4.5) is the image used for the original ONT pileups."
        methylation_docker: "Image containing scripts/methylation/methylation_qc.py."
        runtime_attr_downsample: "Optional RuntimeAttr override. Defaults: 2 CPU, 4 GiB, 50 GB disk, preemptible."
        runtime_attr_pbcpg: "Optional RuntimeAttr override. Defaults: 4 CPU, 16 GiB, preemptible."
        runtime_attr_modkit: "Optional RuntimeAttr override. Defaults: 4 CPU, 8 GiB, preemptible."
    }

    String contig = sub(region, ":.*", "")

    scatter (s in range(length(sample_ids))) {
        call SubsetAndDownsample {
            input:
                bam = bams[s],
                sample_id = sample_ids[s],
                platform = platforms[s],
                region = region,
                target_depth = target_depth,
                depth_min_mapq = depth_min_mapq,
                seed = downsample_seed,
                requester_pays_project = requester_pays_project,
                docker = samtools_docker,
                runtime_attr_override = runtime_attr_downsample
        }

        if (SubsetAndDownsample.meets_target) {
            if (platforms[s] == "PacBio") {
                call PbCpgPileup {
                    input:
                        bam = SubsetAndDownsample.downsampled_bam,
                        bai = SubsetAndDownsample.downsampled_bai,
                        ref_fa = pacbio_ref_fa,
                        ref_fai = pacbio_ref_fai,
                        sample_id = sample_ids[s],
                        min_mapq = pbcpg_min_mapq,
                        docker = pbcpg_docker,
                        runtime_attr_override = runtime_attr_pbcpg
                }
            }
            if (platforms[s] == "ONT") {
                call ModkitPileup {
                    input:
                        bam = SubsetAndDownsample.downsampled_bam,
                        bai = SubsetAndDownsample.downsampled_bai,
                        ref_fa = ont_ref_fa,
                        ref_fai = ont_ref_fai,
                        region = region,
                        sample_id = sample_ids[s],
                        modkit_args = modkit_args,
                        docker = modkit_docker,
                        runtime_attr_override = runtime_attr_modkit
                }
            }
            File pileup_bed = select_first([PbCpgPileup.bed, ModkitPileup.bed])
        }

        if (SubsetAndDownsample.meets_target && platforms[s] == "PacBio") {
            File pacbio_bed = select_first([pileup_bed])
            String pacbio_id = sample_ids[s]
        }
        if (SubsetAndDownsample.meets_target && platforms[s] == "ONT") {
            File ont_bed = select_first([pileup_bed])
            String ont_id = sample_ids[s]
        }
    }

    Array[File] pacbio_beds = select_all(pacbio_bed)
    Array[File] ont_beds = select_all(ont_bed)

    call Helpers.ConcatTsvs as ConcatDownsampleSummary {
        input:
            tsvs = SubsetAndDownsample.summary_tsv,
            prefix = "~{prefix}.downsample_summary",
            preserve_header = true,
            docker = methylation_docker
    }

    if (length(pacbio_beds) > 0) {
        call MQC.MethylationQC as PacBioQC {
            input:
                methylation_beds = pacbio_beds,
                sample_ids = select_all(pacbio_id),
                platform = "PacBio",
                prefix = "~{prefix}.PacBio",
                pacbio_beta_source = pacbio_beta_source,
                ref_fa = qc_ref_fa,
                ref_fai = qc_ref_fai,
                cpg_islands_bed = cpg_islands_bed,
                contigs = [contig],
                min_coverage = min_coverage,
                matched_coverage = matched_coverage,
                min_samples_per_site = min_samples_per_site,
                max_neighbor_distance = max_neighbor_distance,
                pca_site_stride = pca_site_stride,
                methylation_docker = methylation_docker
        }
    }

    if (length(ont_beds) > 0) {
        call MQC.MethylationQC as OntQC {
            input:
                methylation_beds = ont_beds,
                sample_ids = select_all(ont_id),
                platform = "ONT",
                prefix = "~{prefix}.ONT",
                pacbio_beta_source = pacbio_beta_source,
                ref_fa = qc_ref_fa,
                ref_fai = qc_ref_fai,
                cpg_islands_bed = cpg_islands_bed,
                contigs = [contig],
                min_coverage = min_coverage,
                matched_coverage = matched_coverage,
                min_samples_per_site = min_samples_per_site,
                max_neighbor_distance = max_neighbor_distance,
                pca_site_stride = pca_site_stride,
                methylation_docker = methylation_docker
        }
    }

    if (length(pacbio_beds) > 0 && length(ont_beds) > 0) {
        call CMQC.CompareMethylationQC as ComparePlatforms {
            input:
                stats_tar_a = select_first([PacBioQC.stats_tar]),
                stats_tar_b = select_first([OntQC.stats_tar]),
                prefix = "~{prefix}.PacBio_vs_ONT",
                min_samples_per_site = min_samples_per_site,
                methylation_docker = methylation_docker
        }
    }

    if (defined(subgroups) && length(subgroup_names) > 0) {
        Array[String] subgroup_labels = select_first([subgroups])

        scatter (g in subgroup_names) {
            scatter (t in range(length(sample_ids))) {
                if (subgroup_labels[t] == g && SubsetAndDownsample.meets_target[t]) {
                    File member_bed = select_first([pileup_bed[t]])
                    String member_id = sample_ids[t]
                    String member_platform = platforms[t]
                }
            }
            Array[File] member_beds = select_all(member_bed)

            if (length(member_beds) > 0) {
                call MQC.MethylationQC as SubgroupQC {
                    input:
                        methylation_beds = member_beds,
                        sample_ids = select_all(member_id),
                        platform = select_all(member_platform)[0],
                        prefix = "~{prefix}.~{g}",
                        pacbio_beta_source = pacbio_beta_source,
                        ref_fa = qc_ref_fa,
                        ref_fai = qc_ref_fai,
                        cpg_islands_bed = cpg_islands_bed,
                        contigs = [contig],
                        min_coverage = min_coverage,
                        matched_coverage = matched_coverage,
                        min_samples_per_site = min_samples_per_site,
                        max_neighbor_distance = max_neighbor_distance,
                        pca_site_stride = pca_site_stride,
                        methylation_docker = methylation_docker
                }
            }
        }
    }

    output {
        File downsample_summary = ConcatDownsampleSummary.concatenated_tsv
        Array[File?] pileup_beds = pileup_bed
        File? pacbio_stats_tar = PacBioQC.stats_tar
        File? pacbio_figures_tar = PacBioQC.figures_tar
        File? ont_stats_tar = OntQC.stats_tar
        File? ont_figures_tar = OntQC.figures_tar
        File? comparison_stats_tar = ComparePlatforms.stats_tar
        File? comparison_figures_tar = ComparePlatforms.figures_tar
        Array[File?]? subgroup_stats_tars = SubgroupQC.stats_tar
        Array[File?]? subgroup_figures_tars = SubgroupQC.figures_tar
    }
}

task SubsetAndDownsample {
    input {
        String bam
        String sample_id
        String platform
        String region
        Float target_depth
        Int depth_min_mapq
        Int seed
        String? requester_pays_project
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String out_bam = "~{sample_id}.downsampled.bam"

    command <<<
        set -euo pipefail
        case '~{platform}' in
            PacBio|ONT) ;;
            *) echo "platform must be PacBio or ONT, got '~{platform}' for ~{sample_id}" >&2; exit 1 ;;
        esac

        # htslib reads gs:// directly with this token and fetches only the region's blocks.
        GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        export GCS_OAUTH_TOKEN
        ~{"export GCS_REQUESTER_PAYS_PROJECT=" + requester_pays_project}

        samtools view -b -o region.bam '~{bam}' '~{region}'
        samtools index region.bam

        # Median depth over covered bases. A mean is inflated, by up to 2x on chr21 and
        # differently per sample, by very deep collapsed repeats; skipping zero-depth bases
        # leaves out assembly gaps.
        median_depth() {
            samtools depth -Q ~{depth_min_mapq} -r '~{region}' "$1" \
                | awk '$3 > 0 { h[$3]++; n++ }
                       END {
                           if (n == 0) { printf "0"; exit }
                           c = 0
                           for (d = 1; c < n / 2; d++) if (d in h) c += h[d]
                           printf "%d", d - 1
                       }'
        }

        source_depth=$(median_depth region.bam)
        fraction=$(awk -v t=~{target_depth} -v d="$source_depth" \
            'BEGIN { if (d > 0 && d >= t) printf "%.6f", t / d; else printf "NA" }')

        if [ "$fraction" = "NA" ]; then
            echo "~{sample_id}: ~{region} depth ${source_depth}x is below the ~{target_depth}x target; dropping it" >&2
            meets_target=false
            # Header-only placeholder so the output exists; nothing downstream reads it.
            samtools view -H -b -o '~{out_bam}' region.bam
            samtools index '~{out_bam}'
            downsampled_depth=NA
        else
            meets_target=true
            # Reads are kept by a hash of the read name, so all alignments of a read stay together.
            samtools view -b --subsample "$fraction" --subsample-seed ~{seed} -o '~{out_bam}' region.bam
            samtools index '~{out_bam}'
            downsampled_depth=$(median_depth '~{out_bam}')
        fi

        echo "$meets_target" > meets_target.txt
        printf 'sample_id\tplatform\tregion\tsource_median_depth\ttarget_depth\tfraction_kept\tdownsampled_median_depth\tmeets_target\n' \
            > '~{sample_id}.downsample.tsv'
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' '~{sample_id}' '~{platform}' '~{region}' \
            "$source_depth" '~{target_depth}' "$fraction" "$downsampled_depth" "$meets_target" \
            >> '~{sample_id}.downsample.tsv'
    >>>

    output {
        File downsampled_bam = out_bam
        File downsampled_bai = "~{out_bam}.bai"
        File summary_tsv = "~{sample_id}.downsample.tsv"
        Boolean meets_target = read_boolean("meets_target.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task PbCpgPileup {
    input {
        File bam
        File bai
        File ref_fa
        File ref_fai
        String sample_id
        Int min_mapq
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        aligned_bam_to_cpg_scores --version
        aligned_bam_to_cpg_scores \
            --threads "$(nproc)" \
            --bam ~{bam} \
            --ref ~{ref_fa} \
            --output-prefix '~{sample_id}' \
            --min-mapq ~{min_mapq}
    >>>

    output {
        File bed = "~{sample_id}.combined.bed.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: 2 * ceil(size(bam, "GB") + size(ref_fa, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ModkitPileup {
    input {
        File bam
        File bai
        File ref_fa
        File ref_fai
        String region
        String sample_id
        String modkit_args
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        modkit --version
        modkit pileup \
            --threads "$(nproc)" \
            --ref ~{ref_fa} \
            --region '~{region}' \
            --log-filepath '~{sample_id}.modkit.log' \
            ~{modkit_args} \
            ~{bam} '~{sample_id}.bed'
        gzip '~{sample_id}.bed'
    >>>

    output {
        File bed = "~{sample_id}.bed.gz"
        File modkit_log = "~{sample_id}.modkit.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(bam, "GB") + size(ref_fa, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
