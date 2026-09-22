version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow HiPhase {
    meta {
        description: [
            "This tool runs PacBio HiPhase (https://github.com/PacificBiosciences/HiPhase) to jointly phase a sample's small-variant, SV and (optionally) TRGT VCFs against its aligned reads. It preprocesses and synchronizes the input VCFs per contig, phases them together and optionally haplotags the BAM. It outputs the phased VCF, per-contig phasing statistics and an optional haplotagged BAM."
        ]
    }

    parameter_meta {
        bam: "Aligned reads for the sample."
        bai: "Index for `bam`."
        small_vcf: "Small-variant (SNV/indel) VCF to phase."
        small_vcf_idx: "Index for `small_vcf`."
        sv_vcf: "SV VCF to phase."
        sv_vcf_idx: "Index for `sv_vcf`."
        trgt_vcf: "TRGT tandem-repeat VCF to additionally phase."
        trgt_vcf_idx: "Index for `trgt_vcf`."
        ref_fa: "From references."
        ref_fai: "From references."
        contigs: "Contigs to phase."
        trgt_min_repeat_unit: "Minimum repeat-unit length retained when filtering the TRGT VCF."
        trgt_normalize: "Whether to normalize the TRGT VCF before phasing."
        trgt_min_length_diff: "Minimum length difference retained when filtering the TRGT VCF."
        trgt_max_catalog_length: "Maximum catalog length retained when filtering the TRGT VCF."
        hiphase_extra_args: "Additional arguments passed to HiPhase."
        run_haplotagging: "Whether to also haplotag the BAM."
        hiphase_vcf: "Phased VCF."
        hiphase_vcf_idx: "Index for the phased VCF."
        hiphase_haplotag_files: "Per-contig haplotag read assignments."
        hiphase_stats: "Per-contig phasing statistics."
        hiphase_blocks: "Per-contig phase blocks."
        hiphase_summary: "Per-contig phasing summaries."
        hiphase_haplotagged_bam: "Haplotagged BAM (only when `run_haplotagging`)."
        hiphase_haplotagged_bam_idx: "Index for the haplotagged BAM (only when `run_haplotagging`)."
    }

    input {
        File bam
        File bai
        File small_vcf
        File small_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        File? trgt_vcf
        File? trgt_vcf_idx
        File ref_fa
        File ref_fai
        Array[String] contigs
        String prefix

        Int? trgt_min_repeat_unit
        Boolean? trgt_normalize
        Int? trgt_min_length_diff
        Int? trgt_max_catalog_length
        String? hiphase_extra_args
        Boolean run_haplotagging = false

        String hiphase_docker
        String hiphase_preprocess_docker
        String utils_docker

        RuntimeAttr? runtime_attr_preprocess_vcf
        RuntimeAttr? runtime_attr_subset_vcf_short
        RuntimeAttr? runtime_attr_subset_vcf_sv
        RuntimeAttr? runtime_attr_subset_vcf_trgt
        RuntimeAttr? runtime_attr_process_trgt_vcf
        RuntimeAttr? runtime_attr_sync_contigs
        RuntimeAttr? runtime_attr_hiphase
        RuntimeAttr? runtime_attr_concat_per_contig
        RuntimeAttr? runtime_attr_concat_across_contigs
        RuntimeAttr? runtime_attr_merge_bams
    }

    call PreprocessVCF {
        input:
            vcf = sv_vcf,
            vcf_idx = sv_vcf_idx,
            prefix = "~{prefix}.preprocessed",
            docker = hiphase_preprocess_docker,
            runtime_attr_override = runtime_attr_preprocess_vcf
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetVcfShort {
            input:
                vcf = small_vcf,
                vcf_idx = small_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.short",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf_short
        }

        call Helpers.SubsetVcfToContig as SubsetVcfSv {
            input:
                vcf = PreprocessVCF.preprocessed_vcf,
                vcf_idx = PreprocessVCF.preprocessed_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.sv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf_sv
        }

        call SyncContigs {
            input:
                small_vcf = SubsetVcfShort.subset_vcf,
                small_vcf_idx = SubsetVcfShort.subset_vcf_idx,
                sv_vcf = SubsetVcfSv.subset_vcf,
                sv_vcf_idx = SubsetVcfSv.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.synced",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_sync_contigs
        }

        if (defined(trgt_vcf) && defined(trgt_vcf_idx)) {
            call Helpers.SubsetVcfToContig as SubsetVcfTRGT {
                input:
                    vcf = select_first([trgt_vcf]),
                    vcf_idx = select_first([trgt_vcf_idx]),
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.trgt",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_vcf_trgt
            }

            call Helpers.FilterTRGTVcf {
                input:
                    vcf = SubsetVcfTRGT.subset_vcf,
                    vcf_idx = SubsetVcfTRGT.subset_vcf_idx,
                    min_repeat_unit = trgt_min_repeat_unit,
                    min_length_diff = trgt_min_length_diff,
                    max_catalog_length = trgt_max_catalog_length,
                    ref_fa = ref_fa,
                    normalize = select_first([trgt_normalize, false]),
                    prefix = "~{prefix}.~{contig}.trgt.processed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_process_trgt_vcf
            }
        }

        call RunHiPhase {
            input:
                bam = bam,
                bai = bai,
                unphased_snp_vcf = SubsetVcfShort.subset_vcf,
                unphased_snp_vcf_idx = SubsetVcfShort.subset_vcf_idx,
                unphased_sv_vcf = SyncContigs.synced_vcf,
                unphased_sv_vcf_idx = SyncContigs.synced_vcf_idx,
                unphased_trgt_vcf = FilterTRGTVcf.processed_vcf,
                unphased_trgt_vcf_idx = FilterTRGTVcf.processed_vcf_idx,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                extra_args = hiphase_extra_args,
                run_haplotagging = run_haplotagging,
                prefix = "~{prefix}.~{contig}.phased",
                docker = hiphase_docker,
                runtime_attr_override = runtime_attr_hiphase
        }

        call Helpers.ConcatVcfsLR as ConcatPerContig {
            input:
                vcfs = select_all([RunHiPhase.phased_snp_vcf, RunHiPhase.phased_sv_vcf, RunHiPhase.phased_trgt_vcf]),
                vcf_idxs = select_all([RunHiPhase.phased_snp_vcf_idx, RunHiPhase.phased_sv_vcf_idx, RunHiPhase.phased_trgt_vcf_idx]),
                prefix = "~{prefix}.~{contig}.concat",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_per_contig
        }
    }

    call Helpers.ConcatVcfsLR {
        input:
            vcfs = ConcatPerContig.concat_vcf,
            vcf_idxs = ConcatPerContig.concat_vcf_idx,
            prefix = "~{prefix}.phased",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_across_contigs
    }

    if (run_haplotagging) {
        call Helpers.MergeBams {
            input:
                bams = select_all(RunHiPhase.haplotagged_bam),
                bais = select_all(RunHiPhase.haplotagged_bam_idx),
                prefix = "~{prefix}.haplotagged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_bams
        }
    }

    output {
        File hiphase_vcf = ConcatVcfsLR.concat_vcf
        File hiphase_vcf_idx = ConcatVcfsLR.concat_vcf_idx
        Array[File] hiphase_haplotag_files = RunHiPhase.haplotag_file
        Array[File] hiphase_stats = RunHiPhase.hiphase_stats
        Array[File] hiphase_blocks = RunHiPhase.hiphase_blocks
        Array[File] hiphase_summary = RunHiPhase.hiphase_summary
        File? hiphase_haplotagged_bam = MergeBams.merged_bam
        File? hiphase_haplotagged_bam_idx = MergeBams.merged_bam_idx
    }
}

task PreprocessVCF {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String docker_dir = "/truvari_intrasample"

    command <<<
        set -euo pipefail

        cp ~{docker_dir}/convert_lower_case.py .

        python convert_lower_case.py -i ~{vcf} -o ~{prefix}.uppercase.vcf
        bgzip ~{prefix}.uppercase.vcf ~{prefix}.uppercase.vcf.gz
        tabix -p vcf ~{prefix}.uppercase.vcf.gz

        bcftools +setGT ~{prefix}.uppercase.vcf.gz --no-version -Oz -o ~{prefix}.vcf.gz -- --target-gt a --new-gt u
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File preprocessed_vcf = "~{prefix}.vcf.gz"
        File preprocessed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task SyncContigs {
    input {
        File small_vcf
        File small_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view --header-only ~{sv_vcf} | grep -v "^##contig" > sv_header_no_contigs.txt
        bcftools view --header-only ~{small_vcf} | grep "^##contig" > small_vcf_contigs.txt

        grep "^##fileformat" sv_header_no_contigs.txt > new_header.txt
        cat small_vcf_contigs.txt >> new_header.txt
        grep -v "^##fileformat" sv_header_no_contigs.txt >> new_header.txt

        bcftools reheader --header new_header.txt --output ~{prefix}.vcf.gz ~{sv_vcf}
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File synced_vcf = "~{prefix}.vcf.gz"
        File synced_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(small_vcf, "GB") + size(sv_vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task RunHiPhase {
    input {
        File bam
        File bai
        File unphased_snp_vcf
        File unphased_snp_vcf_idx
        File unphased_sv_vcf
        File unphased_sv_vcf_idx
        File? unphased_trgt_vcf
        File? unphased_trgt_vcf_idx
        File ref_fa
        File ref_fai
        String? extra_args
        Boolean run_haplotagging
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        hiphase \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --bam ~{bam} \
            ~{if run_haplotagging then "--output-bam " + prefix + "_haplotagged.bam" else ""} \
            --reference ~{ref_fa} \
            --global-realignment-cputime 300 \
            --vcf ~{unphased_snp_vcf} \
            --output-vcf ~{prefix}_phased_snp.vcf.gz \
            --vcf ~{unphased_sv_vcf} \
            --output-vcf ~{prefix}_phased_sv.vcf.gz \
            ~{if defined(unphased_trgt_vcf) then "--vcf " + unphased_trgt_vcf + " --output-vcf " + prefix + "_phased_trgt.vcf.gz" else ""} \
            --haplotag-file ~{prefix}_phased_sv_haplotag.tsv \
            --stats-file ~{prefix}.stats.csv \
            --blocks-file ~{prefix}.blocks.tsv \
            --summary-file ~{prefix}.summary.tsv \
            --verbose \
            ~{if defined(extra_args) then extra_args else ""}

        ~{if run_haplotagging then "samtools index " + prefix + "_haplotagged.bam" else ""}
    >>>

    output {
        File phased_snp_vcf = "~{prefix}_phased_snp.vcf.gz"
        File phased_snp_vcf_idx = "~{prefix}_phased_snp.vcf.gz.tbi"
        File phased_sv_vcf   = "~{prefix}_phased_sv.vcf.gz"
        File phased_sv_vcf_idx = "~{prefix}_phased_sv.vcf.gz.tbi"
        File? phased_trgt_vcf = select_first(["~{prefix}_phased_trgt.vcf.gz"])
        File? phased_trgt_vcf_idx =select_first(["~{prefix}_phased_trgt.vcf.gz.tbi"])
        File haplotag_file = "~{prefix}_phased_sv_haplotag.tsv"
        File hiphase_stats = "~{prefix}.stats.csv"
        File hiphase_blocks = "~{prefix}.blocks.tsv"
        File hiphase_summary = "~{prefix}.summary.tsv"
        File? haplotagged_bam = select_first(["~{prefix}_haplotagged.bam"])
        File? haplotagged_bam_idx = select_first(["~{prefix}_haplotagged.bam.bai"])
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(bam, "GB")) + 50,
        boot_disk_gb: 100,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
