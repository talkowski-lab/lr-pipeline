version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow HiPhase {
    meta {
        description: [
            "This tool runs PacBio HiPhase (https://github.com/PacificBiosciences/HiPhase) to jointly phase a sample's small-variant, SV and (optionally) TRGT VCFs against its aligned reads. It preprocesses and synchronizes the input VCFs per contig and phases them together. When a TRGT VCF is given, each contig is phased a second time without it, and the small variants and SVs that HiPhase leaves unphased inside a TRGT call (`PF=TR_OVERLAP`) take their phased genotype and phase set from that second pass; the TRGT records excluded from phasing by the TRGT filters are kept in the output unphased. It outputs the phased VCF, per-contig counts of the genotypes filled from the second pass and per-contig phasing statistics."
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
        trgt_min_repeat_unit: "Minimum repeat-unit length a TRGT record needs to be phased; shorter ones are kept in the output unphased."
        trgt_normalize: "Whether to normalize the TRGT VCF before phasing."
        trgt_min_length_diff: "Minimum length difference retained when filtering the TRGT VCF."
        trgt_max_catalog_length: "Maximum catalog length retained when filtering the TRGT VCF."
        hiphase_extra_args: "Additional arguments passed to HiPhase."
        hiphase_vcf: "Phased VCF."
        hiphase_vcf_idx: "Index for the phased VCF."
        hiphase_haplotag_files: "Per-contig haplotag read assignments."
        hiphase_stats: "Per-contig phasing statistics."
        hiphase_blocks: "Per-contig phase blocks."
        hiphase_summary: "Per-contig phasing summaries."
        hiphase_fill_stats: "Per-contig counts of heterozygous, unphased and second-pass-filled genotypes (only when `trgt_vcf` is given)."
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

        String hiphase_docker
        String hiphase_preprocess_docker
        String utils_docker

        RuntimeAttr? runtime_attr_preprocess_vcf
        RuntimeAttr? runtime_attr_subset_vcf_short
        RuntimeAttr? runtime_attr_subset_vcf_sv
        RuntimeAttr? runtime_attr_subset_vcf_trgt
        RuntimeAttr? runtime_attr_process_trgt_vcf
        RuntimeAttr? runtime_attr_sync_contigs
        RuntimeAttr? runtime_attr_subset_bam
        RuntimeAttr? runtime_attr_hiphase
        RuntimeAttr? runtime_attr_hiphase_notrgt
        RuntimeAttr? runtime_attr_fill_tr_overlap_phase
        RuntimeAttr? runtime_attr_concat_per_contig
        RuntimeAttr? runtime_attr_concat_across_contigs
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
                vcf = SubsetVcfSv.subset_vcf,
                vcf_idx = SubsetVcfSv.subset_vcf_idx,
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

            # HiPhase requires every VCF to declare the same contigs, and the small-variant VCF can carry ones TRGT's lacks, e.g. chrEBV
            call SyncContigs as SyncTRGTContigs {
                input:
                    small_vcf = SubsetVcfShort.subset_vcf,
                    small_vcf_idx = SubsetVcfShort.subset_vcf_idx,
                    vcf = FilterTRGTVcf.processed_vcf,
                    vcf_idx = FilterTRGTVcf.processed_vcf_idx,
                    prefix = "~{prefix}.~{contig}.trgt.synced",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_sync_contigs
            }
        }

        # Each shard streams only its contig's reads, so the full BAM is never localized once per contig
        call Helpers.SubsetBamToRegions as SubsetBam {
            input:
                bam = bam,
                bai = bai,
                region_file = write_lines([contig]),
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_bam
        }

        call RunHiPhase {
            input:
                bam = SubsetBam.subset_bam,
                bai = SubsetBam.subset_bai,
                unphased_snp_vcf = SubsetVcfShort.subset_vcf,
                unphased_snp_vcf_idx = SubsetVcfShort.subset_vcf_idx,
                unphased_sv_vcf = SyncContigs.synced_vcf,
                unphased_sv_vcf_idx = SyncContigs.synced_vcf_idx,
                unphased_trgt_vcf = SyncTRGTContigs.synced_vcf,
                unphased_trgt_vcf_idx = SyncTRGTContigs.synced_vcf_idx,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                extra_args = hiphase_extra_args,
                prefix = "~{prefix}.~{contig}.phased",
                docker = hiphase_docker,
                runtime_attr_override = runtime_attr_hiphase
        }

        call Helpers.ConcatVcfsLR as ConcatPerContig {
            input:
                vcfs = select_all([RunHiPhase.phased_snp_vcf, RunHiPhase.phased_sv_vcf, RunHiPhase.phased_trgt_vcf, FilterTRGTVcf.excluded_vcf]),
                vcf_idxs = select_all([RunHiPhase.phased_snp_vcf_idx, RunHiPhase.phased_sv_vcf_idx, RunHiPhase.phased_trgt_vcf_idx, FilterTRGTVcf.excluded_vcf_idx]),
                prefix = "~{prefix}.~{contig}.concat",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_per_contig
        }

        # HiPhase leaves variants inside a TRGT call unphased, so a second pass without TRGT supplies their phase
        if (defined(trgt_vcf) && defined(trgt_vcf_idx)) {
            call RunHiPhase as RunHiPhaseNoTRGT {
                input:
                    bam = SubsetBam.subset_bam,
                    bai = SubsetBam.subset_bai,
                    unphased_snp_vcf = SubsetVcfShort.subset_vcf,
                    unphased_snp_vcf_idx = SubsetVcfShort.subset_vcf_idx,
                    unphased_sv_vcf = SyncContigs.synced_vcf,
                    unphased_sv_vcf_idx = SyncContigs.synced_vcf_idx,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    extra_args = hiphase_extra_args,
                    prefix = "~{prefix}.~{contig}.notrgt.phased",
                    docker = hiphase_docker,
                    runtime_attr_override = runtime_attr_hiphase_notrgt
            }

            call FillTROverlapPhase {
                input:
                    vcf = ConcatPerContig.concat_vcf,
                    vcf_idx = ConcatPerContig.concat_vcf_idx,
                    notrgt_vcfs = [RunHiPhaseNoTRGT.phased_snp_vcf, RunHiPhaseNoTRGT.phased_sv_vcf],
                    notrgt_vcf_idxs = [RunHiPhaseNoTRGT.phased_snp_vcf_idx, RunHiPhaseNoTRGT.phased_sv_vcf_idx],
                    prefix = "~{prefix}.~{contig}.filled",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_fill_tr_overlap_phase
            }
        }

        File contig_vcf = select_first([FillTROverlapPhase.filled_vcf, ConcatPerContig.concat_vcf])
        File contig_vcf_idx = select_first([FillTROverlapPhase.filled_vcf_idx, ConcatPerContig.concat_vcf_idx])
    }

    call Helpers.ConcatVcfsLR {
        input:
            vcfs = contig_vcf,
            vcf_idxs = contig_vcf_idx,
            prefix = "~{prefix}.phased",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_across_contigs
    }

    output {
        File hiphase_vcf = ConcatVcfsLR.concat_vcf
        File hiphase_vcf_idx = ConcatVcfsLR.concat_vcf_idx
        Array[File] hiphase_haplotag_files = RunHiPhase.haplotag_file
        Array[File] hiphase_stats = RunHiPhase.hiphase_stats
        Array[File] hiphase_blocks = RunHiPhase.hiphase_blocks
        Array[File] hiphase_summary = RunHiPhase.hiphase_summary
        Array[File] hiphase_fill_stats = select_all(FillTROverlapPhase.counts_tsv)
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
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view --header-only ~{vcf} | grep -v "^##contig" > header_no_contigs.txt
        bcftools view --header-only ~{small_vcf} | grep "^##contig" > small_vcf_contigs.txt

        grep "^##fileformat" header_no_contigs.txt > new_header.txt
        cat small_vcf_contigs.txt >> new_header.txt
        grep -v "^##fileformat" header_no_contigs.txt >> new_header.txt

        bcftools reheader --header new_header.txt --output ~{prefix}.vcf.gz ~{vcf}
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File synced_vcf = "~{prefix}.vcf.gz"
        File synced_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(small_vcf, "GB") + size(vcf, "GB")) + 10,
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
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        hiphase \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --bam ~{bam} \
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

task FillTROverlapPhase {
    input {
        File vcf
        File vcf_idx
        Array[File] notrgt_vcfs
        Array[File] notrgt_vcf_idxs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

notrgt = {}
for path in ["~{sep='", "' notrgt_vcfs}"]:
    with pysam.VariantFile(path) as vcf_in:
        for record in vcf_in:
            call = record.samples[0]
            gt = call.get("GT")
            if record.id and call.phased and gt and None not in gt and len(set(gt)) > 1 and call.get("PS") is not None:
                notrgt[record.id] = (gt, call["PS"])

n_het = n_unphased = n_filled = 0
with pysam.VariantFile("~{vcf}") as vcf_in:
    sample = list(vcf_in.header.samples)[0]
    with pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=vcf_in.header) as vcf_out:
        for record in vcf_in:
            call = record.samples[0]
            gt = call.get("GT")
            if gt and None not in gt and len(set(gt)) > 1:
                n_het += 1
                if not call.phased:
                    n_unphased += 1
                    match = notrgt.get(record.id) if record.id else None
                    if match is not None:
                        call["GT"] = match[0]
                        call.phased = True
                        call["PS"] = match[1]
                        n_filled += 1
            vcf_out.write(record)

with open("~{prefix}.tsv", "w") as handle:
    handle.write("sample\tn_het\tn_het_unphased\tn_filled\tn_het_unphased_after_fill\n")
    handle.write(f"{sample}\t{n_het}\t{n_unphased}\t{n_filled}\t{n_unphased - n_filled}\n")
CODE

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File filled_vcf = "~{prefix}.vcf.gz"
        File filled_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File counts_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(vcf, "GB") + size(notrgt_vcfs, "GB")) + 10,
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
