version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergeHiPhaseCallsets {
    input {
        Array[File] phased_vcfs
        Array[File] phased_vcf_idxs
        File ref_fa
        File ref_fai
        Array[String] contigs
        String prefix

        Boolean merge_trgt
        String merge_args = "--merge id"

        String trgt_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_trgt
        RuntimeAttr? runtime_attr_drop_fields
        RuntimeAttr? runtime_attr_subset_integrated
        RuntimeAttr? runtime_attr_add_trgt_end
        RuntimeAttr? runtime_attr_extract_trgt_ps
        RuntimeAttr? runtime_attr_fix_al_header
        RuntimeAttr? runtime_attr_merge_trgt
        RuntimeAttr? runtime_attr_annotate_trgt_ps
        RuntimeAttr? runtime_attr_merge_integrated
        RuntimeAttr? runtime_attr_concat_trgt
        RuntimeAttr? runtime_attr_concat_integrated
    }

    if (merge_trgt) {
        scatter (i in range(length(phased_vcfs))) {
            call Helpers.SubsetVcfByArgs as SubsetIntegrated {
                input:
                    vcf = phased_vcfs[i],
                    vcf_idx = phased_vcf_idxs[i],
                    exclude_args = 'INFO/allele_type = "trv"',
                    prefix = "~{prefix}.~{i}.integrated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_integrated
            }

            call Helpers.DropVcfFields {
                input:
                    vcf = SubsetIntegrated.subset_vcf,
                    vcf_idx = SubsetIntegrated.subset_vcf_idx,
                    drop_fields = "FORMAT/AD,FORMAT/PL",
                    prefix = "~{prefix}.~{i}.integrated.no_ad",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_drop_fields
            }

            call Helpers.SubsetVcfByArgs as SubsetTRGT {
                input:
                    vcf = phased_vcfs[i],
                    vcf_idx = phased_vcf_idxs[i],
                    include_args = 'INFO/allele_type = "trv"',
                    prefix = "~{prefix}.~{i}.trgt",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_trgt
            }

            call AddTRGTEndTag {
                input:
                    vcf = SubsetTRGT.subset_vcf,
                    vcf_idx = SubsetTRGT.subset_vcf_idx,
                    prefix = "~{prefix}.~{i}.trgt.with_end",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_add_trgt_end
            }

            call ExtractTRGTPSTags {
                input:
                    vcf = AddTRGTEndTag.vcf_with_end,
                    vcf_idx = AddTRGTEndTag.vcf_with_end_idx,
                    prefix = "~{prefix}.~{i}.trgt",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_trgt_ps
            }

            call FixALHeader {
                input:
                    vcf = AddTRGTEndTag.vcf_with_end,
                    vcf_idx = AddTRGTEndTag.vcf_with_end_idx,
                    prefix = "~{prefix}.~{i}.trgt.fixed_al",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_fix_al_header
            }
        }
    }

    scatter (contig in contigs) {
        call Helpers.MergeVcfs as MergeIntegratedVcfs {
            input:
                vcfs = select_first([DropVcfFields.dropped_vcf, phased_vcfs]),
                vcf_idxs = select_first([DropVcfFields.dropped_vcf_idx, phased_vcf_idxs]),
                prefix = "~{prefix}.~{contig}.integrated",
                contig = contig,
                extra_args = merge_args,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_integrated
        }

        if (merge_trgt) {
            call TRGTMergeContig {
                input:
                    vcfs = select_first([FixALHeader.fixed_vcf]),
                    vcf_idxs = select_first([FixALHeader.fixed_vcf_idx]),
                    contig = contig,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    prefix = "~{prefix}.~{contig}.trgt",
                    docker = trgt_docker,
                    runtime_attr_override = runtime_attr_merge_trgt
            }

            call AnnotateTRGTPS {
                input:
                    merged_vcf = TRGTMergeContig.merged_vcf,
                    merged_vcf_idx = TRGTMergeContig.merged_vcf_idx,
                    ps_tsvs = select_first([ExtractTRGTPSTags.ps_tsv]),
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.trgt",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_annotate_trgt_ps
            }
        }
    }

    call Helpers.ConcatVcfs as ConcatIntegratedVcfs {
        input:
            vcfs = MergeIntegratedVcfs.merged_vcf,
            vcf_idxs = MergeIntegratedVcfs.merged_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.integrated",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_integrated
    }

    if (merge_trgt) {
        call Helpers.ConcatVcfs as ConcatTRGTVcfs {
            input:
                vcfs = select_all(AnnotateTRGTPS.annotated_vcf),
                vcf_idxs = select_all(AnnotateTRGTPS.annotated_vcf_idx),
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.trgt",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_trgt
        }
    }

    output {
        File hiphase_merged_integrated_vcf = ConcatIntegratedVcfs.concat_vcf
        File hiphase_merged_integrated_vcf_idx = ConcatIntegratedVcfs.concat_vcf_idx
        File? hiphase_merged_trgt_vcf = ConcatTRGTVcfs.concat_vcf
        File? hiphase_merged_trgt_vcf_idx = ConcatTRGTVcfs.concat_vcf_idx
    }
}

task ExtractTRGTPSTags {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

with pysam.VariantFile("~{vcf}") as vcf_in:
    sample = next(iter(vcf_in.header.samples))
    with open("~{prefix}.ps.tsv", "w") as out:
        out.write(f"#SAMPLE\t{sample}\n")
        for rec in vcf_in:
            ps_val = rec.samples[sample].get("PS")
            if ps_val is not None:
                out.write(f"{rec.chrom}\t{rec.pos}\t{ps_val}\n")
CODE
    >>>

    output {
        File ps_tsv = "~{prefix}.ps.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 2,
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

task AddTRGTEndTag {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} \
        | awk 'BEGIN { FS = OFS = "\t"; has_end_header = 0 }
            /^##INFO=<ID=END,/ {
                has_end_header = 1
                print
                next
            }
            /^#CHROM/ {
                if (!has_end_header) {
                    print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"
                }
                print
                next
            }
            /^#/ {
                print
                next
            }
            {
                end_val = $2 + length($4) - 1

                if ($8 == "." || $8 == "") {
                    $8 = "END=" end_val
                } else if ($8 ~ /(^|;)END=/) {
                    n = split($8, info_parts, ";")
                    for (i = 1; i <= n; i++) {
                        if (info_parts[i] ~ /^END=/) {
                            info_parts[i] = "END=" end_val
                        }
                    }
                    $8 = info_parts[1]
                    for (i = 2; i <= n; i++) {
                        $8 = $8 ";" info_parts[i]
                    }
                } else {
                    $8 = $8 ";END=" end_val
                }

                print
            }' \
        | bgzip -c > ~{prefix}.vcf.gz

        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File vcf_with_end = "~{prefix}.vcf.gz"
        File vcf_with_end_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size([vcf, vcf_idx], "GB")) + 5,
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

task FixALHeader {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -h ~{vcf} \
        | sed 's/^##FORMAT=<ID=AL,.*$/##FORMAT=<ID=AL,Number=.,Type=Integer,Description="Length of each allele">/' \
        > new_header.txt

        bcftools reheader -h new_header.txt ~{vcf} \
        | bcftools view -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File fixed_vcf = "~{prefix}.vcf.gz"
        File fixed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, vcf_idx], "GB")) + 10,
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

task TRGTMergeContig {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String contig
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo "~{sep='\n' vcfs}" > vcf_list.txt

        trgt merge \
            --vcf-list vcf_list.txt \
            --contig ~{contig} \
            --genome ~{ref_fa} \
            --output-type z \
            --output ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 10,
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

task AnnotateTRGTPS {
    input {
        File merged_vcf
        File merged_vcf_idx
        Array[File] ps_tsvs
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo "~{sep='\n' ps_tsvs}" > ps_tsv_list.txt

        while IFS= read -r tsv_file; do
            awk -F '\t' -v contig="~{contig}" \
                'NR==1{sample=$2} NR>1 && $1==contig{print sample"\t"$2"\t"$3}' \
                "$tsv_file"
        done < ps_tsv_list.txt > contig_ps.tsv

        python3 <<CODE
import pysam
from collections import defaultdict

ps_lookup = defaultdict(dict)

with open("contig_ps.tsv") as f:
    for line in f:
        sample_name, pos_str, ps_val = line.rstrip("\n").split("\t")
        ps_lookup[int(pos_str)][sample_name] = int(ps_val)

with pysam.VariantFile("~{merged_vcf}") as vcf_in:
    header = vcf_in.header
    with pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=header) as vcf_out:
        for rec in vcf_in:
            if rec.pos in ps_lookup:
                for sample_name, ps_val in ps_lookup[rec.pos].items():
                    if sample_name in rec.samples:
                        rec.samples[sample_name]["PS"] = ps_val
            vcf_out.write(rec)
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(merged_vcf, "GB")) + ceil(size(ps_tsvs, "GB")) + 5,
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
