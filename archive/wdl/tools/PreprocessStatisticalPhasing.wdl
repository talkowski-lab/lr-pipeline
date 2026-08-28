version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow PreprocessStatisticalPhasing {
    input {
        File vcf
        File vcf_idx
        String prefix

        Boolean remove_annotations
        String? annotations_to_remove

        String docker

        RuntimeAttr? runtime_attr_add_score_annotations
        RuntimeAttr? runtime_attr_remove_annotations
        RuntimeAttr? runtime_attr_remove_trgt_overlapping_calls
    }

    call AddScoreAnnotations {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = "~{prefix}.annotated",
            docker = docker,
            runtime_attr_override = runtime_attr_add_score_annotations
    }

    if (remove_annotations) {
        call RemoveAnnotations {
            input:
                vcf = AddScoreAnnotations.annotated_vcf,
                vcf_idx = AddScoreAnnotations.annotated_vcf_idx,
                annotations_to_remove = select_first([annotations_to_remove, ""]),
                prefix = "~{prefix}.annotated.filtered",
                docker = docker,
                runtime_attr_override = runtime_attr_remove_annotations
        }
    }

    call RemoveTRGTOverlappingCalls {
        input:
            vcf = select_first([RemoveAnnotations.filtered_vcf, AddScoreAnnotations.annotated_vcf]),
            vcf_idx = select_first([RemoveAnnotations.filtered_vcf_idx, AddScoreAnnotations.annotated_vcf_idx]),
            prefix = "~{prefix}.annotated.no_trgt_ol",
            docker = docker,
            runtime_attr_override = runtime_attr_remove_trgt_overlapping_calls
    }

    output {
        File annotated_vcf = AddScoreAnnotations.annotated_vcf
        File annotated_vcf_idx = AddScoreAnnotations.annotated_vcf_idx
        File filtered_vcf = RemoveTRGTOverlappingCalls.filtered_vcf
        File filtered_vcf_idx = RemoveTRGTOverlappingCalls.filtered_vcf_idx
    }
}

task AddScoreAnnotations {
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

        def add_score_annotations(input_vcf, output_vcf):
            vcf_in = pysam.VariantFile(input_vcf, "r")
            vcf_in.header.add_meta("INFO", items=[("ID", "SCORE_TRGT_PRIO"),
                                                  ("Number", "1"),
                                                  ("Type", 'Integer'),
                                                  ("Description", "FixVariantCollisions score with TRGT prioritization")])
            vcf_in.header.add_meta("INFO", items=[("ID", "SCORE_DVSV_PRIO"),
                                                  ("Number", "1"),
                                                  ("Type", 'Integer'),
                                                  ("Description", "FixVariantCollisions score with DV+SV integration prioritization")])
            vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)
            for v in vcf_in:
                source = v.info["SOURCE"]
                if source == "DeepVariant":
                    trgt_score = 1
                    dvsv_score = 10
                elif source == "HPRC_SV_Integration":
                    trgt_score = 10
                    dvsv_score = 100
                elif source == "TRGT":
                    trgt_score = 100
                    dvsv_score = 1
                else:
                    print("Error: unrecognized source: " + source + " at " + v.contig + ":" + str(v.pos) + " (ID: " + v.id + ")")
                new_v = v.copy()
                new_v.info["SCORE_TRGT_PRIO"] = trgt_score
                new_v.info["SCORE_DVSV_PRIO"] = dvsv_score
                vcf_out.write(new_v)
            vcf_in.close()
            vcf_out.close()
            print(f"Annotated VCF written to: {output_vcf}")

        input_vcf = "~{vcf}"
        output_vcf = "~{prefix}.vcf"
        add_score_annotations(input_vcf, output_vcf)

        CODE

        bgzip ~{prefix}.vcf
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf, "GiB") * 3),
        disk_gb: 15 + ceil(size(vcf, "GiB") * 3),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
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

task RemoveAnnotations {
    input {
        File vcf
        File vcf_idx
        String annotations_to_remove
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate -x ~{annotations_to_remove} ~{vcf} -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf, "GiB") * 3),
        disk_gb: 15 + ceil(size(vcf, "GiB") * 3),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
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

task RemoveTRGTOverlappingCalls {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools filter -i INFO/TR_ENVELOPED=0 ~{vcf} -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(vcf, "GiB") * 3),
        disk_gb: 15 + ceil(size(vcf, "GiB") * 3),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
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
