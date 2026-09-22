version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MinimapAlignment {
    meta {
        description: [
            "This workflow leverages Minimap2 (https://github.com/lh3/minimap2) in order to align a sample's maternal and paternal assemblies to a reference."
        ]
    }

    parameter_meta {
        assembly_mat: "Maternal assembly."
        assembly_pat: "Paternal assembly."
        ref_fa: "From references."
        ref_fai: "From references."
        sample_id: "ID of the sample being aligned."
        minimap_flags: "Parameters to use when running Minimap2."
        minimap_threads: "Number of alignment threads."
        minimap_assembled_bam_mat: "Aligned maternal-assembly BAM."
        minimap_assembled_bai_mat: "Index for the maternal BAM."
        minimap_assembled_paf_mat: "Maternal-assembly PAF alignment."
        minimap_assembled_bam_pat: "Aligned paternal-assembly BAM."
        minimap_assembled_bai_pat: "Index for the paternal BAM."
        minimap_assembled_paf_pat: "Paternal-assembly PAF alignment."
    }

    input {
        File assembly_mat
        File assembly_pat
        File ref_fa
        File ref_fai
        String prefix

        String sample_id
        String minimap_flags = "-a -x asm20 --cs --eqx"
        Int minimap_threads = 32

        String minimap_docker

        RuntimeAttr? runtime_attr_align_asm2ref
    }

    call AlignAssembly as AlignMat {
        input:
            assembly_fa = assembly_mat,
            hap = 1,
            flags = minimap_flags,
            threads = minimap_threads,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            sample_id = sample_id,
            prefix = prefix,
            docker = minimap_docker,
            runtime_attr_override = runtime_attr_align_asm2ref
    }

    call AlignAssembly as AlignPat {
        input:
            assembly_fa = assembly_pat,
            hap = 2,
            flags = minimap_flags,
            threads = minimap_threads,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            sample_id = sample_id,
            prefix = prefix,
            docker = minimap_docker,
            runtime_attr_override = runtime_attr_align_asm2ref
    }

    output {
        File minimap_assembled_bam_mat = AlignMat.bamOut
        File minimap_assembled_bai_mat = AlignMat.baiOut
        File minimap_assembled_paf_mat = AlignMat.pafOut
        File minimap_assembled_bam_pat = AlignPat.bamOut
        File minimap_assembled_bai_pat = AlignPat.baiOut
        File minimap_assembled_paf_pat = AlignPat.pafOut
    }
}

task AlignAssembly {
    input {
        File assembly_fa
        Int hap
        String flags
        Int threads
        File ref_fa
        File ref_fai
        String sample_id
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String out_prefix = "~{prefix}.~{sample_id}-asm_h~{hap}.minimap2"
    Int mm2_threads = threads - 4

    command <<<
        set -euo pipefail

        minimap2 \
            -t ~{mm2_threads} \
            ~{flags} \
            ~{ref_fa} \
            ~{assembly_fa} \
        | samtools sort -@4 -o "~{out_prefix}.bam"

        samtools index -@3 "~{out_prefix}.bam"

        samtools view -h "~{out_prefix}.bam" \
        | k8 $(which paftools.js) sam2paf \
            -L \
            - \
        > "~{out_prefix}.paf"
    >>>

    output {
        File bamOut = "~{out_prefix}.bam"
        File baiOut = "~{out_prefix}.bam.bai"
        File pafOut = "~{out_prefix}.paf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: threads,
        mem_gb: threads * 2,
        disk_gb: ceil(size(assembly_fa, "GB") + size(ref_fa, "GB")) + 5,
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
