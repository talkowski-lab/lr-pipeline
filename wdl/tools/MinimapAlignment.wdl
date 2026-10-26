version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MinimapAlignment {
    meta {
        description: [
            "This workflow leverages Minimap2 (https://github.com/lh3/minimap2) in order to align a sample's two haplotype assemblies to a reference."
        ]
    }

    parameter_meta {
        assembly_hap1: "Haplotype 1 assembly."
        assembly_hap2: "Haplotype 2 assembly."
        ref_fa: "From references."
        ref_fai: "From references."
        sample_id: "ID of the sample being aligned."
        minimap_flags: "Parameters to use when running Minimap2."
        minimap_threads: "Number of alignment threads."
        minimap_assembled_bam_hap1: "Aligned haplotype 1 assembly BAM."
        minimap_assembled_bai_hap1: "Index for the haplotype 1 BAM."
        minimap_assembled_paf_hap1: "Haplotype 1 assembly PAF alignment."
        minimap_assembled_bam_hap2: "Aligned haplotype 2 assembly BAM."
        minimap_assembled_bai_hap2: "Index for the haplotype 2 BAM."
        minimap_assembled_paf_hap2: "Haplotype 2 assembly PAF alignment."
    }

    input {
        File assembly_hap1
        File assembly_hap2
        File ref_fa
        File ref_fai
        String prefix

        String sample_id
        String minimap_flags = "-a -x asm20 --cs --eqx"
        Int minimap_threads = 32

        String minimap2_docker

        RuntimeAttr? runtime_attr_align_asm2ref
    }

    call AlignAssembly as AlignHap1 {
        input:
            assembly_fa = assembly_hap1,
            hap = 1,
            flags = minimap_flags,
            threads = minimap_threads,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            sample_id = sample_id,
            prefix = prefix,
            docker = minimap2_docker,
            runtime_attr_override = runtime_attr_align_asm2ref
    }

    call AlignAssembly as AlignHap2 {
        input:
            assembly_fa = assembly_hap2,
            hap = 2,
            flags = minimap_flags,
            threads = minimap_threads,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            sample_id = sample_id,
            prefix = prefix,
            docker = minimap2_docker,
            runtime_attr_override = runtime_attr_align_asm2ref
    }

    output {
        File minimap_assembled_bam_hap1 = AlignHap1.bam
        File minimap_assembled_bai_hap1 = AlignHap1.bai
        File minimap_assembled_paf_hap1 = AlignHap1.paf
        File minimap_assembled_bam_hap2 = AlignHap2.bam
        File minimap_assembled_bai_hap2 = AlignHap2.bai
        File minimap_assembled_paf_hap2 = AlignHap2.paf
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

    String out_prefix = "~{prefix}.~{sample_id}-asm_hap~{hap}.minimap2"
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
        File bam = "~{out_prefix}.bam"
        File bai = "~{out_prefix}.bam.bai"
        File paf = "~{out_prefix}.paf"
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
