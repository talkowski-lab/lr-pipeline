# Derived from broadinstitute/long-read-pipelines:
# https://github.com/broadinstitute/long-read-pipelines/blob/main/wdl/tasks/Assembly/Hifiasm.wdl

version 1.0

import "../utils/Structs.wdl"

workflow Hifiasm {
    meta {
        description: [
            "This tool assembles a sample's long reads into a haplotype-resolved de novo assembly using hifiasm (https://github.com/chhylp123/hifiasm). Reads are converted to FASTQ, assembled in bubble-phasing mode and the resulting assembly graphs are converted to bgzipped FASTA.",
            "Without parental or Hi-C data the two haplotype assignments are arbitrary and switch between bubbles, so 'hap1' and 'hap2' do not correspond to the maternal and paternal haplotypes. Downstream callers that assume parental phase should not rely on which output a contig came from."
        ]
    }

    parameter_meta {
        bams: "Unaligned BAMs for the sample, one per SMRT cell."
        hifiasm_hap1_fa: "Bgzipped FASTA of the first haplotype assembly."
        hifiasm_hap2_fa: "Bgzipped FASTA of the second haplotype assembly."
        hifiasm_primary_fa: "Bgzipped FASTA of the primary contig assembly."
        hifiasm_hap1_gfa: "Assembly graph for the first haplotype assembly."
        hifiasm_hap2_gfa: "Assembly graph for the second haplotype assembly."
        hifiasm_primary_gfa: "Assembly graph for the primary contig assembly."
        hifiasm_log: "Console log from the hifiasm run, including the inferred coverage histogram."
    }

    input {
        Array[File] bams
        String prefix

        String hifiasm_docker

        RuntimeAttr? runtime_attr_assemble_haplotigs
    }

    call AssembleHaplotigs {
        input:
            bams = bams,
            prefix = prefix,
            docker = hifiasm_docker,
            runtime_attr_override = runtime_attr_assemble_haplotigs
    }

    output {
        File hifiasm_hap1_fa = AssembleHaplotigs.hap1_fa
        File hifiasm_hap2_fa = AssembleHaplotigs.hap2_fa
        File hifiasm_primary_fa = AssembleHaplotigs.primary_fa
        File hifiasm_hap1_gfa = AssembleHaplotigs.hap1_gfa
        File hifiasm_hap2_gfa = AssembleHaplotigs.hap2_gfa
        File hifiasm_primary_gfa = AssembleHaplotigs.primary_gfa
        File hifiasm_log = AssembleHaplotigs.log
    }
}

task AssembleHaplotigs {
    input {
        Array[File] bams
        String prefix
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        for bam in ~{sep=' ' bams}; do
            samtools fastq -@ 4 "$bam"
        done > reads.fastq

        hifiasm \
            -o ~{prefix} \
            -t ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            reads.fastq \
            2>&1 | tee ~{prefix}.hifiasm.log

        rm reads.fastq

        # Convert each assembly graph to bgzipped FASTA, which is what the assembly-based callers consume
        for gfa in ~{prefix}.bp.p_ctg.gfa ~{prefix}.bp.hap1.p_ctg.gfa ~{prefix}.bp.hap2.p_ctg.gfa; do
            awk '/^S/{print ">"$2; print $3}' "$gfa" \
                | bgzip -@ 4 -c > "${gfa%.gfa}.fa.gz"
        done
    >>>

    output {
        File hap1_fa = "~{prefix}.bp.hap1.p_ctg.fa.gz"
        File hap2_fa = "~{prefix}.bp.hap2.p_ctg.fa.gz"
        File primary_fa = "~{prefix}.bp.p_ctg.fa.gz"
        File hap1_gfa = "~{prefix}.bp.hap1.p_ctg.gfa"
        File hap2_gfa = "~{prefix}.bp.hap2.p_ctg.gfa"
        File primary_gfa = "~{prefix}.bp.p_ctg.gfa"
        File log = "~{prefix}.hifiasm.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 48,
        mem_gb: 192,
        disk_gb: 10 * ceil(size(bams, "GB")) + 50,
        boot_disk_gb: 10,
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
