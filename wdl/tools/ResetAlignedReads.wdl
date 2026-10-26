version 1.0

import "../utils/Structs.wdl"

workflow ResetAlignedReads {
    meta {
        description: [
            "This tool turns a sample's aligned long-read BAM back into an unaligned BAM downsampled to a target coverage, for cohorts whose reads are only available aligned. `samtools reset` strips the alignment and the aligner and haplotagging tags and drops secondary and supplementary records, and `samtools view --subsample` keeps a random fraction of read names in the same stream, so the full-depth unaligned BAM is never written. A sample already at or below the target coverage is kept whole.",
            "Coverage is the sequenced bases of primary mapped reads divided by the summed contig lengths in the BAM header. Unlike a mosdepth mean depth it counts soft-clipped and inserted bases and skips deletions, so it can run a few percent higher. Read kinetics, `MM`/`ML` methylation tags and read groups are kept. It runs in the minimap2 image, since `samtools reset` needs samtools 1.17 or later."
        ]
    }

    parameter_meta {
        bam: "Aligned reads for the sample."
        target_coverage: "Coverage to downsample to."
        subsample_seed: "Seed for the read-name subsampling."
        tags_to_remove: "Read tags to strip, besides the alignment itself."
        unaligned_bam: "Unaligned reads downsampled to `target_coverage`."
        coverage: "Coverage of `bam` before downsampling."
        subsample_fraction: "Fraction of reads kept."
    }

    input {
        File bam
        String prefix

        Float target_coverage = 30
        Int subsample_seed = 42
        Array[String] tags_to_remove = ["MD", "NM", "AS", "SA", "XA", "PC", "PS", "HP", "cg", "cm", "cs", "de", "dv", "ms", "nn", "rl", "s1", "s2", "tp", "ts", "mc", "mg", "mi", "rm"]

        String minimap2_docker

        RuntimeAttr? runtime_attr_reset_and_subsample_reads
    }

    call ResetAndSubsampleReads {
        input:
            bam = bam,
            target_coverage = target_coverage,
            subsample_seed = subsample_seed,
            tags_to_remove = tags_to_remove,
            prefix = prefix,
            docker = minimap2_docker,
            runtime_attr_override = runtime_attr_reset_and_subsample_reads
    }

    output {
        File unaligned_bam = ResetAndSubsampleReads.unaligned_bam
        Float coverage = ResetAndSubsampleReads.coverage
        Float subsample_fraction = ResetAndSubsampleReads.subsample_fraction
    }
}

task ResetAndSubsampleReads {
    input {
        File bam
        Float target_coverage
        Int subsample_seed
        Array[String] tags_to_remove
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int threads = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    command <<<
        set -euo pipefail

        # samtools fasta already skips secondary and supplementary records, so -F 0x904 leaves primary mapped reads
        bases=$(samtools fasta -@ ~{threads} -F 0x904 ~{bam} | awk 'NR % 2 == 0 { bases += length($0) } END { printf "%.0f\n", bases }')
        ref_length=$(samtools view -H ~{bam} | awk -F '\t' '$1 == "@SQ" { for (i = 2; i <= NF; i++) if ($i ~ /^LN:/) total += substr($i, 4) } END { printf "%.0f\n", total }')
        awk -v bases="${bases}" -v ref_length="${ref_length}" 'BEGIN { printf "%.2f\n", bases / ref_length }' > coverage.txt
        awk -v coverage="$(cat coverage.txt)" -v target=~{target_coverage} 'BEGIN { f = target / coverage; printf "%.2f\n", f < 1 ? f : 1 }' > fraction.txt

        # Reset and subsample in one stream, so the full-depth unaligned BAM is never written
        samtools reset -@ ~{threads} --remove-tag ~{sep=',' tags_to_remove} ~{bam} \
            | samtools view -@ ~{threads} -b --subsample "$(cat fraction.txt)" --subsample-seed ~{subsample_seed} -o ~{prefix}.unaligned.bam -
    >>>

    output {
        File unaligned_bam = "~{prefix}.unaligned.bam"
        Float coverage = read_float("coverage.txt")
        Float subsample_fraction = read_float("fraction.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(bam, "GB")) + 20,
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
