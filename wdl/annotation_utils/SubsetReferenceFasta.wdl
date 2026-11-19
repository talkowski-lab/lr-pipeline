version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SubsetReferenceFasta {
    input {
        File ref_fa
        File ref_fai
        String prefix

        Array[String] contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
        Array[String] minimap_presets = ["asm20", "sr"]

        String ref_docker

        RuntimeAttr? runtime_attr_subset_fasta
        RuntimeAttr? runtime_attr_bwa_index
        RuntimeAttr? runtime_attr_minimap_index
    }

    call SubsetFasta {
        input:
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            contigs = contigs,
            prefix = prefix,
            docker = ref_docker,
            runtime_attr_override = runtime_attr_subset_fasta
    }

    call BwaIndex {
        input:
            ref_fa = SubsetFasta.subset_fa,
            prefix = prefix,
            docker = ref_docker,
            runtime_attr_override = runtime_attr_bwa_index
    }

    scatter (preset in minimap_presets) {
        call MinimapIndex {
            input:
                ref_fa = SubsetFasta.subset_fa,
                preset = preset,
                prefix = prefix,
                docker = ref_docker,
                runtime_attr_override = runtime_attr_minimap_index
        }
    }

    output {
        File subset_fa = SubsetFasta.subset_fa
        File subset_fai = SubsetFasta.subset_fai
        File subset_dict = SubsetFasta.subset_dict
        Array[File] subset_bwa_idx = BwaIndex.bwa_idx
        Array[File] subset_minimap_idx = MinimapIndex.mmi
    }
}

task SubsetFasta {
    input {
        File ref_fa
        File ref_fai
        Array[String] contigs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        ln -s ~{ref_fa} ref.fa
        ln -s ~{ref_fai} ref.fa.fai
        samtools faidx ref.fa ~{sep=' ' contigs} > ~{prefix}.fa
        samtools faidx ~{prefix}.fa
        samtools dict ~{prefix}.fa > ~{prefix}.dict
        cut -f1,2 ~{prefix}.fa.fai >&2
    >>>

    output {
        File subset_fa = "~{prefix}.fa"
        File subset_fai = "~{prefix}.fa.fai"
        File subset_dict = "~{prefix}.dict"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(ref_fa, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task BwaIndex {
    input {
        File ref_fa
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cp ~{ref_fa} ~{prefix}.fa
        bwa index -a bwtsw ~{prefix}.fa
    >>>

    output {
        Array[File] bwa_idx = ["~{prefix}.fa.amb", "~{prefix}.fa.ann", "~{prefix}.fa.bwt", "~{prefix}.fa.pac", "~{prefix}.fa.sa"]
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: 6 * ceil(size(ref_fa, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task MinimapIndex {
    input {
        File ref_fa
        String preset
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        minimap2 -x ~{preset} -t "$(nproc)" -d ~{prefix}.~{preset}.mmi ~{ref_fa}
    >>>

    output {
        File mmi = "~{prefix}.~{preset}.mmi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 32,
        disk_gb: 6 * ceil(size(ref_fa, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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
