version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CreateSampleReadCounts {
    input {
        Array[File] mosdepth_bed_files
        Array[File] mosdepth_bed_indices
        File ref_dict
        Array[String] contigs
        String prefix

        Int bin_size
        String sample_id

        String utils_docker

        RuntimeAttr? runtime_attr_bin
        RuntimeAttr? runtime_attr_merge
    }

    scatter (i in range(length(contigs))) {
        call BinMosDepthCounts {
            input:
                mosdepth_bed = mosdepth_bed_files[i],
                mosdepth_bed_idx = mosdepth_bed_indices[i],
                contig = contigs[i],
                bin_size = bin_size,
                prefix = "~{prefix}.~{contigs[i]}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_bin
        }
    }

    call MergeBinnedCounts {
        input:
            binned_counts = BinMosDepthCounts.binned_counts,
            binned_counts_indices = BinMosDepthCounts.binned_counts_idx,
            ref_dict = ref_dict,
            sample_id = sample_id,
            prefix = "~{prefix}.counts",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File binned_read_counts = MergeBinnedCounts.merged_counts
    }
}

task BinMosDepthCounts {
    input {
        File mosdepth_bed
        File mosdepth_bed_idx
        String contig
        Int bin_size
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/scripts/helper/bin_mosdepth.py \
            --input ~{mosdepth_bed} \
            --output ~{prefix}.tsv \
            --bin-size ~{bin_size} \
            --coordinate-system one-based \
            --contig ~{contig} \
            --output-contig ~{contig} \
            --truncate-depth

        bgzip ~{prefix}.tsv
        tabix -s 1 -b 2 -e 3 ~{prefix}.tsv.gz
    >>>

    output {
        File binned_counts = "~{prefix}.tsv.gz"
        File binned_counts_idx = "~{prefix}.tsv.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(mosdepth_bed, "GB")) + 10,
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

task MergeBinnedCounts {
    input {
        Array[File] binned_counts
        Array[File] binned_counts_indices
        File ref_dict
        String sample_id
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        grep "^@" ~{ref_dict} > ~{prefix}.tsv
        echo -e "@RG\tID:GATKCopyNumber\tSM:~{sample_id}" >> ~{prefix}.tsv
        echo -e "CONTIG\tSTART\tEND\tCOUNT" >> ~{prefix}.tsv

        export LC_ALL=C
        zcat ~{sep=' ' binned_counts} | sort -s -k1,1 >> ~{prefix}.tsv

        bgzip ~{prefix}.tsv
    >>>

    output {
        File merged_counts = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(binned_counts, "GB")) + 10,
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
