version 1.0

import "../utils/Structs.wdl"

workflow CreateCramIndex {
    input {
        File cram
        File ref_fa
        File ref_fai
        String gcs_output_dir
        
        String utils_docker
        
        RuntimeAttr? runtime_attr_index
    }

    call RunIndexCram {
        input:
            cram = cram,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            gcs_output_dir = gcs_output_dir,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_index
    }

    output {
        String crai_gcs_path = RunIndexCram.crai_gcs_path
    }
}

task RunIndexCram {
    input {
        File cram
        File ref_fa
        File ref_fai
        String gcs_output_dir
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String base_name = basename(cram)
    String output_dir = sub(gcs_output_dir, "/+$", "")
    String gcs_output_file = output_dir + "/" + base_name + ".crai"

    command <<<
        set -euo pipefail

        samtools --reference ~{ref_fa} index ~{cram} "~{base_name}.crai"

        gsutil -m cp "~{base_name}.crai" "~{gcs_output_file}"
    >>>

    output {
        String crai_gcs_path = gcs_output_file
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(cram, "GB") + size(ref_fa, "GB")) + 10,
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
