version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CreateFastqFromS3Reads {
    input {
        Array[String] addresses
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_download_convert
        RuntimeAttr? runtime_attr_merge_fastq
    }

    scatter (i in range(length(addresses))) {
        call DownloadConvert {
            input:
                address = addresses[i],
                prefix = "~{prefix}.~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_download_convert
        }
    }

    call Helpers.MergeFastq {
        input:
            fastq_gz_files = DownloadConvert.fastq_gz,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_fastq
    }

    output {
        File merged_fastq_gz = MergeFastq.merged_fastq_gz
    }
}

task DownloadConvert {
    input {
        String address
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        aws s3 --no-sign-request cp ~{address} .

        FILE_NAME=$(basename ~{address})
        if [[ ${FILE_NAME} == *.bam ]]; then
            set +o pipefail
            SAMPLE_READ=$(samtools view "${FILE_NAME}" | head -1)
            set -o pipefail
            if echo "${SAMPLE_READ}" | grep -q "MM:Z:"; then
                MOD_TAGS="MM,ML"
            elif echo "${SAMPLE_READ}" | grep -q "Mm:Z:"; then
                MOD_TAGS="Mm,Ml"
            else
                echo "ERROR: No base modification tags (MM/ML or Mm/Ml) found in ${FILE_NAME}" >&2
                exit 1
            fi
            samtools fastq -@ 8 -T ${MOD_TAGS} -n ${FILE_NAME} | gzip > ~{prefix}.fastq.gz
        elif [[ ${FILE_NAME} == *.fastq.gz ]]; then
            mv ${FILE_NAME} ~{prefix}.fastq.gz
        elif [[ ${FILE_NAME} == *.fastq ]]; then
            gzip -c ${FILE_NAME} > ~{prefix}.fastq.gz
        fi
    >>>

    output {
        File fastq_gz = "~{prefix}.fastq.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 8,
        disk_gb: 200,
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
