version 1.0

# Helpers tasks that have since been removed from wdl/utils/Helpers.wdl, restored from the
# last commit that carried each one so the archived workflows calling them still validate

import "../../../wdl/utils/Structs.wdl"

task CheckSampleMatch {
    input {
        File vcf_a
        File vcf_a_idx
        File vcf_b
        File vcf_b_idx
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf_a} | sort > samples_a.txt
        bcftools query -l ~{vcf_b} | sort > samples_b.txt

        if ! diff -q samples_a.txt samples_b.txt > /dev/null; then
            echo "ERROR: Sample mismatch between VCFs"
            echo "Samples only in VCF A:"
            comm -23 samples_a.txt samples_b.txt
            echo "Samples only in VCF B:"
            comm -13 samples_a.txt samples_b.txt
            exit 1
        fi

        echo "SUCCESS: Samples match between VCFs"
    >>>

    output {
        String status = "success"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf_a, vcf_b], "GB")) + 10,
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

task CreateShardsFromVcfIndex {
    input {
        File vcf_idx
        File ref_fai
        Int shard_bin_size
        Boolean use_ssd = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        touch input.vcf.gz
        ln -sf ~{vcf_idx} input.vcf.gz.tbi
        tabix -l input.vcf.gz > contigs.txt

        python3 <<CODE
import math

contig_lengths = {}
with open("~{ref_fai}", "r") as f:
    for line in f:
        parts = line.strip().split("\t")
        contig_lengths[parts[0]] = int(parts[1])

shard_bin_size = ~{shard_bin_size}

with open("contigs.txt", "r") as f:
    contigs = [line.strip() for line in f if line.strip()]

with open("~{prefix}.txt", "w") as out:
    for contig in contigs:
        if contig not in contig_lengths:
            continue
        max_pos = contig_lengths[contig]
        shard_count = int(math.ceil(max_pos / shard_bin_size))
        for shard_index in range(shard_count):
            start = shard_index * shard_bin_size + 1
            end = min((shard_index + 1) * shard_bin_size, max_pos)
            out.write(f"{contig}:{start}-{end}\n")
CODE
    >>>

    output {
        Array[String] shard_regions = read_lines("~{prefix}.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeFastq {
    input {
        Array[File] fastq_gz_files
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{sep=" " fastq_gz_files} > ~{prefix}.merged.fastq.gz
    >>>

    output {
        File merged_fastq_gz = "~{prefix}.merged.fastq.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(fastq_gz_files, "GB")) + 50,
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

task SubsetBamToContig {
    input {
        File bam
        File bai
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: { localization_optional: true }
        bai: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            -h \
            -b \
            -o ~{prefix}.bam \
            ~{bam} \
            ~{contig}

        samtools index \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{prefix}.bam
    >>>

    output {
        File contig_bam = "~{prefix}.bam"
        File contig_bai = "~{prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 2,
        disk_gb: 20,
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

task SubsetVcfToSampleList {
    input {
        File vcf
        File vcf_idx
        Array[String] samples
        String? extra_args
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat > samples.txt <<EOF
~{sep='\n' samples}
EOF

        bcftools annotate \
            -x INFO/AC,INFO/AN,INFO/AF,INFO/AQ,FORMAT/AD,FORMAT/PL \
            ~{vcf} \
        | bcftools view \
            --samples-file samples.txt \
            --min-ac 1 \
            ~{if defined(extra_args) then extra_args else ""} \
            -Oz -o ~{prefix}.vcf.gz
        
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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

task TransferAWSToGCS {
    input {
        String aws_path
        String output_gcs_path
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        filename=$(basename "~{aws_path}")

        aws s3 cp --no-sign-request "~{aws_path}" "$filename"

        gsutil cp "$filename" "~{output_gcs_path}"
    >>>

    output {
        String gcs_path = output_gcs_path
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 600,
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
