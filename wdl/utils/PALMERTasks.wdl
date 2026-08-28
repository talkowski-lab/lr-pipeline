version 1.0

import "Structs.wdl"

task RunPALMERShard {
    input {
        File bam
        File bai
        String mode
        String mei_type
        File ref_fa
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        dir=$(pwd)

        mv ~{bam} ./
        mv ~{bai} ./
        bam_base=$(basename ~{bam})
        chrom=$(echo $bam_base | sed 's/\.bam$//' | rev | cut -f1 -d '_' | rev)

        mkdir -p "${chrom}"
        /PALMER/PALMER \
            --input ${bam_base} \
            --ref_fa ~{ref_fa} \
            --ref_ver GRCh38 \
            --type ~{mei_type} \
            --mode ~{mode} \
            --output "~{prefix}" \
            --chr $chrom \
            --workdir "${dir}/${chrom}/"

        sed -i "s/$/\t~{mei_type}/" ${chrom}/~{prefix}_calls.txt
        sed -i "s/$/\t~{mei_type}/" ${chrom}/~{prefix}_TSD_reads.txt
        mv ${chrom}/~{prefix}_calls.txt ~{prefix}_calls_shard.txt
        mv ${chrom}/~{prefix}_TSD_reads.txt ~{prefix}_tsd_reads_shard.txt
    >>>

    output {
        File calls_shard = "~{prefix}_calls_shard.txt"
        File tsd_reads_shard = "~{prefix}_tsd_reads_shard.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4,
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

task MergePALMEROutputs {
    input {
        Array[File] calls_shards
        Array[File] tsd_reads_shards
        String mei_type
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        head -n1 ~{calls_shards[0]} > ~{prefix}_~{mei_type}_calls.txt
        for f in ~{sep=' ' calls_shards}; do
            grep -v '^cluster_id' $f >> ~{prefix}_~{mei_type}_calls.txt || true
        done

        head -n1 ~{tsd_reads_shards[0]} > ~{prefix}_~{mei_type}_tsd_reads.txt
        for f in ~{sep=' ' tsd_reads_shards}; do
            grep -v '^cluster_id' $f >> ~{prefix}_~{mei_type}_tsd_reads.txt || true
        done
    >>>

    output {
        File calls = "~{prefix}_~{mei_type}_calls.txt"
        File tsd_reads = "~{prefix}_~{mei_type}_tsd_reads.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(calls_shards, "GB") + size(tsd_reads_shards, "GB")) + 10,
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

task ConvertPALMERToVcf {
    input {
        File palmer_calls
        File palmer_tsd_reads
        String mei_type
        String sample
        File ref_fa
        File ref_fai
        String haplotype
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python /opt/scripts/mei/PALMER_to_vcf.py \
            --palmer_calls ~{palmer_calls} \
            --palmer_tsd_reads ~{palmer_tsd_reads} \
            --mei_type ~{mei_type} \
            --sample ~{sample} \
            --ref_fa ~{ref_fa} \
            --ref_fai ~{ref_fai} \
            --haplotype "~{haplotype}" \
        | bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o ~{prefix}.palmer_calls.vcf.gz

        tabix -p vcf ~{prefix}.palmer_calls.vcf.gz
    >>>

    output {
        File vcf = "~{prefix}.palmer_calls.vcf.gz"
        File vcf_idx = "~{prefix}.palmer_calls.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(palmer_calls, "GB") + size(palmer_tsd_reads, "GB") + size(ref_fa, "GB")) + 10,
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
