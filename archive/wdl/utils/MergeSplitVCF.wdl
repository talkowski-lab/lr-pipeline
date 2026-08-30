
version 1.0

import "Structs.wdl"

task SplitFile {
    input {
        File file
        Int shards_per_chunk
        String cohort_prefix
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 0,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
        set -euo pipefail

        split -l ~{shards_per_chunk} ~{file} -a 4 -d "~{cohort_prefix}.shard."
    >>>

    output {
        Array[File] chunks = glob("~{cohort_prefix}.*")
    }
}

task CombineVCFs {
    input {
        Array[File] vcf_files
        String sv_base_mini_docker
        String cohort_prefix
        Boolean sort_after_merge
        Boolean naive=true
        Boolean allow_overlaps = false
        Array[File]? vcf_idx
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_files, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = if !sort_after_merge then 5.0 else 20.0
    
    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 0,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    String merged_vcf_name="~{cohort_prefix}.merged.vcf.gz"
    String sorted_vcf_name="~{cohort_prefix}.merged.sorted.vcf.gz"
    String naive_str = if naive then '-n' else ''
    String overlap_str = if allow_overlaps then '--allow-overlaps' else ''

    command <<<
        set -euo pipefail

        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools concat ~{naive_str} ~{overlap_str} --no-version -Oz --file-list vcfs_sorted.list --output ~{merged_vcf_name}
        
        if [ "~{sort_after_merge}" = "true" ]; then
            mkdir -p tmp
            bcftools sort ~{merged_vcf_name} -Oz --output ~{sorted_vcf_name} -T tmp/
            tabix ~{sorted_vcf_name}
        else 
            tabix ~{merged_vcf_name}
        fi
    >>>

    output {
        File combined_vcf = if sort_after_merge then sorted_vcf_name else merged_vcf_name
        File combined_vcf_idx = if sort_after_merge then sorted_vcf_name + ".tbi" else merged_vcf_name + ".tbi"
    }
}

task CombineVCFsCombineSamples {
    input {
        Array[File] vcf_files
        String merged_filename
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_files, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 0,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, default_attr])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, default_attr.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, default_attr.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, default_attr.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, default_attr.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, default_attr.boot_disk_gb])
    }

    command <<<
        set -euo pipefail

        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        
        for vcf in $(cat vcfs_sorted.list);
        do
            bcftools annotate -x ^FORMAT/GT,FORMAT/AD,FORMAT/DP,FORMAT/GQ,FORMAT/PL -Oz -o "$vcf"_stripped.vcf.gz $vcf
            tabix "$vcf"_stripped.vcf.gz
            echo "$vcf"_stripped.vcf.gz >> vcfs_sorted_stripped.list
        done

        bcftools merge -m none --force-samples --no-version -Oz --file-list vcfs_sorted_stripped.list --output ~{merged_filename}_merged.vcf.gz
        tabix ~{merged_filename}_merged.vcf.gz
    >>>

    output {
        File combined_vcf = "~{merged_filename}_merged.vcf.gz"
        File combined_vcf_idx = "~{merged_filename}_merged.vcf.gz.tbi"
    }
}
