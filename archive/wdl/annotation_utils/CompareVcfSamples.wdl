version 1.0

import "../utils/Structs.wdl"

workflow CompareVcfSamples {
    input {
        File vcf
        File vcf_idx
        String prefix

        Array[String] sample_ids

        String utils_docker

        RuntimeAttr? runtime_attr_extract_samples
    }

    call ExtractSamples {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            sample_ids = sample_ids,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_extract_samples
    }

    output {
        File samples_summary_counts = ExtractSamples.summary
        File samples_common = ExtractSamples.common_samples
        File samples_vcf_only = ExtractSamples.vcf_only_samples
        File samples_sample_list_only = ExtractSamples.sample_list_only_samples
    }
}

task ExtractSamples {
    input {
        File vcf
        File vcf_idx
        Array[String] sample_ids
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        # Extract samples from VCF
        bcftools query -l ~{vcf} | grep -v '^$' > ~{prefix}.samples.txt
        
        # Write input sample_ids to a file
        cat > input_samples.txt <<EOF
~{sep='\n' sample_ids}
EOF
        
        # Sort both files for comparison
        sort ~{prefix}.samples.txt > vcf_samples_sorted.txt
        sort input_samples.txt > input_samples_sorted.txt
        
        # Find common samples (in both)
        comm -12 vcf_samples_sorted.txt input_samples_sorted.txt > ~{prefix}.common_samples.txt
        
        # Find samples only in VCF
        comm -23 vcf_samples_sorted.txt input_samples_sorted.txt > ~{prefix}.vcf_only_samples.txt
        
        # Find samples only in input sample_ids
        comm -13 vcf_samples_sorted.txt input_samples_sorted.txt > ~{prefix}.sample_list_only_samples.txt
        
        # Create summary file with counts
        cat > ~{prefix}.summary.txt <<SUMMARY
common_samples: $(wc -l < ~{prefix}.common_samples.txt)
vcf_only_samples: $(wc -l < ~{prefix}.vcf_only_samples.txt)
sample_list_only_samples: $(wc -l < ~{prefix}.sample_list_only_samples.txt)
SUMMARY
    >>>

    output {
        File summary = "~{prefix}.summary.txt"
        File common_samples = "~{prefix}.common_samples.txt"
        File vcf_only_samples = "~{prefix}.vcf_only_samples.txt"
        File sample_list_only_samples = "~{prefix}.sample_list_only_samples.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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
