version 1.0

import "utils/Structs.wdl"

workflow MergeWithAgglovar {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs

        String prefix
        String run_agglovar_merge_script = "https://raw.githubusercontent.com/talkowski-lab/gnomad-lr/main/scripts/agglovar/run_agglovar_merge.py"
        
        Float? ro_min
        Float? size_ro_min
        Int? offset_max
        Float? offset_prop_max
        Boolean match_ref = false
        Boolean match_alt = false
        Float? match_prop_min

        String agglovar_docker

        RuntimeAttr? runtime_attr_run_agglovar_merge
    }

    call RunAgglovarMerge {
        input:
            vcfs = vcfs,
            vcf_idxs = vcf_idxs,
            prefix = prefix,
            run_agglovar_merge_script = run_agglovar_merge_script,
            ro_min = ro_min,
            size_ro_min = size_ro_min,
            offset_max = offset_max,
            offset_prop_max = offset_prop_max,
            match_ref = match_ref,
            match_alt = match_alt,
            match_prop_min = match_prop_min,
            docker = agglovar_docker,
            runtime_attr_override = runtime_attr_run_agglovar_merge
    }

    output {
        File merged_vcf = RunAgglovarMerge.merged_vcf
        File merged_vcf_index = RunAgglovarMerge.merged_vcf_index
    }
}

task RunAgglovarMerge {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String run_agglovar_merge_script
        String docker

        Float? ro_min
        Float? size_ro_min
        Int? offset_max
        Float? offset_prop_max
        Boolean match_ref = false
        Boolean match_alt = false
        Float? match_prop_min

        RuntimeAttr? runtime_attr_override
    }

    String ro_min_arg = if defined(ro_min) then "--ro_min ~{ro_min}" else ""
    String size_ro_min_arg = if defined(size_ro_min) then "--size_ro_min ~{size_ro_min}" else ""
    String offset_max_arg = if defined(offset_max) then "--offset_max ~{offset_max}" else ""
    String offset_prop_max_arg = if defined(offset_prop_max) then "--offset_prop_max ~{offset_prop_max}" else ""
    String match_ref_arg = if match_ref then "--match_ref" else ""
    String match_alt_arg = if match_alt then "--match_alt" else ""
    String match_prop_min_arg = if defined(match_prop_min) then "--match_prop_min ~{match_prop_min}" else ""

    command <<<
        set -euo pipefail

        curl -s ~{run_agglovar_merge_script} > run_agglovar_merge.py
        chmod +x run_agglovar_merge.py

        for vcf in ~{sep=' ' vcfs}; do
            if [ ! -f "${vcf}.tbi" ] && [ ! -f "${vcf}.csi" ]; then
                bcftools index -t "${vcf}"
            fi
        done

        python run_agglovar_merge.py \
            --vcfs ~{sep=' ' vcfs} \
            --out_vcf "~{prefix}.agglovar.merged.vcf" \
            ~{ro_min_arg} \
            ~{size_ro_min_arg} \
            ~{offset_max_arg} \
            ~{offset_prop_max_arg} \
            ~{match_ref_arg} \
            ~{match_alt_arg} \
            ~{match_prop_min_arg}

        bgzip "~{prefix}.agglovar.merged.vcf"
        bcftools index -t "~{prefix}.agglovar.merged.vcf.gz"
    >>>

    output {
        File agglovar_merged_vcf = "~{prefix}.agglovar.merged.vcf.gz"
        File agglovar_merged_vcf_idx = "~{prefix}.agglovar.merged.vcf.gz.tbi"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: 3*ceil(size(vcfs, "GB")) + 10,
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
