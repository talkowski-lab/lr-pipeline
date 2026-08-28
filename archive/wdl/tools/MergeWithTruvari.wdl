version 1.0

import "utils/Structs.wdl"

workflow MergeWithTruvari {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs

        String prefix
        String? truvari_params
        String? bcftools_merge_params
        File? preprocess_script

        File ref_fa
        File ref_fai

        String merge_docker
        String truvari_docker
        
        RuntimeAttr? runtime_attr_override
    }
    
    call BcftoolsMerge {
        input:
            vcfs = vcfs,
            vcf_idxs = vcf_idxs,
            prefix = prefix,
            params = bcftools_merge_params,
            preprocess_script = preprocess_script,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            docker = merge_docker
    }

    call Truvari {
        input:
            prefix = prefix,
            vcf = BcftoolsMerge.merged_vcf,
            params = truvari_params,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            docker = truvari_docker
    }

    output {
        File truvari_collapsed_vcf = Truvari.truvari_collapsed_vcf
        File truvari_collapsed_vcf_idx = Truvari.truvari_collapsed_vcf_idx
    }
}

task BcftoolsMerge {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String params=""
        File? preprocess_script
        File ref_fa
        File ref_fai
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
      
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        for vcf in ~{sep=" " vcfs}; do
            if [[ $vcf == *vcf.gz ]]; then suffix=".vcf.gz"; else suffix=".vcf"; fi
            bcftools norm --multiallelics - --output-type z $vcf > $(basename $vcf $suffix).norm.vcf.gz
            tabix $(basename $vcf $suffix).norm.vcf.gz
        done

        if [ -z ~{preprocess_script} ]; then
            bcftools merge ~{params} --threads ${N_THREADS} -m none *norm.vcf.gz | bgzip > tmp.vcf.gz
        else
            for vcf in *norm.vcf.gz; do
                python ~{preprocess_script} $vcf ~{ref_fa} | bgzip > $(basename $vcf ".vcf.gz").preprocessed.vcf.gz
                tabix $(basename $vcf ".vcf.gz").preprocessed.vcf.gz
            done
            bcftools merge ~{params} --threads ${N_THREADS} -m none *preprocessed.vcf.gz | bgzip > tmp.vcf.gz
        fi
        
        bcftools norm --multiallelics - --output-type z tmp.vcf.gz > ~{prefix}.bcftools_merge.vcf.gz
        tabix ~{prefix}.bcftools_merge.vcf.gz
        rm -f tmp.vcf.gz*
    >>>

    output {
        File merged_vcf = "~{prefix}.bcftools_merge.vcf.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 8*ceil(size(vcfs, 'GB')),
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

task Truvari {
    input {
        String prefix
        String params="-k common --refdist 100 --pctseq 0 --pctsize 0"
        File vcf
        File ref_fa
        File ref_fai
        String docker
        RuntimeAttr? runtime_attr_override
    }
    
    command <<<
        set -euo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        MEM=$(grep '^MemTotal' /proc/meminfo | awk '{ print int($2/1000000) }')

        tabix ~{vcf}

        truvari collapse \
            -i ~{vcf} \
            -c removed.vcf.gz \
            ~{params} \
            -f ~{ref_fa} \
            | bcftools sort --max-mem ${MEM}G -Oz -o ~{prefix}.truvari_collapsed.vcf.gz
        tabix ~{prefix}.truvari_collapsed.vcf.gz
    >>>
    
    output {
        File truvari_collapsed_vcf = "~{prefix}.truvari_collapsed.vcf.gz"
        File truvari_collapsed_vcf_idx = "~{prefix}.truvari_collapsed.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf,"GB") + size(ref_fa,"GB")) + 10,
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
