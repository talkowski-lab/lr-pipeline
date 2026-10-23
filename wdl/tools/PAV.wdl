version 1.0

import "../utils/Structs.wdl"

workflow PAV {
    meta {
        description: [
            "This tool runs PAV (https://github.com/EichlerLab/pav) on a sample's phased haplotype assemblies to call variants against the reference. It outputs the sample's tabix-indexed VCF, along with tarballs of the full PAV results and log directories."
        ]
    }

    parameter_meta {
        mat_haplotype: "Maternal haplotype assembly. Without parental data, either hifiasm haplotype."
        pat_haplotype: "Paternal haplotype assembly. Without parental data, the other hifiasm haplotype."
        sample_id: "Sample ID."
        ref_fa: "From references."
        ref_fai: "From references."
        pav_results_tarball: "Tarball of the full PAV results directory."
        pav_log_tarball: "Tarball of the full PAV log directory."
        pav_vcf: "Called VCF."
        pav_vcf_idx: "Tabix index for `pav_vcf`."
    }

    input {
        File mat_haplotype
        File pat_haplotype
        String sample_id
        String prefix

        File ref_fa
        File ref_fai

        String pav_docker

        RuntimeAttr? runtime_attr_run_pav
    }

    call CallPAV {
        input:
            mat_haplotype = mat_haplotype,
            pat_haplotype = pat_haplotype,
            sample_id = sample_id,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            prefix = prefix,
            docker = pav_docker,
            runtime_attr_override = runtime_attr_run_pav
    }

    output {
        File pav_results_tarball = CallPAV.results_tar
        File pav_log_tarball = CallPAV.log_tar
        File pav_vcf = CallPAV.vcf
        File pav_vcf_idx = CallPAV.vcf_idx
    }
}

task CallPAV {
    input {
        File mat_haplotype
        File pat_haplotype
        String sample_id
        String prefix

        File ref_fa
        File ref_fai

        String docker

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 16,
        mem_gb: 10,
        disk_gb: 3 * ceil(size(mat_haplotype, "GB") + size(pat_haplotype, "GB") + size(ref_fa, "GB")) + 10,
        boot_disk_gb: 20,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int effective_cpu = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    command <<<
        set -euo pipefail

        if [[ "~{ref_fa}" == *.gz ]]; then
            ln -s ~{ref_fa} ref.fa.gz
            ln -s ~{ref_fai} ref.fa.gz.fai
        else
            echo "Compressing reference with bgzip..."
            bgzip -c ~{ref_fa} > ref.fa.gz
            samtools faidx ref.fa.gz
        fi

        python3 <<'CODE'
import os
import json

ref_fa_abs = os.path.abspath("ref.fa.gz")
config = {"reference": ref_fa_abs}
with open("pav.json", "w") as f:
    json.dump(config, f)

os.makedirs("asms", exist_ok=True)
mat_link = os.path.abspath("asms/~{sample_id}_mat.fa.gz")
pat_link = os.path.abspath("asms/~{sample_id}_pat.fa.gz")
os.symlink("~{mat_haplotype}", mat_link)
os.symlink("~{pat_haplotype}", pat_link)

with open("assemblies.tsv", "w") as f:
    f.write("NAME\tHAP_mat\tHAP_pat\n")
    f.write(f"~{sample_id}\t{mat_link}\t{pat_link}\n")

CODE

        python3 -m pav3 batch --cores ~{effective_cpu} || [ $? -eq 1 ]
        bcftools index --tbi ~{sample_id}.vcf.gz

        tar -zcf ~{prefix}.pav_results.tar.gz results
        tar -zcf ~{prefix}.pav_log.tar.gz log
    >>>

    output {
        File results_tar = "~{prefix}.pav_results.tar.gz"
        File log_tar = "~{prefix}.pav_log.tar.gz"
        File vcf = "~{sample_id}.vcf.gz"
        File vcf_idx = "~{sample_id}.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        continueOnReturnCode: [0, 1]
    }
}
