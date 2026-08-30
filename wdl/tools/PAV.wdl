version 1.0

import "../utils/Structs.wdl"

workflow PAV {
    input {
        Array[File] mat_haplotypes
        Array[File] pat_haplotypes
        Array[String] sample_ids
        String prefix

        File ref_fa
        File ref_fai

        String pav_docker

        RuntimeAttr? runtime_attr_run_pav
    }

    call CallPAV {
        input:
            mat_haplotypes = mat_haplotypes,
            pat_haplotypes = pat_haplotypes,
            sample_ids = sample_ids,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            prefix = prefix,
            docker = pav_docker,
            runtime_attr_override = runtime_attr_run_pav
    }

    output {
        File pav_results_tarball = CallPAV.results_tar
        File pav_log_tarball = CallPAV.log_tar
        Array[File] pav_vcfs = CallPAV.vcfs
        Array[File] pav_vcf_idx = CallPAV.vcf_idx

        File? debug_sam = CallPAV.debug_sam
        Array[File]? debug_temp = CallPAV.debug_temp
    }
}

task CallPAV {
    input {
        Array[File] mat_haplotypes
        Array[File] pat_haplotypes
        Array[String] sample_ids
        String prefix

        File ref_fa
        File ref_fai

        String docker

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 16,
        mem_gb: 10,
        disk_gb: 3 * ceil(size(mat_haplotypes, "GB") + size(pat_haplotypes, "GB") + size(ref_fa, "GB")) + 10,
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
sample_ids = "~{sep=' ' sample_ids}".split(' ')
mat_files = "~{sep=' ' mat_haplotypes}".split(' ')
pat_files = "~{sep=' ' pat_haplotypes}".split(' ')

if len(sample_ids) != len(mat_files) or len(sample_ids) != len(pat_files):
    raise ValueError(f"Input array lengths must match: {len(sample_ids)} samples, {len(mat_files)} maternal, {len(pat_files)} paternal")

with open("assemblies.tsv", "w") as f:
    f.write("NAME\tHAP_mat\tHAP_pat\n")
    for i, sample_id in enumerate(sample_ids):
        mat_link = f"asms/{sample_id}_mat.fa.gz"
        pat_link = f"asms/{sample_id}_pat.fa.gz"
        os.symlink(os.path.abspath(mat_files[i]), mat_link)
        os.symlink(os.path.abspath(pat_files[i]), pat_link)
        
        mat_link_abs = os.path.abspath(mat_link)
        pat_link_abs = os.path.abspath(pat_link)
        f.write(f"{sample_id}\t{mat_link_abs}\t{pat_link_abs}\n")

CODE

        python3 -m pav3 batch --cores ~{effective_cpu}

        tar -zcf ~{prefix}.pav_results.tar.gz results
        tar -zcf ~{prefix}.pav_log.tar.gz log
    >>>

    output {
        File results_tar = "~{prefix}.pav_results.tar.gz"
        File log_tar = "~{prefix}.pav_log.tar.gz"
        Array[File] vcfs = glob("*.vcf.gz")
        Array[File] vcf_idx = glob("*.vcf.gz.csi")

        File? debug_sam = "temp/HG00733/align/trim-none/align_qry_mat.sam.gz"
        Array[File]? debug_temp = glob("temp/**/*")
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
