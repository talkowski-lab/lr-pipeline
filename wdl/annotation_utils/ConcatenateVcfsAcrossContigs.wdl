version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ConcatenateVcfsAcrossContigs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        String prefix

        Boolean drop_genotypes = false

        String utils_docker

        RuntimeAttr? runtime_attr_validate_contigs
        RuntimeAttr? runtime_attr_strip_genotypes
        RuntimeAttr? runtime_attr_concat
    }

    call ValidateContigOrder {
        input:
            n_vcfs = length(vcfs),
            n_vcf_idxs = length(vcf_idxs),
            contigs = contigs,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_validate_contigs
    }

    if (drop_genotypes) {
        scatter (i in range(length(vcfs))) {
            call Helpers.StripGenotypes {
                input:
                    vcf = vcfs[i],
                    vcf_idx = vcf_idxs[i],
                    prefix = "~{prefix}.~{contigs[i]}.no_geno",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_strip_genotypes
            }
        }
    }

    Array[File] vcfs_to_concat = select_first([StripGenotypes.stripped_vcf, vcfs])
    Array[File] vcf_idxs_to_concat = select_first([StripGenotypes.stripped_vcf_idx, vcf_idxs])

    call Helpers.ConcatVcfs {
        input:
            vcfs = vcfs_to_concat,
            vcf_idxs = vcf_idxs_to_concat,
            allow_overlaps = false,
            naive = true,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File concat_vcf = ConcatVcfs.concat_vcf
        File concat_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task ValidateContigOrder {
    input {
        Int n_vcfs
        Int n_vcf_idxs
        Array[String] contigs
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        n_contigs=$(wc -l < ~{write_lines(contigs)})

        if [[ ~{n_vcfs} -ne ~{n_vcf_idxs} || ~{n_vcfs} -ne "$n_contigs" ]]; then
            echo "ERROR: vcfs (~{n_vcfs}), vcf_idxs (~{n_vcf_idxs}), and contigs ($n_contigs) must have the same length." >&2
            exit 1
        fi

        duplicate_contigs=$(sort ~{write_lines(contigs)} | uniq -d)
        if [[ -n "$duplicate_contigs" ]]; then
            echo "ERROR: ConcatenateVcfsAcrossContigs expects exactly one VCF per contig. Duplicate contigs:" >&2
            echo "$duplicate_contigs" >&2
            exit 1
        fi
    >>>

    output {
        String status = "success"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: 10,
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
