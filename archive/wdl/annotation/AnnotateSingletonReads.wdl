version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow AnnotateSingletonReads {
    meta {
        description: [
            "This utility flags variants that look like single-read artifacts. Working one contig at a time it recomputes `AC`, then adds a `SINGLE_READ_SUPPORT` FILTER to any variant whose allele count is at or below two and whose alternate allele is supported by exactly one read in exactly one sample, and concatenates the per-contig results."
        ]
    }

    parameter_meta {
        vcf: "Cohort VCF to flag."
        vcf_idx: "Index for `vcf`."
        contigs: "Contigs to process."
        singleton_filtered_vcf: "VCF whose single-read-supported variants carry the `SINGLE_READ_SUPPORT` FILTER."
        singleton_filtered_vcf_idx: "Index for `singleton_filtered_vcf`."
    }

    input {
        File vcf
        File vcf_idx

        String prefix
        Array[String] contigs

        String utils_docker

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_find_singletons
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset
            }
        }

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])

        call FindSingletonReads {
            input:
                vcf = contig_vcf,
                vcf_idx = contig_vcf_idx,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_find_singletons
        }
    }

    if (!single_contig) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = FindSingletonReads.filtered_vcf,
                vcf_idxs = FindSingletonReads.filtered_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.singleton_reads_filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File singleton_filtered_vcf = select_first([ConcatVcfs.concat_vcf, FindSingletonReads.filtered_vcf[0]])
        File singleton_filtered_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, FindSingletonReads.filtered_vcf_idx[0]])
    }
}

task FindSingletonReads {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools +fill-tags ~{vcf} -Oz -o tagged.vcf.gz -- -t AC
        tabix -p vcf tagged.vcf.gz

        python3 <<'EOF'
from pysam import VariantFile

vcf_in = VariantFile("tagged.vcf.gz")

header = vcf_in.header
header.filters.add("SINGLE_READ_SUPPORT", None, None, "Variant supported by only one read in a single sample")

vcf_out = VariantFile("~{prefix}.filtered.vcf", 'w', header=header)

samples = list(vcf_in.header.samples)

for rec in vcf_in:
    is_suspicious = False

    if rec.info['AC'][0] <= 2:
        alt_AD_in_called_samples = []
        for sample in samples:
            if 'AD' in rec.samples[sample]:
                if len(rec.samples[sample]['AD']) == 2:
                    alt_AD = rec.samples[sample]['AD'][1]
                    if alt_AD > 0:
                        alt_AD_in_called_samples.append(alt_AD)
        if len(alt_AD_in_called_samples) == 1 and alt_AD_in_called_samples[0] == 1:
            is_suspicious = True

    if is_suspicious:
        rec.filter.add("SINGLE_READ_SUPPORT")

    vcf_out.write(rec)

vcf_in.close()
vcf_out.close()
EOF

        bgzip "~{prefix}.filtered.vcf"
        tabix "~{prefix}.filtered.vcf.gz"
    >>>

    output {
        File filtered_vcf = "~{prefix}.filtered.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.filtered.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: ceil(size(vcf, "GB")) + 5,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 20,
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
