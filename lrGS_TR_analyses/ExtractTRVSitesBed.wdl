version 1.0

import "../wdl/utils/Structs.wdl"

# Extract tandem-repeat variant (TRV) records (variant ID containing "-TRV-") from per-contig VCFs into one BED.
#
# `vcfs` / `vcf_idxs` / `contigs` are parallel per-contig arrays. Per contig (scattered), ExtractTRVBed streams the VCF
# with bcftools and writes one row per TRV record:
#   #chrom start end ID MOTIFS n_alt AF_total AC_total AN ref_size alt_sizes AF_list AC_list FILTER
#   start / end   repeat span in BED coordinates: start = POS (0-based start after the VCF padding base),
#                 end = POS + len(REF) - 1 (equals the TRID span)
#   n_alt         number of ALT alleles
#   AF_total      sum(AC) / AN; AC_total = sum(AC)
#   ref_size      len(REF) - 1 (padding base removed); alt_sizes = len(ALT) - 1 per ALT allele, in ALT order
#   AF_list / AC_list   per-ALT INFO/AF and INFO/AC as in the VCF
# ConcatTRVBeds then merges the contigs into one sorted, bgzipped and tabix-indexed BED.

workflow ExtractTRVSitesBed {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        String prefix

        String bcftools_docker

        RuntimeAttr? runtime_attr_extract_trv_bed
        RuntimeAttr? runtime_attr_concat_trv_beds
    }

    scatter (i in range(length(contigs))) {
        call ExtractTRVBed {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                prefix = "~{prefix}.~{contigs[i]}",
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_extract_trv_bed
        }
    }

    call ConcatTRVBeds {
        input:
            beds = ExtractTRVBed.trv_bed,
            prefix = "~{prefix}.TRV_sites",
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat_trv_beds
    }

    output {
        File trv_bed = ConcatTRVBeds.trv_bed
        File trv_bed_idx = ConcatTRVBeds.trv_bed_idx
        Array[File] per_contig_trv_beds = ExtractTRVBed.trv_bed
    }
}

task ExtractTRVBed {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -i 'ID~"-TRV-"' -Ou ~{vcf} \
        | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/MOTIFS\t%INFO/AC\t%INFO/AF\t%INFO/AN\t%FILTER\n' \
        | awk -F'\t' 'BEGIN {
                OFS = "\t"
                print "#chrom", "start", "end", "ID", "MOTIFS", "n_alt", "AF_total", "AC_total", "AN", "ref_size", "alt_sizes", "AF_list", "AC_list", "FILTER"
            }
            {
                n = split($5, alt, ",")
                sizes = ""
                for (i = 1; i <= n; i++) sizes = sizes (i > 1 ? "," : "") (length(alt[i]) - 1)
                k = split($7, ac, ",")
                ac_total = 0
                for (i = 1; i <= k; i++) if (ac[i] != ".") ac_total += ac[i]
                af_total = ($9 != "." && $9 > 0) ? sprintf("%.6g", ac_total / $9) : "."
                print $1, $2, $2 + length($4) - 1, $3, $6, n, af_total, ac_total, $9, length($4) - 1, sizes, $8, $7, $10
            }' \
        | gzip -c > ~{prefix}.TRV_sites.bed.gz
    >>>

    output {
        File trv_bed = "~{prefix}.TRV_sites.bed.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 10,
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

task ConcatTRVBeds {
    input {
        Array[File] beds
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        first=1
        for f in ~{sep=" " beds}; do
            if [ "$first" -eq 1 ]; then
                gzip -dc "$f" | sed -n '1p' > header.txt
                first=0
            fi
            gzip -dc "$f" | tail -n +2
        done | sort -k1,1 -k2,2n > body.txt
        cat header.txt body.txt | bgzip -c > ~{prefix}.bed.gz
        tabix -p bed -S 1 ~{prefix}.bed.gz
    >>>

    output {
        File trv_bed = "~{prefix}.bed.gz"
        File trv_bed_idx = "~{prefix}.bed.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size(beds, "GB")) + 10,
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
