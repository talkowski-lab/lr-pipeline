version 1.0

import "Structs.wdl"

# Table of gnomAD v4.1 genome PASS sites with allele frequencies.
#
# Streams the public gnomAD v4.1 genome sites VCFs over HTTPS (tabix-indexed,
# so each shard reads only its own region; nothing is fully downloaded),
# split into shard_size-bp shards per contig. Each shard writes PASS records
# as:
#   chrom  pos  end  allele_type  length  AF  AC  AN
# where end = pos + len(REF) - 1, allele_type is gnomAD's INFO/allele_type
# (snv / ins / del / mixed), and length is 1 for SNVs, len(REF) - len(ALT)
# for deletions, len(ALT) - len(REF) for insertions, and
# max(len(REF), len(ALT)) otherwise. Records are biallelic in the release
# VCFs, so AF/AC are single values. A record is assigned to the shard
# containing its POS (--regions-overlap pos), so records spanning a shard
# boundary are not duplicated. Shards are concatenated into one bgzipped
# table per contig, plus a per-contig count summary.

workflow GnomadV4GenomePassSitesAF {
    input {
        Array[String] contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
        String prefix = "gnomad.genomes.v4.1.sites.PASS"

        String vcf_url_template = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.CONTIG.vcf.bgz"
        Int shard_size = 10000000

        String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"

        RuntimeAttr? runtime_attr_make_shards
        RuntimeAttr? runtime_attr_extract_shard
        RuntimeAttr? runtime_attr_concat_contigs
    }

    call MakeShards {
        input:
            contigs = contigs,
            vcf_url_template = vcf_url_template,
            shard_size = shard_size,
            prefix = prefix,
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_make_shards
    }

    scatter (shard in read_tsv(MakeShards.shards)) {
        call ExtractShard {
            input:
                vcf_url = shard[1],
                region = shard[2],
                prefix = "~{prefix}.~{shard[0]}",
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_extract_shard
        }
    }

    call ConcatContigs {
        input:
            shard_tables = ExtractShard.table,
            contigs = contigs,
            prefix = prefix,
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat_contigs
    }

    output {
        Array[File] contig_tables = ConcatContigs.tables
        Array[File] contig_table_idxs = ConcatContigs.table_idxs
        File site_counts = ConcatContigs.counts
    }
}

task MakeShards {
    input {
        Array[String] contigs
        String vcf_url_template
        Int shard_size
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        for contig in ~{sep=" " contigs}; do
            url=$(echo "~{vcf_url_template}" | sed "s/CONTIG/${contig}/")
            len=$(bcftools view -h "${url}" | awk -v c="${contig}" 'match($0, /^##contig=<ID=[^,]+,length=[0-9]+/) { split(substr($0, 14, RLENGTH - 13), a, ",length="); if (a[1] == c) print a[2] }')
            if [ -z "${len}" ]; then
                echo "No ##contig length for ${contig} in ${url}" >&2
                exit 1
            fi
            awk -v c="${contig}" -v u="${url}" -v L="${len}" -v s=~{shard_size} 'BEGIN {
                n = 0
                for (start = 1; start <= L; start += s) {
                    end = start + s - 1; if (end > L) end = L
                    printf "%s.%05d\t%s\t%s:%d-%d\n", c, n, u, c, start, end
                    n++
                }
            }'
        done > ~{prefix}.shards.tsv
    >>>

    output {
        File shards = "~{prefix}.shards.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10,
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

task ExtractShard {
    input {
        String vcf_url
        String region
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query \
            --regions-overlap pos \
            -r ~{region} \
            -i 'FILTER="PASS"' \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/allele_type\t%INFO/AF\t%INFO/AC\t%INFO/AN\n' \
            "~{vcf_url}" \
        | awk -F'\t' -v OFS='\t' '{
            r = length($3); a = length($4)
            if ($5 == "snv") len = 1
            else if ($5 == "del") len = r - a
            else if ($5 == "ins") len = a - r
            else len = (r > a ? r : a)
            print $1, $2, $2 + r - 1, $5, len, $6, $7, $8
        }' \
        | gzip -c > ~{prefix}.tsv.gz
    >>>

    output {
        File table = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10,
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

task ConcatContigs {
    input {
        Array[File] shard_tables
        Array[String] contigs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        ls ~{sep=" " shard_tables} > shard_list.txt
        printf "contig\tn_sites\tn_snv\tn_ins\tn_del\tn_other\n" > ~{prefix}.site_counts.tsv
        for contig in ~{sep=" " contigs}; do
            out=~{prefix}.${contig}.tsv.gz
            # shard files are named <prefix>.<contig>.<5-digit index>.tsv.gz, so a lexical sort is genomic order
            awk -F/ -v p="~{prefix}.${contig}." '{
                b = $NF
                if (substr(b, 1, length(p)) == p && substr(b, length(p) + 1) ~ /^[0-9][0-9][0-9][0-9][0-9]\.tsv\.gz$/) print b "\t" $0
            }' shard_list.txt | sort -k1,1 | cut -f2 > files.${contig}.txt
            if [ ! -s files.${contig}.txt ]; then
                echo "No shard tables for ${contig}" >&2
                exit 1
            fi
            { printf "#chrom\tpos\tend\tallele_type\tlength\tAF\tAC\tAN\n"; xargs cat < files.${contig}.txt | gzip -dc; } \
                | bgzip -c > "${out}"
            tabix -s 1 -b 2 -e 3 "${out}"
            bgzip -dc "${out}" | awk -F'\t' -v c="${contig}" 'NR > 1 { n++; t[$4]++ } END {
                print c "\t" n + 0 "\t" t["snv"] + 0 "\t" t["ins"] + 0 "\t" t["del"] + 0 "\t" n - t["snv"] - t["ins"] - t["del"]
            }' >> ~{prefix}.site_counts.tsv
        done
    >>>

    output {
        Array[File] tables = glob("~{prefix}.chr*.tsv.gz")
        Array[File] table_idxs = glob("~{prefix}.chr*.tsv.gz.tbi")
        File counts = "~{prefix}.site_counts.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(shard_tables, "GB")) + 20,
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
