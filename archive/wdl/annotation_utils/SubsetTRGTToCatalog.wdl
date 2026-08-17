version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SubsetTRGTToCatalog {
    input {
        File trgt_full_merged_vcf
        File trgt_full_merged_vcf_idx
        File trgt_catalog_bed_gz
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_extract_ids
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_filter_to_catalog
        RuntimeAttr? runtime_attr_concat_vcf
    }

    scatter (contig in contigs) {
        call ExtractCatalogForContig {
            input:
                catalog_bed_gz = trgt_catalog_bed_gz,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_extract_ids
        }

        call Helpers.SubsetVcfToContig as SubsetVcf {
            input:
                vcf = trgt_full_merged_vcf,
                vcf_idx = trgt_full_merged_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call FilterToCatalog {
            input:
                vcf = SubsetVcf.subset_vcf,
                vcf_idx = SubsetVcf.subset_vcf_idx,
                catalog_ids = ExtractCatalogForContig.catalog_ids,
                prefix = "~{prefix}.~{contig}.filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter_to_catalog
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = FilterToCatalog.filtered_vcf,
            vcf_idxs = FilterToCatalog.filtered_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.catalog_subset",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcf
    }

    output {
        File trgt_merged_vcf = ConcatVcfs.concat_vcf
        File trgt_merged_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task ExtractCatalogForContig {
    input {
        File catalog_bed_gz
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        zcat ~{catalog_bed_gz} | \
            awk -F'\t' -v contig="~{contig}" '
                $1 == contig {
                    start_pos = $2;
                    n = split($4, fields, ";");
                    found_id = 0;
                    for (i = 1; i <= n; i++) {
                        if (fields[i] ~ /^ID=/) {
                            sub(/^ID=/, "", fields[i]);
                            print fields[i] "\t" start_pos;
                            found_id = 1;
                            break;
                        }
                    }
                    if (!found_id) {
                        print "WARNING: No ID found for " contig ":" start_pos > "/dev/stderr";
                    }
                }
            ' > ~{prefix}.txt

        echo "Extracted $(wc -l < ~{prefix}.txt) catalog entries for ~{contig}" >&2
    >>>

    output {
        File catalog_ids = "~{prefix}.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(catalog_bed_gz, "GB")) + 10,
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

task FilterToCatalog {
    input {
        File vcf
        File vcf_idx
        File catalog_ids
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

# Parse catalog
catalog = set()
with open('~{catalog_ids}', 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            catalog.add((str(parts[0]), str(parts[1])))

# Process VCF
seen = set()
vcf_in = pysam.VariantFile('~{vcf}')
vcf_out = pysam.VariantFile('~{prefix}.vcf.gz', 'wz', header=vcf_in.header)
for record in vcf_in:
    pos = record.pos
    trid = record.info.get('TRID')
    if isinstance(trid, tuple):
        trid = ','.join(trid)
    key = (str(trid), str(pos))

    if key in catalog and key not in seen:
        seen.add(key)
        vcf_out.write(record)

vcf_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
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
