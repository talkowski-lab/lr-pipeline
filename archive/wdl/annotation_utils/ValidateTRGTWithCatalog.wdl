version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ValidateTRGTWithCatalog {
    meta {
        description: [
            "This utility finds TRGT calls whose coordinates or motifs disagree with the catalog they were genotyped against. Each contig is sharded, checked against the catalog, and the incongruent records are concatenated into one VCF."
        ]
    }

    parameter_meta {
        trgt_vcf: "TRGT callset to validate."
        trgt_vcf_idx: "Index for `trgt_vcf`."
        contigs: "Contigs to process."
        trgt_catalog_bed_gz: "TRGT catalog BED (gzipped) the calls are validated against."
        records_per_shard: "Number of variants to keep within a single shard during validation."
        incongruent_vcf: "VCF holding the calls that disagree with the catalog."
        incongruent_vcf_idx: "Index for `incongruent_vcf`."
    }

    input {
        File trgt_vcf
        File trgt_vcf_idx
        Array[String] contigs
        String prefix

        File trgt_catalog_bed_gz

        Int? records_per_shard

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_catalog
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_validate
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat_contigs
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = trgt_vcf,
                vcf_idx = trgt_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Helpers.SubsetTsvToContig {
            input:
                tsv = trgt_catalog_bed_gz,
                contig = contig,
                compressed_tsv = true,
                prefix = "~{prefix}.~{contig}.catalog",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_catalog
        }

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = SubsetVcfToContig.subset_vcf,
                    vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] shard_vcfs = select_first([ShardVcfByRecords.shards, [SubsetVcfToContig.subset_vcf]])
        Array[File] shard_vcf_idxs = select_first([ShardVcfByRecords.shard_idxs, [SubsetVcfToContig.subset_vcf_idx]])

        scatter (shard in range(length(shard_vcfs))) {
            call FindIncongruentTRGTVariants {
                input:
                    vcf = shard_vcfs[shard],
                    vcf_idx = shard_vcf_idxs[shard],
                    catalog_bed = SubsetTsvToContig.subset_tsv,
                    prefix = "~{prefix}.~{contig}.incongruent.shard_~{shard}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_validate
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatVcfs as ConcatShards {
                input:
                    vcfs = FindIncongruentTRGTVariants.incongruent_vcf,
                    vcf_idxs = FindIncongruentTRGTVariants.incongruent_vcf_idx,
                    allow_overlaps = false,
                    naive = true,
                    prefix = "~{prefix}.~{contig}.incongruent",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File contig_incongruent_vcf = select_first([ConcatShards.concat_vcf, FindIncongruentTRGTVariants.incongruent_vcf[0]])
        File contig_incongruent_vcf_idx = select_first([ConcatShards.concat_vcf_idx, FindIncongruentTRGTVariants.incongruent_vcf_idx[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatVcfs as ConcatContigs {
            input:
                vcfs = contig_incongruent_vcf,
                vcf_idxs = contig_incongruent_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.incongruent",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_contigs
        }
    }

    output {
        File incongruent_vcf = select_first([ConcatContigs.concat_vcf, contig_incongruent_vcf[0]])
        File incongruent_vcf_idx = select_first([ConcatContigs.concat_vcf_idx, contig_incongruent_vcf_idx[0]])
    }
}

task FindIncongruentTRGTVariants {
    input {
        File vcf
        File vcf_idx
        File catalog_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
import gzip

catalog_positions = set()
with open("~{catalog_bed}") as catalog_file:
    for line in catalog_file:
        if not line.strip() or line.startswith("#"):
            continue
        chrom, pos, end = line.rstrip("\n").split("\t")[:3]
        catalog_positions.add((chrom, int(pos), int(end)))

with gzip.open("~{vcf}", "rt") as vcf_in, open("incongruent.vcf", "w") as vcf_out:
    for line in vcf_in:
        if line.startswith("#"):
            vcf_out.write(line)
            continue
        chrom, pos, _, ref, *_ = line.rstrip("\n").split("\t")
        end = int(pos) + len(ref) - 1
        if (chrom, int(pos), end) not in catalog_positions:
            vcf_out.write(line)
PYCODE

        bgzip -c incongruent.vcf > ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File incongruent_vcf = "~{prefix}.vcf.gz"
        File incongruent_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, catalog_bed], "GB")) + 5,
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
