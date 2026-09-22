version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateDbSNP {
    meta {
        description: [
            "This workflow annotates each variant in the input VCF with its dbSNP reference SNP identifier (rsID). It matches variants against a per-contig dbSNP VCF on CHROM, POS, REF and ALT, emitting a TSV mapping each matched variant to its `dbSNP_ID`."
        ]
    }

    parameter_meta {
        vcf: "VCF to annotate."
        vcf_idx: "Index for VCF to annotate."
        dbsnp_vcf: "From references."
        dbsnp_vcf_idx: "From references."
        contigs: "Contigs to annotate within the input VCF."
        shard_bin_size: "If set, shards each contig into regions each containing roughly this many combined VCF + dbSNP records, run in parallel."
        annotations_tsv_dbsnp: "TSV mapping variants to their dbSNP identifiers."
    }

    input {
        File vcf
        File vcf_idx
        File dbsnp_vcf
        File dbsnp_vcf_idx
        Array[String] contigs
        String prefix

        Int? shard_bin_size

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_dbsnp_vcf
        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig as SubsetVcf {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_vcf
            }

            call Helpers.SubsetVcfToContig as SubsetDbSNPVcf {
                input:
                    vcf = dbsnp_vcf,
                    vcf_idx = dbsnp_vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.dbsnp",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_dbsnp_vcf
            }
        }

        File contig_vcf = select_first([SubsetVcf.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcf.subset_vcf_idx, vcf_idx])
        File contig_dbsnp_vcf = select_first([SubsetDbSNPVcf.subset_vcf, dbsnp_vcf])
        File contig_dbsnp_vcf_idx = select_first([SubsetDbSNPVcf.subset_vcf_idx, dbsnp_vcf_idx])

        if (defined(shard_bin_size)) {
            call Helpers.CreateContigShards {
                input:
                    vcfs = [contig_vcf, contig_dbsnp_vcf],
                    vcf_idxs = [contig_vcf_idx, contig_dbsnp_vcf_idx],
                    contig = contig,
                    shard_bin_size = select_first([shard_bin_size]),
                    prefix = "~{prefix}.~{contig}.shards",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_create_shards
            }

            scatter (i in range(length(CreateContigShards.shard_regions))) {
                String shard_region = CreateContigShards.shard_regions[i]

                call Helpers.SubsetVcfToRegion as SubsetVcfShard {
                    input:
                        vcf = contig_vcf,
                        vcf_idx = contig_vcf_idx,
                        region = shard_region,
                        prefix = "~{prefix}.~{contig}.shard_~{i}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_vcf
                }

                call Helpers.SubsetVcfToRegion as SubsetDbSNPVcfShard {
                    input:
                        vcf = contig_dbsnp_vcf,
                        vcf_idx = contig_dbsnp_vcf_idx,
                        region = shard_region,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.dbsnp",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_dbsnp_vcf
                }

                call AnnotateDbSNPIds as AnnotateDbSNPIdsShard {
                    input:
                        vcf = SubsetVcfShard.subset_vcf,
                        vcf_idx = SubsetVcfShard.subset_vcf_idx,
                        dbsnp_vcf = SubsetDbSNPVcfShard.subset_vcf,
                        dbsnp_vcf_idx = SubsetDbSNPVcfShard.subset_vcf_idx,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.annotated",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_annotate
                }
            }

            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = AnnotateDbSNPIdsShard.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.dbsnp_annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        if (!defined(shard_bin_size)) {
            call AnnotateDbSNPIds {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    dbsnp_vcf = contig_dbsnp_vcf,
                    dbsnp_vcf_idx = contig_dbsnp_vcf_idx,
                    prefix = "~{prefix}.~{contig}.annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_annotate
            }
        }

        File contig_annotations_tsv = select_first([ConcatShards.concatenated_tsv, AnnotateDbSNPIds.annotations_tsv])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = contig_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.dbsnp_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotations_tsv_dbsnp = select_first([ConcatTsvs.concatenated_tsv, contig_annotations_tsv[0]])
    }
}

task AnnotateDbSNPIds {
    input {
        File vcf
        File vcf_idx
        File dbsnp_vcf
        File dbsnp_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Extract variants from input VCF
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            ~{vcf} \
            > input_variants.txt

        # Extract all variants from dbSNP VCF
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            ~{dbsnp_vcf} \
            > dbsnp_variants.txt

        # Match variants with dbSNP variants on CHROM, POS, REF and any matching ALT
        awk 'BEGIN{OFS="\t"}
        NR==FNR {
            chrom = $1
            pos = $2
            ref = $3
            alt_string = $4
            id = $5
            n = split(alt_string, alts, ",")
            for (i = 1; i <= n; i++) {
                key = chrom "\t" pos "\t" ref "\t" alts[i]
                dbsnp[key] = id
            }
            next
        }
        {
            key = $1 "\t" $2 "\t" $3 "\t" $4
            if (key in dbsnp) {
                print $1, $2, $3, $4, $5, dbsnp[key]
            }
        }' dbsnp_variants.txt input_variants.txt \
            > ~{prefix}.annotations.tsv
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, dbsnp_vcf], "GB")) + 10,
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
