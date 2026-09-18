version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateGnomADSTR {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        File gnomad_tr_json
        Float trv_reciprocal_overlap = 0.7

        String utils_docker

        RuntimeAttr? runtime_attr_subset_contigs
        RuntimeAttr? runtime_attr_subset_trv
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_annotate_gnomad_str
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat_vcf
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
                    runtime_attr_override = runtime_attr_subset_contigs
            }
        }

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])

        call Helpers.SubsetVcfByArgs {
            input:
                vcf = contig_vcf,
                vcf_idx = contig_vcf_idx,
                include_args = "INFO/allele_type=\"trv\"",
                prefix = "~{prefix}.~{contig}.trv",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_trv
        }

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = SubsetVcfByArgs.subset_vcf,
                    vcf_idx = SubsetVcfByArgs.subset_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.trv",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [SubsetVcfByArgs.subset_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [SubsetVcfByArgs.subset_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call AnnotateGnomADSTRLoci {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    gnomad_tr_json = gnomad_tr_json,
                    trv_reciprocal_overlap = trv_reciprocal_overlap,
                    prefix = "~{prefix}.~{contig}.gnomad_str.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_annotate_gnomad_str
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = AnnotateGnomADSTRLoci.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.gnomad_str",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, AnnotateGnomADSTRLoci.annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.gnomad_str_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcf
        }
    }

    output {
        File annotations_tsv_gnomad_str = select_first([ConcatTsvs.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task AnnotateGnomADSTRLoci {
    input {
        File vcf
        File vcf_idx
        File gnomad_tr_json
        Float trv_reciprocal_overlap
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Convert JSON catalog to BED
        python3 <<PYCODE
import json

with open("~{gnomad_tr_json}") as f:
    catalog = json.load(f)

with open("catalog.bed", "w") as out:
    for entry in catalog:
        if not entry or not entry.get("LocusId"):
            continue
        locus_id = entry["LocusId"]
        region_str = entry.get("MainReferenceRegion")
        if not region_str:
            continue
        repeat_unit = entry.get("RepeatUnit", "")
        chrom, coords = region_str.rsplit(":", 1)
        start_str, end_str = coords.split("-")
        out.write(f"{chrom}\t{start_str}\t{end_str}\t{locus_id}\t{repeat_unit}\n")
PYCODE

        sort -k1,1 -k2,2n catalog.bed > catalog.sorted.bed

        # Extract VCF variants to BED
        bcftools query \
            -f '%CHROM\t%POS0\t%END\t%ID\t%INFO/MOTIFS\n' \
            ~{vcf} \
        | sort -k1,1 -k2,2n > variants.bed

        # Condition 1: variant fully enveloped within a reference region
        bedtools intersect \
            -a variants.bed \
            -b catalog.sorted.bed \
            -f 1.0 \
            -wo \
        | awk 'BEGIN{OFS="\t"} {print $4, $9}' \
            > cond1_matches.tsv

        # Condition 2: reciprocal overlap meets the threshold and at least one motif matches
        bedtools intersect \
            -a variants.bed \
            -b catalog.sorted.bed \
            -f ~{trv_reciprocal_overlap} \
            -r \
            -wo \
        | awk 'BEGIN{OFS="\t"} {
            n = split($5, m, ",")
            for (i = 1; i <= n; i++) {
                if (toupper(m[i]) == toupper($10)) {
                    print $4, $9
                    break
                }
            }
        }' \
        > cond2_matches.tsv

        # Merge the two condition sets, giving condition 1 priority per variant ID
        cat cond1_matches.tsv cond2_matches.tsv \
            | sort -k1,1 -u \
            > all_matches.tsv

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ~{vcf} \
            | awk 'BEGIN{OFS="\t"} NR==FNR{locus[$1]=$2; next} ($5 in locus){print $1,$2,$3,$4,$5,locus[$5]}' \
            all_matches.tsv - \
            | sort -k1,1 -k2,2n \
            > ~{prefix}.annotations.tsv
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
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
