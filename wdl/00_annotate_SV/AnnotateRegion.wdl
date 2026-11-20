version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateRegion {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        File simple_repeats_bed
        File seg_dup_bed
        File repeat_masked_bed

        String annotate_region_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_annotate_region
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
                    runtime_attr_override = runtime_attr_subset_vcf
            }
        }

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.region_annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [contig_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [contig_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call AnnotateGenomicContext {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    simple_repeats_bed = simple_repeats_bed,
                    seg_dup_bed = seg_dup_bed,
                    repeat_masked_bed = repeat_masked_bed,
                    prefix = "~{prefix}.~{contig}.region_annotated.shard_~{i}",
                    docker = annotate_region_docker,
                    runtime_attr_override = runtime_attr_annotate_region
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = AnnotateGenomicContext.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.region_annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, AnnotateGenomicContext.annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.region_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcf
        }
    }

    output {
        File annotations_tsv_region = select_first([ConcatTsvs.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task AnnotateGenomicContext {
    input {
        File vcf
        File vcf_idx
        File simple_repeats_bed
        File seg_dup_bed
        File repeat_masked_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        TMPPATH=$(mktemp -d)

        # Extract variant sites: CHROM, POS0, END, ID, allele_type, allele_length, POS, REF, ALT
        # Normalize type and length to match GATK-SV representation
        bcftools query \
            -f '%CHROM\t%POS0\t%END\t%ID\t%INFO/allele_type\t%INFO/allele_length\t%POS\t%REF\t%ALT\n' \
            ~{vcf} \
        | awk 'BEGIN{OFS="\t"} {
            if (tolower($5) == "trv") {
                $5 = "SNV"
                $6 = length($8)
            } else {
                t = toupper($5)
                if (index(t, "DEL") > 0) $5 = "DEL"
                else if (index(t, "INS") > 0) $5 = "INS"
                else if (index(t, "DUP") > 0 || index(t, "NUMT") > 0) $5 = "DUP"
                else $5 = "SNV"
                if ($6 < 0) $6 = -$6
            }
            print
        }' > ${TMPPATH}/tmp.sites

        # Create left breakpoint BED
        awk '{print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$6}' ${TMPPATH}/tmp.sites \
            | grep -v "^#" | sort -k1,1 -k2,2n > ${TMPPATH}/le_bp

        # Create right breakpoint BED
        awk '{print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6}' ${TMPPATH}/tmp.sites \
            | grep -v "^#" | sort -k1,1 -k2,2n > ${TMPPATH}/ri_bp

        # Extract large DEL/DUP (> 5kb by span)
        awk '{if ($5=="DEL" || $5=="DUP") print}' ${TMPPATH}/tmp.sites \
            | awk '{if ($3-$2 > 5000) print}' \
            | cut -f1-5 \
            | sort -k1,1 -k2,2n > ${TMPPATH}/lg_cnv

        # Breakpoint coverage against each region type
        bedtools coverage -a ${TMPPATH}/le_bp -b ~{simple_repeats_bed} | awk '{if ($NF>0) print}' > ${TMPPATH}/le_bp.vs.SR
        bedtools coverage -a ${TMPPATH}/le_bp -b ~{seg_dup_bed}        | awk '{if ($NF>0) print}' > ${TMPPATH}/le_bp.vs.SD
        bedtools coverage -a ${TMPPATH}/le_bp -b ~{repeat_masked_bed}  | awk '{if ($NF>0) print}' > ${TMPPATH}/le_bp.vs.RM

        bedtools coverage -a ${TMPPATH}/ri_bp -b ~{simple_repeats_bed} | awk '{if ($NF>0) print}' > ${TMPPATH}/ri_bp.vs.SR
        bedtools coverage -a ${TMPPATH}/ri_bp -b ~{seg_dup_bed}        | awk '{if ($NF>0) print}' > ${TMPPATH}/ri_bp.vs.SD
        bedtools coverage -a ${TMPPATH}/ri_bp -b ~{repeat_masked_bed}  | awk '{if ($NF>0) print}' > ${TMPPATH}/ri_bp.vs.RM

        # Large CNV body coverage
        if [ "$(wc -l < ${TMPPATH}/lg_cnv)" -gt 0 ]; then
            bedtools coverage -a ${TMPPATH}/lg_cnv -b ~{simple_repeats_bed} > ${TMPPATH}/lg_cnv.vs.SR
            bedtools coverage -a ${TMPPATH}/lg_cnv -b ~{seg_dup_bed}        > ${TMPPATH}/lg_cnv.vs.SD
            bedtools coverage -a ${TMPPATH}/lg_cnv -b ~{repeat_masked_bed}  > ${TMPPATH}/lg_cnv.vs.RM
        fi

        # Run script to assign genomic context
        Rscript /opt/gnomad-lr/scripts/annotation/annotate_genomic_context.R \
            -i ${TMPPATH}/tmp.sites \
            -o ${TMPPATH}/annotations.unsorted.tsv \
            -p ${TMPPATH}

        sort -k1,1 -k2,2n ${TMPPATH}/annotations.unsorted.tsv > ~{prefix}.annotations.tsv
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(simple_repeats_bed, "GB") + size(seg_dup_bed, "GB") + size(repeat_masked_bed, "GB")) + 10,
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
