version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateMEDs {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        Int del_breakpoint_window = 500
        Float del_reciprocal_overlap = 0.0
        Float del_sequence_similarity = 0.7
        Float del_size_similarity = 0.7

        File mei_catalog

        String utils_docker

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_bedtools
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat_shards
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

        call Helpers.SubsetVcfByArgs as SubsetDeletions {
            input:
                vcf = contig_vcf,
                vcf_idx = contig_vcf_idx,
                include_args = 'INFO/allele_type="del"',
                prefix = "~{prefix}.~{contig}.del_subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = SubsetDeletions.subset_vcf,
                    vcf_idx = SubsetDeletions.subset_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.med_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [SubsetDeletions.subset_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [SubsetDeletions.subset_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call ExtractDeletionsToBed {
                input:
                    vcf = vcfs_to_process[i],
                    prefix = "~{prefix}.~{contig}.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_bedtools
            }

            call IntersectMED {
                input:
                    bed_a = ExtractDeletionsToBed.del_bed,
                    bed_b = mei_catalog,
                    del_reciprocal_overlap = del_reciprocal_overlap,
                    prefix = "~{prefix}.~{contig}.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_bedtools
            }

            call GenerateMedAnnotationTable {
                input:
                    intersect_bed = IntersectMED.intersect_bed,
                    del_breakpoint_window = del_breakpoint_window,
                    del_reciprocal_overlap = del_reciprocal_overlap,
                    del_sequence_similarity = del_sequence_similarity,
                    del_size_similarity = del_size_similarity,
                    prefix = "~{prefix}.~{contig}.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_annotate
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = GenerateMedAnnotationTable.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.med_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, GenerateMedAnnotationTable.annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs as MergeAnnotations {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.med_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotations_tsv_meds = select_first([MergeAnnotations.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task ExtractDeletionsToBed {
    input {
        File vcf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -i 'INFO/allele_type=="del"' ~{vcf} \
            | bcftools query -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%ID\n' \
            > ~{prefix}.del.bed
    >>>

    output {
        File del_bed = "~{prefix}.del.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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

task IntersectMED {
    input {
        File bed_a
        File bed_b
        Float del_reciprocal_overlap
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bedtools intersect \
            -a ~{bed_a} \
            -b ~{bed_b} \
            -wa \
            -wb \
            -r \
            -f ~{del_reciprocal_overlap} \
            > ~{prefix}.intersect.bed
    >>>

    output {
        File intersect_bed = "~{prefix}.intersect.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bed_a, "GB") + size(bed_b, "GB")) + 5,
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

task GenerateMedAnnotationTable {
    input {
        File intersect_bed
        Int del_breakpoint_window
        Float del_reciprocal_overlap
        Float del_sequence_similarity
        Float del_size_similarity
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import edlib
import pandas as pd

def get_me_type(designation):
    d = str(designation)
    if d.startswith('SINE') or d.startswith('Alu'):
        return 'ALU'
    elif d.startswith('LINE') or d.startswith('L1'):
        return 'LINE'
    elif d.startswith('Retroposon') or d.startswith('SVA'):
        return 'SVA'
    return None

def compute_sequence_similarity(seq_a, seq_b):
    seq_a = seq_a.upper()
    seq_b = seq_b.upper()
    if not seq_a and not seq_b:
        return 1.0
    if not seq_a or not seq_b:
        return 0.0
    if seq_a == seq_b:
        return 1.0
    result = edlib.align(seq_a, seq_b, task="distance")
    max_len = max(len(seq_a), len(seq_b))
    if max_len == 0:
        return 1.0
    return 1.0 - (result["editDistance"] / max_len)

def passes_criteria(del_start, del_end, del_seq, med_start, med_end, med_seq, size_similarity, reciprocal_overlap, sequence_similarity, breakpoint_window):
    del_len = abs(del_end - del_start)
    med_len = abs(med_end - med_start)

    if del_len == 0 or med_len == 0:
        return False

    size_sim = min(del_len, med_len) / max(del_len, med_len)
    overlap = max(0, min(del_end, med_end) - max(del_start, med_start))
    rec_ovl = overlap / min(del_len, med_len)

    if size_sim < size_similarity:
        return False
    if rec_ovl < reciprocal_overlap:
        return False
    if del_seq and med_seq:
        seq_sim = compute_sequence_similarity(del_seq, med_seq)
        if seq_sim < sequence_similarity:
            return False
    if abs(del_start - med_start) + abs(del_end - med_end) > breakpoint_window:
        return False

    return True

bed = pd.read_csv("~{intersect_bed}", sep='\t', header=None)

size_similarity = float(~{del_size_similarity})
reciprocal_overlap = float(~{del_reciprocal_overlap})
sequence_similarity = float(~{del_sequence_similarity})
breakpoint_window = int(~{del_breakpoint_window})

annotations = []

for _, row in bed.iterrows():
    # bedtools intersect -wa -wb outputs:
    # A (del.bed): 0=chrom, 1=pos, 2=end, 3=ref, 4=alt, 5=id
    # B (MED): 6=chrom, 7=start, 8=end, 9=ID, 10=seq, 11=designation, 12=sub_family
    del_chrom = row[0]
    del_start = int(row[1])
    del_end = int(row[2])
    del_ref = str(row[3])
    del_alt = row[4]
    del_id = row[5]

    med_start = int(row[7])
    med_end = int(row[8])
    med_seq = str(row[10])
    designation = row[11]
    sub_family = row[12]

    # strip VCF anchor base; empty for symbolic deletions (REF = "N")
    del_seq = del_ref[1:] if len(del_ref) > 1 else ""

    me_type = get_me_type(designation)

    if not me_type:
        continue

    if passes_criteria(del_start, del_end, del_seq, med_start, med_end, med_seq, size_similarity, reciprocal_overlap, sequence_similarity, breakpoint_window):
        annotations.append([
            del_chrom,
            del_start,
            del_ref,
            del_alt,
            del_id,
            me_type,
            sub_family
        ])

out_df = pd.DataFrame(annotations, columns=['CHROM', 'POS', 'REF', 'ALT', 'ID', 'ME_TYPE', 'SUB_FAMILY'])
out_df.drop_duplicates(subset=['ID'], inplace=True)
out_df.sort_values(['CHROM', 'POS'], inplace=True)
out_df.to_csv("~{prefix}.annotations.tsv", sep='\t', index=False, header=False)

EOF
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(intersect_bed, "GB")) + 5,
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
