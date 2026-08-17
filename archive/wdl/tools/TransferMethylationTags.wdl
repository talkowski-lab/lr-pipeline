version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow TransferMethylationTags {
    input {
        File aligned_bam
        File aligned_bai
        Array[String] contigs
        String prefix

        Array[String] unaligned_bam_paths
        Boolean gcs_paths = false
        Boolean recreate_bam = false
        String mm_tag = "MM"
        String ml_tag = "ML"

        String utils_docker

        RuntimeAttr? runtime_attr_extract_tags
        RuntimeAttr? runtime_attr_merge_tags
        RuntimeAttr? runtime_attr_extract_contig
        RuntimeAttr? runtime_attr_transfer_tags
        RuntimeAttr? runtime_attr_sort_contig
        RuntimeAttr? runtime_attr_merge_bams
        RuntimeAttr? runtime_attr_merge_unaligned_bams
    }

    if (recreate_bam) {
        call MergeUnalignedBams {
            input:
                unaligned_bam_paths = unaligned_bam_paths,
                gcs_paths = gcs_paths,
                prefix = "~{prefix}.tagged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_unaligned_bams
        }
    }

    if (!recreate_bam) {
        scatter (path in unaligned_bam_paths) {
            call ExtractMethylationTags {
                input:
                    unaligned_bam_path = path,
                    gcs_paths = gcs_paths,
                    mm_tag = mm_tag,
                    ml_tag = ml_tag,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_tags
            }
        }

        call Helpers.ConcatTsvs as MergeTagsTsvs {
            input:
                tsvs = ExtractMethylationTags.tags_tsv,
                sort_output = false,
                compressed_tsvs = true,
                compressed_output = true,
                prefix = "~{prefix}.all.tags",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_tags
        }

        scatter (contig in contigs) {
            call Helpers.SubsetBamToContig {
                input:
                    bam = aligned_bam,
                    bai = aligned_bai,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_contig
            }

            call TransferTagsToContig {
                input:
                    contig_bam = SubsetBamToContig.contig_bam,
                    tags_tsv = MergeTagsTsvs.concatenated_tsv,
                    prefix = "~{prefix}.~{contig}.tagged",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_transfer_tags
            }

            call SortIndexBam {
                input:
                    unsorted_bam = TransferTagsToContig.tagged_unsorted_bam,
                    prefix = "~{prefix}.~{contig}.sorted",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_sort_contig
            }
        }

        call Helpers.MergeBams {
            input:
                bams = SortIndexBam.tagged_bam,
                bais = SortIndexBam.tagged_bai,
                prefix = "~{prefix}.tagged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_bams
        }
    }

    output {
        File methylation_tagged_bam = select_first([MergeUnalignedBams.merged_bam, MergeBams.merged_bam])
        File methylation_tagged_bai = select_first([MergeUnalignedBams.merged_bam_idx, MergeBams.merged_bam_idx])
        File? methylation_tags = MergeTagsTsvs.concatenated_tsv
    }
}

task MergeUnalignedBams {
    input {
        Array[String] unaligned_bam_paths
        Boolean gcs_paths
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        while IFS= read -r path; do
            bam=$(basename "$path")
            ~{if gcs_paths then "gsutil cp" else "aws s3 --no-sign-request cp"} "$path" "$bam"
            echo "$bam" >> bam_list.txt
        done < ~{write_lines(unaligned_bam_paths)}

        samtools merge \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            -f \
            -o ~{prefix}.bam \
            -b bam_list.txt

        samtools index \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{prefix}.bam
    >>>

    output {
        File merged_bam = "~{prefix}.bam"
        File merged_bam_idx = "~{prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 2000,
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

task ExtractMethylationTags {
    input {
        String unaligned_bam_path
        Boolean gcs_paths
        String mm_tag
        String ml_tag
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String bam_basename = basename(unaligned_bam_path, ".bam")

    command <<<
        set -euo pipefail

        ~{if gcs_paths then "gsutil cp" else "aws s3 --no-sign-request cp"} \
            ~{unaligned_bam_path} \
            ~{bam_basename}.bam

        python3 <<CODE
import gzip
import pysam

with pysam.AlignmentFile("~{bam_basename}.bam", "rb", check_sq=False) as ubam:
    with gzip.open("~{bam_basename}.tags.tsv.gz", "wt") as out:
        for read in ubam.fetch(until_eof=True):
            try:
                mm = read.get_tag('~{mm_tag}')
            except KeyError:
                mm = ''

            try:
                ml = ','.join(str(v) for v in read.get_tag('~{ml_tag}'))
            except KeyError:
                ml = ''

            out.write(f"{read.query_name}\t{mm}\t{ml}\n")
CODE

        rm ~{bam_basename}.bam
    >>>

    output {
        File tags_tsv = "~{bam_basename}.tags.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 800,
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

task TransferTagsToContig {
    input {
        File contig_bam
        File tags_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import array
import gzip
import pysam

# First pass: collect the read names present in this contig to only load relevant tags
contig_read_names = set()
with pysam.AlignmentFile("~{contig_bam}", "rb") as abam:
    for read in abam:
        contig_read_names.add(read.query_name)

# Second pass: stream through tags to only keep contig reads
tags_dict = {}
with gzip.open("~{tags_tsv}", 'rt') as f:
    for line in f:
        read_name, mm, ml_str = line.rstrip('\n').split('\t')
        if read_name in contig_read_names:
            tags_dict[read_name] = (mm, [int(v) for v in ml_str.split(',') if v])

# Third pass: stream through contig BAM to add tags for matched reads
with pysam.AlignmentFile("~{contig_bam}", "rb") as abam:
    with pysam.AlignmentFile("~{prefix}.unsorted.bam", "wb", header=abam.header) as outbam:
        for read in abam:
            if read.query_name in tags_dict:
                mm, ml = tags_dict[read.query_name]
                read.set_tag('MM', mm)
                read.set_tag('ML', array.array('B', ml))
            outbam.write(read)
CODE
    >>>

    output {
        File tagged_unsorted_bam = "~{prefix}.unsorted.bam"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(contig_bam, "GB")) + ceil(size(tags_tsv, "GB")) + 10,
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

task SortIndexBam {
    input {
        File unsorted_bam
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        samtools sort \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            -o ~{prefix}.bam \
            ~{unsorted_bam}

        samtools index \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{prefix}.bam
    >>>

    output {
        File tagged_bam = "~{prefix}.bam"
        File tagged_bai = "~{prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 5 * ceil(size(unsorted_bam, "GB")) + 10,
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
