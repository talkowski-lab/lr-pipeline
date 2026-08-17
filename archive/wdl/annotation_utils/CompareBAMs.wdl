version 1.0

import "../utils/Structs.wdl"

workflow CompareBAMs {
    input {
        File bam1
        File bam2
        String bam1_name
        String bam2_name
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_compare_bams
    }

    call CompareBamReadCounts {
        input:
            bam1 = bam1,
            bam2 = bam2,
            bam1_name = bam1_name,
            bam2_name = bam2_name,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare_bams
    }

    output {
        File comparison_tsv = CompareBamReadCounts.comparison_tsv
        File per_read_tsv = CompareBamReadCounts.per_read_tsv
    }
}

task CompareBamReadCounts {
    input {
        File bam1
        File bam2
        String bam1_name
        String bam2_name
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Stream each BAM, extract read_id + seq MD5 + seq_len, sort by read_id in parallel
        samtools view \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{bam1} \
            | python3 -c '
import sys, hashlib
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    seq = f[9]
    h = hashlib.md5(seq.encode()).hexdigest() if seq != "*" else ""
    length = 0 if seq == "*" else len(seq)
    print(f"{f[0]}\t{h}\t{length}")
' | LC_ALL=C sort -k1,1 -T . -S 4G > bam1_reads.tsv &
        PID1=$!

        samtools view \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{bam2} \
            | python3 -c '
import sys, hashlib
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    seq = f[9]
    h = hashlib.md5(seq.encode()).hexdigest() if seq != "*" else ""
    length = 0 if seq == "*" else len(seq)
    print(f"{f[0]}\t{h}\t{length}")
' | LC_ALL=C sort -k1,1 -T . -S 4G > bam2_reads.tsv &
        PID2=$!

        wait $PID1
        wait $PID2

        BAM1_COUNT=$(wc -l < bam1_reads.tsv)
        BAM2_COUNT=$(wc -l < bam2_reads.tsv)

        # Inner join: read_id\tbam1_md5\tbam1_len\tbam2_md5\tbam2_len
        LC_ALL=C join -t $'\t' -1 1 -2 1 bam1_reads.tsv bam2_reads.tsv > matched_reads.tsv

        MATCHED_COUNT=$(wc -l < matched_reads.tsv)
        SAME_LEN_COUNT=$(awk -F'\t' '$3 == $5 {n++} END {print n+0}' matched_reads.tsv)
        SAME_SEQ_COUNT=$(awk -F'\t' '$2 == $4 {n++} END {print n+0}' matched_reads.tsv)

        # Full outer join for per-read length TSV
        printf "read_id\t~{bam1_name}_len\t~{bam2_name}_len\n" > ~{prefix}.per_read_lengths.tsv
        LC_ALL=C join -t $'\t' -1 1 -2 1 -a 1 -a 2 -e "" -o "0,1.3,2.3" \
            bam1_reads.tsv bam2_reads.tsv >> ~{prefix}.per_read_lengths.tsv

        # Summary TSV
        printf "metric\tvalue\n" > ~{prefix}.bam_comparison.tsv
        printf "~{bam1_name}_total_reads\t%s\n" "$BAM1_COUNT" >> ~{prefix}.bam_comparison.tsv
        printf "~{bam2_name}_total_reads\t%s\n" "$BAM2_COUNT" >> ~{prefix}.bam_comparison.tsv
        printf "matched_id_reads\t%s\n" "$MATCHED_COUNT" >> ~{prefix}.bam_comparison.tsv
        printf "matched_id_reads_same_sequence_length\t%s\n" "$SAME_LEN_COUNT" >> ~{prefix}.bam_comparison.tsv
        printf "matched_id_reads_same_sequence\t%s\n" "$SAME_SEQ_COUNT" >> ~{prefix}.bam_comparison.tsv
    >>>

    output {
        File comparison_tsv = "~{prefix}.bam_comparison.tsv"
        File per_read_tsv = "~{prefix}.per_read_lengths.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: ceil(size(bam1, "GB") + size(bam2, "GB")) + 20,
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
