version 1.0

import "../utils/BedtoolsClosestSV.wdl"
import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "../utils/TruvariMatch.wdl"

task CollectMatchedIDsAndINFO {
    input {
        File annotation_tsv
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth_snv
        File vcf_truth_snv_idx
        File vcf_truth_sv
        File vcf_truth_sv_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import subprocess

annotation_tsv = "~{annotation_tsv}"
vcf_eval = "~{vcf_eval}"
vcf_truth_snv = "~{vcf_truth_snv}"
vcf_truth_sv = "~{vcf_truth_sv}"
prefix = "~{prefix}"

eval_to_truth = {}
eval_ids = set()
truth_ids = set()

with open(annotation_tsv) as f:
    for line in f:
        fields = line.strip().split('\t')
        eval_id = fields[4]
        truth_id = fields[6]
        eval_to_truth[eval_id] = truth_id
        eval_ids.add(eval_id)
        truth_ids.add(truth_id)

eval_info = {}
cmd = f"bcftools query -f '%ID\\t%INFO\\n' {vcf_eval}"
proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
for line in proc.stdout:
    parts = line.strip().split('\t', 1)
    if len(parts) == 2 and parts[0] in eval_ids:
        eval_info[parts[0]] = parts[1]
proc.wait()

truth_info = {}
for vcf in [vcf_truth_snv, vcf_truth_sv]:
    cmd = f"bcftools query -f '%ID\\t%INFO\\n' {vcf}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        parts = line.strip().split('\t', 1)
        if len(parts) == 2 and parts[0] in truth_ids:
            truth_info[parts[0]] = parts[1]
    proc.wait()

with open(f"{prefix}.matched_with_info.tsv", 'w') as out:
    for eval_id, truth_id in eval_to_truth.items():
        eval_inf = eval_info.get(eval_id, '.')
        truth_inf = truth_info.get(truth_id, '.')
        out.write(f"{eval_id}\t{truth_id}\t{eval_inf}\t{truth_inf}\n")

EOF
    >>>

    output {
        File matched_with_info_tsv = "~{prefix}.matched_with_info.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf_eval, "GB") + size(vcf_truth_snv, "GB") + size(vcf_truth_sv, "GB")) + 10,
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

task ExtractVepHeader {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -h \
            ~{vcf} \
        | awk 'BEGIN{IGNORECASE=1} /^##INFO=<ID=(vep|csq),/ {print}' > ~{prefix}_vep_header.txt
    >>>

    output {
        File vep_header_txt = "~{prefix}_vep_header.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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

task ShardedMatchedVariants {
    input {
        File matched_with_info_tsv
        Int records_per_shard
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p shards

        cat ~{matched_with_info_tsv} \
            | awk 'BEGIN{c=0;f=0} {print > sprintf("shards/matched.%06d.tsv", int(c/~{records_per_shard})) ; c++} END{ }'
    >>>

    output {
        Array[File] shard_tsvs = glob("shards/*.tsv")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(matched_with_info_tsv, "GB")) + 5,
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

task ComputeShardBenchmarks {
    input {
        File matched_shard_tsv
        File eval_vep_header
        File truth_vep_header
        String skip_vep_categories
        String contig
        String shard_label
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/gnomad-lr/scripts/benchmark/compute_benchmarks_shard.py \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --matched_shard_tsv ~{matched_shard_tsv} \
            --eval_vep_header ~{eval_vep_header} \
            --truth_vep_header ~{truth_vep_header} \
            --shard_label ~{shard_label} \
            ~{if skip_vep_categories != "" then "--skip_vep_categories " + skip_vep_categories else ""}
    >>>

    output {
        File af_pairs_tsv = "~{prefix}.shard_~{shard_label}.af_pairs.tsv"
        File vep_pairs_tsv = "~{prefix}.shard_~{shard_label}.vep_pairs.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(matched_shard_tsv, "GB")) + 5,
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

task MergeShardBenchmarks {
    input {
        Array[File] af_pair_tsvs
        Array[File] vep_pair_tsvs
        String skip_vep_categories
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/gnomad-lr/scripts/benchmark/merge_benchmarks_from_pairs.py \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --af_pair_tsvs ~{sep=',' af_pair_tsvs} \
            --vep_pair_tsvs ~{sep=',' vep_pair_tsvs} \
            ~{if skip_vep_categories != "" then "--skip_vep_categories " + skip_vep_categories else ""}
    >>>

    output {
        File plot_tarball = "~{prefix}.benchmarks.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 50 * ceil(size(af_pair_tsvs, "GB")) + 5,
        boot_disk_gb: 50,
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

task ComputeSummaryForContig {
    input {
        File eval_vcf
        File eval_vcf_idx
        File annotation_tsv
        File matched_with_info_tsv
        File eval_vep_header
        File truth_vep_header
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/gnomad-lr/scripts/benchmark/compute_summary_for_contig.py \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --eval_vcf ~{eval_vcf} \
            --annotation_tsv ~{annotation_tsv} \
            --matched_with_info_tsv ~{matched_with_info_tsv} \
            --eval_vep_header ~{eval_vep_header} \
            --truth_vep_header ~{truth_vep_header}
    >>>

    output {
        File benchmark_summary_tsv = "~{prefix}.benchmark_summary.tsv"
        File summary_stats_tsv = "~{prefix}.summary_stats.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 50 * ceil(size(eval_vcf, "GB")) + 5,
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

task MergePlotTarballs {
    input {
        Array[File] tarballs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p final_results/AF_plots
        mkdir -p final_results/VEP_plots

        for tarball in ~{sep=' ' tarballs}; do
            tar -xvf $tarball --strip-components=1 -C final_results
        done

        tar -czf ~{prefix}.plots.tar.gz final_results/
    >>>

    output {
        File merged_tarball = "~{prefix}.plots.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10 * ceil(size(tarballs, "GB")) + 5,
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
