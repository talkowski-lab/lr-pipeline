version 1.0

# CADD-SV v2 (https://github.com/kircherlab/CADD-SV) coordinate-based scores for an SV site table.
# 1. PrepareShards: converts the site table (columns #chrom start end name svtype, as in svtk BEDs such as
#    gnomAD_SV_v3.PASS.key_info.gz) to CADD-SV BED records (chrom start end type) for the SV types in sv_types; DEL, DUP and
#    INV keep their interval, INS (any INS subtype, e.g. INS:ME:ALU) becomes a 1 bp interval at start; other types and
#    non-standard contigs are dropped. Records are split into shards of shard_size SVs and an SV-to-record key is written.
# 2. RunCaddSV (one job per shard): extracts the CADD-SV annotation bundle (dependencies.tar.gz, staged once to GCS) and runs
#    caddsv run with the Conda environments pre-built in the docker image.
# 3. MergeScores: joins the shard scores back to SV names by (chrom, start, end, type).
# Outputs: per-SV table (name svtype chrom start end type CADD-SV_PHRED CADD-SV_score), the full feature tables of all
# shards concatenated, and the list of SVs not scored.

import "Structs.wdl"

workflow CaddSV {
    input {
        File sv_table
        String prefix

        String sv_types = "DEL,DUP,INS,INV"
        Int shard_size = 50000
        Int threads = 8

        File annotations_tar

        String utils_docker
        String caddsv_docker

        RuntimeAttr? runtime_attr_prepare_shards
        RuntimeAttr? runtime_attr_run_caddsv
        RuntimeAttr? runtime_attr_merge_scores
    }

    parameter_meta {
        sv_table: "SV site table (gzipped TSV or BED) with #chrom, start, end, name, svtype columns"
        prefix: "Prefix for output files"
        sv_types: "Comma-separated svtype parent classes to score (CADD-SV supports DEL, DUP, INS, INV)"
        shard_size: "SVs per CADD-SV job"
        threads: "Snakemake cores per job"
        annotations_tar: "CADD-SV v2.0 annotation bundle (dependencies.tar.gz from kircherlab.bihealth.org), staged in GCS"
        utils_docker: "Image with python3 (e.g. gd-breakpoint-profiles)"
        caddsv_docker: "Image built from dockerfiles/Dockerfile.CaddSV"
    }

    call PrepareShards {
        input:
            sv_table = sv_table,
            sv_types = sv_types,
            shard_size = shard_size,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_prepare_shards
    }

    scatter (shard in PrepareShards.shards) {
        call RunCaddSV {
            input:
                shard_bed = shard,
                annotations_tar = annotations_tar,
                threads = threads,
                prefix = basename(shard, ".bed"),
                docker = caddsv_docker,
                runtime_attr_override = runtime_attr_run_caddsv
        }
    }

    call MergeScores {
        input:
            scores = RunCaddSV.scores,
            sv_key = PrepareShards.sv_key,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_scores
    }

    output {
        File caddsv_per_sv = MergeScores.per_sv
        File caddsv_features = MergeScores.features
        File caddsv_unscored = MergeScores.unscored
    }
}

task PrepareShards {
    input {
        File sv_table
        String sv_types
        Int shard_size
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import gzip

types = set("~{sv_types}".split(","))
contigs = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
opener = gzip.open if "~{sv_table}".endswith(".gz") else open
records, key, n_in = [], [], 0
with opener("~{sv_table}", "rt") as f:
    header = f.readline().rstrip("\n").lstrip("#").split("\t")
    i_c, i_s, i_e = header.index("chrom"), header.index("start"), header.index("end")
    i_n, i_t = header.index("name"), header.index("svtype")
    for line in f:
        r = line.rstrip("\n").split("\t")
        n_in += 1
        chrom, svtype = r[i_c], r[i_t]
        t = svtype.split(":")[0]
        if t not in types or chrom not in contigs:
            continue
        start = int(r[i_s])
        end = start + 1 if t == "INS" else int(r[i_e])
        if end <= start:
            continue
        records.append((chrom, start, end, t))
        key.append((r[i_n], svtype, chrom, start, end, t))
uniq = sorted(set(records), key=lambda x: (x[0], x[1], x[2], x[3]))
for k in range(0, len(uniq), ~{shard_size}):
    with open(f"~{prefix}.shard_{k // ~{shard_size}:04d}.bed", "w") as out:
        for rec in uniq[k:k + ~{shard_size}]:
            out.write("\t".join(map(str, rec)) + "\n")
with gzip.open("~{prefix}.caddsv_key.tsv.gz", "wt") as out:
    out.write("name\tsvtype\tchrom\tstart\tend\ttype\n")
    for rec in key:
        out.write("\t".join(map(str, rec)) + "\n")
print(f"{n_in} SVs read; {len(key)} of types {sorted(types)} kept; {len(uniq)} unique records in "
      f"{(len(uniq) + ~{shard_size} - 1) // ~{shard_size}} shards")
EOF
    >>>

    output {
        Array[File] shards = glob("~{prefix}.shard_*.bed")
        File sv_key = "~{prefix}.caddsv_key.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(sv_table, "GB")) + 10,
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

task RunCaddSV {
    input {
        File shard_bed
        File annotations_tar
        Int threads
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p annotations
        tar -I pigz -xf ~{annotations_tar} --strip-components=1 -C annotations
        rm -f ~{annotations_tar}

        cp ~{shard_bed} ~{prefix}.bed
        caddsv run ~{prefix}.bed \
            --annotations-dir annotations \
            --output-dir caddsv_out \
            --conda-prefix "${CADD_SV_CONDA_PREFIX}" \
            --threads ~{threads}

        gzip -c caddsv_out/scored/~{prefix}_score.tsv > ~{prefix}.caddsv_score.tsv.gz
    >>>

    output {
        File scores = "~{prefix}.caddsv_score.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 32,
        disk_gb: 4 * ceil(size(annotations_tar, "GB")) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeScores {
    input {
        Array[File] scores
        File sv_key
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import gzip

files = "~{sep=',' scores}".split(",")
score, header = {}, None
with gzip.open("~{prefix}.caddsv_features.tsv.gz", "wt") as feat:
    for fn in files:
        with gzip.open(fn, "rt") as f:
            h = f.readline().rstrip("\n").split("\t")
            if header is None:
                header = h
                feat.write("\t".join(h) + "\n")
                i_p, i_s = h.index("CADD-SV_PHRED"), h.index("CADD-SV_score")
            for line in f:
                feat.write(line)
                r = line.rstrip("\n").split("\t")
                score[(r[0], r[1], r[2], r[3])] = (r[i_p], r[i_s])
n_ok = n_miss = 0
with gzip.open("~{sv_key}", "rt") as f, gzip.open("~{prefix}.caddsv_per_sv.tsv.gz", "wt") as out, \
        open("~{prefix}.caddsv_unscored.tsv", "w") as miss:
    kh = f.readline().rstrip("\n")
    out.write(kh + "\tCADD-SV_PHRED\tCADD-SV_score\n")
    miss.write(kh + "\n")
    for line in f:
        r = line.rstrip("\n").split("\t")
        s = score.get((r[2], r[3], r[4], r[5]))
        if s is None:
            miss.write(line)
            n_miss += 1
        else:
            out.write("\t".join(r + list(s)) + "\n")
            n_ok += 1
print(f"{n_ok} SVs scored, {n_miss} not scored")
EOF
    >>>

    output {
        File per_sv = "~{prefix}.caddsv_per_sv.tsv.gz"
        File features = "~{prefix}.caddsv_features.tsv.gz"
        File unscored = "~{prefix}.caddsv_unscored.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 16,
        disk_gb: 3 * ceil(size(scores, "GB") + size(sv_key, "GB")) + 10,
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
