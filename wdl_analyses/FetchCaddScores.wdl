version 1.0

# CADD scores (GRCh38 v1.7, precomputed) for every variant inside a set of regions, fetched remotely: the bgzipped,
# tabix-indexed CADD files are read with pysam, region by region, without downloading them.
# 1. SplitBed: the input BED (chrom start end ID, 0-based start; header and comment lines skipped) is split into shards
#    of regions_per_shard regions, in input order.
# 2. FetchCadd (one job per shard): for each region, every score in the CADD SNV file (all possible SNVs) and, with
#    include_indels, the CADD indel file (indels seen in gnomAD v4 genomes); each region is retried on a read error.
# 3. MergeCadd: shard outputs concatenated in input order, sorted, bgzipped and tabix-indexed.
# Output BED columns: chrom start end region_id ref alt cadd_raw cadd_phred cadd_source (start = pos - 1, end = pos;
# cadd_source snv or indel). The chromosome is written as in the input BED ('chr' dropped only for the lookup).
# cadd_snv_tsv / cadd_indel_tsv: https:// (default: the CADD server) or gs:// copies, each with its .tbi next to it. Runs
# whose workers have no internet access (e.g. inside a service perimeter) need gs:// copies.

import "Structs.wdl"

workflow FetchCaddScores {
    input {
        File bed
        String prefix

        Int regions_per_shard = 500
        Boolean include_indels = true
        String cadd_snv_tsv = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
        String cadd_indel_tsv = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz"
        String? requester_pays_project

        String cadd_docker

        RuntimeAttr? runtime_attr_split_bed
        RuntimeAttr? runtime_attr_fetch_cadd
        RuntimeAttr? runtime_attr_merge_cadd
    }

    parameter_meta {
        bed: "Regions: chrom start end ID (BED, 0-based start), tab-separated; header/comment lines are skipped"
        prefix: "Prefix for output files"
        regions_per_shard: "Regions per shard (one job each)"
        include_indels: "Also fetch CADD scores of gnomAD v4 indels"
        cadd_snv_tsv: "CADD whole-genome SNV score file (https:// or gs://, .tbi next to it)"
        cadd_indel_tsv: "CADD indel score file (https:// or gs://, .tbi next to it)"
        requester_pays_project: "Google project billed when the CADD copies are in a requester-pays bucket"
        cadd_docker: "Image with python3 and pysam (e.g. the read-recurrence image)"
    }

    call SplitBed {
        input:
            bed = bed,
            regions_per_shard = regions_per_shard,
            prefix = prefix,
            docker = cadd_docker,
            runtime_attr_override = runtime_attr_split_bed
    }

    scatter (shard in SplitBed.shards) {
        call FetchCadd {
            input:
                regions = shard,
                include_indels = include_indels,
                cadd_snv_tsv = cadd_snv_tsv,
                cadd_indel_tsv = cadd_indel_tsv,
                requester_pays_project = requester_pays_project,
                prefix = basename(shard, ".bed"),
                docker = cadd_docker,
                runtime_attr_override = runtime_attr_fetch_cadd
        }
    }

    call MergeCadd {
        input:
            shard_beds = FetchCadd.scores,
            shard_summaries = FetchCadd.summary,
            prefix = prefix,
            docker = cadd_docker,
            runtime_attr_override = runtime_attr_merge_cadd
    }

    output {
        File cadd_bed = MergeCadd.cadd_bed
        File cadd_bed_idx = MergeCadd.cadd_bed_idx
        File region_summary = MergeCadd.summary
    }
}

task SplitBed {
    input {
        File bed
        Int regions_per_shard
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # regions only (no header, comment, track or browser lines), split in input order; zero-padded shard numbers
        mkdir -p shards
        grep -v -E '^(#|track|browser)' ~{bed} | awk -F'\t' 'NF >= 3 && $2 ~ /^[0-9]+$/' > regions.bed
        n=$(wc -l < regions.bed)
        if [[ ${n} -eq 0 ]]; then
            echo "no regions in ~{bed}" >&2
            exit 1
        fi
        split -l ~{regions_per_shard} -d -a 5 --additional-suffix .bed regions.bed shards/~{prefix}.shard_
        echo "${n} regions -> $(ls shards | wc -l) shards of <= ~{regions_per_shard}" >&2
    >>>

    output {
        Array[File] shards = glob("shards/*.bed")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 3 * ceil(size(bed, "GB")) + 10,
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

task FetchCadd {
    input {
        File regions
        Boolean include_indels
        String cadd_snv_tsv
        String cadd_indel_tsv
        String? requester_pays_project
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # gs:// CADD copies: OAuth token from the VM metadata server (htslib reads GCS_OAUTH_TOKEN)
        if [[ "~{cadd_snv_tsv}" == gs://* || "~{cadd_indel_tsv}" == gs://* ]]; then
            export GCS_OAUTH_TOKEN=$(curl -sS -H "Metadata-Flavor: Google" \
                http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
                | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')
        fi
        ~{"export GCS_REQUESTER_PAYS_PROJECT=" + requester_pays_project}

        python3 - ~{regions} ~{prefix}.cadd.bed ~{prefix}.summary.tsv ~{cadd_snv_tsv} \
            ~{if include_indels then cadd_indel_tsv else "none"} <<'PY'
import sys
import time

import pysam

regions_path, out_path, summary_path, snv_url, indel_url = sys.argv[1:6]
sources = [("snv", snv_url)] + ([("indel", indel_url)] if indel_url != "none" else [])
tabs = {k: pysam.TabixFile(u) for k, u in sources}


def fetch(kind, chrom, start, end, tries=5):
    """All lines of a region; on a read error the file is reopened and the region read again."""
    for attempt in range(1, tries + 1):
        try:
            return list(tabs[kind].fetch(chrom, start, end))
        except ValueError:
            return []  # contig not in the CADD file (alt contigs, chrM, ...)
        except OSError as e:
            print(f"  {kind} {chrom}:{start}-{end}: read error ({e}); retry {attempt}/{tries}", file=sys.stderr)
            if attempt == tries:
                raise
            time.sleep(10 * attempt)
            tabs[kind] = pysam.TabixFile(dict(sources)[kind])


n_regions = n_scores = 0
with open(regions_path) as fh, open(out_path, "w") as out, open(summary_path, "w") as summ:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        chrom, start, end = f[0], int(f[1]), int(f[2])
        rid = f[3] if len(f) > 3 and f[3] else f"{chrom}:{start + 1}-{end}"
        cc = chrom[3:] if chrom.startswith("chr") else chrom
        counts = {}
        for kind, _ in sources:
            rows = fetch(kind, cc, start, end)
            counts[kind] = len(rows)
            for r in rows:
                c = r.split("\t")
                pos = int(c[1])
                out.write(f"{chrom}\t{pos - 1}\t{pos}\t{rid}\t{c[2]}\t{c[3]}\t{c[4]}\t{c[5]}\t{kind}\n")
        n_regions += 1
        n_scores += sum(counts.values())
        summ.write(f"{rid}\t{chrom}\t{start}\t{end}\t{counts.get('snv', 0)}\t{counts.get('indel', 0)}\n")
print(f"{n_regions} regions, {n_scores} scores -> {out_path}", file=sys.stderr)
PY
        gzip ~{prefix}.cadd.bed
    >>>

    output {
        File scores = "~{prefix}.cadd.bed.gz"
        File summary = "~{prefix}.summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 20,
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

task MergeCadd {
    input {
        Array[File] shard_beds
        Array[File] shard_summaries
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # concatenate, sort by chromosome and position (version order: chr2 before chr10), bgzip and index
        while read -r f; do gzip -cd "${f}"; done < ~{write_lines(shard_beds)} \
            | LC_ALL=C sort -t $'\t' -k1,1V -k2,2n -k5,5 -k6,6 -S 2G -T . > ~{prefix}.cadd.bed
        n=$(wc -l < ~{prefix}.cadd.bed)
        python3 - ~{prefix}.cadd.bed <<'PY'
import sys

import pysam

path = sys.argv[1]
pysam.tabix_compress(path, path + ".gz", force=True)
pysam.tabix_index(path + ".gz", seq_col=0, start_col=1, end_col=2, zerobased=True, force=True)
PY
        rm ~{prefix}.cadd.bed
        # per region: scores found in each file
        printf 'region_id\tchrom\tstart\tend\tn_snv_scores\tn_indel_scores\n' > ~{prefix}.cadd.region_summary.tsv
        while read -r f; do cat "${f}"; done < ~{write_lines(shard_summaries)} >> ~{prefix}.cadd.region_summary.tsv
        echo "${n} scores in ~{prefix}.cadd.bed.gz" >&2
    >>>

    output {
        File cadd_bed = "~{prefix}.cadd.bed.gz"
        File cadd_bed_idx = "~{prefix}.cadd.bed.gz.tbi"
        File summary = "~{prefix}.cadd.region_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 6 * ceil(size(shard_beds, "GB")) + 20,
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
