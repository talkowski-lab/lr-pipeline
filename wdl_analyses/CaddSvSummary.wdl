version 1.0

# CADD (GRCh38 v1.7, precomputed) summaries per SV, from a BED of SV intervals (chrom start end name svtype, e.g. from
# scripts/miscellaneous/sv_bed_for_cadd.py: insertions as POS +/- 10 bp windows). Two summaries per SV:
#   sv_summary          over the whole interval
#   breakpoint_summary  over +/- breakpoint_flank bp windows at the two breakpoints (first base start+1 and last base
#                       end; the union when they overlap); insertions (svtype INS / INS:*) use their interval as is
# Columns of both: chrom start end name svtype n_bases n_snv mean_phred max_phred n_ge10 n_ge15 n_ge20 n_ge30 n_indel
#   max_indel_phred (breakpoint_summary adds left_max_phred right_max_phred). n_bases = positions summarized; n_snv =
#   scored SNVs (all 3 alternative alleles per base); indels = CADD-scored gnomAD v4 indels starting in the interval;
#   NA when there is nothing to summarize.
# 1. SplitByChrom: one BED per chromosome (input order kept within a chromosome).
# 2. ScoreChrom (one job per chromosome): the chromosome's CADD SNV and indel scores are streamed once (pysam; https or
#    gs://, resuming after read errors) into per-base arrays, from which every SV of the chromosome is summarized.
# 3. MergeSummaries: one sorted, bgzipped, tabix-indexed table per summary.
# cadd_snv_tsv / cadd_indel_tsv: https:// (default: the CADD server) or gs:// copies with .tbi; runs whose workers have no
# internet access need gs:// copies.

import "Structs.wdl"

workflow CaddSvSummary {
    input {
        File sv_bed
        String prefix

        Int breakpoint_flank = 10
        String cadd_snv_tsv = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz"
        String cadd_indel_tsv = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz"
        String? requester_pays_project

        String cadd_docker

        RuntimeAttr? runtime_attr_split_by_chrom
        RuntimeAttr? runtime_attr_score_chrom
        RuntimeAttr? runtime_attr_merge_summaries
    }

    parameter_meta {
        sv_bed: "SV intervals: chrom start end name svtype (0-based start; header/comment lines skipped)"
        prefix: "Prefix for output files"
        breakpoint_flank: "bp on each side of each breakpoint for the breakpoint summary"
        cadd_snv_tsv: "CADD whole-genome SNV score file (https:// or gs://, .tbi next to it)"
        cadd_indel_tsv: "CADD gnomAD v4 indel score file (https:// or gs://, .tbi next to it)"
        requester_pays_project: "Google project billed when the CADD copies are in a requester-pays bucket"
        cadd_docker: "Image with python3, pysam, numpy and pandas (e.g. gd-breakpoint-profiles)"
    }

    call SplitByChrom {
        input:
            sv_bed = sv_bed,
            prefix = prefix,
            docker = cadd_docker,
            runtime_attr_override = runtime_attr_split_by_chrom
    }

    scatter (chrom_bed in SplitByChrom.chrom_beds) {
        call ScoreChrom {
            input:
                chrom_bed = chrom_bed,
                breakpoint_flank = breakpoint_flank,
                cadd_snv_tsv = cadd_snv_tsv,
                cadd_indel_tsv = cadd_indel_tsv,
                requester_pays_project = requester_pays_project,
                prefix = basename(chrom_bed, ".bed"),
                docker = cadd_docker,
                runtime_attr_override = runtime_attr_score_chrom
        }
    }

    call MergeSummaries {
        input:
            sv_summaries = ScoreChrom.sv_summary,
            breakpoint_summaries = ScoreChrom.breakpoint_summary,
            prefix = prefix,
            docker = cadd_docker,
            runtime_attr_override = runtime_attr_merge_summaries
    }

    output {
        File sv_summary = MergeSummaries.sv_summary
        File sv_summary_idx = MergeSummaries.sv_summary_idx
        File breakpoint_summary = MergeSummaries.breakpoint_summary
        File breakpoint_summary_idx = MergeSummaries.breakpoint_summary_idx
    }
}

task SplitByChrom {
    input {
        File sv_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p chroms
        if [[ "~{sv_bed}" == *.gz ]]; then gzip -cd ~{sv_bed}; else cat ~{sv_bed}; fi \
            | grep -v -E '^(#|track|browser)' \
            | awk -F'\t' -v OFS='\t' -v p="chroms/~{prefix}" 'NF >= 5 && $2 ~ /^[0-9]+$/ {print > (p "." $1 ".bed")}'
        for f in chroms/*.bed; do echo "$(basename "${f}"): $(wc -l < "${f}") SVs" >&2; done
    >>>

    output {
        Array[File] chrom_beds = glob("chroms/*.bed")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 4 * ceil(size(sv_bed, "GB")) + 10,
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

task ScoreChrom {
    input {
        File chrom_bed
        Int breakpoint_flank
        String cadd_snv_tsv
        String cadd_indel_tsv
        String? requester_pays_project
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [[ "~{cadd_snv_tsv}" == gs://* || "~{cadd_indel_tsv}" == gs://* ]]; then
            export GCS_OAUTH_TOKEN=$(curl -sS -H "Metadata-Flavor: Google" \
                http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
                | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')
        fi
        ~{"export GCS_REQUESTER_PAYS_PROJECT=" + requester_pays_project}

        python3 - ~{chrom_bed} ~{breakpoint_flank} ~{cadd_snv_tsv} ~{cadd_indel_tsv} ~{prefix} <<'PY'
import io
import sys
import time

import numpy as np
import pandas as pd
import pysam

bed_path, flank, snv_url, indel_url, prefix = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
sv = pd.read_csv(bed_path, sep="\t", header=None, usecols=range(5), names=["chrom", "start", "end", "name", "svtype"],
                 dtype={"chrom": str, "start": np.int64, "end": np.int64, "name": str, "svtype": str})
chrom = sv.chrom.iloc[0]
cc = chrom[3:] if chrom.startswith("chr") else chrom
L = int(sv.end.max()) + flank + 2  # array index = 1-based position
slot = {"A": 0, "C": 1, "G": 2, "T": 3}


def stream(url, contig, handle, batch=1_000_000):
    """Stream all lines of contig (positions < L) from a tabix-indexed file (local, https or gs:// via pysam), parsed in
    batches with pandas into handle(chunk); after a read error the file is reopened and the stream resumes at the last
    position read (handle must be idempotent per position)."""
    last, fails = None, 0

    def flush(buf):
        nonlocal last
        if not buf:
            return
        types = {"pos": np.int64, "ref": str, "alt": str, "phred": np.float32}
        c = pd.read_csv(io.StringIO("\n".join(buf)), sep="\t", header=None, usecols=[1, 2, 3, 5],
                        names=["pos", "ref", "alt", "phred"], dtype=types)
        c = c[c.pos < L]
        if len(c):
            handle(c)
            last = int(c.pos.iloc[-1])

    while True:
        try:
            tf = pysam.TabixFile(url)
            buf = []
            for line in tf.fetch(contig, 0 if last is None else last - 1, L):
                buf.append(line)
                if len(buf) >= batch:
                    flush(buf)
                    buf = []
            flush(buf)
            return
        except ValueError:
            return  # contig not in the file
        except OSError as e:
            fails += 1
            print(f"  {url} {contig}: read error after position {last} ({e}); retry {fails}/8", file=sys.stderr)
            if fails > 8:
                raise
            time.sleep(20 * fails)


# SNVs: per base and alternative allele, PHRED x 100 + 1 (0 = no score)
snv = np.zeros((L, 4), dtype=np.uint16)


def put_snv(c):
    s = c.alt.map(slot)
    ok = s.notna().values
    snv[c.pos.values[ok], s.values[ok].astype(np.int64)] = np.round(c.phred.values[ok] * 100).astype(np.uint16) + 1


t0 = time.time()
stream(snv_url, cc, put_snv)
print(f"{chrom}: SNV scores read in {time.time() - t0:.0f} s", file=sys.stderr)

# per-base summaries, in blocks
cnt = np.zeros(L, dtype=np.uint8)
psum = np.zeros(L, dtype=np.float32)
pmax = np.zeros(L, dtype=np.uint16)
thr = [10, 15, 20, 30]
ge = np.zeros((len(thr), L), dtype=np.uint8)
for a in range(0, L, 10_000_000):
    b = min(L, a + 10_000_000)
    x = snv[a:b]
    v = x > 0
    cnt[a:b] = v.sum(1)
    psum[a:b] = np.where(v, (x.astype(np.float32) - 1) / 100, 0).sum(1)
    pmax[a:b] = x.max(1)
    for k, t in enumerate(thr):
        ge[k, a:b] = (x >= t * 100 + 1).sum(1)
snv = None  # free the per-allele array

# indels: collected (few per chromosome), duplicates from a resumed stream dropped, then per start position the max
# PHRED x 100 + 1 and the count
indels = []
t0 = time.time()
stream(indel_url, cc, lambda c: indels.append(c))
imax = np.zeros(L, dtype=np.uint16)
icnt = np.zeros(L, dtype=np.uint8)
if indels:
    ind = pd.concat(indels).drop_duplicates(["pos", "ref", "alt"])
    p = ind.pos.values
    np.maximum.at(imax, p, (np.round(ind.phred.values * 100) + 1).astype(np.uint16))
    up, n = np.unique(p, return_counts=True)
    icnt[up] = np.minimum(n, 255)
print(f"{chrom}: {sum(len(c) for c in indels)} indel scores read in {time.time() - t0:.0f} s", file=sys.stderr)


def stats(pos_slices):
    """Summary over the union of array slices [(a, b), ...] (array index = 1-based position)."""
    if len(pos_slices) == 1:
        a, b = pos_slices[0]
        take = [arr[a:b] for arr in (cnt, psum, pmax, icnt, imax)] + [ge[k, a:b] for k in range(len(thr))]
    else:
        idx = np.unique(np.concatenate([np.arange(a, b) for a, b in pos_slices]))
        take = [arr[idx] for arr in (cnt, psum, pmax, icnt, imax)] + [ge[k, idx] for k in range(len(thr))]
    c, ps, pm, ic, im = take[:5]
    n_bases, n = len(c), int(c.sum())
    mx = int(pm.max()) if n_bases else 0
    imx = int(im.max()) if n_bases else 0
    return ([n_bases, n, f"{ps.sum() / n:.3f}" if n else "NA", f"{(mx - 1) / 100:.2f}" if mx else "NA"]
            + [int(g.sum()) for g in take[5:]] + [int(ic.sum()), f"{(imx - 1) / 100:.2f}" if imx else "NA"])


def fmt_max(arr_slice):
    m = int(arr_slice.max()) if len(arr_slice) else 0
    return f"{(m - 1) / 100:.2f}" if m else "NA"


with open(f"{prefix}.sv_summary.tsv", "w") as o1, open(f"{prefix}.breakpoint_summary.tsv", "w") as o2:
    for r in sv.itertuples(index=False):
        a, b = max(1, r.start + 1), min(L, r.end + 1)  # positions start+1..end
        base = [r.chrom, r.start, r.end, r.name, r.svtype]
        o1.write("\t".join(map(str, base + stats([(a, max(a, b))]))) + "\n")
        if r.svtype == "INS" or r.svtype.startswith("INS:"):
            wins = [(a, max(a, b))]
            lmax = rmax = fmt_max(pmax[a:b])
        else:
            l0, l1 = max(1, r.start + 1 - flank), min(L, r.start + 1 + flank + 1)
            r0, r1 = max(1, r.end - flank), min(L, r.end + flank + 1)
            wins = [(l0, max(l1, r1))] if r0 <= l1 else [(l0, l1), (r0, r1)]
            lmax, rmax = fmt_max(pmax[l0:l1]), fmt_max(pmax[r0:r1])
        o2.write("\t".join(map(str, base + stats(wins) + [lmax, rmax])) + "\n")
print(f"{chrom}: {len(sv)} SVs summarized", file=sys.stderr)
PY
    >>>

    output {
        File sv_summary = "~{prefix}.sv_summary.tsv"
        File breakpoint_summary = "~{prefix}.breakpoint_summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
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

task MergeSummaries {
    input {
        Array[File] sv_summaries
        Array[File] breakpoint_summaries
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cols="chrom\tstart\tend\tname\tsvtype\tn_bases\tn_snv\tmean_phred\tmax_phred\tn_ge10\tn_ge15\tn_ge20\tn_ge30\tn_indel\tmax_indel_phred"
        merge() {
            local list=$1 out=$2 header=$3
            { printf "#%b\n" "${header}"; while read -r f; do cat "${f}"; done < "${list}" \
                | LC_ALL=C sort -t $'\t' -k1,1V -k2,2n -k3,3n -S 2G -T .; } | bgzip -c > "${out}"
            tabix -f -s 1 -b 2 -e 3 -0 -c '#' "${out}"
            echo "$(basename "${out}"): $(( $(zcat "${out}" | wc -l) - 1 )) SVs" >&2
        }
        merge ~{write_lines(sv_summaries)} ~{prefix}.cadd.sv_summary.bed.gz "${cols}"
        merge ~{write_lines(breakpoint_summaries)} ~{prefix}.cadd.breakpoint_summary.bed.gz "${cols}\tleft_max_phred\tright_max_phred"
    >>>

    output {
        File sv_summary = "~{prefix}.cadd.sv_summary.bed.gz"
        File sv_summary_idx = "~{prefix}.cadd.sv_summary.bed.gz.tbi"
        File breakpoint_summary = "~{prefix}.cadd.breakpoint_summary.bed.gz"
        File breakpoint_summary_idx = "~{prefix}.cadd.breakpoint_summary.bed.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 4 * ceil(size(sv_summaries, "GB") + size(breakpoint_summaries, "GB")) + 20,
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
