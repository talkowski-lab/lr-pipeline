version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AlignQueriesIterative {
    input {
        File queries_fa
        String prefix

        Int min_piece_len = 30
        Int max_rounds = 20
        Int sr_max_len = 500

        File ref_fa

        String minimap_docker

        RuntimeAttr? runtime_attr_iterative_align
    }

    call IterativeAlign {
        input:
            queries_fa = queries_fa,
            min_piece_len = min_piece_len,
            max_rounds = max_rounds,
            sr_max_len = sr_max_len,
            ref_fa = ref_fa,
            prefix = prefix,
            docker = minimap_docker,
            runtime_attr_override = runtime_attr_iterative_align
    }

    output {
        File segments_tsv = IterativeAlign.segments_tsv
    }
}

task IterativeAlign {
    input {
        File queries_fa
        Int min_piece_len
        Int max_rounds
        Int sr_max_len
        File ref_fa
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        minimap2 -x asm20 -t "$(nproc)" -d ref.asm20.mmi ~{ref_fa}
        minimap2 -x sr -t "$(nproc)" -d ref.sr.mmi ~{ref_fa}

        cat > iterate.py <<'EOF'
import re
import subprocess
import sys

queries_fa, min_len, max_rounds, sr_max, threads, out = sys.argv[1:7]
min_len, max_rounds, sr_max = int(min_len), int(max_rounds), int(sr_max)
CIG = re.compile(r"(\d+)([MIDNSHP=X])")

seqs, name = {}, None
for line in open(queries_fa):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]
        seqs[name] = []
    elif line:
        seqs[name].append(line.upper())
seqs = {k: "".join(v) for k, v in seqs.items()}


def align(pieces, preset, index):
    best = {}
    if not pieces:
        return best
    with open("pieces.fa", "w") as o:
        for q, s, e in pieces:
            o.write(f">{q}|{s}|{e}\n{seqs[q][s:e]}\n")
    sam = subprocess.run(["minimap2", "-a", "-x", preset, "-t", threads, index, "pieces.fa"],
                         check=True, capture_output=True, text=True).stdout
    for line in sam.splitlines():
        if line.startswith("@"):
            continue
        c = line.split("\t")
        flag = int(c[1])
        if flag & 4 or flag & 256 or flag & 2048:
            continue
        ops = [(int(n), op) for n, op in CIG.findall(c[5])]
        lead = ops[0][0] if ops[0][1] in "SH" else 0
        tail = ops[-1][0] if ops[-1][1] in "SH" else 0
        m = sum(n for n, op in ops if op in "M=X")
        ins = sum(n for n, op in ops if op == "I")
        dl = sum(n for n, op in ops if op in "DN")
        qlen = sum(n for n, op in ops if op in "MIS=XH")
        nm = next((int(t[5:]) for t in c[11:] if t.startswith("NM:i:")), 0)
        rev = bool(flag & 16)
        qs, qe = (tail, qlen - lead) if rev else (lead, qlen - tail)
        rs = int(c[3])
        best[c[0]] = dict(qs=qs, qe=qe, chrom=c[2], rs=rs, re=rs + m + dl - 1, strand="-" if rev else "+",
                          ident=(m - (nm - ins - dl)) / max(1, m + ins + dl), mapq=int(c[4]), preset=preset)
    return best


# pieces still to align: (query ID, start, end), 0-based half-open on the inserted sequence
todo = [(q, 0, len(s)) for q, s in seqs.items() if len(s) >= min_len]
segments = [(q, 0, len(s), "", "", "", "", "", "", 0, "shorter than min_piece_len")
            for q, s in seqs.items() if len(s) < min_len]
for rnd in range(1, max_rounds + 1):
    if not todo:
        break
    best = align([p for p in todo if p[2] - p[1] > sr_max], "asm20", "ref.asm20.mmi")
    best.update(align([p for p in todo if p[2] - p[1] <= sr_max], "sr", "ref.sr.mmi"))
    nxt = []
    for q, s, e in todo:
        b = best.get(f"{q}|{s}|{e}")
        if b is None:
            segments.append((q, s, e, "", "", "", "", "", "", rnd, "unaligned"))
            continue
        a0, a1 = s + b["qs"], s + b["qe"]
        segments.append((q, a0, a1, f"{b['chrom']}:{b['rs']}-{b['re']}", b["strand"], f"{b['ident']:.3f}", b["mapq"],
                         a1 - a0, b["preset"], rnd, "aligned"))
        for p0, p1 in ((s, a0), (a1, e)):
            if p1 - p0 >= min_len:
                nxt.append((q, p0, p1))
            elif p1 > p0:
                segments.append((q, p0, p1, "", "", "", "", "", "", rnd, "shorter than min_piece_len"))
    todo = nxt
for q, s, e in todo:
    segments.append((q, s, e, "", "", "", "", "", "", max_rounds, "not aligned within max_rounds"))
with open(out, "w") as o:
    o.write("ID\tq_start\tq_end\tref\tstrand\tidentity\tmapq\taligned_len\tpreset\tround\tstatus\n")
    for r in sorted(segments, key=lambda r: (r[0], r[1])):
        o.write(f"{r[0]}\t{r[1] + 1}\t{r[2]}\t" + "\t".join(str(v) for v in r[3:]) + "\n")
EOF

        python3 iterate.py ~{queries_fa} ~{min_piece_len} ~{max_rounds} ~{sr_max_len} "$(nproc)" ~{prefix}.iterative_segments.tsv
    >>>

    output {
        File segments_tsv = "~{prefix}.iterative_segments.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 32,
        disk_gb: 6 * ceil(size(ref_fa, "GB")) + ceil(size(queries_fa, "GB")) + 10,
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
