version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow SummarizeGDCallDepth {
    input {
        Array[File] gd_vcfs
        Array[File] gd_vcf_idxs
        Array[File] rd_matrices
        Array[File] rd_matrix_idxs
        Array[String] batches
        String prefix

        Int flank_bins = 25
        Int event_bins = 50
        Int min_flank_bp = 10000
        Int traces_per_class = 5
        Int ref_traces = 0
        File? sd_pairs_bed

        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_summarize_batch_depth
        RuntimeAttr? runtime_attr_combine_depth_summaries
    }

    scatter (i in range(length(batches))) {
        call SummarizeBatchDepth {
            input:
                gd_vcfs = gd_vcfs,
                gd_vcf_idxs = gd_vcf_idxs,
                rd_matrix = rd_matrices[i],
                rd_matrix_idx = rd_matrix_idxs[i],
                batch = batches[i],
                flank_bins = flank_bins,
                event_bins = event_bins,
                min_flank_bp = min_flank_bp,
                traces_per_class = traces_per_class,
                ref_traces = ref_traces,
                sd_pairs_bed = sd_pairs_bed,
                prefix = "~{prefix}.~{batches[i]}",
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_summarize_batch_depth
        }
    }

    call CombineDepthSummaries {
        input:
            bin_stats = SummarizeBatchDepth.bin_stats,
            traces = SummarizeBatchDepth.traces,
            windows = SummarizeBatchDepth.windows[0],
            prefix = prefix,
            docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_combine_depth_summaries
    }

    output {
        File depth_bin_stats = CombineDepthSummaries.combined_bin_stats
        File depth_traces = CombineDepthSummaries.combined_traces
        File call_windows = CombineDepthSummaries.combined_windows
    }
}

task SummarizeBatchDepth {
    input {
        Array[File] gd_vcfs
        Array[File] gd_vcf_idxs
        File rd_matrix
        File rd_matrix_idx
        String batch
        Int flank_bins
        Int event_bins
        Int min_flank_bp
        Int traces_per_class
        Int ref_traces
        File? sd_pairs_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
        import gzip
        import random
        import numpy as np
        import pysam

        batch = "~{batch}"
        flank_bins, event_bins = ~{flank_bins}, ~{event_bins}
        min_flank = ~{min_flank_bp}
        n_traces = ~{traces_per_class}
        n_ref_traces = ~{ref_traces}
        random.seed(batch)

        sd_pairs = {}
        pairs_bed = "~{default='' sd_pairs_bed}"
        if pairs_bed:
            with (gzip.open(pairs_bed, "rt") if pairs_bed.endswith(".gz") else open(pairs_bed)) as f:
                cols = f.readline().lstrip("#").rstrip("\n").split("\t")
                for line in f:
                    r = dict(zip(cols, line.rstrip("\n").split("\t")))
                    sd_pairs[(r["chrom"], int(r["start"]), int(r["end"]))] = (
                        int(r["left_sd_start"]), int(r["left_sd_end"]), int(r["right_sd_start"]), int(r["right_sd_end"]))

        with gzip.open("~{rd_matrix}", "rt") as f:
            samples = f.readline().rstrip("\n").split("\t")[3:]
        tbx = pysam.TabixFile("~{rd_matrix}", index="~{rd_matrix_idx}")

        calls = []
        for path in "~{sep=',' gd_vcfs}".split(","):
            vcf = pysam.VariantFile(path)
            vcf.subset_samples([s for s in samples if s in vcf.header.samples])
            for rec in vcf:
                cn = {s: rec.samples[s]["RD_CN"] for s in rec.samples}
                ecn = {s: rec.samples[s]["ECN"] for s in rec.samples}
                calls.append((rec.id, rec.chrom, rec.start, rec.stop, rec.info["SVTYPE"], cn, ecn))

        n_bins = 2 * flank_bins + event_bins
        profiles = []
        cache = {}
        for cid, chrom, start, end, svtype, cn, ecn in calls:
            if (chrom, start, end) not in cache:
                length = end - start
                flank = max(length, min_flank)
                lo, hi = start - flank, end + flank
                if (chrom, start, end) in sd_pairs:
                    ls, le, rs, re_ = sd_pairs[(chrom, start, end)]
                    if le - ls > length:
                        lo = min(lo, ls - length)
                    if re_ - rs > length:
                        hi = max(hi, re_ + length)
                edges = np.concatenate([
                    np.linspace(lo, start, flank_bins + 1)[:-1],
                    np.linspace(start, end, event_bins + 1)[:-1],
                    np.linspace(end, hi, flank_bins + 1)])
                depth = np.zeros((n_bins, len(samples)))
                width = np.zeros(n_bins)
                for line in tbx.fetch(chrom, max(0, int(edges[0])), int(edges[-1])):
                    row = line.split("\t")
                    s, e = int(row[1]), int(row[2])
                    b = np.searchsorted(edges, (s + e) / 2, side="right") - 1
                    if 0 <= b < n_bins:
                        depth[b] += np.array(row[3:], dtype=float)
                        width[b] += e - s
                cache[(chrom, start, end)] = (edges, depth, width)
            edges, depth, width = cache[(chrom, start, end)]
            profiles.append((cid, chrom, start, end, svtype, cn, ecn, edges, depth, width))

        flank_idx = np.r_[0:flank_bins, flank_bins + event_bins:n_bins]
        per_bp = [p[8][flank_idx] / np.maximum(p[9][flank_idx], 1)[:, None] for p in profiles]
        scale = np.median(np.vstack(per_bp), axis=0)
        scale[scale <= 0] = np.nan

        with open("~{prefix}.windows.tsv", "w") as w, gzip.open("~{prefix}.bin_stats.tsv.gz", "wt") as bs, \
                gzip.open("~{prefix}.traces.tsv.gz", "wt") as tr:
            w.write("call_id\tchrom\tstart\tend\tsvtype\tbin\tbin_start\tbin_end\tregion\n")
            bs.write("call_id\tbatch\tclass\tbin\tn\tsum\tsumsq\n")
            tr.write("call_id\tbatch\tsample\tclass\t" + "\t".join(f"b{i}" for i in range(n_bins)) + "\n")
            for cid, chrom, start, end, svtype, cn, ecn, edges, depth, width in profiles:
                for b in range(n_bins):
                    region = "left" if b < flank_bins else ("right" if b >= flank_bins + event_bins else "event")
                    w.write(f"{cid}\t{chrom}\t{start}\t{end}\t{svtype}\t{b}\t{int(edges[b])}\t{int(edges[b + 1])}\t{region}\n")
                x = depth / np.maximum(width, 1)[:, None] / scale[None, :]
                ref_cols = [j for j, s in enumerate(samples) if cn.get(s) is not None and cn.get(s) == ecn.get(s)]
                ref = np.nanmedian(x[:, ref_cols], axis=1) if ref_cols else np.nanmedian(x, axis=1)
                ref[ref <= 0] = np.nan
                cn_est = 2 * x / ref[:, None]
                classes = {}
                for j, s in enumerate(samples):
                    c = cn.get(s)
                    if c is None or ecn.get(s) is None:
                        continue
                    label = "ref" if c == ecn[s] else f"CN{c}"
                    classes.setdefault(label, []).append(j)
                for label, cols in classes.items():
                    v = cn_est[:, cols]
                    ok = np.isfinite(v)
                    n = ok.sum(axis=1)
                    sm = np.where(ok, v, 0).sum(axis=1)
                    sq = np.where(ok, v * v, 0).sum(axis=1)
                    for b in range(n_bins):
                        bs.write(f"{cid}\t{batch}\t{label}\t{b}\t{n[b]}\t{sm[b]:.4f}\t{sq[b]:.4f}\n")
                    k = n_ref_traces if label == "ref" else n_traces
                    if k > 0:
                        for j in random.sample(cols, min(k, len(cols))):
                            vals = "\t".join("NA" if not np.isfinite(t) else f"{t:.3f}" for t in cn_est[:, j])
                            tr.write(f"{cid}\t{batch}\t{samples[j]}\t{label}\t{vals}\n")
        CODE
    >>>

    output {
        File bin_stats = "~{prefix}.bin_stats.tsv.gz"
        File traces = "~{prefix}.traces.tsv.gz"
        File windows = "~{prefix}.windows.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(size(rd_matrix, "GB") + size(gd_vcfs, "GB")) + 20,
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

task CombineDepthSummaries {
    input {
        Array[File] bin_stats
        Array[File] traces
        File windows
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
        import gzip

        for name, paths in (("bin_stats", "~{sep=',' bin_stats}"), ("traces", "~{sep=',' traces}")):
            with gzip.open(f"~{prefix}.{name}.tsv.gz", "wt") as out:
                for i, path in enumerate(paths.split(",")):
                    with gzip.open(path, "rt") as f:
                        header = f.readline()
                        if i == 0:
                            out.write(header)
                        for line in f:
                            out.write(line)
        CODE
        gzip -c ~{windows} > ~{prefix}.windows.tsv.gz
    >>>

    output {
        File combined_bin_stats = "~{prefix}.bin_stats.tsv.gz"
        File combined_traces = "~{prefix}.traces.tsv.gz"
        File combined_windows = "~{prefix}.windows.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bin_stats, "GB") + size(traces, "GB")) + 10,
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
