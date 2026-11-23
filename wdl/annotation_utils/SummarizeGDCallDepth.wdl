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
        File? segdup_bed
        Int raster_px = 2000

        String sv_pipeline_docker
        String? plot_docker

        RuntimeAttr? runtime_attr_summarize_batch_depth
        RuntimeAttr? runtime_attr_combine_depth_summaries
        RuntimeAttr? runtime_attr_plot_region_depth
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

    if (defined(plot_docker)) {
        call PlotRegionDepth {
            input:
                traces = CombineDepthSummaries.combined_traces,
                windows = CombineDepthSummaries.combined_windows,
                bin_stats = CombineDepthSummaries.combined_bin_stats,
                sd_pairs_bed = sd_pairs_bed,
                segdup_bed = segdup_bed,
                raster_px = raster_px,
                prefix = prefix,
                docker = select_first([plot_docker]),
                runtime_attr_override = runtime_attr_plot_region_depth
        }
    }

    output {
        File depth_bin_stats = CombineDepthSummaries.combined_bin_stats
        File depth_traces = CombineDepthSummaries.combined_traces
        File call_windows = CombineDepthSummaries.combined_windows
        File? region_depth_pdf = PlotRegionDepth.region_depth_pdf
        File? region_counts = PlotRegionDepth.region_counts
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

task PlotRegionDepth {
    input {
        File traces
        File windows
        File bin_stats
        File? sd_pairs_bed
        File? segdup_bed
        Int raster_px
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat > plot_gd_region_carrier_depth.R <<'RSCRIPT'
        # Read-depth profiles of every carrier of the GD calls at each SD-pair region, from lr-annotation SummarizeGDCallDepth run with
        # traces_per_class large enough to keep all carriers. One page per region (the DEL and DUP records of a region together):
        # every carrier is a thin line (estimated copy number per bin, 2 = diploid), coloured DEL vs DUP; carriers are not averaged.
        # Non-carriers (samples at their expected copy number in both the DEL and the DUP record): if the traces include sampled
        # "ref" profiles (SummarizeGDCallDepth ref_traces > 0), they are drawn as grey thin lines (from the DEL record when present,
        # excluding the region's carriers); otherwise, with --bin-stats, as a grey band, mean +/- 2 SD per bin. --bin-stats also
        # gives the non-carrier count in the header.
        # Header: number of samples with het DEL (RD_CN 1), hom DEL (RD_CN 0), het DUP (RD_CN 3) and hom DUP (RD_CN >= 4).
        # Tracks above the profiles: the two SDs of the pair flanking the region (--pairs-bed) and all hg38 SDs in the window
        # (--segdup). The carrier lines are drawn to a PNG (ragg) and placed in the panel, so pages with many carriers stay small.
        # Multi-page PDF, Arial (Quartz on macOS; cairo_pdf elsewhere, where fontconfig maps Arial to Liberation Sans).
        suppressPackageStartupMessages({
            library(optparse)
            library(data.table)
            library(ggplot2)
            library(patchwork)
            library(ragg)
            library(png)
        })

        opt <- parse_args(OptionParser(
            usage = "%prog --traces X.traces.tsv.gz --windows X.windows.tsv.gz -o out.pdf [--pairs-bed P.bed.gz --segdup SD.bed.gz]",
            description = "One page per GD region with every DEL and DUP carrier's depth profile (SummarizeGDCallDepth outputs).",
            option_list = list(
                make_option("--traces", type = "character", help = "traces.tsv.gz (call_id batch sample class b0..bN) [required]"),
                make_option("--windows", type = "character", help = "windows.tsv.gz (call_id chrom start end svtype bin ...) [required]"),
                make_option("--bin-stats", type = "character", default = NULL,
                            help = "bin_stats.tsv.gz (call_id batch class bin n sum sumsq) for the non-carrier band"),
                make_option("--pairs-bed", type = "character", default = NULL,
                            help = "SD-pair BED (#chrom start end ... left_sd_start left_sd_end right_sd_start right_sd_end)"),
                make_option("--segdup", type = "character", default = NULL, help = "hg38 SD BED (chrom start end ...), .gz ok"),
                make_option("--regions", type = "character", default = NULL, help = "optional file of regions chrom:start-end to plot"),
                make_option("--raster-px", type = "integer", default = 2000, help = "width in pixels of the carrier-line image [%default]"),
                make_option(c("-o", "--out"), type = "character", help = "output multi-page PDF [required]"),
                make_option("--counts", type = "character", default = NULL, help = "optional TSV of carrier counts per region")
            )
        ))
        for (k in c("traces", "windows", "out")) if (is.null(opt[[k]])) stop("--", k, " is required (use -h for help)")

        font <- "Arial"
        sv_cols <- c(`non-carrier` = "grey70", DEL = "#C0392B", DUP = "#2166AC")

        win <- fread(opt$windows)
        win[, mid := (bin_start + bin_end) / 2]
        win[, region := paste0(chrom, ":", start, "-", end)]
        bins <- unique(win[, .(region, bin, mid, bin_start, bin_end)])
        regions <- unique(win[, .(region, chrom, start, end)])
        regions <- regions[order(factor(chrom, levels = paste0("chr", c(1:22, "X", "Y"))), start)]
        if (!is.null(opt$regions)) regions <- regions[region %in% readLines(opt$regions)]
        call_region <- unique(win[, .(call_id, region, svtype)])

        tr <- fread(opt$traces)
        tr <- merge(tr, call_region, by = "call_id")
        ref_tr <- tr[class == "ref"]
        tr <- tr[class != "ref"]
        if (nrow(ref_tr) > 0) {
            ref_tr[, keep := if ("DEL" %in% svtype) "DEL" else "DUP", by = region]
            ref_tr <- ref_tr[svtype == keep][!unique(tr[, .(region, batch, sample)]), on = .(region, batch, sample)]
        }
        tr[, cn := as.integer(sub("^CN", "", class))]
        tr[, group := fifelse(svtype == "DEL", fifelse(cn == 0, "hom DEL", "het DEL"), fifelse(cn >= 4, "hom DUP", "het DUP"))]
        counts <- tr[, .(n = uniqueN(sample)), by = .(region, group)]
        counts <- dcast(counts, region ~ group, value.var = "n", fill = 0)
        for (g in c("het DEL", "hom DEL", "het DUP", "hom DUP")) if (!g %in% names(counts)) counts[, (g) := 0L]
        if (!is.null(opt$counts)) {
            fwrite(merge(regions, counts, by = "region", all.x = TRUE)[, c("region", "het DEL", "hom DEL", "het DUP", "hom DUP")],
                   opt$counts, sep = "\t")
        }

        ref_band <- NULL
        if (!is.null(opt$`bin-stats`)) {
            st <- fread(opt$`bin-stats`)
            st <- merge(st[, .(n = sum(n), sum = sum(sum), sumsq = sum(sumsq)), by = .(call_id, class, bin)], call_region, by = "call_id")
            # DUP carriers are "ref" in the DEL record (and DEL carriers in the DUP record): remove the other type's carriers
            ref <- st[class == "ref", .(region, svtype, bin, n, sum, sumsq)]
            car <- st[class != "ref", .(n = sum(n), sum = sum(sum), sumsq = sum(sumsq)), by = .(region, svtype, bin)]
            car[, svtype := fifelse(svtype == "DEL", "DUP", "DEL")]
            ref <- merge(ref, car, by = c("region", "svtype", "bin"), all.x = TRUE, suffixes = c("", ".o"))
            for (k in c("n.o", "sum.o", "sumsq.o")) ref[is.na(get(k)), (k) := 0]
            ref[, `:=`(n = n - n.o, sum = sum - sum.o, sumsq = sumsq - sumsq.o)]
            ref <- ref[order(region, bin, svtype)][, .SD[1], by = .(region, bin)]
            ref[, mean := sum / n]
            ref[, sd := sqrt(pmax(sumsq / n - mean^2, 0))]
            ref_band <- ref[, .(region, bin, n, mean, sd)]
        }

        flank_sd <- NULL
        if (!is.null(opt$`pairs-bed`)) {
            pb <- fread(opt$`pairs-bed`)
            setnames(pb, 1:3, c("chrom", "start", "end"))
            flank_sd <- rbind(pb[, .(chrom, start, end, sd_start = left_sd_start, sd_end = left_sd_end)],
                              pb[, .(chrom, start, end, sd_start = right_sd_start, sd_end = right_sd_end)])
        }
        segdup <- NULL
        if (!is.null(opt$segdup)) {
            segdup <- fread(opt$segdup, header = FALSE, select = 1:3, col.names = c("chrom", "sd_start", "sd_end"))
        }

        bin_cols <- grep("^b[0-9]+$", names(tr), value = TRUE)
        fmt <- function(x) format(x, big.mark = ",")

        dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
        if (capabilities("aqua")) quartz(type = "pdf", file = opt$out, width = 16, height = 10) else
            cairo_pdf(opt$out, width = 16, height = 10, onefile = TRUE)
        for (i in seq_len(nrow(regions))) {
            rg <- regions[i]
            b <- bins[region == rg$region]
            xlim <- c(min(b$bin_start), max(b$bin_end)) / 1e6
            t <- tr[region == rg$region]
            rt <- ref_tr[region == rg$region]
            cnt <- counts[region == rg$region]
            get_n <- function(g) if (nrow(cnt) == 0) 0 else cnt[[g]]
            header <- paste0("het DEL: ", fmt(get_n("het DEL")), "    hom DEL: ", fmt(get_n("hom DEL")),
                             "    het DUP: ", fmt(get_n("het DUP")), "    hom DUP: ", fmt(get_n("hom DUP")))
            p <- ggplot() +
                annotate("rect", xmin = rg$start / 1e6, xmax = rg$end / 1e6, ymin = -Inf, ymax = Inf, fill = "#FBE9E7") +
                geom_hline(yintercept = c(1, 2, 3, 4, 5), colour = "grey80", linetype = "dashed", linewidth = 0.5)
            rb <- if (!is.null(ref_band)) merge(ref_band[region == rg$region], b[, .(bin, mid)], by = "bin") else NULL
            if (!is.null(rb) && nrow(rb) > 0) header <- paste0(header, "    non-carriers: ", fmt(max(rb$n)))
            if (nrow(t) > 0 || nrow(rt) > 0) {
                if (nrow(rt) > 0) rt[, svtype := "non-carrier"]
                m <- melt(rbind(rt, t, fill = TRUE)[, c("batch", "sample", "svtype", bin_cols), with = FALSE],
                          id.vars = c("batch", "sample", "svtype"), variable.name = "bin", value.name = "cn")
                m[, bin := as.integer(sub("^b", "", bin))]
                m <- merge(m, b[, .(bin, mid)], by = "bin")
                setorder(m, svtype, batch, sample, bin)
                n_type <- m[, .(n = uniqueN(paste(batch, sample))), by = svtype]
                lp <- ggplot() +
                    annotate("rect", xmin = rg$start / 1e6, xmax = rg$end / 1e6, ymin = -Inf, ymax = Inf, fill = "#FBE9E7") +
                    geom_hline(yintercept = c(1, 2, 3, 4, 5), colour = "grey80", linetype = "dashed", linewidth = 0.5)
                if (nrow(rt) == 0 && !is.null(rb) && nrow(rb) > 0) {
                    lp <- lp + geom_ribbon(data = rb, aes(x = mid / 1e6, ymin = mean - 2 * sd, ymax = mean + 2 * sd), fill = "grey70",
                                           alpha = 0.8)
                }
                for (sv in intersect(c("non-carrier", "DEL", "DUP"), n_type$svtype)) {
                    a <- min(0.6, max(if (sv == "non-carrier") 0.06 else 0.02, 15 / n_type[svtype == sv]$n))
                    lp <- lp + geom_line(data = m[svtype == sv], aes(x = mid / 1e6, y = cn, group = paste(batch, sample)),
                                         colour = sv_cols[[sv]], linewidth = 0.5, alpha = a)
                }
                lp <- lp + coord_cartesian(xlim = xlim, ylim = c(0, 6), expand = FALSE) + theme_void() +
                    theme(plot.margin = margin(0, 0, 0, 0))
                png_file <- tempfile(fileext = ".png")
                cur <- dev.cur()
                agg_png(png_file, width = opt$`raster-px`, height = round(opt$`raster-px` / 2.4), background = "white", res = 150)
                print(lp)
                dev.off()
                dev.set(cur)
                img <- readPNG(png_file)[, , 1:3]
                unlink(png_file)
                p <- p + annotation_raster(img, xmin = xlim[1], xmax = xlim[2], ymin = 0, ymax = 6) +
                    geom_point(data = data.frame(x = xlim[1], y = -1, svtype = c("non-carrier", "DEL", "DUP")), aes(x = x, y = y, colour = svtype))
            }
            p <- p +
                scale_colour_manual(values = sv_cols, name = NULL, drop = FALSE, limits = c("non-carrier", "DEL", "DUP"),
                                    labels = c(if (nrow(rt) > 0) paste0("non-carriers (random ", fmt(uniqueN(rt[, .(batch, sample)])), ")")
                                               else "non-carriers (mean +/- 2 SD)", "DEL carriers", "DUP carriers")) +
                guides(colour = guide_legend(override.aes = list(shape = 15, size = 6))) +
                coord_cartesian(xlim = xlim, ylim = c(0, 6), expand = FALSE) +
                scale_y_continuous(breaks = 0:6) +
                labs(x = paste0(rg$chrom, " position (Mb)"), y = "Estimated copy number") +
                theme_classic(base_size = 22, base_family = font) +
                theme(legend.position = "top")
            fs <- if (!is.null(flank_sd)) flank_sd[chrom == rg$chrom & start == rg$start & end == rg$end] else NULL
            sd <- if (!is.null(segdup)) segdup[chrom == rg$chrom & sd_end > xlim[1] * 1e6 & sd_start < xlim[2] * 1e6] else NULL
            trk <- ggplot() +
                annotate("rect", xmin = rg$start / 1e6, xmax = rg$end / 1e6, ymin = -Inf, ymax = Inf, fill = "#FBE9E7") +
                scale_y_continuous(breaks = c(2, 1), labels = c("Flanking SDs", "hg38 SDs"), limits = c(0.5, 2.5)) +
                coord_cartesian(xlim = xlim, expand = FALSE) +
                labs(x = NULL, y = NULL) +
                theme_classic(base_size = 22, base_family = font) +
                theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
            if (!is.null(fs) && nrow(fs) > 0) {
                trk <- trk + geom_rect(data = fs, aes(xmin = sd_start / 1e6, xmax = sd_end / 1e6, ymin = 1.65, ymax = 2.35),
                                       fill = "#E08214", colour = "#8C510A", linewidth = 0.4)
            }
            if (!is.null(sd) && nrow(sd) > 0) {
                trk <- trk + geom_rect(data = sd, aes(xmin = pmax(sd_start / 1e6, xlim[1]), xmax = pmin(sd_end / 1e6, xlim[2]),
                                                      ymin = 0.65, ymax = 1.35), fill = "grey35", alpha = 0.25, colour = NA)
            }
            title <- paste0(rg$chrom, ":", fmt(rg$start), "-", fmt(rg$end), "  (", format(round((rg$end - rg$start) / 1000, 1), nsmall = 1),
                            " kb)")
            print(trk / p + plot_layout(heights = c(1, 4)) +
                      plot_annotation(title = title, subtitle = header,
                                      theme = theme(plot.title = element_text(size = 22, face = "bold", family = font),
                                                    plot.subtitle = element_text(size = 20, family = font))))
        }
        invisible(dev.off())
        cat("wrote", opt$out, "(", nrow(regions), "regions )\n")
        RSCRIPT

        Rscript plot_gd_region_carrier_depth.R \
            --traces ~{traces} \
            --windows ~{windows} \
            --bin-stats ~{bin_stats} \
            ~{"--pairs-bed " + sd_pairs_bed} \
            --segdup ~{default="/opt/gd_scripts/references/hg38.SegDup.sorted.merged.bed.gz" segdup_bed} \
            --raster-px ~{raster_px} \
            --counts ~{prefix}.region_counts.tsv \
            -o ~{prefix}.region_depth.pdf
    >>>

    output {
        File region_depth_pdf = "~{prefix}.region_depth.pdf"
        File region_counts = "~{prefix}.region_counts.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 32,
        disk_gb: 10 * ceil(size(traces, "GB") + size(bin_stats, "GB")) + 20,
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
