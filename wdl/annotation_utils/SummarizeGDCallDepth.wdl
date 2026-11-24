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
        File? genes_bed
        File? gene_info
        File? manual_labels
        File? keep_samples
        Boolean nahr_check = false
        Int outside_bins = 5
        Float max_dev = 0.5
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
                genes_bed = genes_bed,
                gene_info = gene_info,
                manual_labels = manual_labels,
                keep_samples = keep_samples,
                nahr_check = nahr_check,
                outside_bins = outside_bins,
                max_dev = max_dev,
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
        File? carrier_flags = PlotRegionDepth.carrier_flags
        File? carrier_labels = PlotRegionDepth.carrier_labels
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
        File? genes_bed
        File? gene_info
        File? manual_labels
        File? keep_samples
        Boolean nahr_check
        Int outside_bins
        Float max_dev
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
        # --nahr-check: carriers whose median CN in the non-SD bins just outside the flanking SDs has not returned to normal are
        # drawn in green as 'CNV, not NAHR' (and listed with --flags). --manual: carriers labelled non_NAHR by hand are drawn in orange
        # (taking precedence over green); regions labelled false get a note in the header; regions labelled blacklist are not plotted.
        # A manual label's sample may be an ID, * (all carriers of that region / svtype) or top<N>_<left|right|any>: the N carriers
        # whose CN just outside the left (right, either) flanking SD is shifted most in their CNV's direction (needs --nahr-check).
        # Region-wide (*) rows are applied before rows naming samples, so the specific label wins.
        # --carrier-labels: per-carrier table of the automatic and manual calls and the final NAHR / non_NAHR status.
        # --genes / --gene-info: gene track between the SD tracks and the depth panel (names grey; green LOEUF_v4 < 0.6; blue
        # pTriplo >= 0.94; red pHaplo >= 0.86). --review: category added to each page title.
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
                make_option("--counts", type = "character", default = NULL, help = "optional TSV of carrier counts per region"),
                make_option("--nahr-check", action = "store_true", default = FALSE,
                            help = paste("flag carriers whose depth does not return to normal just outside the flanking SDs as",
                                         "'CNV carrier, not NAHR' (green); needs --pairs-bed and --segdup")),
                make_option("--outside-bins", type = "integer", default = 5,
                            help = "non-SD bins next to each flanking SD (outer side) used by --nahr-check [%default]"),
                make_option("--max-dev", type = "double", default = 0.5,
                            help = "--nahr-check: carrier flagged if the median CN of those bins is shifted by more than this in the CNV's direction [%default]"),
                make_option("--flags", type = "character", default = NULL, help = "optional TSV of the per-carrier --nahr-check results"),
                make_option("--keep-samples", type = "character", default = NULL,
                            help = "optional file of sample IDs (one per line): only these samples are plotted and counted"),
                make_option("--genes", type = "character", default = NULL,
                            help = "gene BED (#chrom start end gene_name ...), .gz ok; adds a gene track between the SD tracks and the depth panel"),
                make_option("--gene-info", type = "character", default = NULL,
                            help = "gene metrics TSV (gene_name LOEUF_v4 pHaplo pTriplo): names in red if pHaplo >= --phaplo, blue if pTriplo >= --ptriplo, green if LOEUF_v4 < --loeuf, else grey"),
                make_option("--loeuf", type = "double", default = 0.6, help = "LOEUF_v4 cutoff for constrained genes [%default]"),
                make_option("--phaplo", type = "double", default = 0.86, help = "pHaplo cutoff for haploinsufficient genes [%default]"),
                make_option("--ptriplo", type = "double", default = 0.94, help = "pTriplo cutoff for triplosensitive genes [%default]"),
                make_option("--review", type = "character", default = NULL,
                            help = "optional review TSV (region, category_name, ...); the category is added to each page title"),
                make_option("--manual", type = "character", default = NULL,
                            help = paste("optional TSV of manual labels (region chrom:start-end, svtype DEL/DUP/*, sample or *, label, note);",
                                         "carriers labelled non_NAHR are drawn in orange; NAHR overrides the --nahr-check flag (no green);",
                                         "regions labelled false get a note in the header; blacklist regions are not plotted;",
                                         "sample may be top<N>_<left|right|any> (N carriers most shifted just outside that SD; needs --nahr-check)")),
                make_option("--carrier-labels", type = "character", default = NULL,
                            help = "optional TSV of every plotted carrier with its automatic and manual calls and final status (NAHR / non_NAHR)")
            )
        ))
        for (k in c("traces", "windows", "out")) if (is.null(opt[[k]])) stop("--", k, " is required (use -h for help)")

        font <- "Arial"
        sv_cols <- c(`non-carrier` = "grey70", DEL = "#C0392B", DUP = "#2166AC", `non-NAHR` = "#2CA02C", `manual non-NAHR` = "#FF7F00")

        win <- fread(opt$windows)
        win[, mid := (bin_start + bin_end) / 2]
        win[, region := paste0(chrom, ":", start, "-", end)]
        bins <- unique(win[, .(region, bin, mid, bin_start, bin_end)])
        regions <- unique(win[, .(region, chrom, start, end)])
        regions <- regions[order(factor(chrom, levels = paste0("chr", c(1:22, "X", "Y"))), start)]
        if (!is.null(opt$regions)) regions <- regions[region %in% readLines(opt$regions)]
        manual <- NULL
        if (!is.null(opt$manual)) {
            manual <- fread(opt$manual, sep = "\t", colClasses = "character")
            regions <- regions[!region %in% manual[label == "blacklist"]$region]
        }
        call_region <- unique(win[, .(call_id, region, svtype)])

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

        tr <- fread(opt$traces)
        tr <- merge(tr, call_region, by = "call_id")
        keep <- NULL
        if (!is.null(opt$`keep-samples`)) {
            keep <- unique(readLines(opt$`keep-samples`))
            tr <- tr[sample %in% keep]
        }
        ref_tr <- tr[class == "ref"]
        tr <- tr[class != "ref"]
        if (nrow(ref_tr) > 0) {
            ref_tr[, keep := if ("DEL" %in% svtype) "DEL" else "DUP", by = region]
            ref_tr <- ref_tr[svtype == keep][!unique(tr[, .(region, batch, sample)]), on = .(region, batch, sample)]
        }
        tr[, cn := as.integer(sub("^CN", "", class))]
        tr[, group := fifelse(svtype == "DEL", fifelse(cn == 0, "hom DEL", "het DEL"), fifelse(cn >= 4, "hom DUP", "het DUP"))]

        tr[, `:=`(manual = NA_character_, score_left = -Inf, score_right = -Inf)]

        # NAHR check: median CN of the --outside-bins non-SD bins just outside each flanking SD; a carrier whose CN there is shifted by
        # more than --max-dev in its CNV's direction (DUP up, DEL down) on either side is a CNV carrier but not NAHR between the SD pair
        tr[, nahr := NA]
        if (opt$`nahr-check`) {
            if (is.null(flank_sd) || is.null(segdup)) stop("--nahr-check needs --pairs-bed and --segdup")
            bcols0 <- grep("^b[0-9]+$", names(tr), value = TRUE)
            for (rg in intersect(regions$region, unique(tr$region))) {
                b <- bins[region == rg][order(bin)]
                r <- regions[region == rg]
                p <- pb[chrom == r$chrom & start == r$start & end == r$end]
                if (nrow(p) == 0) next
                o <- segdup[chrom == r$chrom]
                in_sd <- vapply(seq_len(nrow(b)), function(k) any(o$sd_start < b$bin_end[k] & o$sd_end > b$bin_start[k]), logical(1))
                lb <- tail(b[bin_end <= p$left_sd_start[1] & !in_sd]$bin, opt$`outside-bins`)
                rb <- head(b[bin_start >= p$right_sd_end[1] & !in_sd]$bin, opt$`outside-bins`)
                lb3 <- tail(b[bin_end <= p$left_sd_start[1]]$bin, 3)
                rb3 <- head(b[bin_start >= p$right_sd_end[1]]$bin, 3)
                idx <- which(tr$region == rg)
                X <- as.matrix(tr[idx, ..bcols0])
                med <- function(bs) if (length(bs) == 0) rep(NA_real_, length(idx)) else apply(X[, paste0("b", bs), drop = FALSE], 1, median, na.rm = TRUE)
                d <- sign(tr$cn[idx] - 2)
                devL <- (med(lb) - 2) * d
                devR <- (med(rb) - 2) * d
                # rank scores for top<N>_<side> labels: largest shift in the CNV direction over the outside bins or the 3 bins next to the SD
                scL <- pmax(fcoalesce(devL, -Inf), fcoalesce((med(lb3) - 2) * d, -Inf))
                scR <- pmax(fcoalesce(devR, -Inf), fcoalesce((med(rb3) - 2) * d, -Inf))
                tr[idx, `:=`(score_left = scL, score_right = scR)]
                tr[idx, `:=`(outside_left_cn = round(med(lb), 2), outside_right_cn = round(med(rb), 2), n_left_bins = length(lb),
                             n_right_bins = length(rb),
                             nahr = fifelse(is.na(devL) & is.na(devR), NA, !(fcoalesce(devL, 0) > opt$`max-dev` | fcoalesce(devR, 0) > opt$`max-dev`)))]
            }
            if (!is.null(opt$flags)) {
                fwrite(tr[region %in% regions$region, .(region, svtype, call_id, batch, sample, cn, outside_left_cn, outside_right_cn,
                                                         n_left_bins, n_right_bins, nahr)], opt$flags, sep = "\t")
            }
        }
        # manual labels: rows with sample * first, then rows naming samples (ID or top<N>_<left|right|any>); within each, file order
        # (later rows win), so a specific sample's label overrides a region-wide one
        if (!is.null(manual)) {
            manual <- manual[order(sample != "*")]
            for (k in seq_len(nrow(manual))) {
                mk <- manual[k]
                if (mk$label == "blacklist") next
                hit <- tr$region == mk$region & (mk$svtype == "*" | tr$svtype == mk$svtype)
                m <- regmatches(mk$sample, regexec("^top([0-9]+)_(left|right|any)$", mk$sample))[[1]]
                if (length(m) == 3) {
                    if (!opt$`nahr-check`) stop("manual top<N>_<side> labels need --nahr-check")
                    sc <- switch(m[3], left = tr$score_left, right = tr$score_right, any = pmax(tr$score_left, tr$score_right))
                    idx <- which(hit)
                    idx <- head(idx[order(-sc[idx])], as.integer(m[2]))
                    hit <- seq_len(nrow(tr)) %in% idx
                } else {
                    hit <- hit & (mk$sample == "*" | tr$sample == mk$sample)
                }
                tr[hit, manual := mk$label]
            }
        }
        if (!is.null(opt$`carrier-labels`)) {
            cl <- tr[region %in% regions$region]
            if (!"nahr" %in% names(cl)) cl[, nahr := NA]
            for (k in c("outside_left_cn", "outside_right_cn")) if (!k %in% names(cl)) cl[, (k) := NA_real_]
            cl[, final := fifelse(manual %in% "non_NAHR", "non_NAHR", fifelse(manual %in% "NAHR", "NAHR",
                          fifelse(nahr %in% FALSE, "non_NAHR", "NAHR")))]
            cl <- merge(cl, regions[, .(region, chrom, start, end)], by = "region")
            fwrite(cl[order(chrom, start, svtype, sample), .(region, chrom, start, end, svtype, call_id, batch, sample, cn, outside_left_cn,
                                                             outside_right_cn, auto_nahr = nahr, manual, final)], opt$`carrier-labels`, sep = "\t")
        }
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

        bin_cols <- grep("^b[0-9]+$", names(tr), value = TRUE)
        fmt <- function(x) format(x, big.mark = ",")

        dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
        genes <- NULL
        if (!is.null(opt$genes)) {
            genes <- fread(opt$genes, select = 1:4, col.names = c("chrom", "start", "end", "gene_name"))
            genes[, col := "grey45"]
            if (!is.null(opt$`gene-info`)) {
                gi <- unique(fread(opt$`gene-info`)[, .(gene_name, LOEUF_v4, pHaplo, pTriplo)], by = "gene_name")
                genes <- merge(genes, gi, by = "gene_name", all.x = TRUE)
                genes[, col := fifelse(pHaplo %in% pHaplo[pHaplo >= opt$phaplo], "#C0392B",
                               fifelse(pTriplo %in% pTriplo[pTriplo >= opt$ptriplo], "#2166AC",
                               fifelse(LOEUF_v4 %in% LOEUF_v4[LOEUF_v4 < opt$loeuf], "#2CA02C", "grey45")))]
            }
        }
        review <- if (!is.null(opt$review)) fread(opt$review) else NULL
        page_h <- if (is.null(genes)) 10 else 12.5
        # gene track: genes stacked in rows so labels do not overlap; label size shrinks for dense windows
        gene_track <- function(ch, xlim, rg) {
            g <- genes[chrom == ch & end > xlim[1] * 1e6 & start < xlim[2] * 1e6]
            tp <- ggplot() + annotate("rect", xmin = rg$start / 1e6, xmax = rg$end / 1e6, ymin = -Inf, ymax = Inf, fill = "#FBE9E7")
            if (nrow(g) == 0) {
                return(tp + annotate("text", x = mean(xlim), y = 0, label = "no protein-coding genes", colour = "grey45", size = 6, family = font) +
                           coord_cartesian(xlim = xlim, ylim = c(-1, 1), expand = FALSE) + labs(x = NULL, y = "Genes") +
                           theme_classic(base_size = 22, base_family = font) +
                           theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line.x = element_blank()))
            }
            g <- g[order(start)]
            g[, `:=`(s = pmax(start / 1e6, xlim[1]), e = pmin(end / 1e6, xlim[2]))]
            g[, mid := (s + e) / 2]
            for (sz in c(7.76, 6, 4.5, 3.5)) {
                half <- nchar(g$gene_name) * 0.55 * sz / 25.4 / 13 * diff(xlim) / 2
                lo <- pmin(g$s, g$mid - half); hi <- pmax(g$e, g$mid + half)
                row <- integer(nrow(g)); ends <- numeric(0)
                for (k in seq_len(nrow(g))) {
                    j <- which(ends < lo[k])[1]
                    if (is.na(j)) { ends <- c(ends, hi[k]); j <- length(ends) } else ends[j] <- hi[k]
                    row[k] <- j
                }
                if (max(row) <= 4 || sz == 3.5) break
            }
            g[, row := row]
            n_row <- max(g$row)
            tp + geom_segment(data = g, aes(x = s, xend = e, y = -row, yend = -row), colour = "grey35", linewidth = 2) +
                geom_text(data = g, aes(x = mid, y = -row + 0.42, label = gene_name, colour = col), size = sz, family = font,
                          fontface = "italic") +
                scale_colour_identity() +
                coord_cartesian(xlim = xlim, ylim = c(-n_row - 0.4, 0.2), expand = FALSE, clip = "off") +
                labs(x = NULL, y = "Genes") +
                theme_classic(base_size = 22, base_family = font) +
                theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line.x = element_blank())
        }
        if (capabilities("aqua")) quartz(type = "pdf", file = opt$out, width = 16, height = page_h) else
            cairo_pdf(opt$out, width = 16, height = page_h, onefile = TRUE)
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
            if (!is.null(manual)) t <- copy(t)[manual %in% "NAHR", nahr := TRUE]
            if (opt$`nahr-check`) {
                header <- paste0(header, "    not NAHR: ", fmt(t[nahr == FALSE, uniqueN(paste(batch, sample))]))
                t <- copy(t)[nahr == FALSE, svtype := "non-NAHR"]
            }
            if (!is.null(manual)) {
                n_man <- t[manual %in% "non_NAHR", uniqueN(paste(batch, sample))]
                if (n_man > 0) header <- paste0(header, "    manual non-NAHR: ", fmt(n_man))
                if (any(t$manual %in% "false")) header <- paste0(header, "    manual: calls look false")
                t <- copy(t)[manual %in% "non_NAHR", svtype := "manual non-NAHR"]
            }
            p <- ggplot() +
                annotate("rect", xmin = rg$start / 1e6, xmax = rg$end / 1e6, ymin = -Inf, ymax = Inf, fill = "#FBE9E7") +
                geom_hline(yintercept = c(1, 2, 3, 4, 5), colour = "grey80", linetype = "dashed", linewidth = 0.5)
            rb <- if (!is.null(ref_band)) merge(ref_band[region == rg$region], b[, .(bin, mid)], by = "bin") else NULL
            if (!is.null(keep)) {
                header <- paste0(header, "    non-carriers: ", fmt(length(keep) - uniqueN(t[, .(sample)])))
            } else if (!is.null(rb) && nrow(rb) > 0) header <- paste0(header, "    non-carriers: ", fmt(max(rb$n)))
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
                for (sv in intersect(c("non-carrier", "DEL", "DUP", "non-NAHR", "manual non-NAHR"), n_type$svtype)) {
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
                    geom_point(data = data.frame(x = xlim[1], y = -1, svtype = names(sv_cols)), aes(x = x, y = y, colour = svtype))
            }
            p <- p +
                scale_colour_manual(values = sv_cols, name = NULL, drop = FALSE,
                                    limits = c("non-carrier", "DEL", "DUP", if (opt$`nahr-check`) "non-NAHR", if (!is.null(manual)) "manual non-NAHR"),
                                    labels = c(if (nrow(rt) > 0) paste0("non-carriers (random ", fmt(uniqueN(rt[, .(batch, sample)])), ")")
                                               else "non-carriers (mean +/- 2 SD)", "DEL carriers", "DUP carriers",
                                               if (opt$`nahr-check`) "CNV, not NAHR", if (!is.null(manual)) "non-NAHR (manual)")) +
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
            if (!is.null(review)) {
                cat_name <- review[region == rg$region]$category_name
                if (length(cat_name) == 1) title <- paste0(title, "   [", gsub("_", " ", sub("^[0-9]_", "", cat_name)), "]")
            }
            panels <- if (is.null(genes)) trk / p + plot_layout(heights = c(1, 4)) else
                trk / gene_track(rg$chrom, xlim, rg) / p + plot_layout(heights = c(1, 1.3, 4))
            caption <- if (is.null(genes)) NULL else
                sprintf("Gene names: grey; green LOEUF_v4 < %s; blue pTriplo >= %s; red pHaplo >= %s", opt$loeuf, opt$ptriplo, opt$phaplo)
            print(panels +
                      plot_annotation(title = title, subtitle = header, caption = caption,
                                      theme = theme(plot.title = element_text(size = 22, face = "bold", family = font),
                                                    plot.subtitle = element_text(size = 20, family = font),
                                                    plot.caption = element_text(size = 16, family = font, colour = "grey30"))))
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
            ~{"--genes " + genes_bed} \
            ~{"--gene-info " + gene_info} \
            ~{"--manual " + manual_labels} \
            ~{"--keep-samples " + keep_samples} \
            ~{if nahr_check then "--nahr-check --flags ~{prefix}.carrier_flags.tsv" else ""} \
            --carrier-labels ~{prefix}.carrier_labels.tsv \
            --outside-bins ~{outside_bins} \
            --max-dev ~{max_dev} \
            --raster-px ~{raster_px} \
            --counts ~{prefix}.region_counts.tsv \
            -o ~{prefix}.region_depth.pdf

        bgzip ~{prefix}.region_counts.tsv
        bgzip ~{prefix}.carrier_labels.tsv
        if [ -f ~{prefix}.carrier_flags.tsv ]; then
            bgzip ~{prefix}.carrier_flags.tsv
        fi
    >>>

    output {
        File region_depth_pdf = "~{prefix}.region_depth.pdf"
        File region_counts = "~{prefix}.region_counts.tsv.gz"
        File? carrier_flags = "~{prefix}.carrier_flags.tsv.gz"
        File carrier_labels = "~{prefix}.carrier_labels.tsv.gz"
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
