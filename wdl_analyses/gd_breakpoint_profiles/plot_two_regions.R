#!/usr/bin/env Rscript
# Two regions side by side (e.g. the two breakpoint SDs of a genomic disorder: left = proximal, right = distal), each as a
# column of SD tracks, copy number (from depth), folded BAF and LD, one line per sample; samples whose ID matches
# --highlight-pattern are coloured (same colour in both columns), the rest grey; black = cohort median.
#
# Per region (1 = left, 2 = right):
#   --baf1/--baf2   BAF table (CHROM POS ID REF ALT + samples), best the clean-site table of find_baf_terminal_loh.py
#   --rd1/--rd2     depth matrix (#Chr Start End + samples), optional; a missing one leaves an empty depth panel
#   --ld1/--ld2     compute_ld.py output prefix (uses <prefix>.binned.tsv), optional
#   --title1/2, --mark1/2 (dashed lines, Mb, comma-separated)
# Shared: --medians (2 rows: sample IDs / genome-wide median depth; needed with any --rd), --segdup (NAME=bed[,...]),
#   --highlight-pattern, --bin-size (CN), --baf-bin, --ld-stat (mean_r2|max_r2), --na-missing, --title.
#   --bin-size / --baf-bin 0 (default) = automatic per region: CN bin = window / 150, BAF bin = window / 40 (rounded to
#   1 kb, at least 1 kb / 2 kb), so 30 kb and 5 Mb windows both plot with a sensible number of points; the bin sizes
#   are shown in each column title. --win1/--win2 "start,end" (bp) fix the plotted window (else taken from the data),
#   so a region without usable sites still gets its panels.
# Copy number = RD / sample median x 2, divided by the cohort median of each bin x 2. Folded BAF = min(BAF, 1 - BAF),
# NA counted as 0 unless --na-missing. Both columns use the same y scales.
# Colours: highlighted samples matching --del-pattern get red-to-orange tones, --dup-pattern blue-to-green tones, others
#   a mixed palette.
# --het-baf adds a het-only panel under the BAF panel: folded BAF of heterozygous sites only (--het-lo < BAF < --het-hi),
#   each het site of a highlighted sample as a point in its colour, the cohort (non-highlighted samples) median and
#   5-95% band per --baf-bin window in black/grey; guides at 1/3 (one extra copy) and 1/2. Deletions leave no hets here.
# --overview-rd adds a full-width panel on top: copy number of every sample across the whole locus (depth matrix spanning
#   both breakpoint SDs and the GD interval; needs --medians), with the SD and SD breakpoint tracks above it, the two
#   zoomed windows shaded and dashed lines at --mark1/--mark2. --overview-window "start,end" (bp) fixes its range,
#   --overview-bin its bin size (0 = automatic, window / 300).
# --page-size N: one page per N DEL + N DUP carriers (other carriers grey on that page); 0 = all on one page.
# --order-by-edges: carriers are ordered (and so grouped into pages) by their copy-number edges: per carrier, the overview
#   copy number is smoothed with a running median over ~25 kb, and the longest run below 1.5 (DEL) or above 2.5 (DUP)
#   gives its left and right edge (DEL carriers: BAF breakpoint estimates when available). Carriers of each type are then
#   clustered (single linkage, both edges within --cluster-tol bp) and whole clusters are packed onto pages (a cluster is
#   only split when it is larger than --page-size). --edges-out writes the edges and cluster IDs.
# Breakpoint arrows: --bp1/--bp2 = ranked.tsv of find_baf_terminal_loh.py (columns sample and --bp-col, default bp_est,
#   in bp); an arrow in the sample's colour marks it in the depth and BAF panels, for highlighted samples matching
#   --bp-pattern (default: the --del-pattern; BAF loss only places deletion breakpoints).
#
# Usage:
#   Rscript plot_two_regions.R --baf1 aou_II/chr16_bp4.renamed.right_loh.clean_sites.baf.tsv.gz \
#     --rd1 aou_II/chr16_bp4.renamed.RD.txt.gz --ld1 aou_II/chr16_bp4.renamed.ld --title1 BP4 \
#     --baf2 aou_II/chr16_bp5.renamed.left_loh.clean_sites.baf.tsv.gz --ld2 aou_II/chr16_bp5.renamed.ld --title2 BP5 \
#     --medians aou_II/GD_samples.RD_median --highlight-pattern "GD_16p11.2-BP4-5_DEL" \
#     --segdup "SD=references/hg38.SegDup.sorted.merged.bed.gz,SD breakpoint=references/GD_NAHR.breakpoint_SDs.bed.gz" \
#     --out aou_II/figures/chr16_bp4_bp5.BP4-5_DEL.pdf --title "16p11.2 BP4-5 DEL carriers"
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(optparse); library(patchwork)})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--baf1"), make_option("--rd1"), make_option("--ld1"), make_option("--title1", default = "left"),
  make_option("--mark1"),
  make_option("--baf2"), make_option("--rd2"), make_option("--ld2"), make_option("--title2", default = "right"),
  make_option("--mark2"),
  make_option("--medians"), make_option("--segdup"), make_option("--out"), make_option("--title", default = ""),
  make_option("--highlight-pattern", default = "GD_16p11.2-BP4-5"),
  make_option("--bin-size", type = "integer", default = 0), make_option("--baf-bin", type = "integer", default = 0),
  make_option("--win1"), make_option("--win2"),
  make_option("--ld-stat", default = "mean_r2"), make_option("--na-missing", action = "store_true", default = FALSE),
  make_option("--width", type = "double", default = 30), make_option("--height", type = "double", default = 22),
  make_option("--del-pattern", default = "_DEL"), make_option("--dup-pattern", default = "_DUP"),
  make_option("--bp1"), make_option("--bp2"), make_option("--bp-col", default = "bp_est"),
  make_option("--bp-pattern", help = "regex of samples that get arrows (default --del-pattern)"),
  make_option("--het-baf", action = "store_true", default = FALSE, help = "add a het-only BAF panel"),
  make_option("--het-lo", type = "double", default = 0.1), make_option("--het-hi", type = "double", default = 0.9),
  make_option("--overview-rd"), make_option("--overview-window"),
  make_option("--overview-bin", type = "integer", default = 0),
  make_option("--order-by-edges", action = "store_true", default = FALSE,
              help = "order carriers by their copy-number edges from the overview depth (needs --overview-rd)"),
  make_option("--edges-out", help = "TSV of the per-carrier copy-number edges (with --order-by-edges)"),
  make_option("--cluster-tol", type = "double", default = 30000,
              help = "with --order-by-edges: carriers whose left and right edges are both within this (bp) are clustered"),
  make_option("--page-size", type = "integer", default = 0,
              help = "with > 0: one page per group of this many DEL plus this many DUP carriers (multi-page PDF)"),
  make_option("--max-legend", type = "integer", default = 16,
              help = "list coloured samples in the legend only up to this many (title gives the DEL/DUP counts)"))))
for (x in c("baf1", "baf2", "out")) if (is.null(opt[[x]])) stop("need --baf1, --baf2 and --out (see -h)")
if ((!is.null(opt$rd1) || !is.null(opt$rd2)) && is.null(opt$medians)) stop("--rd1/--rd2 need --medians")
font_family <- if (capabilities("aqua") || capabilities("cairo")) "Arial" else "sans"
fixed <- c("CHROM", "POS", "ID", "REF", "ALT")

read_sd <- function(spec, region_chrom, xr) {
  if (is.null(spec) || !nzchar(spec)) return(NULL)
  parts <- strsplit(strsplit(spec, ",")[[1]], "=", fixed = TRUE)
  rbindlist(lapply(parts, function(p) {
    nm <- if (length(p) > 1) p[1] else basename(p[1])
    b <- fread(p[length(p)], header = "auto")
    setnames(b, 1:3, c("chrom", "start", "end"))
    lab <- intersect(c("breakpoint", "name"), names(b))
    b[, label := if (length(lab)) as.character(get(lab[1])) else ""]
    b <- b[chrom == region_chrom & end / 1e6 > xr[1] & start / 1e6 < xr[2]]
    b[, .(track = nm, start = pmax(start / 1e6, xr[1]), end = pmin(end / 1e6, xr[2]), label)]
  }))
}

med <- if (!is.null(opt$medians)) {
  m <- fread(opt$medians, header = FALSE); setNames(as.numeric(unlist(m[2])), unlist(m[1]))
} else NULL

auto_bin <- function(span, div, lo) max(lo, round(span / div / 1000) * 1000)
# Add NA rows for empty bins between the first and last bin, so lines break across gaps instead of bridging them
fill_bins <- function(d, b, v) {
  if (!nrow(d)) return(d)
  grid <- CJ(sample = unique(d$sample), bin = seq(min(d$bin), max(d$bin), by = b))
  merge(grid, d[, c("sample", "bin", v), with = FALSE], by = c("sample", "bin"), all.x = TRUE)
}
load_region <- function(baf, rd, ld, title, win) {
  b <- fread(baf, na.strings = c("NA", "."))
  if (!"CHROM" %in% names(b)) stop(baf, ": no CHROM column")
  b[, CHROM := as.character(CHROM)]
  samples <- setdiff(names(b), fixed)
  r <- NULL
  if (!is.null(rd)) {
    r <- fread(rd)
    setnames(r, 1:3, c("chrom", "start", "end"))
  }
  win <- if (!is.null(win)) as.numeric(strsplit(win, ",")[[1]]) else
    range(c(b$POS, if (!is.null(r)) c(r$start, r$end)))
  chrom <- if (nrow(b)) b$CHROM[1] else if (!is.null(r) && nrow(r)) r$chrom[1] else NA
  span <- diff(win)
  cnb <- if (opt$`bin-size` > 0) opt$`bin-size` else auto_bin(span, 150, 1000)
  bfb <- if (opt$`baf-bin` > 0) opt$`baf-bin` else auto_bin(span, 40, 2000)
  bl <- melt(b[, c("POS", samples), with = FALSE], id.vars = "POS", variable.name = "sample", value.name = "BAF",
             variable.factor = FALSE)
  bl[, mBAF := pmin(BAF, 1 - BAF)]
  if (!opt$`na-missing`) bl[is.na(mBAF), mBAF := 0]
  bl[, bin := floor(POS / bfb) * bfb + bfb / 2]
  mb <- bl[, .(mBAF = mean(mBAF, na.rm = TRUE), n = sum(!is.na(mBAF))), by = .(sample, bin)][n > 0]
  hs <- bl[!is.na(BAF) & BAF > opt$`het-lo` & BAF < opt$`het-hi`, .(sample, POS, bin, hBAF = pmin(BAF, 1 - BAF))]
  cn <- NULL
  if (!is.null(r) && nrow(r)) {
    if (nrow(b) && r$chrom[1] != chrom) warning(title, ": depth chromosome ", r$chrom[1], " differs from BAF ", chrom)
    rs <- intersect(names(r)[-(1:3)], names(med))
    if (!length(rs)) stop("no sample of the depth matrix has a median in --medians (check that both use the same, renamed, sample IDs)")
    r[, bin := floor(start / cnb) * cnb + cnb / 2]
    rl <- melt(r[, c("bin", rs), with = FALSE], id.vars = "bin", variable.name = "sample", value.name = "RD",
               variable.factor = FALSE)
    cn <- rl[, .(RD = sum(RD), n = .N), by = .(sample, bin)]
    cn[, CN := RD / (med[sample] * n) * 2]
    cn[, CN := CN / median(CN) * 2, by = bin]
    cn <- fill_bins(cn[is.finite(CN)], cnb, "CN")   # bins with zero cohort depth (gaps) become breaks
  }
  poly <- NULL; w <- NA
  if (!is.null(ld) && file.exists(paste0(ld, ".binned.tsv")) && nrow(bn <- fread(paste0(ld, ".binned.tsv")))) {
    w <- min(diff(sort(unique(c(bn$bin1, bn$bin2)))))
    bn[, v := get(opt$`ld-stat`)]
    poly <- bn[, .(x = c(bin1 + bin2, bin1 + bin2 + w, bin1 + bin2 + 2 * w, bin1 + bin2 + w) / 2,
                   y = c(bin2 - bin1, bin2 - bin1 + w, bin2 - bin1, bin2 - bin1 - w) / 2),
               by = .(id = seq_len(nrow(bn)), v)]
  }
  xr <- win / 1e6
  samples <- unique(c(samples, if (!is.null(cn)) cn$sample))
  title <- sprintf("%s\n(CN %g kb, BAF %g kb%s bins)", title, cnb / 1000, bfb / 1000,
                   if (is.na(w)) "" else sprintf(", LD %g kb", w / 1000))
  list(chrom = chrom, samples = samples, mb = mb, hs = hs, cn = cn, poly = poly, w = w, xr = xr, title = title,
       sd = read_sd(opt$segdup, chrom, xr), cnb = cnb, bfb = bfb)
}
R <- list(load_region(opt$baf1, opt$rd1, opt$ld1, opt$title1, opt$win1),
          load_region(opt$baf2, opt$rd2, opt$ld2, opt$title2, opt$win2))
marks <- lapply(list(opt$mark1, opt$mark2), function(m) if (is.null(m)) numeric(0) else as.numeric(strsplit(m, ",")[[1]]))

# Whole-locus copy number for the overview panel
ov <- NULL
if (!is.null(opt$`overview-rd`)) {
  if (is.null(med)) stop("--overview-rd needs --medians")
  r <- fread(opt$`overview-rd`)
  setnames(r, 1:3, c("chrom", "start", "end"))
  ow <- if (!is.null(opt$`overview-window`)) as.numeric(strsplit(opt$`overview-window`, ",")[[1]]) else
    range(c(r$start, r$end))
  ob <- if (opt$`overview-bin` > 0) opt$`overview-bin` else auto_bin(diff(ow), 300, 1000)
  rs <- intersect(names(r)[-(1:3)], names(med))
  if (!length(rs)) stop("no sample of the depth matrix has a median in --medians (check that both use the same, renamed, sample IDs)")
  r[, bin := floor(start / ob) * ob + ob / 2]
  rl <- melt(r[, c("bin", rs), with = FALSE], id.vars = "bin", variable.name = "sample", value.name = "RD",
             variable.factor = FALSE)
  ocn <- rl[, .(RD = sum(RD), n = .N), by = .(sample, bin)]
  ocn[, CN := RD / (med[sample] * n) * 2]
  ocn[, CN := CN / median(CN) * 2, by = bin]
  ocn <- fill_bins(ocn[is.finite(CN)], ob, "CN")
  ov <- list(cn = ocn, xr = ow / 1e6, chrom = r$chrom[1], bin = ob, samples = rs)
}

hl <- sort(grep(opt$`highlight-pattern`, unique(c(R[[1]]$samples, R[[2]]$samples, ov$samples)), value = TRUE))
cat("highlighted samples (", length(hl), "):\n", paste(hl, collapse = "\n"), "\n", sep = "")
tones <- function(lo, hi, n) if (n == 1) hi else grDevices::colorRampPalette(c(lo, hi))(n)
is_del <- grepl(opt$`del-pattern`, hl); is_dup <- grepl(opt$`dup-pattern`, hl) & !is_del
pal <- setNames(rep(c("#1B5E20", "#FF8C00", "#6A3D9A", "#00897B", "#8D6E63", "#C0CA33"), length.out = length(hl)), hl)
if (any(is_del)) pal[is_del] <- tones("#7F0000", "#FF8C00", sum(is_del))   # dark red -> orange
if (any(is_dup)) pal[is_dup] <- tones("#08306B", "#2E9E44", sum(is_dup))   # dark blue -> green
show_leg <- length(hl) <= opt$`max-legend`
bp_pat <- if (is.null(opt$`bp-pattern`)) opt$`del-pattern` else opt$`bp-pattern`
read_bp <- function(f) {
  if (is.null(f)) return(NULL)
  b <- fread(f)[sample %in% hl & grepl(bp_pat, sample)]
  b <- b[!is.na(get(opt$`bp-col`)), .(sample, x = get(opt$`bp-col`) / 1e6)]
  b[, k := match(sample, hl)]
  cat(basename(f), "breakpoints (Mb):", paste(sprintf("%s %.3f", sub(".*_", "", b$sample), b$x), collapse = "; "), "\n")
  b
}
bps <- list(read_bp(opt$bp1), read_bp(opt$bp2))
cn_all <- rbindlist(lapply(R, `[[`, "cn"))
cn_max <- if (nrow(cn_all)) max(4, quantile(cn_all$CN, 0.999, na.rm = TRUE)) else 4
theme_set(theme_bw(base_size = 22, base_family = font_family) +
          theme(panel.grid.minor = element_blank(), strip.background = element_blank()))

band <- function(d, v) d[, .(med = median(get(v), na.rm = TRUE)), by = bin]
lines_panel <- function(d, v, ylim, xr, ylab, mk, extra = NULL) {
  coh <- band(d, v)
  ggplot() + extra +
    geom_line(data = d[!sample %in% hl], aes(bin / 1e6, .data[[v]], group = sample), colour = "grey70", alpha = 0.25,
              linewidth = 0.4) +
    geom_line(data = coh, aes(bin / 1e6, med), colour = "black", linewidth = 1.4) +
    geom_line(data = d[sample %in% hl], aes(bin / 1e6, .data[[v]], colour = sample), linewidth = 0.9) +
    { if (length(mk)) geom_vline(xintercept = mk, linetype = "dashed", colour = "#1F5FA8") } +
    scale_colour_manual(values = pal, name = NULL, drop = FALSE, guide = if (show_leg) "legend" else "none") +
    coord_cartesian(xlim = xr, ylim = ylim) + labs(x = NULL, y = ylab)
}
arrows_layer <- function(bp, ylim) {
  if (is.null(bp) || !nrow(bp)) return(NULL)
  top <- ylim[2]; step <- diff(ylim) * 0.045
  b <- copy(bp)[, y0 := top - (rank(x, ties.method = "first") - 1) %% 4 * step]
  geom_segment(data = b, aes(x = x, xend = x, y = y0, yend = y0 - step * 1.6, colour = sample), linewidth = 1.3,
               arrow = arrow(length = unit(0.18, "inches"), type = "closed"), inherit.aes = FALSE, show.legend = FALSE)
}
empty_panel <- function(xr, ylim, ylab, msg) {
  ggplot() + annotate("text", x = mean(xr), y = mean(ylim), label = msg, size = 7.76, family = font_family,
                      colour = "grey40") + coord_cartesian(xlim = xr, ylim = ylim) + labs(x = NULL, y = ylab)
}
column <- function(reg, mk, first, bp) {
  ylab <- function(s) if (first) s else NULL
  sd <- if (!is.null(reg$sd) && nrow(reg$sd)) {
    s <- copy(reg$sd); s[, track := factor(track, levels = rev(unique(track)))]; s[label %in% c("NA", NA), label := ""]
    ggplot(s) + geom_rect(aes(xmin = start, xmax = end, ymin = as.numeric(track) - 0.4, ymax = as.numeric(track) + 0.4),
                          fill = "#E8A33D", colour = "grey30", linewidth = 0.2) +
      geom_text(aes(x = (start + end) / 2, y = as.numeric(track), label = label), size = 7.76 * 0.7,
                family = font_family, fontface = "bold") +
      { if (length(mk)) geom_vline(xintercept = mk, linetype = "dashed", colour = "#1F5FA8") } +
      scale_y_continuous(breaks = seq_along(levels(s$track)), labels = if (first) levels(s$track) else NULL,
                         expand = c(0, 0.3)) +
      coord_cartesian(xlim = reg$xr) + labs(x = NULL, y = NULL, title = reg$title) +
      theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"))
  } else ggplot() + labs(title = reg$title) + theme_void()
  p_cn <- if (!is.null(reg$cn)) lines_panel(reg$cn, "CN", c(0, cn_max), reg$xr, ylab("copy number"), mk,
                                           geom_hline(yintercept = c(1, 2, 3), linetype = "dotted", colour = "grey40")) +
    arrows_layer(bp, c(0, cn_max)) else
    empty_panel(reg$xr, c(0, cn_max), ylab("copy number"), "depth not available")
  p_b <- lines_panel(reg$mb, "mBAF", c(0, 0.5), reg$xr, ylab("folded BAF"), mk) +
    arrows_layer(bp, c(0, 0.5))
  p_ld <- if (!is.null(reg$poly)) {
    ggplot(reg$poly, aes(x / 1e6, -y / 1e6, group = id, fill = v)) + geom_polygon(colour = NA) +
      { if (length(mk)) geom_vline(xintercept = mk, linetype = "dashed", colour = "#1F5FA8") } +
      scale_fill_gradientn(colours = c("white", "#FEE0D2", "#FC9272", "#DE2D26", "#67000D"), limits = c(0, 1),
                           name = if (opt$`ld-stat` == "mean_r2") expression("mean r"^2) else expression("max r"^2)) +
      coord_cartesian(xlim = reg$xr) + labs(x = paste(reg$chrom, "position (Mb)"), y = ylab("LD")) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())
  } else empty_panel(reg$xr, c(0, 1), ylab("LD"), "LD not available") + labs(x = paste(reg$chrom, "position (Mb)"))
  p_h <- NULL
  if (opt$`het-baf`) {
    coh <- reg$hs[!sample %in% hl, .(med = median(hBAF), lo = quantile(hBAF, 0.05), hi = quantile(hBAF, 0.95)), by = bin]
    p_h <- ggplot() + geom_ribbon(data = coh, aes(bin / 1e6, ymin = lo, ymax = hi), fill = "grey85") +
      geom_line(data = coh, aes(bin / 1e6, med), colour = "black", linewidth = 1.2) +
      geom_hline(yintercept = 1 / 3, linetype = "dotted", colour = "grey40") +
      geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey40") +
      geom_point(data = reg$hs[sample %in% hl], aes(POS / 1e6, hBAF, colour = sample), size = 1.8, alpha = 0.8,
                 show.legend = FALSE) +
      { if (length(mk)) geom_vline(xintercept = mk, linetype = "dashed", colour = "#1F5FA8") } +
      scale_colour_manual(values = pal, name = NULL, drop = FALSE, guide = if (show_leg) "legend" else "none") +
      coord_cartesian(xlim = reg$xr, ylim = c(0, 0.5)) + labs(x = NULL, y = ylab("het-only folded BAF"))
  }
  Filter(Negate(is.null), list(sd, p_cn, p_b, p_h, p_ld))
}
make_fig <- function(page_label) {
  c1 <- column(R[[1]], marks[[1]], TRUE, bps[[1]])
  c2 <- column(R[[2]], marks[[2]], FALSE, bps[[2]])
  grid_hts <- c(0.18, 1.2, 1, if (opt$`het-baf`) 1, 1.4)
  fig <- wrap_plots(c(rbind(c1, c2)), ncol = 2, byrow = TRUE) + plot_layout(heights = grid_hts, guides = "collect")
  if (!is.null(ov)) {
    mk <- c(marks[[1]], marks[[2]])
    zoom <- data.table(xmin = c(R[[1]]$xr[1], R[[2]]$xr[1]), xmax = c(R[[1]]$xr[2], R[[2]]$xr[2]))
    osd <- read_sd(opt$segdup, ov$chrom, ov$xr)
    p_osd <- if (!is.null(osd) && nrow(osd)) {
      osd[, track := factor(track, levels = rev(unique(track)))]; osd[label %in% c("NA", NA), label := ""]
      ggplot(osd) + geom_rect(aes(xmin = start, xmax = end, ymin = as.numeric(track) - 0.4, ymax = as.numeric(track) + 0.4),
                              fill = "#E8A33D", colour = "grey30", linewidth = 0.2) +
        geom_text(aes(x = (start + end) / 2, y = as.numeric(track), label = label), size = 7.76 * 0.7,
                  family = font_family, fontface = "bold") +
        { if (length(mk)) geom_vline(xintercept = mk, linetype = "dashed", colour = "#1F5FA8") } +
        scale_y_continuous(breaks = seq_along(levels(osd$track)), labels = levels(osd$track), expand = c(0, 0.3)) +
        coord_cartesian(xlim = ov$xr) +
        labs(x = NULL, y = NULL, title = sprintf("whole locus (CN %g kb bins; shaded = zoomed SD windows below)", ov$bin / 1000)) +
        theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
              plot.title = element_text(hjust = 0.5, face = "bold"))
    } else ggplot() + theme_void()
    ocn_max <- max(4, quantile(ov$cn$CN, 0.999, na.rm = TRUE))
    p_ocn <- lines_panel(ov$cn, "CN", c(0, ocn_max), ov$xr, "copy number", mk,
                         list(geom_rect(data = zoom, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                                        fill = "#9ECAE1", alpha = 0.25, inherit.aes = FALSE),
                              geom_hline(yintercept = c(1, 2, 3), linetype = "dotted", colour = "grey40"))) +
      labs(x = paste(ov$chrom, "position (Mb)")) + guides(colour = "none")
    fig <- wrap_plots(p_osd, p_ocn, fig, ncol = 1, heights = c(0.18, 1.3, sum(grid_hts))) +
      plot_layout(guides = "collect")
  }
  fig <- fig +
    plot_annotation(title = sprintf("%s%s\ncoloured: %d DEL in red to orange, %d DUP in blue to green; arrows = BAF breakpoint estimates of DEL",
                                    opt$title, page_label, sum(is_del), sum(is_dup)),
                    theme = theme(plot.title = element_text(size = 22, family = font_family)))
  fig
}

# Pages: all carriers on one page, or --page-size DEL + --page-size DUP carriers per page (colours re-spread per page)
nat <- function(x) x[order(sub("_[0-9]+$", "", x), suppressWarnings(as.numeric(sub(".*_", "", x))))]
hl_all <- hl; bps_all <- bps
dels <- nat(hl_all[grepl(opt$`del-pattern`, hl_all)]); dups <- nat(hl_all[grepl(opt$`dup-pattern`, hl_all) & !hl_all %in% dels])
edges <- NULL
if (opt$`order-by-edges` && !is.null(ov)) {
  # Depth edges: smoothed copy number (running median over ~25 kb); bins beyond the threshold (< 1.5 DEL, > 2.5 DUP)
  # form runs, gaps shorter than ~20 kb between runs are closed, and the run covering the GD midpoint (mean of
  # --mark1/--mark2; else the longest run) gives the carrier's left and right edge
  k_bins <- max(3, round(25000 / ov$bin)) %/% 2 * 2 + 1
  gap_bins <- max(1, round(20000 / ov$bin))
  gd_mid <- if (length(marks[[1]]) && length(marks[[2]])) mean(c(marks[[1]][1], marks[[2]][1])) * 1e6 else NA
  e <- ov$cn[sample %in% c(dels, dups) & !is.na(CN)][order(sample, bin)]
  edges <- e[, {
    sm <- if (.N >= k_bins) stats::runmed(CN, k_bins, endrule = "median") else CN
    flag <- if (sample[1] %in% dels) sm < 1.5 else sm > 2.5
    r <- rle(flag)
    short_gap <- !r$values & r$lengths <= gap_bins & seq_along(r$values) > 1 & seq_along(r$values) < length(r$values)
    r$values[short_gap] <- TRUE
    flag <- inverse.rle(r); r <- rle(flag); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
    runs <- which(r$values)
    hit <- if (!is.na(gd_mid)) runs[bin[starts[runs]] <= gd_mid & bin[ends[runs]] >= gd_mid] else integer(0)
    best <- if (length(hit)) hit[1] else if (length(runs)) runs[which.max(r$lengths[runs])] else integer(0)
    if (length(best)) .(depth_left = bin[starts[best]], depth_right = bin[ends[best]]) else
      .(depth_left = NA_real_, depth_right = NA_real_)
  }, by = sample]
  edges[, type := ifelse(sample %in% dels, "DEL", "DUP")]
  # BAF breakpoint estimates (DEL carriers): left SD and right SD
  bl <- if (!is.null(bps_all[[1]])) bps_all[[1]][, .(sample, baf_left = x * 1e6)] else data.table(sample = character(0), baf_left = numeric(0))
  br <- if (!is.null(bps_all[[2]])) bps_all[[2]][, .(sample, baf_right = x * 1e6)] else data.table(sample = character(0), baf_right = numeric(0))
  edges <- merge(merge(edges, bl, by = "sample", all.x = TRUE), br, by = "sample", all.x = TRUE)
  edges[, `:=`(left_edge = fifelse(type == "DEL" & !is.na(baf_left), baf_left, depth_left),
               right_edge = fifelse(type == "DEL" & !is.na(baf_right), baf_right, depth_right),
               edge_source = fifelse(type == "DEL" & !is.na(baf_left), "BAF", "depth"))]
  # Breakpoint clusters per type: single linkage on (left, right) with the Chebyshev distance, cut at --cluster-tol
  edges[, cluster := {
    ok <- !is.na(left_edge) & !is.na(right_edge)
    cl <- rep(NA_integer_, .N)
    if (sum(ok) > 1) {
      m <- cbind(left_edge[ok], right_edge[ok])
      cl[ok] <- cutree(hclust(dist(m, method = "maximum"), method = "single"), h = opt$`cluster-tol`)
    } else if (sum(ok) == 1) cl[ok] <- 1L
    cl
  }, by = type]
  edges[, cl_left := median(left_edge, na.rm = TRUE), by = .(type, cluster)]
  edges[, cl_right := median(right_edge, na.rm = TRUE), by = .(type, cluster)]
  setorder(edges, type, cl_left, cl_right, left_edge, right_edge, na.last = TRUE)
  edges[, cluster := match(paste(cl_left, cl_right), unique(paste(cl_left, cl_right))), by = type]
  dels <- edges[type == "DEL", sample]; dups <- edges[type == "DUP", sample]
  if (!is.null(opt$`edges-out`)) fwrite(edges, opt$`edges-out`, sep = "\t")
  print(edges[, .(type, cluster, sample = sub(".*_", "#", sample), left_Mb = round(left_edge / 1e6, 3),
                  right_Mb = round(right_edge / 1e6, 3), edge_source)], nrows = 200)
}
page_note <- function(sel) {
  if (is.null(edges)) return("")
  f <- function(t) {
    x <- edges[sample %in% sel & type == t & !is.na(left_edge)]
    if (!nrow(x)) return(NULL)
    sprintf("%s cluster%s %s: edges %.2f-%.2f / %.2f-%.2f Mb", t, if (uniqueN(x$cluster) > 1) "s" else "",
            paste(unique(x$cluster), collapse = ","), min(x$left_edge) / 1e6, max(x$left_edge) / 1e6,
            min(x$right_edge) / 1e6, max(x$right_edge) / 1e6)
  }
  paste0("; ", paste(c(f("DEL"), f("DUP")), collapse = "; "))
}
ps <- opt$`page-size`
# Pages: consecutive chunks of --page-size, or (with breakpoint clusters) whole clusters packed up to --page-size
pack <- function(x, t) {
  if (ps <= 0) return(list(x))
  if (is.null(edges) || !length(x)) return(split(x, ceiling(seq_along(x) / ps)))
  cl <- edges[type == t][match(x, sample), cluster]
  pages <- list(); cur <- character(0)
  for (g in unique(cl)) {
    m <- x[cl == g]
    if (length(m) > ps) {  # a cluster larger than a page gets pages of its own
      if (length(cur)) { pages[[length(pages) + 1]] <- cur; cur <- character(0) }
      pages <- c(pages, split(m, ceiling(seq_along(m) / ps)))
    } else if (length(cur) + length(m) > ps) {
      pages[[length(pages) + 1]] <- cur; cur <- m
    } else cur <- c(cur, m)
  }
  if (length(cur)) pages[[length(pages) + 1]] <- cur
  pages
}
del_pages <- pack(dels, "DEL"); dup_pages <- pack(dups, "DUP")
np <- max(1, length(del_pages), length(dup_pages))
chunk <- function(x, k) {
  pg <- if (identical(x, dels)) del_pages else dup_pages
  if (k <= length(pg)) pg[[k]] else character(0)
}
ht <- opt$height + if (!is.null(ov)) 6 else 0
if (capabilities("aqua")) quartz(type = "pdf", file = opt$out, width = opt$width, height = ht) else
  if (capabilities("cairo")) cairo_pdf(opt$out, width = opt$width, height = ht, onefile = TRUE) else
    pdf(opt$out, width = opt$width, height = ht)
for (k in seq_len(np)) {
  hl <- c(chunk(dels, k), chunk(dups, k), if (ps == 0) setdiff(hl_all, c(dels, dups)))
  is_del <- hl %in% dels; is_dup <- hl %in% dups
  pal <- setNames(rep(c("#1B5E20", "#FF8C00", "#6A3D9A", "#00897B", "#8D6E63", "#C0CA33"), length.out = length(hl)), hl)
  if (any(is_del)) pal[is_del] <- tones("#7F0000", "#FF8C00", sum(is_del))   # dark red -> orange
  if (any(is_dup)) pal[is_dup] <- tones("#08306B", "#2E9E44", sum(is_dup))   # dark blue -> green
  show_leg <- length(hl) <= max(opt$`max-legend`, 2 * ps)
  bps <- lapply(bps_all, function(b) if (is.null(b)) NULL else b[sample %in% hl])
  lab <- sub("^; ", "", paste0(if (np > 1) sprintf("page %d of %d", k, np) else "", page_note(hl)))
  print(make_fig(if (nzchar(lab)) paste0(" (", lab, ")") else ""))
}
invisible(dev.off())
cat("wrote", opt$out, "\n")
