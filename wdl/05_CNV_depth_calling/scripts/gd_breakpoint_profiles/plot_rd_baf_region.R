#!/usr/bin/env Rscript
# Read depth (copy number) and folded BAF across a region, one line per sample, all samples together; samples whose ID
# matches --highlight-pattern (regex) are coloured, the rest grey; black = cohort median.
#
# Inputs:
#   --rd       depth matrix: #Chr Start End + one column per sample (e.g. 100 bp bins), plain or .gz; optional: without
#              --rd (and --medians) the depth panel and page 2 are left out
#   --medians  2 rows, tab-separated: sample IDs / median depth per sample (same per-bin scale as --rd)
#   --baf      BAF table: CHROM POS ID REF ALT + one column per sample; best the clean-site table of
#              find_baf_terminal_loh.py (high-confidence biallelic SNVs)
# Copy number: CN = RD / sample median x 2, then (default) divided by the cohort median CN of the same bin and x 2, which
#   removes bin-level mappability/GC effects (assumes most samples are CN 2 in each bin; --no-bin-norm skips it).
#   Bins are merged to --bin-size (sum of RD over the merged bins) before normalising.
# Folded BAF: min(BAF, 1 - BAF) (0 = hom, 0.5 = balanced het); NA counts as 0 unless --na-missing; mean per --baf-bin.
#
# Pages:
#   1  CN (top) and folded BAF (bottom), all samples
#   2  CN, one facet per highlighted sample, with the cohort median and 5-95% band
#   3  folded BAF, one facet per highlighted sample, same layout
#
# Usage:
#   Rscript plot_rd_baf_region.R --rd aou_II/chr16_bp4.renamed.RD.txt.gz --medians aou_II/chr16_bp4.renamed.RD_median \
#     --baf aou_II/chr16_bp4.renamed.right_loh.clean_sites.baf.tsv.gz --highlight-pattern "GD_16p11.2-BP4-5" \
#     --out aou_II/figures/chr16_bp4.renamed.rd_baf.pdf [--bin-size 2000] [--baf-bin 10000] [--title "chr16 BP4"]
#     [--no-bin-norm] [--na-missing] [--ncol 4]
#     [--segdup "SD=references/hg38.SegDup.sorted.merged.bed.gz,GD BP=references/GD_NAHR.breakpoint_SDs.bed.gz"]
#   --segdup adds a segmental-duplication track above the depth and BAF panels of page 1 (one row per NAME=bed).
#     [--ld-prefix aou_II/chr16_bp4.renamed.ld [--ld-stat mean_r2|max_r2] [--mark 29.464,29.534]]
#   --ld-prefix adds an LD triangle (compute_ld.py <prefix>.binned.tsv: window-pair r^2) under the BAF panel of page 1,
#   on the same position axis; --mark draws dashed lines (Mb) through all page-1 panels.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(optparse); library(patchwork)})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--rd"), make_option("--medians"), make_option("--baf"), make_option("--out"),
  make_option("--highlight-pattern", default = "GD_16p11.2-BP4-5"), make_option("--title", default = ""),
  make_option("--bin-size", type = "integer", default = 2000), make_option("--baf-bin", type = "integer", default = 10000),
  make_option("--no-bin-norm", action = "store_true", default = FALSE),
  make_option("--na-missing", action = "store_true", default = FALSE),
  make_option("--ncol", type = "integer", default = 4),
  make_option("--segdup", help = "NAME=bed[,NAME=bed] SD tracks drawn above page 1"),
  make_option("--ld-prefix", help = "compute_ld.py output prefix; adds an LD triangle under page 1"),
  make_option("--ld-stat", default = "mean_r2"), make_option("--mark", help = "positions in Mb, comma-separated"))))
for (x in c("baf", "out")) if (is.null(opt[[x]])) stop("need --baf and --out (see -h)")
has_rd <- !is.null(opt$rd)
if (has_rd && is.null(opt$medians)) stop("--rd needs --medians")
font_family <- if (capabilities("aqua") || capabilities("cairo")) "Arial" else "sans"
open_pdf <- function(file, w, h) {
  if (capabilities("aqua")) quartz(type = "pdf", file = file, width = w, height = h) else
    if (capabilities("cairo")) cairo_pdf(file, width = w, height = h, onefile = TRUE) else pdf(file, width = w, height = h)
}

# Copy number from depth
# SD track: --segdup "NAME=bed[,NAME=bed]" (plain or .gz; first 3 columns chrom/start/end; label from a 'breakpoint' or
# 'name' column, or a 4th text column, when present), one row per file, clipped to the plotted window.
sd_track <- function(spec, region_chrom, xr, font_family, expand = TRUE) {
  if (is.null(spec) || !nzchar(spec)) return(NULL)
  parts <- strsplit(strsplit(spec, ",")[[1]], "=", fixed = TRUE)
  tr <- rbindlist(lapply(seq_along(parts), function(i) {
    nm <- if (length(parts[[i]]) > 1) parts[[i]][1] else basename(parts[[i]][1])
    f <- parts[[i]][length(parts[[i]])]
    b <- fread(f, header = "auto")
    setnames(b, 1:3, c("chrom", "start", "end"))
    lab_col <- intersect(c("breakpoint", "name"), names(b))
    b[, label := if (length(lab_col)) as.character(get(lab_col[1])) else
                 if (ncol(b) >= 4 && is.character(b[[4]])) b[[4]] else ""]
    b <- b[chrom == region_chrom & end / 1e6 > xr[1] & start / 1e6 < xr[2]]
    if (!nrow(b)) return(NULL)
    b[, .(track = nm, start = pmax(start / 1e6, xr[1]), end = pmin(end / 1e6, xr[2]), label)]
  }))
  if (!nrow(tr)) return(NULL)
  tr[, track := factor(track, levels = rev(unique(track)))]
  tr[label %in% c("NA", NA), label := ""]
  ggplot(tr) + geom_rect(aes(xmin = start, xmax = end, ymin = as.numeric(track) - 0.4, ymax = as.numeric(track) + 0.4),
                         fill = "#E8A33D", colour = "grey30", linewidth = 0.2) +
    geom_text(aes(x = (start + end) / 2, y = as.numeric(track), label = label), size = 7.76 * 0.7,
              family = font_family, fontface = "bold") +
    scale_y_continuous(breaks = seq_along(levels(tr$track)), labels = levels(tr$track), expand = c(0, 0.3)) +
    coord_cartesian(xlim = xr, expand = expand) + labs(x = NULL, y = NULL) +
    theme_bw(base_size = 22, base_family = font_family) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())
}

b <- fread(opt$baf, na.strings = c("NA", "."))
cn <- data.table(sample = character(0), bin = numeric(0), CN = numeric(0))
if (has_rd) {
  rd <- fread(opt$rd)
  setnames(rd, 1:3, c("chrom", "start", "end"))
  samples <- names(rd)[-(1:3)]
  if (rd$chrom[1] != b$CHROM[1]) warning("depth file chromosome (", rd$chrom[1], ") differs from the BAF table (", b$CHROM[1], ")")
  med <- fread(opt$medians, header = FALSE)
  med <- setNames(as.numeric(unlist(med[2])), unlist(med[1]))
  miss <- setdiff(samples, names(med))
  if (length(miss)) stop(length(miss), " samples have no median, e.g. ", miss[1])
  rd[, bin := floor(start / opt$`bin-size`) * opt$`bin-size` + opt$`bin-size` / 2]
  rdl <- melt(rd[, c("bin", samples), with = FALSE], id.vars = "bin", variable.name = "sample", value.name = "RD",
              variable.factor = FALSE)
  cn <- rdl[, .(RD = sum(RD), n = .N), by = .(sample, bin)]
  cn[, CN := RD / (med[sample] * n) * 2]
  if (!opt$`no-bin-norm`) cn[, CN := CN / median(CN) * 2, by = bin]
} else {
  samples <- setdiff(names(b), c("CHROM", "POS", "ID", "REF", "ALT"))
}
chrom <- b$CHROM[1]

# Folded BAF
bs <- intersect(setdiff(names(b), c("CHROM", "POS", "ID", "REF", "ALT")), samples)
bl <- melt(b[, c("POS", bs), with = FALSE], id.vars = "POS", variable.name = "sample", value.name = "BAF",
           variable.factor = FALSE)
bl[, mBAF := pmin(BAF, 1 - BAF)]
if (!opt$`na-missing`) bl[is.na(mBAF), mBAF := 0]
bl[, bin := floor(POS / opt$`baf-bin`) * opt$`baf-bin` + opt$`baf-bin` / 2]
mb <- bl[, .(mBAF = mean(mBAF, na.rm = TRUE), n = sum(!is.na(mBAF))), by = .(sample, bin)][n > 0]

hl <- grep(opt$`highlight-pattern`, samples, value = TRUE)
cat(length(samples), "samples;", uniqueN(cn$bin), "CN bins of", opt$`bin-size`, "bp;", nrow(b), "BAF sites;",
    length(hl), "highlighted:\n", paste(hl, collapse = "\n "), "\n")
pal <- setNames(rep(c("#B2182B", "#1F5FA8", "#1B5E20", "#FF8C00", "#6A3D9A", "#00897B", "#8D6E63", "#E91E63",
                      "#5E35B1", "#C0CA33", "#039BE5", "#6D4C41"), length.out = max(1, length(hl)))[seq_along(hl)], hl)
band <- function(d, v) d[, .(med = median(get(v)), lo = quantile(get(v), 0.05), hi = quantile(get(v), 0.95)), by = bin]
cn_c <- band(cn, "CN"); mb_c <- band(mb, "mBAF")
cn_max <- if (has_rd) max(4, quantile(cn$CN, 0.999)) else 4

theme_set(theme_bw(base_size = 22, base_family = font_family) +
          theme(panel.grid.minor = element_blank(), strip.background = element_blank()))
xl <- paste(chrom, "position (Mb)")
ycn <- sprintf("copy number (%g kb bins)", opt$`bin-size` / 1000)
ybaf <- sprintf("folded BAF (%g kb mean)", opt$`baf-bin` / 1000)
all_lines <- function(d, coh, v, ylab, ylim) {
  ggplot() +
    geom_line(data = d[!sample %in% hl], aes(bin / 1e6, .data[[v]], group = sample), colour = "grey70", alpha = 0.25,
              linewidth = 0.4) +
    geom_line(data = coh, aes(bin / 1e6, med), colour = "black", linewidth = 1.4) +
    geom_line(data = d[sample %in% hl], aes(bin / 1e6, .data[[v]], colour = sample), linewidth = 0.9) +
    scale_colour_manual(values = pal, name = NULL) + coord_cartesian(xlim = xr, ylim = ylim) + labs(x = xl, y = ylab)
}
xr <- range(c(cn$bin, mb$bin)) / 1e6
p_cn <- if (has_rd) all_lines(cn, cn_c, "CN", ycn, c(0, cn_max)) +
  geom_hline(yintercept = c(1, 2, 3), linetype = "dotted", colour = "grey40") + labs(x = NULL) else NULL
p_b <- all_lines(mb, mb_c, "mBAF", ybaf, c(0, 0.5))
trk <- sd_track(opt$segdup, chrom, xr, font_family)

# LD triangle: each window pair drawn as a diamond centred between the two windows, depth = their distance / 2
p_ld <- NULL
if (!is.null(opt$`ld-prefix`)) {
  bn <- fread(paste0(opt$`ld-prefix`, ".binned.tsv"))
  w <- min(diff(sort(unique(c(bn$bin1, bn$bin2)))))
  bn[, v := get(opt$`ld-stat`)]
  poly <- bn[, .(x = c(bin1 + bin2, bin1 + bin2 + w, bin1 + bin2 + 2 * w, bin1 + bin2 + w) / 2,
                 y = c(bin2 - bin1, bin2 - bin1 + w, bin2 - bin1, bin2 - bin1 - w) / 2), by = .(id = seq_len(nrow(bn)), v)]
  p_ld <- ggplot(poly, aes(x / 1e6, -y / 1e6, group = id, fill = v)) + geom_polygon(colour = NA) +
    scale_fill_gradientn(colours = c("white", "#FEE0D2", "#FC9272", "#DE2D26", "#67000D"), limits = c(0, 1),
                         name = if (opt$`ld-stat` == "mean_r2") expression("mean r"^2) else expression("max r"^2)) +
    coord_cartesian(xlim = xr) + labs(x = xl, y = sprintf("LD (%g kb windows)", w / 1000)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())
  p_b <- p_b + labs(x = NULL)
}
marks <- if (!is.null(opt$mark)) as.numeric(strsplit(opt$mark, ",")[[1]]) else numeric(0)
add_marks <- function(p) if (length(marks)) p + geom_vline(xintercept = marks, linetype = "dashed", colour = "#1F5FA8") else p
panels <- Filter(Negate(is.null), list(trk, if (has_rd) add_marks(p_cn), add_marks(p_b), if (!is.null(p_ld)) add_marks(p_ld)))
hts <- c(if (!is.null(trk)) 0.15, if (has_rd) 1.2, 1, if (!is.null(p_ld)) 1.4)
page1 <- wrap_plots(panels, ncol = 1) + plot_layout(guides = "collect", heights = hts) +
  plot_annotation(title = sprintf("%s  %d samples; coloured: '%s' (%d)", opt$title, length(samples),
                                  opt$`highlight-pattern`, length(hl)),
                  theme = theme(plot.title = element_text(size = 22, family = font_family)))

per_sample <- function(d, coh, v, ylab, ylim, extra = NULL) {
  x <- merge(d[sample %in% hl], coh, by = "bin")
  x[, sample := factor(sample, levels = hl)]
  ggplot(x, aes(bin / 1e6)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey85") + extra +
    geom_line(aes(y = med), colour = "black", linewidth = 0.7) +
    geom_line(aes(y = .data[[v]], colour = sample), linewidth = 1, show.legend = FALSE) +
    scale_colour_manual(values = pal) + coord_cartesian(ylim = ylim) +
    facet_wrap(~sample, ncol = opt$ncol, labeller = label_wrap_gen(30)) + labs(x = xl, y = ylab) +
    theme(strip.text = element_text(size = 14))
}
nr <- ceiling(max(1, length(hl)) / opt$ncol)
open_pdf(opt$out, 22, (if (is.null(p_ld)) 14 else 20) - (if (has_rd) 0 else 5))
print(page1)
if (length(hl)) {
  if (has_rd) print(per_sample(cn, cn_c, "CN", ycn, c(0, cn_max),
                               geom_hline(yintercept = c(1, 2, 3), linetype = "dotted", colour = "grey40")) +
                      ggtitle(paste(opt$title, "copy number, coloured samples vs cohort (grey band = 5-95%)")))
  print(per_sample(mb, mb_c, "mBAF", ybaf, c(0, 0.5)) +
          ggtitle(paste(opt$title, "folded BAF, coloured samples vs cohort (grey band = 5-95%)")))
}
invisible(dev.off())
cat("wrote", opt$out, "\n")
