#!/usr/bin/env Rscript
# Several breakpoint-SD windows of one GD cluster side by side, under a whole-region depth panel:
#   top     SD + SD breakpoint tracks and copy number across the whole region (--overview-rd), windows shaded
#   columns one per window (--windows "NAME=baf,rd,start,end;NAME=...", plotted left to right): SD tracks, copy number,
#           folded BAF
# Carriers matching --highlight-pattern are coloured: --del-pattern red to orange, --dup-pattern blue to green; all other
# samples grey; black = cohort median. Dashed lines at --marks (Mb, comma-separated; e.g. GD edges) in every panel.
# Copy number = RD / genome-wide sample median (--medians) x 2, divided by the cohort median of each bin x 2; folded BAF
# = min(BAF, 1 - BAF) with NA counted as 0. Bin sizes: --bin-size (CN, default 1000) and --baf-bin (0 = window / 40).
#
# Usage:
#   Rscript plot_multi_regions.R --overview-rd 22q11/22q11_A-D.region.renamed.RD.txt.gz --medians GD_samples.RD_median \
#     --windows "A=22q11/22q11_A.loh.clean_sites.baf.tsv.gz,-,18601988,19032295;B=...;C=...;D=..." \
#     --highlight-pattern "GD_22q11\\.21-(AB|AD|BD|CD)_" --marks 18.924718,20.299685 \
#     --segdup "SD=hg38.SegDup.sorted.merged.bed.gz,SD breakpoint=GD_NAHR.breakpoint_SDs.bed.gz" \
#     --out figures/22q11_A-D.breakpoint_profile.pdf --title "22q11.21 A-D"
#   The rd field of a window may be "-" to cut the window out of --overview-rd.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(optparse); library(patchwork)})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--overview-rd"), make_option("--medians"),
  make_option("--windows", help = "NAME=baf,rd,start,end;NAME=... (windows left to right)"),
  make_option("--highlight-pattern", default = "GD_"), make_option("--del-pattern", default = "_DEL_"),
  make_option("--dup-pattern", default = "_DUP_"), make_option("--marks"), make_option("--segdup"),
  make_option("--bin-size", type = "integer", default = 1000), make_option("--baf-bin", type = "integer", default = 0),
  make_option("--max-legend", type = "integer", default = 30),
  make_option("--cn-max", type = "double", default = 6, help = "upper limit of the copy-number axis"),
  make_option("--out"), make_option("--title", default = ""),
  make_option("--width", type = "double", default = 0), make_option("--height", type = "double", default = 20))))
for (x in c("overview-rd", "medians", "out")) if (is.null(opt[[x]])) stop("need --overview-rd, --medians and --out")
if (is.null(opt$windows)) stop("need --windows")
win_specs <- strsplit(opt$windows, ";", fixed = TRUE)[[1]]
font_family <- if (capabilities("aqua") || capabilities("cairo")) "Arial" else "sans"
m <- fread(opt$medians, header = FALSE); med <- setNames(as.numeric(unlist(m[2])), unlist(m[1]))
marks <- if (!is.null(opt$marks)) as.numeric(strsplit(opt$marks, ",")[[1]]) else numeric(0)

read_rd <- function(path) { r <- fread(path); setnames(r, 1:3, c("chrom", "start", "end")); r }
to_cn <- function(r, b, win = NULL) {
  if (!is.null(win)) r <- r[end > win[1] & start < win[2]]
  rs <- intersect(names(r)[-(1:3)], names(med))
  if (!length(rs)) stop("no sample of the depth matrix has a median in --medians (check that both use the same, renamed, sample IDs)")
  r[, bin := floor(start / b) * b + b / 2]
  rl <- melt(r[, c("bin", rs), with = FALSE], id.vars = "bin", variable.name = "sample", value.name = "RD",
             variable.factor = FALSE)
  cn <- rl[, .(RD = sum(RD), n = .N), by = .(sample, bin)]
  cn[, CN := RD / (med[sample] * n) * 2]
  cn[, CN := CN / median(CN) * 2, by = bin]
  cn <- cn[is.finite(CN)]
  if (!nrow(cn)) return(cn)
  grid <- CJ(sample = unique(cn$sample), bin = seq(min(cn$bin), max(cn$bin), by = b))  # NA rows break lines at gaps
  merge(grid, cn[, .(sample, bin, CN)], by = c("sample", "bin"), all.x = TRUE)
}
read_sd <- function(spec, region_chrom, xr) {
  if (is.null(spec)) return(NULL)
  parts <- strsplit(strsplit(spec, ",")[[1]], "=", fixed = TRUE)
  rbindlist(lapply(parts, function(p) {
    b <- fread(p[length(p)], header = "auto"); setnames(b, 1:3, c("chrom", "start", "end"))
    lab <- intersect(c("breakpoint", "name"), names(b))
    b[, label := if (length(lab)) as.character(get(lab[1])) else ""]
    b <- b[chrom == region_chrom & end / 1e6 > xr[1] & start / 1e6 < xr[2]]
    b[, .(track = p[1], start = pmax(start / 1e6, xr[1]), end = pmin(end / 1e6, xr[2]), label)]
  }))
}

ov_rd <- read_rd(opt$`overview-rd`)
chrom <- ov_rd$chrom[1]
ov_xr <- range(c(ov_rd$start, ov_rd$end)) / 1e6
ov <- to_cn(ov_rd, max(opt$`bin-size`, round(diff(ov_xr) * 1e6 / 3000 / 1000) * 1000))
wins <- lapply(win_specs, function(w) {
  nm <- sub("=.*", "", w); f <- strsplit(sub("^[^=]*=", "", w), ",")[[1]]
  win <- as.numeric(f[3:4])
  b <- fread(f[1], na.strings = c("NA", "."))
  s <- setdiff(names(b), c("CHROM", "POS", "ID", "REF", "ALT"))
  bfb <- if (opt$`baf-bin` > 0) opt$`baf-bin` else max(2000, round(diff(win) / 40 / 1000) * 1000)
  bl <- melt(b[, c("POS", s), with = FALSE], id.vars = "POS", variable.name = "sample", value.name = "BAF",
             variable.factor = FALSE)
  bl[, mBAF := pmin(BAF, 1 - BAF)][is.na(mBAF), mBAF := 0]
  bl[, bin := floor(POS / bfb) * bfb + bfb / 2]
  mb <- bl[, .(mBAF = mean(mBAF)), by = .(sample, bin)]
  rd <- if (f[2] == "-") ov_rd else read_rd(f[2])
  list(name = nm, xr = win / 1e6, cn = to_cn(rd, opt$`bin-size`, win), mb = mb, bfb = bfb, samples = s)
})

hl <- sort(grep(opt$`highlight-pattern`, unique(c(ov$sample, unlist(lapply(wins, `[[`, "samples")))), value = TRUE))
nat <- function(x) x[order(sub("_[0-9]+$", "", x), suppressWarnings(as.numeric(sub(".*_", "", x))))]
dels <- nat(hl[grepl(opt$`del-pattern`, hl)]); dups <- nat(hl[grepl(opt$`dup-pattern`, hl) & !hl %in% dels])
hl <- c(dels, dups, setdiff(hl, c(dels, dups)))
tones <- function(lo, hi, n) if (n == 1) hi else grDevices::colorRampPalette(c(lo, hi))(n)
pal <- setNames(rep("#6A3D9A", length(hl)), hl)
if (length(dels)) pal[dels] <- tones("#7F0000", "#FF8C00", length(dels))   # dark red -> orange
if (length(dups)) pal[dups] <- tones("#08306B", "#2E9E44", length(dups))   # dark blue -> green
cat("coloured:", length(dels), "DEL,", length(dups), "DUP\n")
show_leg <- length(hl) <= opt$`max-legend`

theme_set(theme_bw(base_size = 22, base_family = font_family) +
          theme(panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold")))
vmarks <- function(xr) { mk <- marks[marks >= xr[1] & marks <= xr[2]]
  if (length(mk)) geom_vline(xintercept = mk, linetype = "dashed", colour = "#1F5FA8") }
sd_panel <- function(xr, first, title = NULL, shade = NULL) {
  s <- read_sd(opt$segdup, chrom, xr)
  if (is.null(s) || !nrow(s)) return(ggplot() + theme_void() + labs(title = title))
  s[, track := factor(track, levels = rev(unique(track)))][label %in% c("NA", NA), label := ""]
  ggplot(s) + geom_rect(aes(xmin = start, xmax = end, ymin = as.numeric(track) - 0.4, ymax = as.numeric(track) + 0.4),
                        fill = "#E8A33D", colour = "grey30", linewidth = 0.2) +
    geom_text(aes(x = (start + end) / 2, y = as.numeric(track), label = label), size = 7.76 * 0.7,
              family = font_family, fontface = "bold") + vmarks(xr) +
    scale_y_continuous(breaks = seq_along(levels(s$track)), labels = if (first) levels(s$track) else NULL,
                       expand = c(0, 0.3)) +
    coord_cartesian(xlim = xr) + labs(x = NULL, y = NULL, title = title) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())
}
lines_panel <- function(d, v, ylim, xr, ylab, extra = NULL) {
  coh <- d[, .(med = median(get(v), na.rm = TRUE)), by = bin]
  ggplot() + extra +
    geom_line(data = d[!sample %in% hl], aes(bin / 1e6, .data[[v]], group = sample), colour = "grey70", alpha = 0.25,
              linewidth = 0.4) +
    geom_line(data = coh, aes(bin / 1e6, med), colour = "black", linewidth = 1.2) +
    geom_line(data = d[sample %in% hl], aes(bin / 1e6, .data[[v]], colour = sample), linewidth = 0.8) +
    vmarks(xr) +
    scale_colour_manual(values = pal, name = NULL, guide = if (show_leg) "legend" else "none") +
    coord_cartesian(xlim = xr, ylim = ylim) + labs(x = NULL, y = ylab)
}
cn_max <- opt$`cn-max`
hlines <- geom_hline(yintercept = c(1, 2, 3), linetype = "dotted", colour = "grey40")
zoom <- data.table(xmin = sapply(wins, function(w) w$xr[1]), xmax = sapply(wins, function(w) w$xr[2]),
                   name = sapply(wins, `[[`, "name"))
p_ov <- lines_panel(ov, "CN", c(0, cn_max), ov_xr, "copy number",
                    list(geom_rect(data = zoom, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                                   fill = "#9ECAE1", alpha = 0.25, inherit.aes = FALSE),
                         geom_text(data = zoom, aes(x = (xmin + xmax) / 2, y = cn_max * 0.97, label = name),
                                   size = 7.76, family = font_family, inherit.aes = FALSE), hlines)) +
  labs(x = paste(chrom, "position (Mb)")) + guides(colour = "none")
top <- sd_panel(ov_xr, TRUE, sprintf("whole region (CN %g kb bins; shaded = SD windows below)",
                                     (ov$bin[2] - ov$bin[1]) / 1000)) / p_ov + plot_layout(heights = c(0.18, 1.3))

cols <- lapply(seq_along(wins), function(i) {
  w <- wins[[i]]; first <- i == 1
  yl <- function(s) if (first) s else NULL
  list(sd_panel(w$xr, first, sprintf("%s (CN %g kb, BAF %g kb bins)", w$name, opt$`bin-size` / 1000, w$bfb / 1000)),
       lines_panel(w$cn, "CN", c(0, cn_max), w$xr, yl("copy number"), hlines),
       lines_panel(w$mb, "mBAF", c(0, 0.5), w$xr, yl("folded BAF")) + labs(x = paste(chrom, "(Mb)")))
})
grid <- wrap_plots(c(do.call(rbind, cols)), ncol = length(wins), byrow = TRUE) +
  plot_layout(heights = c(0.18, 1.2, 1), guides = "collect")
fig <- wrap_plots(top, grid, ncol = 1, heights = c(1.5, 2.4)) +
  plot_annotation(title = sprintf("%s\ncoloured: %d DEL in red to orange, %d DUP in blue to green", opt$title,
                                  length(dels), length(dups)),
                  theme = theme(plot.title = element_text(size = 22, family = font_family)))
wd <- if (opt$width > 0) opt$width else 10 + 9 * length(wins)
if (capabilities("aqua")) quartz(type = "pdf", file = opt$out, width = wd, height = opt$height) else
  if (capabilities("cairo")) cairo_pdf(opt$out, width = wd, height = opt$height) else pdf(opt$out, width = wd, height = opt$height)
print(fig)
invisible(dev.off())
cat("wrote", opt$out, "\n")
