#!/usr/bin/env Rscript
# Recurrence (dot) plot of one read against its reference region, from read_ref_recurrence.py outputs.
# x = reference position (with SD and gene tracks underneath), y = position in the read; one point per matching window (forward strand blue, reverse
# complement red; shade = identity). The read's own alignment blocks are drawn as grey lines underneath, so repeats,
# duplications, inversions and deletion junctions stand out as off-diagonal or broken diagonals.
#
# Usage:
#   Rscript plot_read_ref_recurrence.R --prefix read1 --out figures/read1.recurrence.pdf [--segdup hg38.SegDup.bed.gz]
#     [--genes gencode.gtf.gz] [--title ...] [--size 12]
#   Rscript plot_read_ref_recurrence.R --hits read1.hits.tsv.gz --chrom chr16 --ref-label hg38 --out ...
#   With --prefix, <prefix>.hits.tsv[.gz] is read and .info.tsv / .aln.tsv[.gz] are used when present; without them
#   the axes span the hits and no alignment blocks are drawn.
#   --segdup / --genes add an SD track and a gene track under the reference (x) axis.
suppressPackageStartupMessages({library(ggplot2); library(data.table); library(optparse); library(patchwork)})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--prefix", help = "output prefix of read_ref_recurrence.py / read_t2t_recurrence.py"),
  make_option("--hits", help = "hits table (.tsv or .tsv.gz) when used without --prefix; info/aln files are then optional"),
  make_option("--chrom", help = "reference chromosome for the y axis (needed without an info file)"),
  make_option("--ref-label", dest = "ref_label", help = "reference name for the y axis, e.g. hg38"),
  make_option("--out"), make_option("--title"), make_option("--segdup"),
  make_option("--genes"), make_option("--size", type = "double", default = 12),
  make_option("--max-copies", dest = "max_copies", type = "integer", default = 5,
              help = "windows with more hits than this (interspersed repeats) are drawn as a faint grey background [%default]"))))
if ((is.null(opt$prefix) && is.null(opt$hits)) || is.null(opt$out)) stop("need --prefix or --hits, and --out (see -h)")
font_family <- if (capabilities("aqua") || capabilities("cairo")) "Arial" else "sans"
find_file <- function(suffix) {
  if (is.null(opt$prefix)) return(NULL)
  f <- paste0(opt$prefix, suffix, c("", ".gz"))
  f <- f[file.exists(f)]
  if (length(f)) f[1] else NULL
}
hits_file <- if (!is.null(opt$hits)) opt$hits else find_file(".hits.tsv")
if (is.null(hits_file)) stop("no hits table found for --prefix ", opt$prefix)
h <- fread(hits_file)
info_file <- find_file(".info.tsv")
I <- if (!is.null(info_file)) { info <- fread(info_file); setNames(as.character(info$value), info$key) } else character()
aln_file <- find_file(".aln.tsv")
aln <- if (!is.null(aln_file)) fread(aln_file) else data.table()
get <- function(key, default) if (key %in% names(I)) I[[key]] else default
win <- get("window", as.character(median(h$read_end - h$read_start + 1)))
region_chrom <- get("chrom", if (is.null(opt$chrom)) "" else opt$chrom)
ref_label <- if (!is.null(opt$ref_label)) opt$ref_label else get("reference", "")
yr <- if ("ref_start" %in% names(I)) as.numeric(c(I[["ref_start"]], I[["ref_end"]])) / 1e6 else
  c(min(h$ref_start), max(h$ref_end)) / 1e6
xr <- if ("read_seg_start" %in% names(I)) as.numeric(c(I[["read_seg_start"]], I[["read_seg_end"]])) / 1000 else
  c(min(h$read_start), max(h$read_end)) / 1000
read_name <- get("read_name", sub("\\.hits\\.tsv(\\.gz)?$", "", basename(hits_file)))
title <- if (is.null(opt$title)) sprintf("%s; %s bp windows%s", read_name, win,
                                         if ("step" %in% names(I)) sprintf(", step %s, edit distance <= %s", I[["step"]],
                                                                           I[["max_edit_distance"]]) else "") else opt$title
if (nrow(aln) && "chrom" %in% names(aln)) aln <- aln[chrom == region_chrom]
h[, n_hits := .N, by = read_start]
h[, rmid := (ref_start + ref_end) / 2e6][, qmid := (read_start + read_end) / 2e3]
theme_set(theme_bw(base_size = 22, base_family = font_family) + theme(panel.grid.minor = element_blank()))
# Main panel: reference on x (Mb), read on y (kb)
p <- ggplot() +
  { if (nrow(aln)) geom_segment(data = aln, aes(x = ifelse(strand == "+", ref_start, ref_end) / 1e6,
                                                xend = ifelse(strand == "+", ref_end, ref_start) / 1e6,
                                                y = read_start / 1e3, yend = read_end / 1e3),
                                colour = "grey75", linewidth = 3, lineend = "round") } +
  geom_point(data = h[n_hits > opt$max_copies], aes(rmid, qmid), colour = "grey70", alpha = 0.15, size = 0.5) +
  geom_point(data = h[n_hits <= opt$max_copies], aes(rmid, qmid, colour = strand, alpha = identity), size = 0.9) +
  scale_colour_manual(values = c("+" = "#1F5FA8", "-" = "#B2182B"), name = "strand",
                      labels = c("+" = "forward", "-" = "reverse complement")) +
  scale_alpha_continuous(range = c(0.25, 1), name = "identity") +
  coord_cartesian(xlim = yr, ylim = xr, expand = FALSE) +
  labs(x = NULL, y = if (get("read_orientation", "original") == "reverse_complemented")
         "position in read, reverse complemented (kb)" else "position in read (kb)", title = title) +
  guides(colour = guide_legend(override.aes = list(size = 4)), alpha = "none")
cap <- paste0(sprintf("faint grey points: windows with > %d hits (interspersed repeats)", opt$max_copies),
                        if (nrow(aln)) sprintf("; grey lines: the read's alignment blocks (primary %s)", get("primary", "")) else "")
p <- p + theme(legend.position = "top", legend.justification = "left", plot.title = element_text(size = 16),
        plot.caption = element_text(size = 13), axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Reference-axis tracks under the main panel (SDs, genes); the bottom track carries the x axis
xlab <- trimws(paste(ref_label, region_chrom, "position (Mb)"))
track_theme <- theme(panel.grid = element_blank(), plot.caption = element_text(size = 13), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     axis.title.y = element_text(angle = 0, vjust = 0.5, size = 18))
tracks <- list()
if (!is.null(opt$segdup)) {
  s <- fread(opt$segdup, header = FALSE)[, 1:3]; setnames(s, c("chrom", "start", "end"))
  s <- s[chrom == region_chrom & end / 1e6 > yr[1] & start / 1e6 < yr[2]]
  tracks$SD <- ggplot(s) + geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0, ymax = 1), fill = "#E8A33D") +
    coord_cartesian(xlim = yr, ylim = c(0, 1), expand = FALSE) + labs(x = NULL, y = "SD") + track_theme
}
if (!is.null(opt$genes)) {
  cmd <- sprintf("%s '%s' | awk -F'\\t' '$1 == \"%s\" && $3 == \"gene\" && $5 > %d && $4 < %d'",
                 if (grepl("\\.gz$", opt$genes)) "gzip -cd" else "cat", opt$genes, region_chrom, round(yr[1] * 1e6), round(yr[2] * 1e6))
  g <- tryCatch(fread(cmd = cmd, header = FALSE, sep = "\t", quote = ""), error = function(e) data.table())
  if (nrow(g)) {
    g <- g[grepl('(gene_type|gene_biotype)[ =]"?protein_coding', V9),
           .(start = V4, end = V5, name = sub('.*gene_name[ =]"?([^";]+).*', "\\1", V9))]
    setorder(g, start)
    # stack genes into lanes so that neither boxes nor labels overlap (label width estimated from its length)
    span <- diff(yr); g[, lab_w := nchar(name) * span * 0.011]
    g[, xs := pmin(start / 1e6, (start + end) / 2e6 - lab_w / 2)][, xe := pmax(end / 1e6, (start + end) / 2e6 + lab_w / 2)]
    lane_end <- numeric(); g[, lane := 0L]
    for (i in seq_len(nrow(g))) {
      l <- which(lane_end < g$xs[i] - span * 0.005)[1]
      if (is.na(l)) { lane_end <- c(lane_end, g$xe[i]); l <- length(lane_end) } else lane_end[l] <- g$xe[i]
      g$lane[i] <- l
    }
    tracks$genes <- ggplot(g) +
      geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -lane + 0.55, ymax = -lane + 0.85), fill = "#3C3C3C") +
      geom_text(aes(x = (start + end) / 2e6, y = -lane + 0.3, label = name), size = 7.76 * 0.6,
                family = font_family, fontface = "italic") +
      coord_cartesian(xlim = yr, ylim = c(-max(g$lane) - 0.1, 1), expand = FALSE, clip = "off") +
      labs(x = NULL, y = "genes") + track_theme
  }
}
panels <- c(list(p), tracks)
heights <- c(1, if (!is.null(tracks$SD)) 0.03, if (!is.null(tracks$genes)) 0.06 * (max(tracks$genes$data$lane) + 1))
for (i in seq_along(panels)) {
  last <- i == length(panels)
  panels[[i]] <- panels[[i]] + labs(x = if (last) xlab else NULL, caption = if (last) cap else NULL) +
    theme(axis.text.x = if (last) element_text() else element_blank(), axis.ticks.x = if (last) element_line() else element_blank())
}
fig <- wrap_plots(panels, ncol = 1, heights = heights)
h_in <- opt$size * sum(heights) / 1
if (capabilities("aqua")) quartz(type = "pdf", file = opt$out, width = opt$size * 1.1, height = h_in) else
  if (capabilities("cairo")) cairo_pdf(opt$out, width = opt$size * 1.1, height = h_in) else
    pdf(opt$out, width = opt$size * 1.1, height = h_in)
print(fig)
invisible(dev.off())
cat("wrote", opt$out, "\n")
