#!/usr/bin/env Rscript
# Recurrence (dot) plot of one read against its reference region(s), from read_ref_recurrence.py,
# read_t2t_recurrence.py or read_split_recurrence.py outputs.
# x = reference position (with SD and gene tracks underneath), y = position in the read; one point per matching window
# (forward strand blue, reverse complement red; shade = identity). The read's own alignment blocks are drawn as grey
# lines underneath, so repeats, duplications, inversions and deletion junctions stand out as off-diagonal or broken
# diagonals. Several reference segments (info key "segments", hits column "segment" = 1-based segment index; e.g. the
# two sides of a split alignment) are drawn side by side, each panel as wide as its segment.
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
info_get <- function(key, default) if (key %in% names(I)) I[[key]] else default
win <- info_get("window", as.character(median(h$read_end - h$read_start + 1)))
region_chrom <- info_get("chrom", if (is.null(opt$chrom)) "" else opt$chrom)
ref_label <- if (!is.null(opt$ref_label)) opt$ref_label else info_get("reference", "")
xr <- if ("read_seg_start" %in% names(I)) as.numeric(c(I[["read_seg_start"]], I[["read_seg_end"]])) / 1000 else
  c(min(h$read_start), max(h$read_end)) / 1000
read_name <- info_get("read_name", sub("\\.hits\\.tsv(\\.gz)?$", "", basename(hits_file)))
title <- if (is.null(opt$title)) sprintf("%s; %s bp windows%s", read_name, win,
                                         if ("step" %in% names(I)) sprintf(", step %s, edit distance <= %s", I[["step"]],
                                                                           I[["max_edit_distance"]]) else "") else opt$title

# Reference segments (one panel each)
segs <- if ("segments" %in% names(I)) {
  r <- strsplit(I[["segments"]], ";")[[1]]
  data.table(chrom = sub(":.*", "", r), start = as.numeric(sub(".*:([0-9]+)-.*", "\\1", r)),
             end = as.numeric(sub(".*-", "", r)))
} else if ("segment" %in% names(h)) {
  # no info file: each segment spans its hits
  h[, .(chrom = region_chrom, start = min(ref_start), end = max(ref_end)), keyby = segment][, segment := NULL][]
} else if ("ref_start" %in% names(I)) {
  data.table(chrom = region_chrom, start = as.numeric(I[["ref_start"]]), end = as.numeric(I[["ref_end"]]))
} else data.table(chrom = region_chrom, start = min(h$ref_start), end = max(h$ref_end))
segs[, label := sprintf("%s:%.3f-%.3f Mb", chrom, start / 1e6, end / 1e6)]
segs[, seg := factor(label, levels = label)]
multi <- nrow(segs) > 1
h[, seg := segs$seg[if ("segment" %in% names(h)) segment else 1L]]
in_seg <- function(d, chrom_col, s_col, e_col) {
  # segment of each feature (first segment it overlaps), features clipped to it; features outside all segments dropped
  if (!nrow(d)) return(d)
  out <- rbindlist(lapply(seq_len(nrow(segs)), function(k) {
    x <- d[get(chrom_col) == segs$chrom[k] & get(e_col) >= segs$start[k] & get(s_col) <= segs$end[k]]
    if (nrow(x)) x[, seg := segs$seg[k]]
    x
  }), fill = TRUE)
  out
}
if (nrow(aln) && "chrom" %in% names(aln)) aln <- in_seg(aln, "chrom", "ref_start", "ref_end") else
  if (nrow(aln)) aln[, seg := segs$seg[1]]
h[, n_hits := .N, by = read_start]
h[, rmid := (ref_start + ref_end) / 2e6][, qmid := (read_start + read_end) / 2e3]
bounds <- segs[, .(x = c(start, end) / 1e6), by = seg]
facet <- facet_grid(cols = vars(seg), scales = "free_x", space = "free_x")
xscale <- scale_x_continuous(expand = c(0, 0))
theme_set(theme_bw(base_size = 22, base_family = font_family) +
            theme(panel.grid.minor = element_blank(), panel.spacing.x = unit(1.2, "lines")))

# Main panel: reference on x (Mb), read on y (kb)
p <- ggplot() +
  geom_blank(data = bounds, aes(x = x, y = xr[1])) +
  { if (nrow(aln)) geom_segment(data = aln, aes(x = ifelse(strand == "+", ref_start, ref_end) / 1e6,
                                                xend = ifelse(strand == "+", ref_end, ref_start) / 1e6,
                                                y = read_start / 1e3, yend = read_end / 1e3),
                                colour = "grey75", linewidth = 3, lineend = "round") } +
  geom_point(data = h[n_hits > opt$max_copies], aes(rmid, qmid), colour = "grey70", alpha = 0.15, size = 0.5) +
  geom_point(data = h[n_hits <= opt$max_copies], aes(rmid, qmid, colour = strand, alpha = identity), size = 0.9) +
  scale_colour_manual(values = c("+" = "#1F5FA8", "-" = "#B2182B"), name = "strand",
                      labels = c("+" = "forward", "-" = "reverse complement")) +
  scale_alpha_continuous(range = c(0.25, 1), name = "identity") +
  facet + xscale + coord_cartesian(ylim = xr, expand = FALSE) +
  labs(x = NULL, y = if (info_get("read_orientation", "original") == "reverse_complemented")
         "position in read, reverse complemented (kb)" else "position in read (kb)", title = title) +
  guides(colour = guide_legend(override.aes = list(size = 4)), alpha = "none") +
  theme(legend.position = "top", legend.justification = "left", plot.title = element_text(size = 16),
        strip.text = if (multi) element_text(size = 16) else element_blank(),
        strip.background = if (multi) element_rect(fill = "grey92") else element_blank())
cap <- paste0(sprintf("faint grey points: windows with > %d hits (interspersed repeats)", opt$max_copies),
              if (nrow(aln)) sprintf("; grey lines: the read's alignment blocks (primary %s)", info_get("primary", "")) else "")

# Reference-axis tracks under the main panel (SDs, genes); the bottom track carries the x axis
chroms <- unique(segs$chrom)
xlab <- trimws(paste(ref_label, if (length(chroms) == 1) chroms else "", "position (Mb)"))
track_theme <- theme(panel.grid = element_blank(), plot.caption = element_text(size = 13),
                     axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     axis.title.y = element_text(angle = 0, vjust = 0.5, size = 18),
                     strip.text = element_blank(), strip.background = element_blank())
tracks <- list()
if (!is.null(opt$segdup)) {
  s <- fread(opt$segdup, header = FALSE)[, 1:3]; setnames(s, c("chrom", "start", "end"))
  s <- in_seg(s, "chrom", "start", "end")
  if (nrow(s)) s <- merge(s, segs[, .(seg, s0 = start, s1 = end)], by = "seg")[, `:=`(start = pmax(start, s0), end = pmin(end, s1))]
  tracks$SD <- ggplot() + geom_blank(data = bounds, aes(x = x, y = 0)) +
    { if (nrow(s)) geom_rect(data = s, aes(xmin = start / 1e6, xmax = end / 1e6, ymin = 0, ymax = 1), fill = "#E8A33D") } +
    facet + xscale + coord_cartesian(ylim = c(0, 1), expand = FALSE) + labs(x = NULL, y = "SD") + track_theme
}
if (!is.null(opt$genes)) {
  awk_cond <- paste(sprintf('($1 == "%s" && $5 > %d && $4 < %d)', segs$chrom, round(segs$start), round(segs$end)), collapse = " || ")
  cmd <- sprintf("%s '%s' | awk -F'\\t' '$3 == \"gene\" && (%s)'",
                 if (grepl("\\.gz$", opt$genes)) "gzip -cd" else "cat", opt$genes, awk_cond)
  g <- tryCatch(fread(cmd = cmd, header = FALSE, sep = "\t", quote = ""), error = function(e) data.table())
  if (nrow(g)) {
    g <- g[grepl('(gene_type|gene_biotype)[ =]"?protein_coding', V9),
           .(chrom = V1, start = V4, end = V5, name = sub('.*gene_name[ =]"?([^";]+).*', "\\1", V9))]
    g <- in_seg(g, "chrom", "start", "end")
  }
  if (nrow(g)) {
    # per segment, stack genes into lanes so that neither boxes nor labels overlap (label width from its length)
    g <- merge(g, segs[, .(seg, span = (end - start) / 1e6)], by = "seg")
    setorder(g, seg, start)
    g[, lab_w := nchar(name) * span * 0.011 * nrow(segs)]
    g[, xs := pmin(start / 1e6, (start + end) / 2e6 - lab_w / 2)][, xe := pmax(end / 1e6, (start + end) / 2e6 + lab_w / 2)]
    g[, lane := {
      lane_end <- numeric(); ln <- integer(.N)
      for (i in seq_len(.N)) {
        l <- which(lane_end < xs[i] - span[i] * 0.005)[1]
        if (is.na(l)) { lane_end <- c(lane_end, xe[i]); l <- length(lane_end) } else lane_end[l] <- xe[i]
        ln[i] <- l
      }
      ln
    }, by = seg]
    tracks$genes <- ggplot(g) + geom_blank(data = bounds, aes(x = x, y = 0)) +
      geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -lane + 0.55, ymax = -lane + 0.85), fill = "#3C3C3C") +
      geom_text(aes(x = (start + end) / 2e6, y = -lane + 0.3, label = name), size = 7.76 * 0.6,
                family = font_family, fontface = "italic") +
      facet + xscale + coord_cartesian(ylim = c(-max(g$lane) - 0.1, 1), expand = FALSE, clip = "off") +
      labs(x = NULL, y = "genes") + track_theme
  }
}
panels <- c(list(p), tracks)
heights <- c(1, if (!is.null(tracks$SD)) 0.03, if (!is.null(tracks$genes)) 0.06 * (max(tracks$genes$data$lane) + 1))
for (i in seq_along(panels)) {
  last <- i == length(panels)
  panels[[i]] <- panels[[i]] + labs(x = if (last) xlab else NULL, caption = if (last) cap else NULL) +
    theme(axis.text.x = if (last) element_text() else element_blank(),
          axis.ticks.x = if (last) element_line() else element_blank())
}
fig <- wrap_plots(panels, ncol = 1, heights = heights)
h_in <- opt$size * sum(heights)
w_in <- opt$size * 1.1 * (if (multi) 1.4 else 1)
if (capabilities("aqua")) quartz(type = "pdf", file = opt$out, width = w_in, height = h_in) else
  if (capabilities("cairo")) cairo_pdf(opt$out, width = w_in, height = h_in) else
    pdf(opt$out, width = w_in, height = h_in)
print(fig)
invisible(dev.off())
cat("wrote", opt$out, "\n")
