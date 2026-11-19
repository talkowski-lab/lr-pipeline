#!/usr/bin/env Rscript
# Insertion SVs by size range, three panels top to bottom:
#   1. number of SVs per size range
#   2. proportion of SVs in the range by repeat class (build_repeat_classes.py): full-length MEIs (solid) and truncated
#      MEIs (transparent) per family, VNTR in MEIs, STR and VNTR (STR_novel / VNTR_novel = not in the current TR catalog, transparent), each
#      satellite type, LCR22 paralogs, duplications / insertions of MEI sequence, MEI nested Alu, other TR;
#      the remaining non-repeat insertions are the unfilled part of the bar
#   3. proportion of the non-repeat insertions in the range that remap aligned, by remap class (bar height = aligned)
# Writes the per-bin numbers. ggplot2 + patchwork, Arial, base_size 22, Quartz PDF.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

opt <- parse_args(OptionParser(
  usage = "%prog --annotations <insertion_annotations.tsv.gz> --classes <repeat_classes.tsv> --out-pdf <pdf> --out-tsv <tsv>",
  description = "SV count, repeat-class composition and remap alignment of non-repeat insertions by size range.",
  option_list = list(
    make_option("--annotations", help = "Per-insertion table from build_insertion_annotations.R (remap_class)"),
    make_option("--classes", help = "Per-insertion table from build_repeat_classes.py"),
    make_option("--breaks", default = "50,100,250,500,1000,2500,5000,10000", help = "Size range lower edges, bp [default %default]"),
    make_option("--min-sat", dest = "min_sat", type = "integer", default = 5,
                help = "Satellite types with fewer insertions are grouped as 'other satellite' [default %default]"),
    make_option("--title", default = "", help = "Title prefix [default none]"),
    make_option("--out-pdf", dest = "out_pdf", default = NULL, help = "Optional output PDF (omit to write the table only)"),
    make_option("--out-tsv", dest = "out_tsv", help = "Output per-bin TSV")
  )
))
for (a in c("annotations", "classes", "out_tsv")) if (is.null(opt[[a]])) stop("--", a, " is required")

x <- merge(fread(opt$classes, na.strings = c("", "NA"))[, .(ID, len, repeat_class)],
           fread(opt$annotations, na.strings = c("", "NA"))[, .(ID, remap_class)], by = "ID")
br <- as.numeric(strsplit(opt$breaks, ",")[[1]])
fmt <- function(v) format(v, big.mark = ",", trim = TRUE, scientific = FALSE)
lab <- c(sprintf("%s-%s", fmt(br[-length(br)]), fmt(br[-1] - 1)), sprintf(">=%s", fmt(br[length(br)])))
x[, size_bin := cut(len, c(br, Inf), right = FALSE, labels = lab)]
x <- x[!is.na(size_bin)]

sat_n <- x[repeat_class %like% "^satellite: ", .N, by = repeat_class]
rare <- sat_n[N < opt$min_sat, repeat_class]
x[, cls := fifelse(repeat_class %in% rare, "satellite: other", repeat_class)]
# MEI-family sequence without insertion hallmarks, grouped for the plot
x[cls %like% "^tandem duplication \\(in ", cls := "tandem duplication of MEI sequence"]
x[cls %like% "^interspersed duplication \\(in ", cls := "interspersed duplication of MEI sequence"]
x[cls %like% " insertion into | insertion, no local source", cls := "Alu inserted into Alu / other TE"]
sat_lv <- c(sat_n[N >= opt$min_sat][order(-N), repeat_class], if (length(rare)) "satellite: other")
lv <- c("MEI full-length Alu", "MEI truncated Alu", "MEI full-length L1", "MEI truncated L1", "MEI full-length SVA", "MEI truncated SVA",
        "VNTR in MEI", "STR", "STR_novel", "VNTR", "VNTR_novel", sat_lv, "segmental duplication paralog", "tandem duplication of MEI sequence",
        "interspersed duplication of MEI sequence", "Alu inserted into Alu / other TE", "MEI nested Alu", "other TR")
base_cols <- c("MEI full-length Alu" = "#E31A1C", "MEI truncated Alu" = "#E31A1C", "MEI full-length L1" = "#33A02C",
               "MEI truncated L1" = "#33A02C", "MEI full-length SVA" = "#6A3D9A", "MEI truncated SVA" = "#6A3D9A",
               "VNTR in MEI" = "#8C564B", STR = "#1F78B4", VNTR = "#A6CEE3", "MEI nested Alu" = "#C49C94", "other TR" = "#BDBDBD", "STR_novel" = "#1F78B4", "VNTR_novel" = "#A6CEE3",
               "segmental duplication paralog" = "#000000", "tandem duplication of MEI sequence" = "#F781BF",
               "interspersed duplication of MEI sequence" = "#984EA3", "Alu inserted into Alu / other TE" = "#FB6A4A")
sat_pal <- c("#FF7F00", "#FFD92F", "#66C2A5", "#1B9E77", "#E7298A", "#A6761D", "#17BECF", "#D95F02", "#636363")
cols <- c(base_cols, setNames(sat_pal[seq_along(sat_lv)], sat_lv))[lv]
alph <- setNames(fifelse(lv %like% "truncated|_novel", 0.4, 1), lv)

n <- x[, .(n = .N), by = size_bin]
r <- x[cls %in% lv, .N, by = .(size_bin, cls)]
r <- merge(r, n, by = "size_bin")[, prop := N / n]
r[, cls := factor(cls, levels = rev(lv))]
rtot <- r[, .(p = sum(prop)), by = size_bin]
nr <- x[repeat_class == "non-repeat (next stage)"]
a <- nr[remap_class != "novel", .N, by = .(size_bin, remap_class)]
nn <- nr[, .(n_nonrep = .N, aligned = sum(remap_class != "novel")), by = size_bin]
a <- merge(a, nn, by = "size_bin")[, prop := N / n_nonrep]

tab <- Reduce(function(p, q) merge(p, q, by = "size_bin", all = TRUE),
              list(n, dcast(r, size_bin ~ cls, value.var = "N", fill = 0), nn,
                   dcast(a, size_bin ~ paste0("aligned_", remap_class), value.var = "N", fill = 0)))
for (k in names(tab)[-1]) set(tab, which(is.na(tab[[k]])), k, 0L)
tab[, repeat_frac := round(rtot[match(tab$size_bin, rtot$size_bin), p], 4)]
tab[, aligned_frac_nonrepeat := round(aligned / n_nonrep, 4)]
tab <- tab[order(size_bin)]
fwrite(tab, opt$out_tsv, sep = "\t")
if (is.null(opt$out_pdf)) quit(save = "no")  # table only (no --out-pdf)

th <- theme_bw(base_size = 22, base_family = "Arial") +
  theme(panel.grid = element_blank(), legend.title = element_blank(), axis.title.x = element_blank())
lbl <- function(d, y, v) geom_text(data = d, aes(x = size_bin, y = {{ y }}, label = {{ v }}), vjust = -0.3, size = 7.76 * 0.75,
                                   family = "Arial", inherit.aes = FALSE)
cls_cols <- c(tandem = "#1F78B4", tandem_inverted = "#A6CEE3", partial = "#FF7F00", interspersed = "#B2DF8A")
a[, remap_class := factor(remap_class, levels = names(cls_cols))]
atot <- nn[, .(size_bin, p = aligned / n_nonrep, lab = sprintf("%.0f%% of %s", 100 * aligned / n_nonrep, fmt(n_nonrep)))]

p1 <- ggplot(n, aes(size_bin, n)) + geom_col(fill = "grey45", width = 0.8) + lbl(n, n, fmt(n)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = fmt) +
  labs(y = "SVs", title = sprintf("%s1. Insertion SVs per size range (n = %s)", opt$title, fmt(nrow(x)))) + th
p2 <- ggplot(r, aes(size_bin, prop, fill = cls, alpha = cls)) + geom_col(width = 0.8) +
  scale_fill_manual(values = cols, breaks = lv, drop = FALSE) + scale_alpha_manual(values = alph, breaks = lv, drop = FALSE) +
  lbl(rtot, p, sprintf("%.0f%%", 100 * p)) + scale_x_discrete(drop = FALSE, limits = levels(x$size_bin)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  guides(fill = guide_legend(ncol = 2, override.aes = list(alpha = alph)), alpha = "none") +
  labs(y = "SVs by repeat class",
       title = sprintf("%s2. Repeat classes (truncated MEIs and _novel = not in the TR catalog transparent; unfilled = non-repeat)", opt$title)) + th
p3 <- ggplot(a, aes(size_bin, prop, fill = remap_class)) + geom_col(width = 0.8) + scale_fill_manual(values = cls_cols, drop = FALSE) +
  lbl(atot, p, lab) + scale_x_discrete(drop = FALSE, limits = levels(x$size_bin)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  labs(x = "Insertion size (bp)", y = "Non-repeat SVs aligned",
       title = sprintf("%s3. Non-repeat insertions aligned by remap (by remap class)", opt$title)) +
  th + theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.title.x = element_text())
quartz(type = "pdf", file = opt$out_pdf, width = 24, height = 24)
print(p1 / p2 / p3 + plot_layout(heights = c(0.8, 1.3, 1)))
invisible(dev.off())
print(tab)
cat("wrote", opt$out_pdf, "\n")
