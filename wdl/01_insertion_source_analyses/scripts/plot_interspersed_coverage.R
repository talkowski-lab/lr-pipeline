#!/usr/bin/env Rscript
# For insertions left after the tandem-dup stage (call_tandem_dups.R, td_call empty), how much of the inserted
# sequence is covered by its interspersed origin = the best-scoring alignment that is NOT adjacent to the site
# (not on the insertion contig within the insertion length of the site). Also the union of all non-adjacent alignments.
# Source location: same chromosome (distance), other chromosome, or unplaced / alt contig.
#   A. histogram of the best non-adjacent alignment's coverage, remap interspersed insertions, coloured by source location
#   B. the same for remap partial insertions
#   C. union of all non-adjacent alignments, both classes
# ggplot2 + patchwork, Arial, base_size 22, Quartz PDF.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})
opt <- parse_args(OptionParser(
  usage = "%prog --tandem <tandem_dups.tsv> --hits <step3_hits.tsv.gz> --out-pdf <pdf> --out-tsv <tsv>",
  description = "Coverage of inserted sequence by the interspersed (non-adjacent) source alignment.",
  option_list = list(
    make_option("--tandem", help = "call_tandem_dups.R per-insertion output"),
    make_option("--hits", help = "Per-hit table from replay_remap_step3.R"),
    make_option("--classes", default = "interspersed,partial,tandem", help = "Remap classes to analyse [default %default]"),
    make_option("--title", default = "", help = "Title prefix"),
    make_option("--adjacent", type = "integer", default = 5, help = "Alignments within this distance of the site are tandem, not interspersed [default %default]"),
    make_option("--out-pdf", dest = "out_pdf", default = NULL, help = "Optional output PDF (omit to write the table only)"),
    make_option("--out-tsv", dest = "out_tsv", help = "Output per-insertion TSV")
  )
))
for (a in c("tandem", "hits", "out_tsv")) if (is.null(opt[[a]])) stop("--", a, " is required")
cls <- strsplit(opt$classes, ",")[[1]]
td <- fread(opt$tandem, na.strings = c("", "NA"))[is.na(td_call) & remap_class %in% cls, .(ID, CHROM, POS, len, remap_class)]
h <- merge(fread(opt$hits)[ID %in% td$ID], td, by = "ID")
h[, `:=`(rchrom = sub(":.*", "", coords), rs = as.numeric(sub("^[^:]+:([0-9]+)-.*", "\\1", coords)), re = as.numeric(sub(".*-", "", coords)))]
h[, site_dist := fifelse(rchrom == CHROM & rs - 1 <= POS & POS <= re, 0, chrom_dist)]
h[, adjacent := rchrom == CHROM & site_dist <= opt$adjacent]
nh <- h[adjacent == FALSE]
union_len <- function(s, e) {
  o <- order(s); s <- s[o]; e <- e[o]; tot <- 0; cs <- s[1]; ce <- e[1]
  for (i in seq_along(s)[-1]) { if (s[i] > ce + 1) { tot <- tot + ce - cs + 1; cs <- s[i]; ce <- e[i] } else ce <- max(ce, e[i]) }
  tot + ce - cs + 1
}
b <- nh[order(ID, -score)][, .SD[1], by = ID][, .(ID, source = coords, source_strand = strand, source_identity = identity, source_mapq = mapq,
                                                    best_cov = pmin(perc, 100) / 100, rchrom, chrom_dist)]
u <- nh[, .(union_cov = min(1, union_len(q_start, q_end) / qlen[1]), n_nonadj = .N), by = ID]
x <- merge(merge(td, b, by = "ID", all.x = TRUE), u, by = "ID", all.x = TRUE)
x[is.na(best_cov), `:=`(best_cov = 0, union_cov = 0, n_nonadj = 0L)]
x[, location := fcase(is.na(rchrom), "no non-adjacent alignment",
                      grepl("_|chrUn|chrM", rchrom), "unplaced / alt contig",
                      rchrom == CHROM, "same chromosome",
                      default = "other chromosome")]
x[, cov_bin := cut(best_cov, c(-0.01, 0.5, 0.8, 0.9, 1.0001), labels = c("< 50%", "50-80%", "80-90%", ">= 90%"))]
fwrite(x[order(remap_class, -best_cov), .(ID, len, remap_class, best_cov = round(best_cov, 3), union_cov = round(union_cov, 3), n_nonadj,
                                          source, source_strand, source_identity, source_mapq, location, chrom_dist)], opt$out_tsv, sep = "\t")
loc_cols <- c("same chromosome" = "#1F78B4", "other chromosome" = "#33A02C", "unplaced / alt contig" = "#FF7F00",
              "no non-adjacent alignment" = "#BDBDBD")
th <- theme_bw(base_size = 22, base_family = "Arial") + theme(panel.grid = element_blank(), legend.title = element_blank())
hist_p <- function(d, v, ttl, xl) ggplot(d, aes({{ v }} * 100, fill = location)) + geom_histogram(breaks = seq(0, 100, 5)) +
  geom_vline(xintercept = c(80, 90), linetype = c("dashed", "solid")) + scale_fill_manual(values = loc_cols, drop = FALSE) +
  labs(x = xl, y = "Insertions", title = ttl) + th
if (is.null(opt$out_pdf)) quit(save = "no")  # table only (no --out-pdf)
pl <- lapply(cls, function(k) hist_p(x[remap_class == k], best_cov, sprintf("%s%s. remap %s (n = %d): best interspersed alignment",
                                                                                opt$title, LETTERS[match(k, cls)], k, x[remap_class == k, .N]),
                                     "% of the inserted sequence covered by the best non-adjacent alignment"))
pl <- c(pl, list(hist_p(x, union_cov, sprintf("%s%s. all classes: all non-adjacent alignments together", opt$title, LETTERS[length(cls) + 1]),
                        "% of the inserted sequence covered by all non-adjacent alignments")))
quartz(type = "pdf", file = opt$out_pdf, width = 20, height = 7 * length(pl))
print(wrap_plots(pl, ncol = 1))
invisible(dev.off())
print(dcast(x[, .N, by = .(remap_class, cov_bin)], remap_class ~ cov_bin, value.var = "N", fill = 0))
print(dcast(x[, .N, by = .(remap_class, location)], remap_class ~ location, value.var = "N", fill = 0))
cat("wrote", opt$out_pdf, "\n")
