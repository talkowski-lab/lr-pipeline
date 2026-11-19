#!/usr/bin/env Rscript
# Compare SVAN insertion annotations (ITYPE_N, FAM_N, DUP_COORD, SOURCE_COORD) with truvari remap
# classes for the same variants, matched by variant ID. Writes the joined table and cross-tabs.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(
  usage = "%prog --remap <remap.tsv> --svan <svan.tsv> --svan-header <header.txt> --joined <summarize_remap joined.tsv.gz> --out-prefix <prefix>",
  description = "Cross-tabulate SVAN insertion types against truvari remap classes.",
  option_list = list(
    make_option("--remap", help = "AnnotateTruvariRemap TSV (no header)"),
    make_option("--svan", help = "AnnotateSVAN TSV (no header: chr pos ref alt name + header fields)"),
    make_option("--svan-header", dest = "svan_header", help = "AnnotateSVAN header file (one field name per line)"),
    make_option("--joined", help = "Joined table from summarize_remap.R (adds len, size_bin, REGION, AF)"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("remap", "svan", "svan_header", "joined", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")

fields <- readLines(opt$svan_header)
svan <- fread(opt$svan, header = FALSE, sep = "\t", quote = "", na.strings = ".",
              col.names = c("chr", "pos", "ref", "alt", "ID", fields),
              select = NULL)
svan <- svan[, .(ID, ITYPE_N, FAM_N, CONFORMATION, STRAND, DUP_COORD, SOURCE_COORD, SRC_GENE, MT_COORD, VNTR_COORD, PERC_RESOLVED)]
stopifnot(!anyDuplicated(svan$ID))

j <- fread(opt$joined)[, .(ID, class, coords, ori, perc, len, size_bin, REGION, AF, SUB_FAMILY)]
d <- merge(j, svan, by = "ID", all.x = TRUE)
d[, in_svan := ID %in% svan$ID]
cat(sprintf("remap insertions: %d; found in SVAN table: %d; missing: %d\n", nrow(d), sum(d$in_svan), sum(!d$in_svan)))

# SVAN insertion type, collapsed to the categories it reports
d[, svan_type := fcase(
  !in_svan, "not in SVAN table",
  is.na(ITYPE_N), "unannotated",
  default = ITYPE_N)]
d[, svan_mei_fam := fifelse(!is.na(FAM_N), FAM_N, "none")]

# Hit chromosome relative to the insertion, from remap_coords
d[, hit_chr := fifelse(is.na(coords) | coords == ".", NA_character_, sub(":.*", "", coords))]
d[, hit_where := fcase(is.na(hit_chr), "no hit",
                      hit_chr == "chr22", "same chrom",
                      grepl("_", hit_chr), "unplaced/random",
                      default = "other chrom")]

# SVAN groups used for the agreement plot
d[, svan_group := fcase(
  svan_type %in% c("solo", "partnered", "orphan", "chimera"), "MEI (solo/partnered/orphan)",
  svan_type == "DUP", "DUP (tandem)",
  svan_type == "DUP_INTERSPERSED", "DUP_INTERSPERSED",
  svan_type == "VNTR", "VNTR",
  svan_type %in% c("NUMT", "PSD", "INV_DUP", "COMPLEX_DUP"), "NUMT/PSD/INV_DUP/COMPLEX_DUP",
  default = svan_type)]

# For SVAN DUP / DUP_INTERSPERSED with a remap hit: does any SVAN DUP_COORD segment start
# within max(len, 1 kb) of the remap hit start on the same contig?
d[, dup_vs_hit := NA_character_]
pos <- function(x) as.numeric(sub("-.*", "", sub("^.*:", "", x)))
idx <- d[, which(svan_type %in% c("DUP", "DUP_INTERSPERSED") & !is.na(DUP_COORD) & hit_where != "no hit")]
d[idx, dup_vs_hit := mapply(function(s, hc, hs, l) {
  s <- strsplit(s, ",")[[1]]
  ch <- sub(":.*", "", s)
  if (any(ch == hc & abs(pos(s) - hs) <= max(l, 1000))) "same locus"
  else if (any(ch == hc)) "same contig, other locus"
  else "different contig"
}, DUP_COORD, hit_chr, pos(coords), len)]

fwrite(d, paste0(opt$out_prefix, ".svan_remap.joined.tsv.gz"), sep = "\t")

xt <- function(rows, cols, name) {
  x <- dcast(d[, .N, by = c(rows, cols)], as.formula(paste(rows, "~", cols)), value.var = "N", fill = 0)
  x[, total := rowSums(.SD), .SDcols = -1]
  fwrite(x, sprintf("%s.xtab_%s.tsv", opt$out_prefix, name), sep = "\t")
  cat("\n==", name, "\n"); print(x)
}
xt("svan_type", "class", "svan_type_by_remap_class")
xt("svan_mei_fam", "class", "svan_family_by_remap_class")
xt("svan_type", "hit_where", "svan_type_by_hit_location")
xt("svan_type", "size_bin", "svan_type_by_size")
d[, big := fifelse(len >= 1000, "ge1kb", "lt1kb")]
for (b in c("ge1kb", "lt1kb")) {
  x <- dcast(d[big == b, .N, by = .(svan_group, class)], svan_group ~ class, value.var = "N", fill = 0)
  x[, total := rowSums(.SD), .SDcols = -1]
  fwrite(x, sprintf("%s.xtab_svan_group_by_remap_class.%s.tsv", opt$out_prefix, b), sep = "\t")
  cat("\n== svan_group_by_remap_class", b, "\n"); print(x)
}
x <- d[!is.na(dup_vs_hit), .N, by = .(svan_type, class, dup_vs_hit)][order(svan_type, class, dup_vs_hit)]
fwrite(x, paste0(opt$out_prefix, ".dup_coord_vs_remap_hit.tsv"), sep = "\t")
cat("\n== SVAN DUP_COORD vs remap hit\n"); print(x)
