#!/usr/bin/env Rscript
# Copy-number annotation of insertions that are extra copies of a nearby reference segment (mCNVs, and tandem copies
# of IGL / other known loci). Per insertion, from the remap step-3 hits on the insertion contig within the insertion
# length + --local-bp of the site (identity >= --min-ident):
#   segment          largest merged block of those hits on the reference = the duplicated unit
#   copies_added     insertion length / segment length
#   total_copies     round(1 + copies_added) = copies of the segment on the carrier haplotype (hg38 has 1)
#   ip_offset / ip_position / ip_where   insertion point relative to the segment, as in call_tandem_dups.R
#   known_locus      the smallest named locus in --loci overlapping the segment (e.g. IGL, C4, CYP2D6/CYP2D7, SMN1/SMN2)
#   cnv_call         CNV_<total_copies>copy[_<known locus>], e.g. CNV_3copy_IGL, when local hits cover >= --min-local-cov
#                    of the insertion; otherwise "partly local (copies uncertain, not called)" or "no local copy"
#   cnv_detail       segment, length, copies added, insertion point, local coverage
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})
opt <- parse_args(OptionParser(
  usage = "%prog --ids <ids.txt> --hits <step3_hits.tsv.gz> --loci <known_cnv_loci.bed> --out <tsv>",
  description = "Annotate copy number, duplicated segment and known locus of multi-copy insertions.",
  option_list = list(
    make_option("--ids", help = "File of insertion IDs"),
    make_option("--hits", help = "Per-hit table from replay_remap_step3.R"),
    make_option("--loci", help = "BED of named loci (fetch_immune_loci.py --cnv-genes)"),
    make_option("--min-ident", dest = "min_ident", type = "double", default = 0.9, help = "Min hit identity [default %default]"),
    make_option("--min-local-cov", dest = "min_local_cov", type = "double", default = 0.7,
                help = "Min fraction of the insertion covered by local hits to call a CNV [default %default]"),
    make_option("--local-bp", dest = "local_bp", type = "integer", default = 10000, help = "Window beyond the insertion length [default %default]"),
    make_option("--out", help = "Output TSV")
  )
))
for (a in c("ids", "hits", "loci", "out")) if (is.null(opt[[a]])) stop("--", a, " is required")
ids <- readLines(opt$ids)
L <- fread(opt$loci, header = FALSE, select = 1:4, col.names = c("chrom", "s", "e", "locus"))
x <- data.table(ID = ids, CHROM = sub("-.*", "", ids), POS = as.numeric(sub("^[^-]+-([0-9]+)-.*", "\\1", ids)),
                len = as.numeric(sub(".*-", "", ids)))
h <- merge(fread(opt$hits)[ID %in% ids & identity >= opt$min_ident], x, by = "ID")
h[, `:=`(rchrom = sub(":.*", "", coords), rs = as.numeric(sub("^[^:]+:([0-9]+)-.*", "\\1", coords)), re = as.numeric(sub(".*-", "", coords)))]
h <- h[rchrom == CHROM & re >= POS - len - opt$local_bp & rs <= POS + len + opt$local_bp]
union_len <- function(s, e) {
  o <- order(s); s <- s[o]; e <- e[o]; tot <- 0; cs <- s[1]; ce <- e[1]
  for (i in seq_along(s)[-1]) { if (s[i] > ce + 1) { tot <- tot + ce - cs + 1; cs <- s[i]; ce <- e[i] } else ce <- max(ce, e[i]) }
  tot + ce - cs + 1
}
seg <- h[order(ID, rs)][, {
  grp <- cumsum(c(TRUE, rs[-1] > cummax(re)[-.N] + 1))
  b <- data.table(rs, re, grp)[, .(bs = min(rs), be = max(re)), by = grp][order(-(be - bs))][1]
  .(seg_start = b$bs, seg_end = b$be, local_cov = round(min(1, union_len(q_start, q_end) / qlen[1]), 3))
}, by = ID]
x <- merge(x, seg, by = "ID", all.x = TRUE)
x[, seg_len := seg_end - seg_start + 1]
x[, copies_added := round(len / seg_len, 2)]
x[, total_copies := round(1 + copies_added)]
x[, ip_offset := POS - seg_start + 1]
x[, ip_position := round(ip_offset / seg_len, 3)]
x[, ip_where := fcase(is.na(seg_len), NA_character_, ip_position <= 0.05, "at / left of the segment's left edge",
                      ip_position >= 0.95, "at / right of the segment's right edge", default = "inside the segment")]
x[, known_locus := vapply(seq_len(.N), function(i) {
  if (is.na(seg_start[i])) return(NA_character_)
  k <- L[chrom == CHROM[i] & s < seg_end[i] & e >= seg_start[i]][order(e - s)]
  if (nrow(k)) gsub("[ /()-]+", "_", sub("_$", "", k$locus[1])) else NA_character_
}, "")]
x[, known_locus := sub("_$", "", known_locus)]
x[, cnv_call := fcase(is.na(seg_len), "no local copy (not a CNV)",
                      local_cov < opt$min_local_cov, "partly local (copies uncertain, not called)",
                      default = paste0("CNV_", total_copies, "copy", fifelse(is.na(known_locus), "", paste0("_", known_locus))))]
x[, cnv_detail := sprintf("segment %s:%s-%s (%s bp), +%s copies (%s on the haplotype), insertion point %s (%.2f), local coverage %.2f",
                          CHROM, format(seg_start, big.mark = ","), format(seg_end, big.mark = ","), format(seg_len, big.mark = ","),
                          copies_added, total_copies, ip_where, ip_position, local_cov)]
setorder(x, CHROM, seg_start)
fwrite(x[, .(ID, CHROM, POS, len, cnv_call, known_locus, seg_start, seg_end, seg_len, copies_added, total_copies, ip_offset,
             ip_position, ip_where, local_cov, cnv_detail)], opt$out, sep = "\t")
print(x[, .(ID, len, cnv_call, seg_len, copies_added, ip_position, local_cov)])
