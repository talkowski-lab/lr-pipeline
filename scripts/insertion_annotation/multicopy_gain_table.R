#!/usr/bin/env Rscript
# For multi-copy local tandem gains (classify_partial_structure.R), report the reference region the copies come
# from and an estimated copy number, and group insertions whose regions overlap into loci (alleles of one site).
#   aligned_region = union of local hits on the reference (identity >= --min-ident; the largest merged block)
#   copies_by_len  = insertion length / aligned region length
#   copies_by_hits = sum of local hit query spans / aligned region length (the classify script's "reuse")
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})
opt <- parse_args(OptionParser(
  usage = "%prog --structure <partial_structure.tsv> --hits <step3_hits.annotated.tsv.gz> --out <tsv>",
  description = "Aligned region and estimated copy number of multi-copy local tandem gains, grouped into loci.",
  option_list = list(
    make_option("--structure", help = "Per-insertion table from classify_partial_structure.R"),
    make_option("--hits", help = "Annotated hit table from annotate_step3_hits.R"),
    make_option("--class", default = "multi-copy local tandem gain", help = "structure value [default %default]"),
    make_option("--min-ident", dest = "min_ident", type = "double", default = 0.9, help = "Min hit identity [default %default]"),
    make_option("--local-bp", dest = "local_bp", type = "integer", default = 10000, help = "Local window beyond the insertion length [default %default]"),
    make_option("--out", help = "Output TSV")
  )
))
for (a in c("structure", "hits", "out")) if (is.null(opt[[a]])) stop("--", a, " is required")
x <- fread(opt$structure, na.strings = c("", "NA"))[structure == opt$class]
x[, `:=`(CHROM = sub("-.*", "", ID), POS = as.numeric(sub("^[^-]+-([0-9]+)-.*", "\\1", ID)))]
h <- fread(opt$hits)[ID %in% x$ID & identity >= opt$min_ident]
h[, `:=`(rchrom = sub(":.*", "", coords), rs = as.numeric(sub("^[^:]+:([0-9]+)-.*", "\\1", coords)), re = as.numeric(sub(".*-", "", coords)))]
h <- merge(h, x[, .(ID, CHROM, POS, len)], by = "ID")
h <- h[rchrom == CHROM & re >= POS - len - opt$local_bp & rs <= POS + len + opt$local_bp]
blocks <- h[order(ID, rs)][, {
  grp <- cumsum(c(TRUE, rs[-1] > cummax(re)[-.N] + 1))
  .SD[, .(bs = min(rs), be = max(re), qsum = sum(q_end - q_start + 1)), by = .(grp = grp)]
}, by = ID]
b <- blocks[order(ID, -(be - bs))][, .SD[1], by = ID]
out <- merge(x[, .(ID, CHROM, POS, len, vcf, svan_type, AF)], b[, .(ID, region_start = bs, region_end = be, qsum)], by = "ID")
out[, region_len := region_end - region_start + 1]
out[, aligned_region := sprintf("%s:%d-%d", CHROM, region_start, region_end)]
out[, copies_by_len := round(len / region_len, 2)]
out[, copies_by_hits := round(qsum / region_len, 2)]
setorder(out, CHROM, region_start)
out[, locus := cumsum(c(TRUE, region_start[-1] > cummax(region_end)[-.N] | CHROM[-1] != CHROM[-.N]))]
out[, `:=`(locus_n_alleles = .N, locus_copy_range = sprintf("%.1f-%.1f", min(copies_by_len), max(copies_by_len))), by = locus]
fwrite(out[, .(locus, locus_n_alleles, locus_copy_range, ID, len, aligned_region, region_len, copies_by_len, copies_by_hits, vcf, svan_type, AF)],
       opt$out, sep = "\t")
cat(sprintf("insertions: %d; loci: %d; loci with >= 2 alleles: %d (%d insertions)\n", nrow(out), uniqueN(out$locus),
            out[locus_n_alleles >= 2, uniqueN(locus)], out[locus_n_alleles >= 2, .N]))
