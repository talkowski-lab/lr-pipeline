#!/usr/bin/env Rscript
# Describe how remap "partial" insertions (no single hit >= 80%) are built, for one ins_call group.
#
# Per insertion, from all its hits (identity >= --min-ident; n_aln = number of such hits):
#   all_cov     fraction of the insertion covered by any hit
#   best_perc   % of the insertion covered by the best-scoring hit
#   local       hits on the insertion contig within --local-bp of the site
#   local_cov   fraction of the insertion covered by local hits
#   local_ref   union bp of the local hits on the reference
#   reuse       (sum of local hit query spans) / local_ref: > 1 means the same reference segment appears
#               several times in the insertion (multi-copy gain)
# structure (first match wins):
#   multi-copy local tandem gain   local_cov >= 0.8 and reuse >= --min-reuse
#   local dup + extra local seq    local_cov >= 0.8 (pieces of the flank, reused < --min-reuse)
#   local dup + non-local seq      local_cov < 0.8 and all_cov >= 0.8
#   local dup + unaligned seq      local_cov >= 0.3 and all_cov < 0.8
#   distal / mostly unaligned      otherwise
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(
  usage = "%prog --annotations <insertion_annotations.tsv.gz> --hits <step3_hits.annotated.tsv.gz> --call 'no repeat' --out-prefix <prefix>",
  description = "Classify the structure of remap partial insertions in one ins_call group.",
  option_list = list(
    make_option("--annotations", help = "Per-insertion table from build_insertion_annotations.R"),
    make_option("--hits", help = "Annotated hit table from annotate_step3_hits.R"),
    make_option("--call", default = "no repeat", help = "ins_call group [default %default]"),
    make_option("--ids", default = NULL, help = "Optional file of insertion IDs to analyse instead of --call / --remap-class"),
    make_option("--remap-class", dest = "remap_class", default = "partial", help = "Remap class [default %default]"),
    make_option("--min-ident", dest = "min_ident", type = "double", default = 0.9, help = "Min hit identity [default %default]"),
    make_option("--local-bp", dest = "local_bp", type = "integer", default = 10000, help = "Local window beyond the insertion length [default %default]"),
    make_option("--min-reuse", dest = "min_reuse", type = "double", default = 1.5, help = "Reuse for a multi-copy gain [default %default]"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("annotations", "hits", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")

x <- fread(opt$annotations, na.strings = c("", "NA"))
x <- if (!is.null(opt$ids)) x[ID %in% readLines(opt$ids)] else x[ins_call == opt$call & remap_class == opt$remap_class]
h <- fread(opt$hits, na.strings = c("", "NA"))[ID %in% x$ID & identity >= opt$min_ident]
h[, `:=`(rchrom = sub(":.*", "", coords), rs = as.numeric(sub("^[^:]+:([0-9]+)-.*", "\\1", coords)), re = as.numeric(sub(".*-", "", coords)))]
h <- merge(h, x[, .(ID, CHROM, POS, alt_len)], by = "ID")
h[, local := rchrom == CHROM & re >= POS - alt_len - opt$local_bp & rs <= POS + alt_len + opt$local_bp]

union_len <- function(s, e) {
  if (!length(s)) return(0)
  o <- order(s); s <- s[o]; e <- e[o]; tot <- 0; cs <- s[1]; ce <- e[1]
  for (i in seq_along(s)[-1]) { if (s[i] > ce + 1) { tot <- tot + ce - cs + 1; cs <- s[i]; ce <- e[i] } else ce <- max(ce, e[i]) }
  tot + ce - cs + 1
}
s <- h[, .(n_aln = .N, n_local = sum(local),
           all_cov = round(union_len(q_start, q_end) / qlen[1], 3),
           best_perc = perc[which.max(score)],
           local_cov = round(union_len(q_start[local], q_end[local]) / qlen[1], 3),
           local_ref = union_len(rs[local], re[local]),
           local_qsum = sum((q_end - q_start + 1)[local]),
           best_hit = coords[which.max(score)]), by = ID]
s[, reuse := round(fifelse(local_ref > 0, local_qsum / local_ref, 0), 2)]
x <- merge(x, s, by = "ID", all.x = TRUE)
x[is.na(n_aln), `:=`(n_aln = 0L, n_local = 0L, all_cov = 0, local_cov = 0, reuse = 0)]
x[, structure := fcase(local_cov >= 0.8 & reuse >= opt$min_reuse, "multi-copy local tandem gain",
                       local_cov >= 0.8, "local dup + extra local seq",
                       all_cov >= 0.8, "local dup + non-local seq",
                       local_cov >= 0.3, "local dup + unaligned seq",
                       default = "distal / mostly unaligned")]
x[, len := alt_len - 1L]
x[, size_bin := cut(len, c(0, 1000, 5000, 20000, Inf), right = FALSE, labels = c("<1k", "1-5k", "5-20k", ">=20k"))]
x[, vcf := fcase(!is.na(allele_subtype), allele_subtype, !is.na(TRID), "in TR", default = "plain ins")]
x[, af05 := !is.na(AF) & AF == 0.5]
setorder(x, structure, -len)
fwrite(x[, .(ID, len, structure, vcf, svan_type, REGION, AF, n_aln, n_local, all_cov, local_cov, best_perc, reuse, best_hit,
             rotated_tandem_dup, complex_insertion)], paste0(opt$out_prefix, ".partial_structure.tsv"), sep = "\t")

o <- function(d, name) { fwrite(d, sprintf("%s.%s.tsv", opt$out_prefix, name), sep = "\t"); cat("\n==", name, "\n"); print(d) }
cat(sprintf("%s: %d insertions\n", if (!is.null(opt$ids)) opt$ids else paste(opt$call, "/", opt$remap_class), nrow(x)))
o(x[, .(n = .N, med_len = as.numeric(median(len)), med_reuse = median(reuse), af_0.5 = sum(af05)), by = structure][order(-n)], "partial_structure_summary")
o(dcast(x[, .N, by = .(structure, vcf)], structure ~ vcf, value.var = "N", fill = 0), "partial_structure_by_vcf")
o(dcast(x[, .N, by = .(structure, size_bin)], structure ~ size_bin, value.var = "N", fill = 0), "partial_structure_by_size")
