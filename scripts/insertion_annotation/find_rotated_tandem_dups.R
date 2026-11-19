#!/usr/bin/env Rscript
# Find "rotated" tandem duplications among remap partial calls: the VCF places the insertion point
# inside the duplicated unit, so the inserted string is [part right of the site][part left of the site]
# and a linear alignment splits it in two at the site. Detected from the AlignRemapQueries PAF:
# on the insertion's contig, same strand (+), one primary segment starts the query and maps starting at
# the site, another ends the query and maps ending at the site, and the two meet on the query.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(
  usage = "%prog --paf <alignments.paf.gz> --info <info.tsv.gz> --sources <sources_per_insertion.tsv.gz> --out-prefix <prefix>",
  description = "Flag remap partial insertions that are complete tandem duplications split at the insertion point.",
  option_list = list(
    make_option("--paf", help = "PAF from AlignRemapQueries.wdl"),
    make_option("--info", help = "INFO table from extract_ins_info.sh (allele_subtype, TRID)"),
    make_option("--sources", help = "Per-insertion table from count_alignment_sources.R (svan_type, structure)"),
    make_option("--tol", type = "integer", default = 50, help = "bp tolerance at query ends and at the site [default %default]"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("paf", "info", "sources", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")
tol <- opt$tol

awk_cmd <- paste("gzip -dc", shQuote(opt$paf),
                 "| awk -F'\\t' -v OFS='\\t' '{tp=\"NA\"; for (i=13;i<=NF;i++) if ($i ~ /^tp:A:/) tp=substr($i,6); print $1,$2,$3,$4,$5,$6,$8,$9,$12,tp}'")
p <- fread(cmd = awk_cmd, header = FALSE, sep = "\t",
           col.names = c("ID", "qlen", "qs", "qe", "strand", "tname", "ts", "te", "mapq", "tp"))[tp == "P"]
p[, ins_chr := sub("-.*", "", ID)]
p[, ins_pos := as.numeric(sub("^[^-]+-([0-9]+)-.*", "\\1", ID))]

r <- p[, {
  a <- which(strand == "+" & tname == ins_chr & qs <= tol & abs(ts - ins_pos) <= tol)
  b <- which(strand == "+" & tname == ins_chr & qe >= qlen - tol & abs(te - ins_pos) <= tol)
  ok <- length(a) > 0 && length(b) > 0 && max(qe[a]) >= min(qs[b]) - tol
  .(rotated = ok, n_seg = .N, qlen = as.numeric(qlen[1]),
    unit_start = if (ok) as.numeric(min(ts[b]) + 1) else NA_real_,
    unit_end = if (ok) as.numeric(max(te[a])) else NA_real_)
}, by = ID]
r[, unit_len := unit_end - unit_start + 1]
r[, len_over_unit := round(qlen / unit_len, 3)]

info <- fread(opt$info, sep = "\t", quote = "", na.strings = ".")[, .(ID, FILTER, AF, REGION, allele_subtype, TRID)]
src <- fread(opt$sources)[, .(ID, svan_type, structure)]
r <- merge(merge(r, info, by = "ID", all.x = TRUE), src, by = "ID", all.x = TRUE)
fwrite(r, paste0(opt$out_prefix, ".rotated_tandem_dups.tsv"), sep = "\t")

cat(sprintf("partials: %d; rotated tandem dups: %d\n", nrow(r), sum(r$rotated)))
out <- function(x, name) { fwrite(x, sprintf("%s.%s.tsv", opt$out_prefix, name), sep = "\t"); cat("\n==", name, "\n"); print(x) }
out(r[rotated == TRUE, .N, by = .(allele_subtype, svan_type)][order(-N)], "rotated_by_vcf_svan")
out(r[rotated == TRUE, .N, by = .(in_TR = !is.na(TRID))], "rotated_by_tr")
out(r[rotated == TRUE, .(n = .N, med_len_over_unit = median(len_over_unit), within_5pct = sum(abs(len_over_unit - 1) <= 0.05))], "rotated_unit_length")
