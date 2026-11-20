#!/usr/bin/env Rscript
# Join truvari anno remap output (AnnotateTruvariRemap TSV) with VCF INFO fields and
# summarize remap classes by insertion size, allele frequency, genomic region, TR and MEI status.
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(
  usage = "%prog --remap <remap.tsv> --info <info.tsv.gz> --out-prefix <prefix>",
  description = "Join remap classes with INFO fields and write summary tables.",
  option_list = list(
    make_option("--remap", help = "AnnotateTruvariRemap TSV (no header: CHROM POS REF ALT ID class coords ori perc)"),
    make_option("--info", help = "TSV with header: ID FILTER SVLEN AF AC SOURCE REGION TR_ENVELOPED SUB_FAMILY ORIGIN"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix (writes <prefix>.joined.tsv.gz and <prefix>.summary_*.tsv)")
  )
))
if (is.null(opt$remap) || is.null(opt$info) || is.null(opt$out_prefix)) stop("--remap, --info and --out-prefix are required")

remap <- fread(opt$remap, header = FALSE, sep = "\t", quote = "",
               col.names = c("CHROM", "POS", "REF", "ALT", "ID", "class", "coords", "ori", "perc"))
info <- fread(opt$info, sep = "\t", quote = "", na.strings = ".")
stopifnot(!anyDuplicated(remap$ID), !anyDuplicated(info$ID))

d <- merge(remap[, .(ID, CHROM, POS, class, coords, ori, perc)], info, by = "ID", all = TRUE)
cat(sprintf("remap rows: %d, info rows: %d, joined: %d, missing class: %d, missing info: %d\n",
            nrow(remap), nrow(info), nrow(d), sum(is.na(d$class)), sum(is.na(d$SVLEN))))

d[, len := abs(SVLEN)]
d[, size_bin := cut(len, c(50, 100, 300, 1000, 5000, 10000, Inf), right = FALSE,
                    labels = c("50-99", "100-299", "300-999", "1-5kb", "5-10kb", ">=10kb"))]
d[, af_bin := cut(AF, c(0, 0.01, 0.05, 0.5, 1.0001), right = FALSE,
                  labels = c("<1%", "1-5%", "5-50%", ">=50%"))]
d[, tr := fifelse(is.na(TR_ENVELOPED), "not TR-enveloped", "TR-enveloped")]
d[, mei := fifelse(is.na(SUB_FAMILY), "no MEI call", sub("^(Alu|L1|SVA).*", "\\1", SUB_FAMILY))]
d[, class := factor(class, c("novel", "partial", "tandem", "tandem_inverted", "tandem_complex", "interspersed"))]

fwrite(d, paste0(opt$out_prefix, ".joined.tsv.gz"), sep = "\t")

tab <- function(by) {
  x <- d[, .N, by = c(by, "class")]
  x[, frac := round(N / sum(N), 4), by = by]
  setorderv(x, c(by, "class"))
  x
}
for (by in c("size_bin", "af_bin", "REGION", "tr", "mei", "FILTER")) {
  fwrite(tab(by), sprintf("%s.summary_%s.tsv", opt$out_prefix, by), sep = "\t")
}
overall <- d[, .N, by = class][, frac := round(N / sum(N), 4)][order(class)]
fwrite(overall, paste0(opt$out_prefix, ".summary_overall.tsv"), sep = "\t")
print(overall)
