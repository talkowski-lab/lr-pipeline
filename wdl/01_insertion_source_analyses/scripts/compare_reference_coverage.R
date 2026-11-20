#!/usr/bin/env Rscript
# Compare how much of each insertion aligns to two references (e.g. hg38 and T2T-CHM13), from two PAFs made
# with the same query sequences (AlignRemapQueries.wdl). Per insertion and reference, using alignments with
# identity >= --min-ident and query span >= --min-bp:
#   cov        fraction of the insertion covered by any alignment
#   best_perc  % of the insertion covered by the best-scoring alignment
#   best_hit   reference interval and strand of that alignment, with identity and mapq
#   contigs    reference contigs hit, by covered query bp
# plus the gain in coverage on the second reference, and the second-reference intervals covering query
# segments that do not align to the first reference ("new" segments >= --min-bp).
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(
  usage = "%prog --paf-a <hg38.paf.gz> --paf-b <t2t.paf.gz> --ids <ids.txt> --out-prefix <prefix>",
  description = "Compare per-insertion alignment coverage between two references.",
  option_list = list(
    make_option("--paf-a", dest = "paf_a", help = "PAF against reference A (e.g. hg38)"),
    make_option("--paf-b", dest = "paf_b", help = "PAF against reference B (e.g. T2T)"),
    make_option("--name-a", dest = "name_a", default = "hg38", help = "Label of reference A [default %default]"),
    make_option("--name-b", dest = "name_b", default = "T2T", help = "Label of reference B [default %default]"),
    make_option("--ids", help = "File with one insertion ID per line"),
    make_option("--labels", default = NULL, help = "Optional per-insertion TSV with SVID or ID column, joined to the output"),
    make_option("--min-ident", dest = "min_ident", type = "double", default = 0.9, help = "Min alignment identity [default %default]"),
    make_option("--min-bp", dest = "min_bp", type = "integer", default = 100, help = "Min alignment query span / new segment, bp [default %default]"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("paf_a", "paf_b", "ids", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")

ids <- readLines(opt$ids)
read_paf <- function(f) {
  cmd <- paste("gzip -dc", shQuote(f), "| awk -F'\\t' -v OFS='\\t' '{as=\"NA\"; for (i=13;i<=NF;i++) if ($i ~ /^AS:i:/) as=substr($i,6); print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,as}'")
  p <- fread(cmd = cmd, header = FALSE, sep = "\t",
             col.names = c("ID", "qlen", "qs", "qe", "strand", "tname", "ts", "te", "nmatch", "alen", "mapq", "score"))
  p[ID %in% ids & tname != "*" & nmatch / alen >= opt$min_ident & qe - qs >= opt$min_bp]
}
segs <- function(s, e) {
  if (!length(s)) return(data.table(s = numeric(), e = numeric()))
  o <- order(s); s <- s[o]; e <- e[o]; out <- list(); cs <- s[1]; ce <- e[1]
  for (i in seq_along(s)[-1]) { if (s[i] > ce) { out[[length(out) + 1]] <- c(cs, ce); cs <- s[i]; ce <- e[i] } else ce <- max(ce, e[i]) }
  out[[length(out) + 1]] <- c(cs, ce)
  data.table(s = sapply(out, `[`, 1), e = sapply(out, `[`, 2))
}
summ <- function(p, nm) {
  p[, .(cov = round(sum(segs(qs, qe)[, e - s]) / qlen[1], 3),
        best_perc = round(100 * (qe - qs)[which.max(score)] / qlen[1], 1),
        best_hit = sprintf("%s:%d-%d(%s,%.1f%%,mq%d)", tname, ts + 1, te, strand, 100 * nmatch / alen, mapq)[which.max(score)],
        n_aln = .N,
        contigs = { b <- p[ID == .BY$ID, .(bp = sum(qe - qs)), by = tname][order(-bp)]; paste(sprintf("%s:%d", b$tname, b$bp), collapse = ",") }),
    by = ID][, setnames(.SD, setdiff(names(.SD), "ID"), paste0(nm, "_", setdiff(names(.SD), "ID")))]
}
a <- read_paf(opt$paf_a); b <- read_paf(opt$paf_b)
x <- merge(merge(data.table(ID = ids), summ(a, opt$name_a), by = "ID", all.x = TRUE), summ(b, opt$name_b), by = "ID", all.x = TRUE)
ca <- paste0(opt$name_a, "_cov"); cb <- paste0(opt$name_b, "_cov")
for (k in c(ca, cb)) set(x, which(is.na(x[[k]])), k, 0)
x[, cov_gain := round(get(cb) - get(ca), 3)]

# reference-B intervals for query segments not covered on reference A
new <- rbindlist(lapply(ids, function(i) {
  sa <- segs(a[ID == i, qs], a[ID == i, qe]); pb <- b[ID == i]
  if (!nrow(pb)) return(NULL)
  pb[, new_bp := mapply(function(s, e) { cov <- if (nrow(sa)) sum(pmax(0, pmin(e, sa$e) - pmax(s, sa$s))) else 0; (e - s) - cov }, qs, qe)]
  pb[new_bp >= opt$min_bp, .(ID, q_start = qs + 1, q_end = qe, new_bp, ref = sprintf("%s:%d-%d", tname, ts + 1, te), strand,
                             identity = round(nmatch / alen, 4), mapq, score)]
}))
if (nrow(new)) {
  x <- merge(x, new[order(-new_bp), .(new_segments_on_b = paste(sprintf("q%d-%d>%s(%s,mq%d)", q_start, q_end, ref, strand, mapq)[seq_len(min(.N, 5))],
                                                              collapse = " | ")), by = ID], by = "ID", all.x = TRUE)
}
if (!is.null(opt$labels)) {
  lab <- fread(opt$labels, na.strings = c("", "NA"))
  if ("SVID" %in% names(lab)) setnames(lab, "SVID", "ID")
  lab <- lab[, intersect(c("ID", "len", "structure", "vcf", "svan_type", "REGION", "AF"), names(lab)), with = FALSE]
  x <- merge(lab, x, by = "ID", all.y = TRUE)
}
setorder(x, -cov_gain)
fwrite(x, paste0(opt$out_prefix, ".reference_coverage.tsv"), sep = "\t")
fwrite(new, paste0(opt$out_prefix, ".new_segments.tsv"), sep = "\t")
cat(sprintf("insertions: %d; mean coverage %s %.3f, %s %.3f; gain >= 0.2: %d; %s coverage >= 0.8: %d (vs %d on %s)\n",
            nrow(x), opt$name_a, mean(x[[ca]]), opt$name_b, mean(x[[cb]]), x[cov_gain >= 0.2, .N], opt$name_b, sum(x[[cb]] >= 0.8),
            sum(x[[ca]] >= 0.8), opt$name_a))
