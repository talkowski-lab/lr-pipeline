#!/usr/bin/env Rscript
# Immune-locus annotation of insertions (fetch_immune_loci.py BED: IGH, IGK, IGL, TRA/TRD, TRB, TRG, MHC, LRC):
#   site_locus     immune locus containing the insertion site
#   source_locus   immune locus of the best-scoring alignment (remap step-3 hits), and immune_aln_frac = fraction of the
#                  inserted sequence covered by alignments that land in any immune locus
#   immune_call    "in <locus> (sequence from <locus>)", "in <locus>", "sequence from <locus>", or empty
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})
opt <- parse_args(OptionParser(
  usage = "%prog --ids <ids.txt> --hits <step3_hits.tsv.gz> --loci <immune_loci.bed> --out <tsv>",
  description = "Annotate insertion sites and alignment sources in immune loci.",
  option_list = list(
    make_option("--ids", help = "File of insertion IDs"),
    make_option("--hits", help = "Per-hit table from replay_remap_step3.R"),
    make_option("--loci", help = "BED of immune loci (chrom start end name ...)"),
    make_option("--min-frac", dest = "min_frac", type = "double", default = 0.5,
                help = "Min fraction of the insertion aligned into immune loci for 'sequence from' [default %default]"),
    make_option("--out", help = "Output TSV")
  )
))
for (a in c("ids", "hits", "loci", "out")) if (is.null(opt[[a]])) stop("--", a, " is required")
ids <- readLines(opt$ids)
L <- fread(opt$loci, header = FALSE, select = 1:4, col.names = c("chrom", "s", "e", "locus"))
locus_of <- function(ch, p) vapply(seq_along(ch), function(i) { k <- L[chrom == ch[i] & s < p[i] & e >= p[i], locus]; if (length(k)) k[1] else NA_character_ }, "")
x <- data.table(ID = ids, CHROM = sub("-.*", "", ids), POS = as.numeric(sub("^[^-]+-([0-9]+)-.*", "\\1", ids)))
x[, site_locus := locus_of(CHROM, POS)]
h <- fread(opt$hits)[ID %in% ids]
h[, `:=`(rchrom = sub(":.*", "", coords), mid = (as.numeric(sub("^[^:]+:([0-9]+)-.*", "\\1", coords)) + as.numeric(sub(".*-", "", coords))) / 2)]
h[, hit_locus := locus_of(rchrom, mid)]
union_len <- function(s, e) {
  if (!length(s)) return(0)
  o <- order(s); s <- s[o]; e <- e[o]; tot <- 0; cs <- s[1]; ce <- e[1]
  for (i in seq_along(s)[-1]) { if (s[i] > ce + 1) { tot <- tot + ce - cs + 1; cs <- s[i]; ce <- e[i] } else ce <- max(ce, e[i]) }
  tot + ce - cs + 1
}
hs <- h[, .(source_locus = hit_locus[which.max(score)], best_hit = coords[which.max(score)],
            immune_aln_frac = round(union_len(q_start[!is.na(hit_locus)], q_end[!is.na(hit_locus)]) / qlen[1], 3),
            immune_loci_hit = paste(sort(unique(na.omit(hit_locus))), collapse = ",")), by = ID]
x <- merge(x, hs, by = "ID", all.x = TRUE)
x[is.na(immune_aln_frac), immune_aln_frac := 0]
# "sequence from" requires the best-scoring alignment itself to land in the locus (secondary hits alone do not count)
x[, from_locus := fifelse(immune_aln_frac >= opt$min_frac & !is.na(source_locus), source_locus, NA_character_)]
x[, immune_call := fcase(!is.na(site_locus) & !is.na(from_locus), sprintf("in %s (sequence from %s)", site_locus, from_locus),
                         !is.na(site_locus), sprintf("in %s", site_locus),
                         !is.na(from_locus), sprintf("sequence from %s", from_locus),
                         default = NA_character_)]
fwrite(x[, .(ID, CHROM, POS, site_locus, source_locus, best_hit, immune_aln_frac, immune_loci_hit, immune_call)], opt$out, sep = "\t")
print(x[, .N, by = immune_call][order(-N)])
