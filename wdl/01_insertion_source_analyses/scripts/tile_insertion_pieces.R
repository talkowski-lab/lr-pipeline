#!/usr/bin/env Rscript
# Find insertions whose different pieces align to different places: tile each insertion with
# non-overlapping alignments and count how many separate reference loci the tiles come from.
#
# Per insertion:
#   pieces    = hits with identity >= --min-ident and query span >= max(--min-piece-bp, --min-piece-frac * len)
#   tiling    = greedy by alignment score: a piece is added unless its query overlap with any piece already
#               chosen exceeds --max-overlap x the shorter of the two
#   loci      = chosen pieces merged on the reference when on the same contig and within --merge-bp
#               (so the two halves of a rotated tandem dup, which meet at the site, are one locus)
#   multi_piece = >= 2 loci, tiles cover >= --min-cov of the insertion, and no single piece covers
#                 >= --max-single of it (one alignment does not explain the insertion)
#   unique_pieces = multi_piece and every chosen piece has mapq >= --min-mapq
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(
  usage = "%prog --hits <step3_hits.annotated.tsv.gz> --out-prefix <prefix>",
  description = "Tile insertions with non-overlapping alignments and flag those built from pieces of several loci.",
  option_list = list(
    make_option("--hits", help = "Annotated hit table from annotate_step3_hits.R"),
    make_option("--min-ident", dest = "min_ident", type = "double", default = 0.9, help = "Minimum piece identity [default %default]"),
    make_option("--min-piece-bp", dest = "min_bp", type = "integer", default = 50, help = "Minimum piece query span, bp [default %default]"),
    make_option("--min-piece-frac", dest = "min_frac", type = "double", default = 0.05, help = "Minimum piece query span, fraction of insertion [default %default]"),
    make_option("--max-overlap", dest = "max_ov", type = "double", default = 0.2, help = "Max query overlap between two pieces, fraction of the shorter [default %default]"),
    make_option("--merge-bp", dest = "merge_bp", type = "integer", default = 1000, help = "Pieces on the same contig within this distance are one locus [default %default]"),
    make_option("--min-cov", dest = "min_cov", type = "double", default = 0.8, help = "Minimum insertion fraction covered by the tiles [default %default]"),
    make_option("--max-single", dest = "max_single", type = "double", default = 0.8, help = "multi_piece requires the largest piece to cover less than this [default %default]"),
    make_option("--min-mapq", dest = "min_mapq", type = "integer", default = 20, help = "mapq for unique_pieces [default %default]"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("hits", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")

h <- fread(opt$hits, na.strings = c("", "NA"))
h[, qspan := q_end - q_start + 1]
h[, `:=`(rchrom = sub(":.*", "", coords), rs = as.numeric(sub("^[^:]+:([0-9]+)-.*", "\\1", coords)), re = as.numeric(sub(".*-", "", coords)))]
p <- h[identity >= opt$min_ident & qspan >= pmax(opt$min_bp, opt$min_frac * qlen)][order(ID, -score, -qspan)]

tile <- function(qs, qe) {
  keep <- logical(length(qs))
  for (i in seq_along(qs)) {
    ok <- TRUE
    for (j in which(keep)) {
      ov <- min(qe[i], qe[j]) - max(qs[i], qs[j]) + 1
      if (ov > opt$max_ov * min(qe[i] - qs[i] + 1, qe[j] - qs[j] + 1)) { ok <- FALSE; break }
    }
    keep[i] <- ok
  }
  keep
}
p[, chosen := tile(q_start, q_end), by = ID]
t <- p[chosen == TRUE][order(ID, q_start)]

locus_id <- function(chr, s, e, gap) {
  o <- order(chr, s); id <- integer(length(s)); n <- 0L; cc <- NA; ce <- -Inf
  for (k in o) {
    if (is.na(cc) || chr[k] != cc || s[k] > ce + gap) { n <- n + 1L; cc <- chr[k]; ce <- e[k] } else ce <- max(ce, e[k])
    id[k] <- n
  }
  id
}
t[, locus := locus_id(rchrom, rs, re, opt$merge_bp), by = ID]
t[, piece := sprintf("q%d-%d>%s:%d-%d(%s,%.1f%%id,mq%d,L%d,%s)", q_start, q_end, rchrom, rs, re, strand, 100 * identity, mapq, locus, hit_label)]

union_len <- function(s, e) {
  o <- order(s); s <- s[o]; e <- e[o]; tot <- 0; cs <- s[1]; ce <- e[1]
  for (i in seq_along(s)[-1]) { if (s[i] > ce + 1) { tot <- tot + ce - cs + 1; cs <- s[i]; ce <- e[i] } else ce <- max(ce, e[i]) }
  tot + ce - cs + 1
}
s <- t[, .(
  n_pieces = .N, n_loci = uniqueN(locus), n_contigs = uniqueN(rchrom),
  tile_cov = round(union_len(q_start, q_end) / qlen[1], 3),
  max_piece_frac = round(max(qspan) / qlen[1], 3),
  min_mapq = min(mapq), min_identity = min(identity),
  piece_labels = paste(sort(unique(hit_label)), collapse = ","),
  pieces = paste(piece, collapse = " | ")
), by = .(ID, len = qlen, remap_class, ins_annotation, allele_subtype, TRID, svan_type, rotated_tandem_dup, REGION, AF, FILTER)]
s[, multi_piece := n_loci >= 2 & tile_cov >= opt$min_cov & max_piece_frac < opt$max_single]
s[, unique_pieces := multi_piece & min_mapq >= opt$min_mapq]
setorder(s, -multi_piece, -unique_pieces, -len)

fwrite(s, paste0(opt$out_prefix, ".tiled_pieces.tsv"), sep = "\t")
fwrite(t[ID %in% s[multi_piece == TRUE, ID], .(ID, piece_order = seq_len(.N), q_start, q_end, qspan, strand, coords, identity, mapq, score, locus, hit_label), by = ID][, ID := NULL],
       paste0(opt$out_prefix, ".multi_piece.pieces.tsv"), sep = "\t")

o <- function(x, name) { fwrite(x, sprintf("%s.%s.tsv", opt$out_prefix, name), sep = "\t"); cat("\n==", name, "\n"); print(x) }
cat(sprintf("insertions with pieces: %d; multi_piece: %d; with all pieces mapq >= %d: %d\n",
            nrow(s), s[multi_piece == TRUE, .N], opt$min_mapq, s[unique_pieces == TRUE, .N]))
o(dcast(s[multi_piece == TRUE, .N, by = .(ins_annotation, remap_class)], ins_annotation ~ remap_class, value.var = "N", fill = 0), "multi_piece_by_annotation")
o(s[multi_piece == TRUE, .N, by = .(n_loci = pmin(n_loci, 5))][order(n_loci)], "multi_piece_n_loci")
