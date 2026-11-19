#!/usr/bin/env Rscript
# Count reference sources per insertion from a minimap2 PAF of the inserted sequences
# (AlignRemapQueries.wdl output). Per insertion:
#   n_aln          all alignments (primary + secondary)
#   n_segments     primary alignments (tp:A:P): distinct pieces of the insertion and their best placement
#   n_loci         distinct reference loci over all alignments (same contig, merged within --merge-bp)
#   n_seg_loci     distinct loci over primary segments only
#   seg_cov        fraction of the insertion covered by the union of primary segments
#   best_cov       largest single-alignment coverage (what truvari remap's perc reports)
#   n_seg_local    primary segments on the insertion's contig within max(len, --local-bp) of its position
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(
  usage = "%prog --paf <alignments.paf.gz> --merged <merged_labels.tsv.gz> --out-prefix <prefix>",
  description = "Count reference sources per insertion from a PAF of inserted sequences.",
  option_list = list(
    make_option("--paf", help = "PAF from AlignRemapQueries.wdl (query name = variant ID)"),
    make_option("--merged", help = "Merged table from merge_svan_remap_mei.R (ID, POS from ID, len, class, svan_type)"),
    make_option("--merge-bp", dest = "merge_bp", type = "integer", default = 10000,
                help = "Alignments on the same contig within this distance are one locus [default %default]"),
    make_option("--local-bp", dest = "local_bp", type = "integer", default = 1000,
                help = "A segment is local if within max(insertion length, this) of the insertion [default %default]"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("paf", "merged", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")

# Keep the 12 standard PAF columns plus the tp tag (P = primary, S = secondary), found by name
awk_cmd <- paste("gzip -dc", shQuote(opt$paf),
                 "| awk -F'\\t' -v OFS='\\t' '{tp=\"NA\"; for (i=13;i<=NF;i++) if ($i ~ /^tp:A:/) tp=substr($i,6); print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,tp}'")
paf <- fread(cmd = awk_cmd, header = FALSE, sep = "\t",
             col.names = c("ID", "qlen", "qs", "qe", "strand", "tname", "tlen", "ts", "te",
                           "nmatch", "alen", "mapq", "tp"))
paf <- paf[tname != "*"]
paf[, ident := round(nmatch / alen, 4)]

m <- fread(opt$merged, na.strings = c("", "NA", "."))[, .(ID, len, class, svan_type, size_bin, REGION)]
m[, ins_chr := sub("-.*", "", ID)]
m[, ins_pos := as.numeric(sub("^[^-]+-([0-9]+)-.*", "\\1", ID))]
paf <- merge(paf, m[, .(ID, ins_chr, ins_pos, len)], by = "ID")
paf[, local := tname == ins_chr & pmin(abs(ts - ins_pos), abs(te - ins_pos)) <= pmax(len, opt$local_bp)]

count_loci <- function(chr, s, e, gap) {
  o <- order(chr, s); chr <- chr[o]; s <- s[o]; e <- e[o]
  n <- 0; cur_chr <- NA; cur_end <- -Inf
  for (i in seq_along(s)) {
    if (is.na(cur_chr) || chr[i] != cur_chr || s[i] > cur_end + gap) { n <- n + 1; cur_chr <- chr[i]; cur_end <- e[i] }
    else cur_end <- max(cur_end, e[i])
  }
  n
}
union_len <- function(s, e) {
  o <- order(s); s <- s[o]; e <- e[o]; tot <- 0; cs <- s[1]; ce <- e[1]
  for (i in seq_along(s)[-1]) { if (s[i] > ce) { tot <- tot + ce - cs; cs <- s[i]; ce <- e[i] } else ce <- max(ce, e[i]) }
  tot + ce - cs
}

per <- paf[, .(
  n_aln = .N,
  n_segments = sum(tp == "P"),
  n_loci = count_loci(tname, ts, te, opt$merge_bp),
  n_seg_loci = count_loci(tname[tp == "P"], ts[tp == "P"], te[tp == "P"], opt$merge_bp),
  seg_cov = round(union_len(qs[tp == "P"], qe[tp == "P"]) / qlen[1], 3),
  best_cov = round(max((qe - qs) / qlen), 3),
  n_seg_local = sum(tp == "P" & local),
  seg_contigs = paste(unique(tname[tp == "P"]), collapse = ","),
  segments = paste(sprintf("%s:%d-%d(%s):q%d-%d,id%.3f", tname[tp == "P"], ts[tp == "P"] + 1, te[tp == "P"], strand[tp == "P"],
                           qs[tp == "P"] + 1, qe[tp == "P"], ident[tp == "P"])[order(qs[tp == "P"])], collapse = "|")
), by = ID]
# Reference reuse among primary segments: summed reference span / union of reference span.
# ~1 = pieces of one copy (split by indels or divergence); >= 1.5 = the same reference stretch appears
# more than once in the insertion (multi-copy tandem gain); ins_copies ~ insertion length / reference union.
reuse <- paf[tp == "P", {
  o <- order(ts); s <- ts[o]; e <- te[o]
  u <- union_len(s, e)
  .(ref_reuse = round(sum(te - ts) / u, 2), ref_union = u, ins_copies = round(qlen[1] / u, 2))
}, by = ID]
per <- merge(per, reuse, by = "ID", all.x = TRUE)
per <- merge(m, per, by = "ID", all.y = TRUE)
per[, structure := fcase(
  n_seg_loci > 1, "segments from >1 locus",
  n_segments == 1, "single segment (< cov threshold)",
  ref_reuse >= 1.5, "one locus, multi-copy (ref reused >= 1.5x)",
  ref_reuse >= 1.1, "one locus, partly reused (1.1-1.5x)",
  default = "one locus, one copy split into pieces")]
cat(sprintf("insertions with alignments: %d\n", nrow(per)))
per[, seg_pattern := fcase(n_seg_local == n_segments, "all segments local",
                           n_seg_local == 0, "all segments distant",
                           default = "local + distant")]
fwrite(per, paste0(opt$out_prefix, ".sources_per_insertion.tsv.gz"), sep = "\t")

bin <- function(x) cut(x, c(0, 1, 2, 3, 5, 10, Inf), labels = c("1", "2", "3", "4-5", "6-10", ">10"))
out <- function(x, name) { fwrite(x, sprintf("%s.%s.tsv", opt$out_prefix, name), sep = "\t"); cat("\n==", name, "\n"); print(x) }
out(per[, .N, by = .(n_segments = bin(n_segments))][order(n_segments)], "n_segments_dist")
out(per[, .N, by = .(n_seg_loci = bin(n_seg_loci))][order(n_seg_loci)], "n_seg_loci_dist")
out(per[, .N, by = .(n_loci_all_alignments = bin(n_loci))][order(n_loci_all_alignments)], "n_loci_dist")
out(per[, .N, by = seg_pattern][order(-N)], "segment_locality")
out(dcast(per[, .N, by = .(structure, seg_pattern)], structure ~ seg_pattern, value.var = "N", fill = 0), "structure_by_locality")
out(per[structure %like% "multi-copy", .(n = .N, med_ref_reuse = median(ref_reuse), med_ins_copies = median(ins_copies),
                                          q25_ins_copies = quantile(ins_copies, .25), q75_ins_copies = quantile(ins_copies, .75))], "multicopy_summary")
out(per[, .(n = .N, med_segments = median(n_segments), med_seg_loci = median(n_seg_loci), med_loci_all = median(n_loci),
            med_seg_cov = median(seg_cov), med_best_cov = median(best_cov)), by = size_bin][order(size_bin)], "by_size")
out(per[, .(n = .N, med_segments = median(n_segments), med_seg_cov = median(seg_cov)), by = svan_type][order(-n)], "by_svan_type")
