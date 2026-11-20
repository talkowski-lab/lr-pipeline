#!/usr/bin/env Rscript
# Replay the classification step of truvari anno remap (Yulia Mostovoy's fork, remap.py remap_entry,
# local rescue off) from a PAF of the full ALT sequences aligned with the same minimap2 command
# (AlignRemapQueries.wdl with keep_anchor=true, classes=all), and export the evidence behind each call.
#
# Fork logic reproduced:
#   every alignment record is a hit; perc = (query_end - query_start + 1) / len(ALT)
#   chrom_dist = min(|hit start - POS|, |hit end - POS|) on the same contig, Inf otherwise
#   full hit: perc >= cov_threshold; partial hit otherwise
#   close_dist / inv_close_dist = smallest chrom_dist among full hits on + / - strand
#   class: novel if no hits; tandem if close_dist <= len(ALT); tandem_inverted if inv_close_dist <= len(ALT);
#          partial if no full hits; interspersed otherwise
#   reported hit: tandem* sorted by (chrom_dist, -perc, -score); others by (-score, -perc, chrom_dist);
#                 with no full hits, the best partial hit by (score, perc, -chrom_dist)
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(
  usage = "%prog --paf <full_alt.paf.gz> --remap <remap.tsv> --out-prefix <prefix> [--cov-threshold 0.8]",
  description = "Replay truvari remap step 3 (classification) and export per-insertion hits and decisions.",
  option_list = list(
    make_option("--paf", help = "PAF of the full ALT sequences (query name = variant ID)"),
    make_option("--remap", help = "Truvari remap TSV (no header: CHROM POS REF ALT ID class coords ori perc)"),
    make_option("--cov-threshold", dest = "cov", type = "double", default = 0.8, help = "Remap --cov-threshold [default %default]"),
    make_option("--out-prefix", dest = "out_prefix", help = "Output prefix")
  )
))
for (a in c("paf", "remap", "out_prefix")) if (is.null(opt[[a]])) stop("--", a, " is required")

remap <- fread(opt$remap, header = FALSE, sep = "\t", quote = "",
               col.names = c("CHROM", "POS", "REF", "ALT", "ID", "remap_class", "remap_coords", "remap_ori", "remap_perc"))
remap[, alt_len := nchar(ALT)]
remap[, c("REF", "ALT") := NULL]

awk_cmd <- paste("gzip -dc", shQuote(opt$paf),
                 "| awk -F'\\t' -v OFS='\\t' '{tp=\"NA\"; as=\"NA\"; for (i=13;i<=NF;i++) { if ($i ~ /^tp:A:/) tp=substr($i,6); if ($i ~ /^AS:i:/) as=substr($i,6) } print $1,$2,$3,$4,$5,$6,$8,$9,$10,$11,$12,tp,as}'")
h <- fread(cmd = awk_cmd, header = FALSE, sep = "\t",
           col.names = c("ID", "qlen", "qs", "qe", "strand", "tname", "ts", "te", "nmatch", "alen", "mapq", "tp", "score"))
h <- h[tname != "*"]
h <- merge(h, remap[, .(ID, CHROM, POS, alt_len)], by = "ID")
h[, `:=`(hit_start = ts + 1, hit_end = te)]
h[, perc := 100 * (qe - qs) / alt_len]
h[, chrom_dist := fifelse(tname == CHROM, pmin(abs(hit_start - POS), abs(hit_end - POS)), Inf)]
h[, full := perc / 100 >= opt$cov]
h[, identity := round(nmatch / alen, 4)]
h[, coords := sprintf("%s:%d-%d", tname, hit_start, hit_end)]

dec <- h[, .(
  n_hits = .N,
  n_full = sum(full),
  n_partial = sum(!full),
  n_full_same_chrom = sum(full & tname == CHROM),
  close_dist = if (any(full & tname == CHROM & strand == "+")) min(chrom_dist[full & tname == CHROM & strand == "+"]) else NA_real_,
  inv_close_dist = if (any(full & tname == CHROM & strand == "-")) min(chrom_dist[full & tname == CHROM & strand == "-"]) else NA_real_,
  max_perc = round(max(perc), 2),
  best_score_hit = coords[order(-score, -perc, chrom_dist)][1],
  best_score_perc = round(perc[order(-score, -perc, chrom_dist)][1], 2),
  best_score_identity = identity[order(-score, -perc, chrom_dist)][1],
  best_score_mapq = mapq[order(-score, -perc, chrom_dist)][1]
), by = ID]

d <- merge(remap, dec, by = "ID", all.x = TRUE)
d[is.na(n_hits), `:=`(n_hits = 0L, n_full = 0L, n_partial = 0L, n_full_same_chrom = 0L)]
d[, replay_class := fcase(
  n_hits == 0, "novel",
  !is.na(close_dist) & close_dist <= alt_len, "tandem",
  !is.na(inv_close_dist) & inv_close_dist <= alt_len, "tandem_inverted",
  n_full == 0, "partial",
  default = "interspersed")]

# reported hit, as the fork orders hit_records
rep <- merge(h, d[, .(ID, replay_class)], by = "ID")
rep <- rep[(replay_class %in% c("tandem", "tandem_inverted") & full) | (!replay_class %in% c("tandem", "tandem_inverted") & full) |
             (replay_class == "partial")]
rep[, ord1 := fifelse(replay_class %in% c("tandem", "tandem_inverted"), chrom_dist, -score)]
rep[, ord2 := fifelse(replay_class %in% c("tandem", "tandem_inverted"), -perc, -perc)]
rep[, ord3 := fifelse(replay_class %in% c("tandem", "tandem_inverted"), -score, chrom_dist)]
rep <- rep[order(ID, ord1, ord2, ord3)][, .SD[1], by = ID][, .(ID, replay_coords = coords, replay_ori = strand, replay_perc = round(perc, 2),
                                                             replay_hit_identity = identity, replay_hit_mapq = mapq, replay_hit_dist = chrom_dist)]
d <- merge(d, rep, by = "ID", all.x = TRUE)
d[, class_match := replay_class == remap_class]
d[, coords_match := fifelse(is.na(remap_coords) | remap_coords == ".", is.na(replay_coords), remap_coords == replay_coords)]

# why each call came out the way it did
d[, decision := fcase(
  replay_class == "novel", "no alignment",
  replay_class == "tandem", sprintf("full + hit %s bp from site <= len %d", close_dist, alt_len),
  replay_class == "tandem_inverted", sprintf("full - hit %s bp from site <= len %d", inv_close_dist, alt_len),
  replay_class == "partial", sprintf("%d hits, none >= %d%% (max %.1f%%)", n_hits, as.integer(100 * opt$cov), max_perc),
  default = sprintf("%d full hits, nearest same-contig + %s / - %s bp, len %d", n_full,
                    fifelse(is.na(close_dist), "none", as.character(close_dist)),
                    fifelse(is.na(inv_close_dist), "none", as.character(inv_close_dist)), alt_len))]
setcolorder(d, c("ID", "CHROM", "POS", "alt_len", "remap_class", "replay_class", "class_match", "decision"))
fwrite(d, paste0(opt$out_prefix, ".step3_decisions.tsv.gz"), sep = "\t")

hits_out <- h[, .(ID, qlen, q_start = qs + 1, q_end = qe, strand, coords, mapq, score, identity, perc = round(perc, 2),
                  chrom_dist, full, type = tp)][order(ID, -score)]
fwrite(hits_out, paste0(opt$out_prefix, ".step3_hits.tsv.gz"), sep = "\t")

cat(sprintf("insertions: %d; hits: %d\n", nrow(d), nrow(h)))
cat("\n== replayed vs remap class\n"); print(dcast(d[, .N, by = .(remap_class, replay_class)], remap_class ~ replay_class, value.var = "N", fill = 0))
cat(sprintf("\nclass agreement: %d / %d; reported hit agreement: %d / %d\n", sum(d$class_match), nrow(d), sum(d$coords_match, na.rm = TRUE), nrow(d)))
s <- d[, .(n = .N, med_hits = median(n_hits), med_full = median(n_full), med_max_perc = median(max_perc, na.rm = TRUE)), by = replay_class]
fwrite(s, paste0(opt$out_prefix, ".step3_summary.tsv"), sep = "\t")
cat("\n== per class\n"); print(s)
