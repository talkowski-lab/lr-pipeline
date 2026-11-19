#!/usr/bin/env bash
# Re-run the insertion classification chain (after RepeatMasker and the remap step-3 replay) for one contig:
# STR/VNTR stage -> MEI stage -> local sources -> repeat classes -> tandem / interspersed duplications ->
# mCNV + immune loci -> CNV calls -> SV-by-size figure.
# Usage: run_classification_chain.sh -d <study_dir> -c <contig prefix, e.g. chr22> [-r <ref.fa>] [-s <sd_blocks.bed>] [-n]
#   -d  study directory holding scripts/, input/, output/, figures/
#   -c  file prefix used in output/ and input/ (e.g. chr22)
#   -r  reference FASTA with .fai for the contig (default <study_dir>/input/hg38/<prefix>.fa)
#   -s  named SD blocks BED for "segmental duplication paralog" (default input/LCR22_blocks.hg38.bed)
#   -n  no PDF figures (tables only; for Linux / pipeline runs without the Quartz device)
# Expects in output/: <c>.insertion_annotations.tsv.gz, <c>.step3_hits(.annotated).tsv.gz, <c>.step3_decisions.tsv.gz,
# <c>.hit_repeats.tsv.gz, <c>.inserted_sequences.fa, <c>.insertion_flanks.tsv, repeatmasker_<c>/<c>.inserted_sequences.fa.out;
# in input/: <c>.<c>.remap_annotations.tsv, REPEATS_hg38.bed.gz, the TRExplorer catalog, LCR22_blocks.hg38.bed,
# known_cnv_loci.hg38.bed, immune_loci.hg38.bed.
set -euo pipefail

usage() { sed -n '2,15p' "$0"; exit 1; }
dir=""; c=""; ref=""; sdb=""; noplot=0
while getopts "d:c:r:s:nh" o; do
  case "$o" in d) dir=$OPTARG ;; c) c=$OPTARG ;; r) ref=$OPTARG ;; s) sdb=$OPTARG ;; n) noplot=1 ;; *) usage ;; esac
done
[[ -z "$dir" || -z "$c" ]] && usage
cd "$dir"
ref=${ref:-input/hg38/$c.fa}
sdb=${sdb:-input/LCR22_blocks.hg38.bed}
O=output; I=input
tr_cat=$I/variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz

python3 scripts/classify_tr_insertions.py --annotations $O/$c.insertion_annotations.tsv.gz --remap $I/$c.$c.remap_annotations.tsv \
  --tr-catalog $tr_cat --ref-fa "$ref" --hit-repeats $O/$c.hit_repeats.tsv.gz --out $O/$c.tr_stage.tsv
python3 scripts/annotate_mei_structure.py --rm-out $O/repeatmasker_$c/$c.inserted_sequences.fa.out --sequences $O/$c.inserted_sequences.fa \
  --flanks $O/$c.insertion_flanks.tsv --annotations $O/$c.insertion_annotations.tsv.gz --tr-stage $O/$c.tr_stage.tsv \
  --repeats $I/REPEATS_hg38.bed.gz --out $O/$c.mei_stage.tsv

awk -F'\t' 'NR>1 && $6 ~ /without insertion hallmarks/ {print $1}' $O/$c.mei_stage.tsv > $I/$c.te_no_hallmark.ids.txt
python3 scripts/locate_local_source.py --ids $I/$c.te_no_hallmark.ids.txt --sequences $O/$c.inserted_sequences.fa --flanks $O/$c.insertion_flanks.tsv \
  --ref-fa "$ref" --annotations $O/$c.insertion_annotations.tsv.gz --out $O/$c.te_no_hallmark.local_source.tsv
python3 scripts/align_to_local_origin.py --local-source $O/$c.te_no_hallmark.local_source.tsv --ids $I/$c.te_no_hallmark.ids.txt \
  --remap $I/$c.$c.remap_annotations.tsv --ref-fa "$ref" --out $O/$c.te_no_hallmark.local_aln.tsv > /dev/null
awk -F'\t' 'NR>1 && $5 == "tandem-repeat sequence, uncatalogued" {print $1}' $O/$c.tr_stage.tsv > $I/$c.uncatalogued_tr.ids.txt
python3 scripts/locate_local_source.py --ids $I/$c.uncatalogued_tr.ids.txt --sequences $O/$c.inserted_sequences.fa --flanks $O/$c.insertion_flanks.tsv \
  --ref-fa "$ref" --annotations $O/$c.insertion_annotations.tsv.gz --out $O/$c.uncatalogued_tr.local_source.tsv

python3 scripts/build_repeat_classes.py --mei-stage $O/$c.mei_stage.tsv --tr-stage $O/$c.tr_stage.tsv \
  --rm-out $O/repeatmasker_$c/$c.inserted_sequences.fa.out --sd-blocks "$sdb" \
  --te-local $O/$c.te_no_hallmark.local_source.tsv --te-local-aln $O/$c.te_no_hallmark.local_aln.tsv --repeats $I/REPEATS_hg38.bed.gz \
  --tr-novel-local $O/$c.uncatalogued_tr.local_source.tsv --out $O/$c.repeat_classes.tsv
cut -f3 $O/$c.repeat_classes.tsv | tail -n +2 | sort | uniq -c | sort -nr > $O/$c.repeat_class_counts.txt

Rscript scripts/call_tandem_dups.R --classes $O/$c.repeat_classes.tsv --hits $O/$c.step3_hits.tsv.gz \
  --decisions $O/$c.step3_decisions.tsv.gz --out-prefix $O/$c.nonrepeat > /dev/null
pdf_opt=(--out-pdf figures/$c.interspersed_dup_coverage.pdf); [[ $noplot == 1 ]] && pdf_opt=()
Rscript scripts/plot_interspersed_coverage.R --tandem $O/$c.nonrepeat.tandem_dups.tsv --hits $O/$c.step3_hits.tsv.gz \
  --title "$c: " ${pdf_opt[@]+"${pdf_opt[@]}"} --out-tsv $O/$c.interspersed_dup_coverage.tsv > /dev/null 2>&1
Rscript scripts/combine_dup_stage.R --tandem $O/$c.nonrepeat.tandem_dups.tsv --interspersed $O/$c.interspersed_dup_coverage.tsv \
  --out-prefix $O/$c.nonrepeat > /dev/null

awk -F'\t' 'NR>1 && $6 == "next stage" {print $1}' $O/$c.nonrepeat.dup_stage.tsv > $I/$c.next_stage.ids.txt
Rscript scripts/classify_partial_structure.R --annotations $O/$c.insertion_annotations.tsv.gz --hits $O/$c.step3_hits.annotated.tsv.gz \
  --ids $I/$c.next_stage.ids.txt --out-prefix $O/$c.next_stage > /dev/null
Rscript scripts/multicopy_gain_table.R --structure $O/$c.next_stage.partial_structure.tsv --hits $O/$c.step3_hits.annotated.tsv.gz \
  --out $O/$c.next_stage.multicopy_gain.tsv > /dev/null
Rscript scripts/annotate_immune_loci.R --ids $I/$c.next_stage.ids.txt --hits $O/$c.step3_hits.tsv.gz --loci $I/immune_loci.hg38.bed \
  --out $O/$c.next_stage.immune.tsv > /dev/null
Rscript scripts/combine_next_stage.R --ids $I/$c.next_stage.ids.txt --structure $O/$c.next_stage.partial_structure.tsv \
  --mcnv $O/$c.next_stage.multicopy_gain.tsv --immune $O/$c.next_stage.immune.tsv --out-prefix $O/$c.next_stage > /dev/null

awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)k[$i]=i} NR>1 && ($k["next_call"] ~ /^mCNV/ || $k["next_call"] ~ /sequence from IGL\)$/) {print $1}' \
  $O/$c.next_stage.next_stage_calls.tsv > $I/$c.cnv_candidates.ids.txt
Rscript scripts/annotate_cnv.R --ids $I/$c.cnv_candidates.ids.txt --hits $O/$c.step3_hits.tsv.gz --loci $I/known_cnv_loci.hg38.bed \
  --out $O/$c.cnv_calls.tsv > /dev/null
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)k[$i]=i} NR>1 && $k["next_call"] == "unresolved (next stage)" {print $1}' \
  $O/$c.next_stage.next_stage_calls.tsv > $I/$c.unresolved.ids.txt
python3 scripts/locate_local_source.py --ids $I/$c.unresolved.ids.txt --sequences $O/$c.inserted_sequences.fa --flanks $O/$c.insertion_flanks.tsv \
  --ref-fa "$ref" --annotations $O/$c.insertion_annotations.tsv.gz --out $O/$c.unresolved.local_source.tsv

pdf_opt=(--out-pdf figures/$c.sv_by_size.pdf); [[ $noplot == 1 ]] && pdf_opt=()
Rscript scripts/plot_sv_by_size.R --annotations $O/$c.insertion_annotations.tsv.gz --classes $O/$c.repeat_classes.tsv --title "$c: " \
  ${pdf_opt[@]+"${pdf_opt[@]}"} --out-tsv $O/$c.sv_by_size.tsv > /dev/null 2>&1
echo "done: $c"
