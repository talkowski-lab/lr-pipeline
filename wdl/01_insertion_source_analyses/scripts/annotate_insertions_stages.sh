#!/usr/bin/env bash
# Insertion annotation pipeline, one contig, in stages (the alignment steps between stages run elsewhere, e.g. as WDL
# tasks). Every stage reads and writes the working directory <W> (input/, output/; scripts/ is linked to this script's
# folder), so the stages can be chained by passing <W> as a tarball.
#
# Usage: annotate_insertions_stages.sh <stage> -w <workdir> -c <contig> [stage options]
#   prepare  -r remap.tsv -v vcf(.gz, indexed) -f ref.fa -s svan.tsv -H svan_header.txt [-m min_len]
#            [-I info.tsv.gz: an existing extract_ins_info.sh table instead of reading -v]
#            -> VCF INFO table, SVAN / remap merge, inserted sequences + flanks, output/<c>.queries.fa (full ALT, for the
#               genome-wide minimap2 alignment)
#   stage1   -p hg38_all.paf.gz -o repeatmasker.out -f ref.fa -R resources_dir
#            -> step-3 replay, repeat annotation of hits, insertion annotations, classification chain (STR / VNTR,
#               satellites, MEI, duplications, CNV / immune), output/<c>.t2t_queries.fa (still unexplained, full ALT)
#   stage2   -t t2t.paf.gz -f ref.fa -R resources_dir
#            -> T2T labels, local-origin duplications, output/<c>.bwa_queries.fa (still unresolved, no anchor)
#   stage3   -b bwa.sam.gz -R resources_dir
#            -> BWA calls (incl. NUMT), output/<c>.iterative_queries.fa (still unresolved, no anchor)
#   stage4   -i hg38_iterative_segments.tsv -j t2t_iterative_segments.tsv -f ref.fa -R resources_dir [-t threads]
#            -> iterative-layout calls, T2T-only labels, repeat-region context, SD-family links, last step (local
#               realignment to the site's SD copy, tandem duplication / CPX_SV / CPX_INS), output/<c>.final_annotation.tsv
# resources_dir holds: REPEATS_hg38.bed.gz, variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz,
# immune_loci.hg38.bed, known_cnv_loci.hg38.bed, hg38.cytoBand.txt, hg38.sd_reference.gd_blocks.bed,
# hg38.sd_reference.loci.bed, hg38.sd_reference.copies.bed, hs1.ncbiRefSeqCurated.txt.gz, genomicSuperDups.txt.gz
# Needs on PATH: python3 (pysam), Rscript (data.table, optparse), bcftools, minimap2 (stage4).
set -euo pipefail

usage() { sed -n '2,25p' "$0"; exit 1; }
[[ $# -lt 1 ]] && usage
stage=$1; shift
W=""; c=""; remap=""; vcf=""; info=""; ref=""; svan=""; svanh=""; minlen=50; paf=""; rm=""; R=""; t2t=""; sam=""; ih=""; it=""; th=2
while getopts "w:c:r:v:I:f:s:H:m:p:o:R:t:b:i:j:h" o; do
  case "$o" in
    w) W=$OPTARG ;; c) c=$OPTARG ;; r) remap=$OPTARG ;; v) vcf=$OPTARG ;; I) info=$OPTARG ;; f) ref=$OPTARG ;; s) svan=$OPTARG ;;
    H) svanh=$OPTARG ;; m) minlen=$OPTARG ;; p) paf=$OPTARG ;; o) rm=$OPTARG ;; R) R=$OPTARG ;; b) sam=$OPTARG ;;
    i) ih=$OPTARG ;; j) it=$OPTARG ;; t) if [[ $stage == stage2 ]]; then t2t=$OPTARG; else th=$OPTARG; fi ;; *) usage ;;
  esac
done
[[ -z "$W" || -z "$c" ]] && usage
S=$(cd "$(dirname "$0")" && pwd)
abs() { [[ -z "$1" ]] && return; [[ "$1" == gs://* ]] && { echo "$1"; return; }; echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }
ref=$(abs "$ref"); remap=$(abs "$remap"); vcf=$(abs "$vcf"); info=$(abs "$info"); svan=$(abs "$svan"); svanh=$(abs "$svanh")
paf=$(abs "$paf"); rm=$(abs "$rm"); t2t=$(abs "$t2t"); sam=$(abs "$sam"); ih=$(abs "$ih"); it=$(abs "$it")
[[ -n "$R" ]] && R=$(cd "$R" && pwd)
mkdir -p "$W/input" "$W/output" "$W/figures"
cd "$W"
[[ -e scripts ]] || ln -s "$S" scripts
I=input; O=output

link_resources() {
  for f in REPEATS_hg38.bed.gz variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz immune_loci.hg38.bed \
           known_cnv_loci.hg38.bed hg38.cytoBand.txt hg38.sd_reference.gd_blocks.bed hg38.sd_reference.loci.bed \
           hg38.sd_reference.copies.bed hs1.ncbiRefSeqCurated.txt.gz genomicSuperDups.txt.gz; do
    [[ -e $I/$f ]] || ln -s "$R/$f" "$I/$f"
  done
}

# final annotation with whatever stage outputs exist so far
final_build() {
  local t2tl=$O/$c.t2t_labels.all.tsv
  if [[ ! -s $t2tl ]]; then printf 'ID\tlen\tt2t_call\tT2T_best_perc\tT2T_best_hit\tT2T_cov\thg38_best_perc\thg38_cov\tt2t_genes\tt2t_label\n' > "$t2tl"; fi
  local opt=()
  [[ -s $O/$c.nearby_dup.origin_alignment.tsv ]] && opt+=(--local-dups "$O/$c.nearby_dup.origin_alignment.tsv")
  [[ -s $O/$c.bwa_calls.tsv ]] && opt+=(--bwa "$O/$c.bwa_calls.tsv")
  [[ -s $O/$c.iterative_classes.tsv ]] && opt+=(--iterative "$O/$c.iterative_classes.tsv")
  [[ -s $O/$c.low_identity.site_context.tsv ]] && opt+=(--site-context "$O/$c.low_identity.site_context.tsv")
  [[ -s $O/$c.sd_links.tsv ]] && opt+=(--sd-links "$O/$c.sd_links.tsv")
  [[ -s $O/$c.last_step_calls.tsv ]] && opt+=(--last-step "$O/$c.last_step_calls.tsv")
  Rscript scripts/build_final_annotation.R --classes $O/$c.repeat_classes.tsv --dups $O/$c.nonrepeat.dup_stage.tsv \
    --cnv $O/$c.cnv_calls.tsv --next $O/$c.next_stage.next_stage_calls.tsv --t2t "$t2tl" ${opt[@]+"${opt[@]}"} \
    --out-prefix $O/$c > /dev/null
  awk -F'\t' 'NR>1 && $3=="unresolved"{print $1}' $O/$c.final_annotation.tsv > $I/$c.final_unresolved.ids.txt
}

# FASTA of the IDs in a list: full ALT from the remap TSV (anchor kept) or the inserted sequence (no anchor)
fasta_of() {  # ids_file out.fa full|insert
  python3 - "$1" "$2" "$3" "$I/$c.$c.remap_annotations.tsv" "$O/$c.inserted_sequences.fa" <<'EOF'
import sys
ids_f, out, mode, remap, ins = sys.argv[1:6]
ids = {x.strip() for x in open(ids_f) if x.strip()}
with open(out, "w") as o:
    if mode == "full":
        for line in open(remap):
            c = line.rstrip("\n").split("\t")
            if c[4] in ids:
                o.write(f">{c[4]}\n{c[3]}\n")
    else:
        keep = False
        for line in open(ins):
            if line.startswith(">"):
                keep = line[1:].split()[0] in ids
            if keep:
                o.write(line)
EOF
}

case $stage in
prepare)
  cp "$remap" $I/$c.$c.remap_annotations.tsv
  if [[ -n "$info" ]]; then cp "$info" $I/$c.V10.ins50.info.tsv.gz; else bash scripts/extract_ins_info.sh "$vcf" "$minlen" $I/$c.V10.ins50.info.tsv.gz; fi
  Rscript scripts/summarize_remap.R --remap $I/$c.$c.remap_annotations.tsv --info $I/$c.V10.ins50.info.tsv.gz --out-prefix $O/$c > /dev/null
  cp "$svan" $I/$c.svan_annotations.tsv; cp "$svanh" $I/$c.svan_annotations.header.txt
  Rscript scripts/compare_svan_remap.R --remap $I/$c.$c.remap_annotations.tsv --svan $I/$c.svan_annotations.tsv \
    --svan-header $I/$c.svan_annotations.header.txt --joined $O/$c.joined.tsv.gz --out-prefix $O/$c > /dev/null
  Rscript scripts/merge_svan_remap_mei.R --joined $O/$c.svan_remap.joined.tsv.gz --len-tol 0.1 --out-prefix $O/$c > /dev/null
  python3 scripts/extract_insertion_sequences.py --remap $I/$c.$c.remap_annotations.tsv --ref-fa "$ref" \
    --out-fa $O/$c.inserted_sequences.fa --out-flanks $O/$c.insertion_flanks.tsv > /dev/null
  cut -f5 $I/$c.$c.remap_annotations.tsv > $I/$c.all.ids.txt
  fasta_of $I/$c.all.ids.txt $O/$c.queries.fa full
  echo "prepare done: $(grep -c '>' $O/$c.queries.fa) insertions"
  ;;
stage1)
  link_resources
  cp "$paf" $I/$c.all.alignments.paf.gz
  mkdir -p $O/repeatmasker_$c; cp "$rm" $O/repeatmasker_$c/$c.inserted_sequences.fa.out
  P=$I/$c.all.alignments.paf.gz
  Rscript scripts/count_alignment_sources.R --paf $P --merged $O/$c.merged_labels.tsv.gz --out-prefix $O/$c.partial > /dev/null
  Rscript scripts/find_rotated_tandem_dups.R --paf $P --info $I/$c.V10.ins50.info.tsv.gz \
    --sources $O/$c.partial.sources_per_insertion.tsv.gz --tol 200 --out-prefix $O/$c.partial.tol200 > /dev/null
  Rscript scripts/replay_remap_step3.R --paf $P --remap $I/$c.$c.remap_annotations.tsv --out-prefix $O/$c > /dev/null
  Rscript scripts/annotate_hit_repeats.R --hits $O/$c.step3_hits.tsv.gz --decisions $O/$c.step3_decisions.tsv.gz \
    --repeats $I/REPEATS_hg38.bed.gz --tr-catalog $I/variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz \
    --info $I/$c.V10.ins50.info.tsv.gz --out-prefix $O/$c > /dev/null
  Rscript scripts/annotate_step3_hits.R --hits $O/$c.step3_hits.tsv.gz --hit-repeats $O/$c.hit_repeats.tsv.gz \
    --decisions $O/$c.step3_decisions.tsv.gz --ins-labels $O/$c.insertion_repeat_labels.tsv.gz \
    --merged $O/$c.merged_labels.tsv.gz --rotated $O/$c.partial.tol200.rotated_tandem_dups.tsv \
    --info $I/$c.V10.ins50.info.tsv.gz --out $O/$c.step3_hits.annotated.tsv.gz > /dev/null
  Rscript scripts/tile_insertion_pieces.R --hits $O/$c.step3_hits.annotated.tsv.gz --out-prefix $O/$c > /dev/null
  Rscript scripts/build_insertion_annotations.R --decisions $O/$c.step3_decisions.tsv.gz --hits $O/$c.step3_hits.annotated.tsv.gz \
    --tiled $O/$c.tiled_pieces.tsv --pieces $O/$c.multi_piece.pieces.tsv --info $I/$c.V10.ins50.info.tsv.gz \
    --out $O/$c.insertion_annotations.tsv.gz > /dev/null
  bash scripts/run_classification_chain.sh -d . -c $c -r "$ref" -s $I/hg38.sd_reference.gd_blocks.bed -n
  final_build
  cp $I/$c.final_unresolved.ids.txt $I/$c.t2t_round2.ids.txt
  fasta_of $I/$c.t2t_round2.ids.txt $O/$c.t2t_queries.fa full
  echo "stage1 done: $(wc -l < $I/$c.t2t_round2.ids.txt) unexplained -> T2T"
  ;;
stage2)
  link_resources
  cp "$t2t" $I/$c.t2t_round2.hs1.alignments.paf.gz
  Rscript scripts/compare_reference_coverage.R --paf-a $I/$c.all.alignments.paf.gz --paf-b $I/$c.t2t_round2.hs1.alignments.paf.gz \
    --ids $I/$c.t2t_round2.ids.txt --min-bp 30 --out-prefix $O/$c.t2t_round2.hg38_vs_T2T > /dev/null
  Rscript scripts/t2t_match_calls.R --coverage $O/$c.t2t_round2.hg38_vs_T2T.reference_coverage.tsv --out-prefix $O/$c.t2t_round2 > /dev/null
  python3 scripts/name_t2t_matches.py --calls $O/$c.t2t_round2.t2t_calls.tsv --genes-table $I/hs1.ncbiRefSeqCurated.txt.gz \
    --out $O/$c.t2t_round2.t2t_labels.tsv > /dev/null
  cp $O/$c.t2t_round2.t2t_labels.tsv $O/$c.t2t_labels.all.tsv
  final_build
  python3 scripts/locate_local_source.py --ids $I/$c.final_unresolved.ids.txt --sequences $O/$c.inserted_sequences.fa \
    --flanks $O/$c.insertion_flanks.tsv --ref-fa "$ref" --annotations $O/$c.insertion_annotations.tsv.gz \
    --out $O/$c.final_unresolved.local_source.tsv > /dev/null
  awk -F'\t' 'NR>1 && $5=="nearby duplication"{print $1}' $O/$c.final_unresolved.local_source.tsv > $I/$c.nearby_dup.ids.txt
  python3 scripts/align_to_local_origin.py --local-source $O/$c.final_unresolved.local_source.tsv --ids $I/$c.nearby_dup.ids.txt \
    --remap $I/$c.$c.remap_annotations.tsv --ref-fa "$ref" --out $O/$c.nearby_dup.origin_alignment.tsv > /dev/null
  final_build
  fasta_of $I/$c.final_unresolved.ids.txt $O/$c.bwa_queries.fa insert
  echo "stage2 done: $(wc -l < $I/$c.final_unresolved.ids.txt) unresolved -> BWA"
  ;;
stage3)
  link_resources
  cp "$sam" $I/$c.unresolved.bwa.sam.gz
  python3 scripts/summarize_bwa_alignments.py --sam $I/$c.unresolved.bwa.sam.gz --queries $O/$c.bwa_queries.fa \
    --cytoband $I/hg38.cytoBand.txt --out $O/$c.bwa_calls.tsv > /dev/null
  final_build
  fasta_of $I/$c.final_unresolved.ids.txt $O/$c.iterative_queries.fa insert
  echo "stage3 done: $(wc -l < $I/$c.final_unresolved.ids.txt) unresolved -> iterative alignment"
  ;;
stage4)
  link_resources
  cp "$ih" $I/$c.unresolved.hg38.iterative_segments.tsv; cp "$it" $I/$c.unresolved.T2T.iterative_segments.tsv
  python3 scripts/classify_iterative_layouts.py --site-segments $I/$c.unresolved.hg38.iterative_segments.tsv \
    --other-segments $I/$c.unresolved.T2T.iterative_segments.tsv --t2t-calls-out $O/$c.iterative.t2t_calls.tsv \
    --out $O/$c.iterative_classes.tsv > /dev/null
  python3 scripts/name_t2t_matches.py --calls $O/$c.iterative.t2t_calls.tsv --genes-table $I/hs1.ncbiRefSeqCurated.txt.gz \
    --out $O/$c.iterative.t2t_labels.tsv > /dev/null
  # one row per ID: round-2 rows, replaced by the iterative label where there is one
  (head -1 $O/$c.t2t_round2.t2t_labels.tsv
   awk -F'\t' 'NR==FNR{if(FNR>1 && $10!="")s[$1]=1; next} FNR>1 && !($1 in s)' $O/$c.iterative.t2t_labels.tsv $O/$c.t2t_round2.t2t_labels.tsv
   awk -F'\t' 'NR>1 && $10!=""' $O/$c.iterative.t2t_labels.tsv) > $O/$c.t2t_labels.all.tsv
  awk -F'\t' 'NR>1 && $6=="low-identity match"{print $1}' $O/$c.iterative_classes.tsv > $I/$c.low_identity.ids.txt
  python3 scripts/site_context.py --ids $I/$c.low_identity.ids.txt --cytoband $I/hg38.cytoBand.txt \
    --sds $I/genomicSuperDups.txt.gz --sd-blocks $I/hg38.sd_reference.gd_blocks.bed \
    --tr-catalog $I/variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz --repeats $I/REPEATS_hg38.bed.gz \
    --out $O/$c.low_identity.site_context.tsv > /dev/null
  final_build
  python3 scripts/link_sd_insertions.py --loci $I/hg38.sd_reference.loci.bed --copies $I/hg38.sd_reference.copies.bed \
    --origin remap_interspersed=$O/$c.interspersed_dup_coverage.tsv,ID,source \
    --origin remap_tandem=$O/$c.nonrepeat.tandem_dups.tsv,ID,unit \
    --origin tr_stage=$O/$c.tr_stage.tsv,ID,best_hit \
    --origin te_local=$O/$c.te_no_hallmark.local_aln.tsv,ID,origin_aligned \
    --origin nearby=$O/$c.nearby_dup.origin_alignment.tsv,ID,origin_aligned \
    --origin bwa=$O/$c.bwa_calls.tsv,ID,aligned_location \
    --origin iterative=$I/$c.unresolved.hg38.iterative_segments.tsv,ID,ref,status:aligned \
    --out-prefix $O/$c.sd_links > /dev/null
  final_build
  python3 scripts/prepare_site_copy_realign.py --ids $I/$c.final_unresolved.ids.txt --segments $I/$c.unresolved.hg38.iterative_segments.tsv \
    --sequences $O/$c.inserted_sequences.fa --loci $I/hg38.sd_reference.loci.bed --ref-fa "$ref" --out-dir $O/site_copy_realign > /dev/null
  (cd $O/site_copy_realign
   tail -n +2 manifest.tsv | cut -f1 | while read -r id; do
     [[ -s $id.insertion.fa ]] || continue
     minimap2 -x asm20 -d $id.asm20.mmi $id.target.fa 2>/dev/null
     minimap2 -x sr -d $id.sr.mmi $id.target.fa 2>/dev/null
     python3 ../../scripts/iterative_align.py --queries $id.insertion.fa --asm20-index $id.asm20.mmi --sr-index $id.sr.mmi \
       --threads "$th" --out $id.local_segments.tsv
     len=$(grep -v '>' $id.insertion.fa | tr -d '\n' | wc -c)
     preset=asm20; [[ $len -le 500 ]] && preset=sr
     for side in right left; do
       minimap2 -c -x $preset -N 50 --secondary=yes -t "$th" $id.$side.fa $id.insertion.fa > $id.$side.paf 2>/dev/null
     done
     rm -f $id.asm20.mmi $id.sr.mmi $id.local_segments.tsv.pieces.tmp.fa
   done)
  python3 scripts/assemble_local_tandem.py --manifest $O/site_copy_realign/manifest.tsv \
    --genome-segments $I/$c.unresolved.hg38.iterative_segments.tsv --loci $I/hg38.sd_reference.loci.bed \
    --out $O/$c.last_step_calls.tsv > /dev/null
  final_build
  echo "stage4 done: $(tail -n +2 $O/$c.final_annotation.tsv | wc -l) insertions, $(wc -l < $I/$c.final_unresolved.ids.txt) unresolved"
  ;;
*) usage ;;
esac
