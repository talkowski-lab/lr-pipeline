# read_recurrence

Recurrence (dot) plots of reads or assembly contigs aligned to hg38 across one region, against hg38 and T2T-CHM13.
Run by `wdl_analyses/ReadRecurrence.wdl` (every read in the region) and `wdl_analyses/SplitReadRecurrence.wdl` (only the
split-aligned reads, each against all its aligned segments side by side) (image built from `dockerfiles/Dockerfile.ReadRecurrence`, scripts in
`/opt/read_recurrence`); every script also runs on its own (`-h` prints usage).

## prepare_region.py
hg38 region +/- flank and the region lifted to T2T (UCSC chain hg38 -> CHM13) +/- flank, each as a small FASTA with one
sequence named `chrom:start-end`.
```
python3 prepare_region.py --region chr16:29324175-29638676 --flank 5000 --hg38-ref hg38.fa --t2t-ref chm13v2.0.fa.bgz \
  --chain grch38-chm13v2.chain --out 16p_bp4 [--t2t-region chr16:...] [--lift-step 1000]
```
Outputs: `<out>.hg38.fa`, `<out>.t2t.fa` (+ .fai), `<out>.t2t_region.txt`, `<out>.lift.tsv` (each lifted point).

## list_region_reads.py
Reads/contigs whose non-secondary alignments overlap the region (>= `--min-overlap` bp; `--spanning-only`: covering the
whole region; at most `--max-reads`, largest overlap first), and one small BAM per read with all its primary and
supplementary alignments (SA tag followed).
```
python3 list_region_reads.py --bam asm_h1.bam --region chr16:29324175-29638676 --out hap1 --outdir reads
```
Outputs: `<out>.reads.tsv` (read_name file_id n_alignments strands overlap_bp aligned_start aligned_end spans_region kept),
`<outdir>/<out>.<file_id>.bam` (+ .bai).

## read_ref_recurrence.py
One read against a reference: the read is cut by its CIGAR to the part aligned to `--region` (+/- `--pad`), split into
`--window` bp windows (`--step` apart; default = window), and every window is matched against `--region` +/- `--pad` with
edit distance <= `--max-ed-frac` x window on both strands (exact seed-and-verify). `--orient-to-ref` reverse-complements
reads that match mostly on the reverse strand so that their diagonal runs forward. `--ref` may be a whole genome or a
`chrom:start-end` segment FASTA.
```
python3 read_ref_recurrence.py --bam hap1.ctg.bam --read-name h1tg003138l --region chr16:29324175-29638676 \
  --ref 16p_bp4.hg38.fa --out hap1__hg38__h1tg003138l --pad 5000 --threads 4 --orient-to-ref
```
Outputs: `<out>.hits.tsv` (read_start read_end strand ref_start ref_end edit_distance identity), `<out>.aln.tsv` (the
read's alignment blocks), `<out>.info.tsv` (read, segments, orientation, settings).

## read_t2t_recurrence.py
As read_ref_recurrence.py, against T2T-CHM13: `--region` is lifted through `--chain` (or `--t2t-region` is used) and the
hg38 alignment blocks are lifted for the guide lines. Extra outputs: `<out>.lift.tsv`.
```
python3 read_t2t_recurrence.py --bam hap1.ctg.bam --read-name h1tg003138l --region chr16:29324175-29638676 \
  --chain grch38-chm13v2.chain --t2t-ref 16p_bp4.t2t.fa --t2t-region $(cat 16p_bp4.t2t_region.txt) \
  --out hap1__T2T__h1tg003138l --pad 5000 --orient-to-ref
```

## list_split_reads.py
Split-aligned reads of a region (>= `--min-alignments` primary/supplementary alignments on the region's chromosome, one
overlapping the region by >= `--min-overlap` bp; secondary alignments ignored) and one bundle per read for
read_split_recurrence.py: the read's alignments (BAM), its reference segments (each alignment +/- `--pad`, merged) as a
small FASTA, and the SD and gene tracks cut to them.
```
python3 list_split_reads.py --bam asm_h1.bam --region chr16:28000000-30500000 --ref hg38.fa --out hap1 --outdir bundles \
  [--pad 20000] [--segdup segdups.bed.gz] [--genes genes.gtf.gz] [--segments chr16:S1-E1,chr16:S2-E2]
```
`--segments` gives every read the same fixed reference segments (e.g. BP4 and BP5 SDs +/- 50 kb).
Outputs: `<out>.split_reads.tsv` (read_name file_id n_alignments alignments segments bundled) and
`<outdir>/<out>.<file_id>.tar.gz` (read.bam, segments.fa, segments.txt, read_name.txt, segdup.bed, genes.gtf).

## read_split_recurrence.py
A read/contig with a split alignment (primary + supplementary, e.g. the two sides of a large deletion) against all its
aligned reference segments together: each alignment's reference span +/- `--pad` is a segment (overlaps merged, or
`--segments`), the read is windowed over the span of those alignments +/- `--pad`, and every window is matched against
every segment. The plot puts the segments side by side, so a deletion shows as one diagonal that jumps between panels.
```
python3 read_split_recurrence.py --bam asm_h1.minimap2.bam --read-name h1tg003138l --region chr16:29268226-29300000 \
  --ref hg38.fa --out 16p_BP4-5.h1tg003138l --pad 20000 --threads 8 --orient-to-ref
Rscript plot_read_ref_recurrence.R --prefix 16p_BP4-5.h1tg003138l --ref-label hg38 --out 16p_BP4-5.h1tg003138l.pdf \
  [--segdup segdups.bed.gz] [--genes genes.gtf.gz]
```
Outputs: as read_ref_recurrence.py, with a `segment` column in `<out>.hits.tsv`, `segments` in `<out>.info.tsv`, and
`<out>.alignments.tsv` (the alignments used, in reference and read coordinates). `--all-chroms` also uses alignments on
other chromosomes.

## plot_read_ref_recurrence.R
Plot from any of the scripts' outputs: reference on x (SD and gene tracks underneath), read on y; blue forward, red
reverse complement, windows with > `--max-copies` hits faint grey, the read's alignment blocks as grey lines; several
reference segments are drawn side by side, each panel as wide as its segment.
```
Rscript plot_read_ref_recurrence.R --prefix hap1__hg38__h1tg003138l --out hap1__hg38__h1tg003138l.pdf \
  [--segdup segdups.bed.gz] [--genes genes.gtf.gz] [--ref-label hg38] [--title ...] [--max-copies 5]
  [--segment-names BP4,BP5] [--regions "BP4 SD=chr16:S-E;BP5 SD=chr16:S-E"] [--sections blocks.tsv]
```

## references/16p11_BP4_BP5_SDs.flank50kb.hg38_blocks.tsv
hg38 block layout of the 16p11.2 BP4 and BP5 SDs +/- 50 kb (segment_structure.py on the SD-vs-SD recurrence): A B1 C1
B2 D1 E | C2 B3 D2 F; copies of one block share the letter (B = the NPIPB unit, C and D = the near-identical BP4/BP5
copies). Use as SplitReadRecurrence `blocks_tsv` with
`ref_segments = chr16:29274176-29688676,chr16:30138532-30393458`, `segment_names = BP4,BP5`,
`region_labels = BP4 SD=chr16:29324176-29638676;BP5 SD=chr16:30188532-30343458`.
