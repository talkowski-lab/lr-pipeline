# lrGS_TR_analyses

Tandem-repeat variant (TRV) analyses of the gnomAD long-read cohort callsets (hprc_hgsvc, AoU_I). A TRV is a record whose
variant ID contains `-TRV-`; its TRID gives the repeat span and motifs. Outputs and findings are described in
`final_vcfs/STR_catalog/SESSION_SUMMARY.md` (local analysis folder).

## Workflows (Terra)

### ExtractTRVSitesBed.wdl
Per contig (scattered over parallel `vcfs` / `vcf_idxs` / `contigs`), TRV records (ID containing `-TRV-`) are written to a
BED; contigs are merged into one sorted, bgzipped, tabix-indexed BED.

| Column | Content |
|---|---|
| `#chrom` `start` `end` | repeat span: start = POS, end = POS + len(REF) − 1 (the TRID span) |
| `ID` `MOTIFS` | variant ID, repeat motifs |
| `n_alt` | number of ALT alleles |
| `AF_total` `AC_total` `AN` | overall AF = ΣAC / AN, summed AC, AN |
| `ref_size` `alt_sizes` | REF size and per-ALT sizes (VCF padding base excluded) |
| `AF_list` `AC_list` | per-ALT AF and AC as in the VCF |
| `FILTER` | VCF FILTER |

Inputs (root entity `LR_contig_set`): `vcfs = this.LR_contigs.<vcf column>`, `vcf_idxs = this.LR_contigs.<idx column>`,
`contigs = this.LR_contigs.contig`, `prefix`, `bcftools_docker`.

### AnalyzeTRVariants.wdl
Per contig:
1. `ExtractTRV` (bcftools): TRV records → `<prefix>.<contig>.TRV.vcf.gz`, sample list, and a `bcftools query` table of site
   fields plus per-sample GT.
2. `AnalyzeTRV` (`analyze_TRV_vcf.py analyze`, standard library only):
   - Site BED: repeat span, REF repeat sequence / length / motif count; per-ALT length, size difference, size difference per
     motif and motif-count difference; AC/AF/AN; genic context (coding > UTR > intronic > intergenic) with gene names, and
     for intergenic TRs the distance to and gene of the closest 5'UTR and 3'UTR.
   - Per-sample size-difference and motif-count-difference matrices (site × sample, `d1,d2` in GT order).
   - Per-sample summary by genic context (PASS only unless `all_filters`).
3. `ConcatSites` → genome-wide `<prefix>.TRV.sites.bed.gz` (+ .tbi); `MergeSummaries` → `<prefix>.TRV.per_sample_summary.tsv`.

Motif count is `INFO/MC_allele` where present (multi-allelic sites), otherwise a greedy exact-match count of MOTIFS
(`mc_source` column). Inputs: `vcfs = this.LR_contigs.hprc_hgsvc_vcf_V10`, `vcf_idxs = this.LR_contigs.hprc_hgsvc_vcf_idx_V10`,
`contigs = this.LR_contigs.contig`, `gtf` (gnomAD-SV r3 GENCODE v39), `analysis_script = analyze_TRV_vcf.py`.
`run_analyze_TRV_vcf.sh <vcf_list> <gencode.gtf.gz> <out_dir>` runs the same steps locally (gs:// inputs need `GCS_OAUTH_TOKEN`).

## Catalog and annotation (from `gnomAD_LR.{cohort}.vep_parsed.annotated.bed.gz`)
- `extract_TRV_sites.sh <input.bed.gz> <output.bed.gz>`: rows whose ID contains `-TRV-`; bgzipped, tabix-indexed, FILTER counts.
- `annotate_TRV_genic_context.sh <TRV.bed.gz> <gencode.gtf.gz> <out.tsv.gz>`: genic context via `bedtools map`, locus length,
  minimum motif length and per-ALT length change.
- `build_trv_locus_table.py --trv-bed B --genic-tsv G --gtf GTF --out loci.tsv.gz`: one row per TRV with TRID span, catalog,
  motif, copies, GC and purity of the REF repeat, ALT alleles (AC > 0), non-ref AF, max |ALT − REF|, in-frame / frameshift
  allele counts, genic context (coding = inside one CDS block; coding_partial otherwise). Input to most plots below.
- `annotate_trv_sites.py --sites <prefix>.TRV.sites.bed.gz --intergenic intergenic_regions_annotated.bed.gz --gene-info
  gene_information_loeufv4_pHaplo_pTriplo...txt.gz --out P.bed.gz`: appends to the AnalyzeTRVariants site BED, for intergenic
  TRVs, the overlapping intergenic region and its category (Tiny / Small / Medium / Desert, or between_regions /
  chrom_not_in_region_file); for coding and intronic TRVs, LOEUF_v4 with its percentile and decile per gene.
- `merge_trv_tables.sh <sites.annotated.bed.gz> <vep_parsed.annotated.TRV.bed.gz> <out.bed.gz>`: one merged TRV table: the
  annotated site BED plus the columns only in the VEP-parsed TRV bed (full REF / ALT, SOURCE, REGION, nhomref / nhet / nhomalt,
  telomere / centromere / N-masked / low-coverage flags). Checks that both inputs list the same IDs in the same order.
- `trv_gene_table.py --trv <merged TRV.bed.gz> [--gene-info G] [--min-motif 1 --max-motif 6] --out genes.tsv`: per-gene STR
  counts (non-ref PASS): coding, UTR, intronic, and intergenic STRs assigned to the gene with the nearest 5'UTR / 3'UTR
  (split by side, with median distance), plus LOEUF_v4.
- `count_variant_overview.sh <annotated.bed.gz> <label> <out.tsv>`: variant counts by class (SNV, del/ins split at 50 bp, TRV,
  others), all and PASS.

## Counts and size distributions
- `plot_TRV_genic_context.py`: counts, locus-length ECDF, allele length-change bins and motif length by genic context.
- `plot_TRV_motif_counts.py`: non-ref site counts vs motif size, cohort × context panels.
- `plot_TRV_motif_length_line.py`: motif-length distribution line plot, cohorts overlaid.
- `plot_coding_TRV_motif_counts.py [--context coding|UTR|intronic|intergenic]`: TRV counts by motif length in one context,
  cohorts overlaid.
- `plot_TRV_motif_counts_by_context.py`: one cohort, counts per exact motif length in 1–10 / 11–100 / 101–1000 bp panels,
  one line per context.
- `plot_str_ref_length_by_context.py --trv <merged TRV.bed.gz> [--size ref_allele|trid_span] --out-prefix P`: STR size
  (REF allele length) across coding, intronic, intergenic non-desert and desert: exact sizes to 150 bp (log y) and an ECDF.
- `plot_TRV_size_distribution.py [--context coding]`: repeat-size distribution (exact 1–100 bp, 10-bp bins to 1 kb, ECDF).

## Mutability and constraint
- `plot_STR_variability.py`: STR variability (distinct ALT alleles; non-ref AF) vs size, copy number and motif, by catalog.
- `plot_STR_constraint_by_context.py --motif M`: variability vs copy number by genic context for one motif length.
- `plot_STR_constraint_allele_class.py --motif M`: the same, split into in-frame / frameshift / no-length-change alleles.
- `plot_STR_constraint_by_loeuf.py --motif M`: coding STRs by gene LOEUF tertile with an intergenic reference.
- `plot_nalt_vs_length_by_context.py --trv <merged TRV.bed.gz> --label L --out-prefix P [--min-copies 2]`: distinct ALT
  alleles per site vs repeat length, one PDF each for coding, intronic, intergenic non-desert (Tiny / Small / Medium) and
  intergenic desert; site density plus median / IQR per length bin by motif class.
- `plot_top_nalt_size_dist.py --trv <merged TRV.bed.gz> --top 10 [--max-motif M] [--contexts ...] --out-prefix P`: for the
  loci with the most distinct ALT alleles, the haplotype distribution of ALT size − REF size (`--mode size_diff`) or of
  allele length with REF included (`--mode length`, y = allele frequency); ALT weighted by AC, REF by AN − ΣAC.
- `str_density_by_context.py --trv <merged TRV.bed.gz> --gtf GTF --intergenic intergenic_regions_annotated.bed.gz --out-prefix P`:
  STRs per Mb (autosomes) in coding (CDS bp), intronic (transcript minus exons), intergenic non-desert and desert
  (regions minus MANE transcripts), overall and by motif size.
- `plot_alt_size_range_by_context.py --trv <merged TRV.bed.gz> --motif 3 --out-prefix P`: per-site ALT allele size range
  (max − min carried ALT length) for one motif size across coding, intronic, intergenic non-desert and desert.
- `classify_coding_TRV_alleles.py`: coding TRV ALT alleles as frameshift / in-frame del / in-frame ins / no length change /
  partial CDS, with allele, site and gene counts by AC.
- `summarize_coding_STR_per_gene.py [--max-motif 6] [--min-ac 1]`: coding STR sites per gene, split into LoF / in-frame /
  no length change.
- `make_trv_deck_figures.py`: deck figures and the coding constraint model (Poisson GLM of ALT allele count on size,
  motif × catalog, GC and purity; obs/exp per gene and per LOEUF decile).
- `plot_coding_oe_motif_subset.py --motifs M ... [--compare-motifs M ...]`: obs/exp by LOEUF decile for a motif-length set
  (imports `make_trv_deck_figures.py`; run with this folder on PYTHONPATH).

## Novelty relative to short-read resources
- `count_variants_in_trv.sh` + `trv_novelty.py`: routine variants inside each TRV (TRID column) and how many lack both
  dbGaP and gnomAD v4 matches; per-TRV BED and proportions of novel TRVs / variants by size, motif and context.
- `extract_indels_in_trvs.sh` + `summarize_coding_trv_indels.py`: indels inside coding TRVs not captured by srGS, binned by
  indel size and by TRV reference size.

## Comparison with gnomAD short-read STRs
- `fetch_gnomad_str.py --dataset gnomad_r4 --out-prefix P`: downloads the gnomAD browser STR loci (GraphQL API; throttled,
  retries on HTTP 429) with allele-size distributions summed over ancestry / sex / quality cohorts.
- `compare_gnomad_str_to_lr.py --loci P.loci.tsv --allele-sizes P.allele_sizes.tsv --trv LABEL=<TRV table> ... --colors C ...
  --out-prefix Q`: one or more long-read datasets (merged TRV table or AnalyzeTRVariants sites bed, tabix-indexed);
  matches each locus to the overlapping long-read TRV (largest overlap), converts long-read alleles to repeat units
  (rounded reference units + size difference / unit length; AC-weighted, REF by AN − ΣAC), and compares distributions
  (per-locus panels, summary scatter of median / p99 / max, KS distance, match_flag for wider / partial TRVs).
- `split_sr_matched_strs.py --loci ... --allele-sizes ... --comparison Q.tsv --long-read-alleles Q.long_read_alleles.tsv
  --trv LABEL=NAME=path ... --prefix X --out-dir D`: splits the comparison into parallel per-dataset files X.sr.* and
  X.NAME.* (loci table, allele distribution in repeat units, matched TRV records) for side-by-side comparison.
- `extract_trv_genotypes.py --trv-bed X.NAME.TRV.bed --vcf-map contigs.tsv [--loci X.sr.loci.tsv] --out-prefix P`: pulls
  the per-sample genotypes of a TRV set from per-contig VCFs (local or gs://, GCS_OAUTH_TOKEN); writes P.vcf.gz and
  P.genotypes.tsv (GT, per-haplotype allele length and gnomAD SR-scale units, TRGT AL / ALLR / SD / MC / AP / PS).
- `flag_str_outliers.py --input LABEL=genotypes=P.genotypes.tsv LABEL2=allele_dist=X.NAME.allele_dist.tsv ... [--loci
  LABEL=X.NAME.loci.tsv ...] [--sr-allele-dist X.sr.allele_dist.tsv] --k 10 --out-prefix O`: per site and cohort, flags
  alleles outside median ± K × IQR (IQR floored at --min-iqr, default 1 unit). Genotype inputs give per-sample outliers
  (with TRGT SD / AP); allele-count inputs (aggregate only, e.g. AoU) give outlier allele sizes with counts. Writes a
  combined per-site table, per-cohort outlier tables and side-by-side per-site strip plots.
- `sr_envelope_outliers.py --trv <merged TRV.bed.gz> --sr-loci X.NAME.loci.tsv --sr-allele-dist X.sr.allele_dist.tsv
  --out-prefix O`: trains a 95% quantile regression of the gnomAD SR reach (SR max − SR median, bp) on motif length at the
  SR-matched loci (LOO-calibrated; LR spread features compared in .model_selection.tsv), then flags rare long-read alleles
  (AC ≤ 5) beyond LR median + predicted reach at every single-motif 2-6 bp STR, classifying each allele by sequence as
  motif_expansion, other_repeat, homopolymer_run or non_repeat_insertion.
- `plot_allele_motif_structure.py --trv <TRV table> --trv-id ID [--k 4 5 6] --out-prefix O`: discovers the repeat units of
  all alleles at one TRV (tandem primitive k-mers, rotations merged) and plots each allele (REF, then ALTs by dominant motif
  and length) painted by motif, with AC on the right; writes per-allele composition and per-motif tables.
- `run_allele_motif_structure_sr_loci.sh <main loci.tsv> <main TRV> <genotypes> <ancestry> <SR allele_sizes> <second loci.tsv>
  <second label> <second TRV> <out_dir> <tag>`: runs plot_allele_motif_structure.py for every SR locus with a matched TRV
  (second cohort in the top panel when matched) and merges the pages into all_SR_loci.<tag>.pdf.
- `str_locus_allele_distribution.py --trv LABEL=path ... --trv-id ID --motif M --sr-sizes P.allele_sizes.tsv --sr-id G
  --out-prefix Q`: one locus in detail: gnomAD SR repeat-unit counts vs long-read alleles measured three ways (allele
  length / motif length, longest pure motif run, motif copies with interruptions allowed); per-allele table, summary, PDF.

## Single-locus tools
- `plot_trv_locus.py --trv <merged TRV.bed.gz> --gtf GTF --catalog TRExplorer=... Vamos=... (--gene G | --region R |
  --gene-table T --count-column n_coding --min-count 10) [--max-motif 6] --out P.pdf [--out-dir D]`: locus plot with one
  track per TR catalog (observed non-ref loci coloured, unobserved grey; total ALT alleles and observed / catalog loci on the
  right); per observed TRV: distinct ALT alleles, total ALT AC, allele length range (shortest to longest carried ALT, REF
  tick) and TR locus size; and gene structure (exons as boxes, CDS taller, introns as lines with strand arrows).
- `vcf_region_allele_lengths.py --vcf V --region chr:start-end --motif-len N --out P.tsv`: per-sample allele lengths of a
  repeat region rebuilt from the variants called inside it (per haplotype: length, estimated unit count, SV size vs REF,
  variant IDs with inferred type and size). Variant types come from REF/ALT (SVLEN for symbolic alleles); no INFO field needed.
- `match_variants_to_trv.py --variants V.tsv --trv-id ID --vcf L=locus.vcf.gz ...`: match external indel calls inside a TR
  locus to gnomAD-LR by applied-haplotype sequence.
- `trv_length_vs_methylation.py --vcf locus.vcf.gz --trv-id ID --meth-bed region_meth.bed.gz`: allele length vs per-haplotype
  CpG methylation at one locus.

## Slide deck
- `build_trv_deck_pptx.py --figures-dir deck/figures --out deck/TRV_STR_deck.pptx`: PowerPoint version of the TRV / STR
  summary deck built from `make_trv_deck_figures.py` output.
