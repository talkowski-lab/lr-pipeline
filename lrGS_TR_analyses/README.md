# lrGS_TR_analyses

Tandem-repeat variant (TRV) analyses of the gnomAD long-read cohort VCFs.

## AnalyzeTRVariants.wdl
Per contig (scattered over parallel `vcfs` / `vcf_idxs` / `contigs` arrays):
1. `ExtractTRV` (bcftools): TRV records (`INFO/allele_type == "trv"`) → `<prefix>.<contig>.TRV.vcf.gz`, sample list, and a `bcftools query` table of site fields plus per-sample GT.
2. `AnalyzeTRV` (`analyze_TRV_vcf.py analyze`, standard library only) writes three outputs:
   - **Site BED.** Columns: repeat span (TRID), REF repeat sequence / length / motif count, and per-ALT length, size difference, size difference per motif and motif-count difference vs REF. It also has AC/AF/AN, genic context (coding > UTR > intronic > intergenic) with gene names, and for intergenic TRs the distance to, and gene of, the closest 5'UTR and 3'UTR.
   - **Per-sample size-difference and motif-count-difference matrices.** Site × sample; each cell is `d1,d2` in GT order.
   - **Per-sample summary by genic context.** PASS sites only unless `all_filters`.
3. `ConcatSites` → genome-wide `<prefix>.TRV.sites.bed.gz` (+ .tbi); `MergeSummaries` → genome-wide `<prefix>.TRV.per_sample_summary.tsv`.

Motif count is `INFO/MC_allele` where present (multi-allelic sites). Otherwise it is a greedy exact-match count of MOTIFS; the `mc_source` column records which. Full definitions are in the `analyze_TRV_vcf.py` docstring.

Inputs on Terra (root entity `LR_contig_set`): `vcfs = this.LR_contigs.hprc_hgsvc_vcf_V10`, `vcf_idxs = this.LR_contigs.hprc_hgsvc_vcf_idx_V10`, `contigs = this.LR_contigs.contig`, `gtf` = gnomAD-SV r3 GENCODE v39 GTF, `analysis_script = analyze_TRV_vcf.py`.

## run_analyze_TRV_vcf.sh
Local driver running the same steps: `run_analyze_TRV_vcf.sh <vcf_list> <gencode.gtf.gz> <out_dir>` (gs:// inputs need `GCS_OAUTH_TOKEN`).

## Bed-table analyses (from the VEP-parsed annotated bed files)
These run on `gnomAD_LR.{cohort}.vep_parsed.annotated.bed.gz`. Outputs and findings are described in `final_vcfs/STR_catalog/SESSION_SUMMARY.md`.
- `extract_TRV_sites.sh <input.bed.gz> <output.bed.gz>`: keeps rows whose ID contains `-TRV-`. Output is bgzipped and tabix-indexed, plus FILTER counts.
- `annotate_TRV_genic_context.sh <TRV.bed.gz> <gencode.gtf.gz> <out.tsv.gz>`: genic context (coding > UTR > intronic > intergenic) via `bedtools map`, plus locus length, minimum motif length and per-ALT length change.
- `plot_TRV_genic_context.py --inputs LABEL=tsv.gz ... --out-prefix P`: counts, locus-length ECDF, allele length-change bins and motif length by genic context.
- `plot_TRV_motif_counts.py --inputs LABEL=GENIC_TSV,TRV_BED ... --out-prefix P`: non-ref site counts vs motif size, cohort × context panels.
- `plot_TRV_motif_length_line.py --inputs LABEL=GENIC_TSV,TRV_BED ... --colors LABEL=HEX ... --context coding --out-prefix P`: motif-length distribution line plot with cohorts overlaid.
- `classify_coding_TRV_alleles.py --genic-tsv G --trv-bed B --gtf GTF --label L --out-prefix P`: classifies coding TRV ALT alleles as frameshift / in-frame del / in-frame ins / no length change / partial CDS, with allele, site and gene counts by AC.
- `summarize_coding_STR_per_gene.py --alleles <classify output alleles.tsv.gz> --out-prefix P [--max-motif 6] [--min-ac 1]`: coding STR sites per gene (motif ≤ 6 bp, locus inside one CDS block), split by site class LoF (any frameshift allele) / in-frame / no length change; writes per-site, per-gene and summary tables.
- `plot_STR_variability.py --trv-bed TRV.bed.gz --label L --out-prefix P [--max-motif 6] [--min-copies 2] [--min-loci 100]`: STR locus variability (distinct ALT alleles with AC > 0; non-ref allele frequency ΣAC/AN) vs STR size, copy number and motif length, split by TR catalog (SOURCE). Writes line plots (vs size/copies, and vs motif at fixed copy number or size), heatmaps and a binned TSV.
