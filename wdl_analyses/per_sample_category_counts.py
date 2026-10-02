#!/usr/bin/env python3
"""Per-sample pLoF/missense/synonymous/intronic/intergenic counts and pLoF gene
attribution, from one chromosome-shard VCF annotated with VEP's raw INFO/vep
CSQ-style string (Allele|Consequence|IMPACT|SYMBOL|Gene|...) plus the existing
INFO/allele_type + INFO/allele_length annotations.

Restricted to FILTER=PASS. A variant can belong to more than one category (the
categories are not mutually exclusive, matching how the existing gnomAD LR
Analyze_* category files were built: e.g. a variant can be both HIGH impact on
one overlapping transcript and missense on another).

pLoF is defined as any transcript annotation carrying one of VEP's HIGH-impact
consequence terms (Ensembl's documented HIGH tier: transcript_ablation,
splice_acceptor_variant, splice_donor_variant, stop_gained, frameshift_variant,
stop_lost, start_lost, transcript_amplification, feature_elongation,
feature_truncation), rather than a direct check against an IMPACT sub-field
position, since the IMPACT sub-field is not independently exposed in this
dataset's bed/VCF outputs. pLoF genes are only credited from the specific
transcript annotations that carry a HIGH-impact term, not all overlapping
transcripts' symbols.

pLoF SNV/indel/SV uses INFO/allele_type + abs(INFO/allele_length), with the
indel/SV boundary at 50bp (matching per_sample_type_counts.sh's size bins).

Usage: per_sample_category_counts.py <in.vcf.gz> <out_prefix>
Writes <out_prefix>.category_counts.tsv with columns:
  sample,
  n_plof_snv, n_plof_indel_del, n_plof_indel_ins, n_plof_sv_del, n_plof_sv_ins,
  n_missense, n_synonymous, n_intronic, n_intergenic,
  plof_snv_genes, plof_indel_genes, plof_sv_genes (comma-joined, sorted, unique),
  n_unique_genes_plof_snv, n_unique_genes_plof_indel, n_unique_genes_plof_sv,
  n_unique_genes_plof_any (union of the three gene sets)
"""
import subprocess
import sys
from collections import defaultdict

HIGH_IMPACT_TERMS = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification", "feature_elongation", "feature_truncation",
}
COUNT_CATEGORIES = [
    "n_plof_snv", "n_plof_indel_del", "n_plof_indel_ins", "n_plof_sv_del", "n_plof_sv_ins",
    "n_missense", "n_synonymous", "n_intronic", "n_intergenic",
]
GENE_CATEGORIES = ["plof_snv_genes", "plof_indel_genes", "plof_sv_genes"]


def parse_vep(vep_field):
    """Return (consequence_term_set, high_impact_gene_symbols) for one variant."""
    all_terms = set()
    high_impact_genes = set()
    for entry in vep_field.split(","):
        fields = entry.split("|")
        if len(fields) < 4:
            continue
        consequence_terms = fields[1].split("&") if fields[1] else []
        symbol = fields[3]
        all_terms.update(consequence_terms)
        if symbol and any(t in HIGH_IMPACT_TERMS for t in consequence_terms):
            high_impact_genes.add(symbol)
    return all_terms, high_impact_genes


def is_non_ref(gt):
    if gt in (".", "./.", ".|."):
        return False
    return any(a not in ("0", ".") for a in gt.replace("|", "/").split("/"))


def main():
    vcf, out_prefix = sys.argv[1], sys.argv[2]

    samples = subprocess.run(
        ["bcftools", "query", "-l", vcf], check=True, capture_output=True, text=True
    ).stdout.split()

    counts = {s: defaultdict(int) for s in samples}
    genes = {s: {g: set() for g in GENE_CATEGORIES} for s in samples}

    query_fmt = "%CHROM\t%POS\t%INFO/allele_type\t%INFO/allele_length\t%INFO/vep[\t%GT]\n"
    proc = subprocess.Popen(
        ["bcftools", "query", "-i", 'FILTER="PASS"', "-f", query_fmt, vcf],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    )

    n_rows = 0
    for line in proc.stdout:
        n_rows += 1
        fields = line.rstrip("\n").split("\t")
        allele_type, allele_length, vep_field = fields[2], fields[3], fields[4]
        gts = fields[5:]

        terms, plof_genes = parse_vep(vep_field)
        is_plof = bool(plof_genes)
        is_missense = "missense_variant" in terms
        is_synonymous = "synonymous_variant" in terms
        is_intronic = "intron_variant" in terms
        is_intergenic = "intergenic_variant" in terms

        if not (is_plof or is_missense or is_synonymous or is_intronic or is_intergenic):
            continue

        plof_count_col = None
        plof_gene_col = None
        if is_plof:
            length = abs(int(allele_length)) if allele_length not in (".", "") else 0
            if allele_type == "snv":
                plof_count_col, plof_gene_col = "n_plof_snv", "plof_snv_genes"
            elif allele_type == "del":
                plof_count_col = "n_plof_indel_del" if length < 50 else "n_plof_sv_del"
                plof_gene_col = "plof_indel_genes" if length < 50 else "plof_sv_genes"
            elif allele_type == "ins":
                plof_count_col = "n_plof_indel_ins" if length < 50 else "n_plof_sv_ins"
                plof_gene_col = "plof_indel_genes" if length < 50 else "plof_sv_genes"

        for sample, gt in zip(samples, gts):
            if not is_non_ref(gt):
                continue
            if plof_count_col is not None:
                counts[sample][plof_count_col] += 1
                genes[sample][plof_gene_col].update(plof_genes)
            if is_missense:
                counts[sample]["n_missense"] += 1
            if is_synonymous:
                counts[sample]["n_synonymous"] += 1
            if is_intronic:
                counts[sample]["n_intronic"] += 1
            if is_intergenic:
                counts[sample]["n_intergenic"] += 1

    returncode = proc.wait()
    if returncode != 0:
        sys.exit(f"bcftools query failed with exit code {returncode}")

    out_path = f"{out_prefix}.category_counts.tsv"
    header = (
        ["sample"] + COUNT_CATEGORIES + GENE_CATEGORIES
        + ["n_unique_genes_plof_snv", "n_unique_genes_plof_indel", "n_unique_genes_plof_sv", "n_unique_genes_plof_any"]
    )
    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for s in samples:
            row = [s] + [str(counts[s][c]) for c in COUNT_CATEGORIES]
            row += [",".join(sorted(genes[s][g])) for g in GENE_CATEGORIES]
            row += [str(len(genes[s][g])) for g in GENE_CATEGORIES]
            row.append(str(len(set.union(*genes[s].values()))))
            out.write("\t".join(row) + "\n")

    print(f"Processed {n_rows} PASS variant rows across {len(samples)} samples; wrote {out_path}")


if __name__ == "__main__":
    main()
