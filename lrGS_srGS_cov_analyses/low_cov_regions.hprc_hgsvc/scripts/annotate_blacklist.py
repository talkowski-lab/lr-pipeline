"""Annotate an hgsvc_hprc.blacklist.p{N}.raw.bed with size + overlap bp against
telomere, centromere, N-masked, segmental-duplication, and simple-repeat
tracks. A region can overlap multiple tracks simultaneously; overlap size is
reported independently per track rather than a single priority label. A
`primary_annotation` column also assigns each region to exactly one category
using priority order centromere > telomere > N-mask > simple_repeat > seg_dup
(first category with nonzero overlap wins; "none" if it matches none of the
5 tracks).

Usage: python3 annotate_blacklist.py 90   (or 70, or 50)
"""
import subprocess
import sys

DIR = "/Users/xzhao/Downloads/analyses/gnomAD_LR_vcf/final_vcfs/low_cov_benchmark/low_cov_regions.hprc_hgsvc"
PARENT = "/Users/xzhao/Downloads/analyses/gnomAD_LR_vcf/final_vcfs/low_cov_benchmark"

TRACKS = {
    "telomere": (f"{PARENT}/hg38.telomere.bed.gz", True),
    "centromere": (f"{PARENT}/hg38.centromere.bed.gz", True),
    "nmasked": (f"{PARENT}/hg38.n_masked_regions.bed.gz", True),
    "segdup": (f"{PARENT}/hg38.SegDup.sorted.merged.bed", False),
    "simprep": (f"{PARENT}/hg38.SimpRep.sorted.merged.bed", False),
}


def primary_annotation(telomere_bp, centromere_bp, n_masked_bp, segdup_bp, simprep_bp):
    if centromere_bp > 0:
        return "centromere"
    if telomere_bp > 0:
        return "telomere"
    if n_masked_bp > 0:
        return "N-mask"
    if simprep_bp > 0:
        return "simple_repeat"
    if segdup_bp > 0:
        return "seg_dup"
    return "none"


def annotate(threshold_pct):
    raw_path = f"{DIR}/hgsvc_hprc.blacklist.p{threshold_pct}.raw.bed"
    indexed_path = f"{DIR}/hgsvc_hprc.blacklist.p{threshold_pct}.raw.indexed.bed"
    out_path = f"{DIR}/hgsvc_hprc.blacklist.p{threshold_pct}.annotated.bed"

    regions = []
    with open(raw_path) as f:
        for line in f:
            chrom, start, end = line.rstrip("\n").split("\t")
            regions.append([chrom, int(start), int(end)])

    with open(indexed_path, "w") as f:
        for i, (chrom, start, end) in enumerate(regions):
            f.write(f"{chrom}\t{start}\t{end}\t{i}\n")

    overlaps = {t: [0] * len(regions) for t in TRACKS}
    for name, (path, is_gz) in TRACKS.items():
        cat_cmd = ["gzcat", path] if is_gz else ["cat", path]
        cat_p = subprocess.run(cat_cmd, capture_output=True, text=True, check=True)
        sort_p = subprocess.run(
            ["sort", "-k1,1", "-k2,2n"], input="\n".join(
                "\t".join(line.split("\t")[:3]) for line in cat_p.stdout.splitlines()
            ), capture_output=True, text=True, check=True,
        )
        intersect_p = subprocess.run(
            ["bedtools", "intersect", "-a", indexed_path, "-b", "-", "-wao"],
            input=sort_p.stdout, capture_output=True, text=True, check=True,
        )
        for line in intersect_p.stdout.splitlines():
            fields = line.split("\t")
            idx = int(fields[3])
            ov = int(fields[-1])
            overlaps[name][idx] += ov

    with open(out_path, "w") as out:
        out.write("chrom\tstart\tend\tsize\ttelomere_bp\tcentromere_bp\tn_masked_bp\tsegdup_bp\tsimprep_bp\tprimary_annotation\n")
        for i, (chrom, start, end) in enumerate(regions):
            ov = [overlaps[t][i] for t in TRACKS]  # order: telomere, centromere, nmasked, segdup, simprep
            telomere_bp, centromere_bp, n_masked_bp, segdup_bp, simprep_bp = ov
            cat = primary_annotation(telomere_bp, centromere_bp, n_masked_bp, segdup_bp, simprep_bp)
            row = [chrom, start, end, end - start] + ov + [cat]
            out.write("\t".join(str(x) for x in row) + "\n")

    print(f"wrote {out_path}  ({len(regions)} regions)")


if __name__ == "__main__":
    thresholds = sys.argv[1:] if len(sys.argv) > 1 else ["90", "70", "50"]
    for t in thresholds:
        annotate(t)
