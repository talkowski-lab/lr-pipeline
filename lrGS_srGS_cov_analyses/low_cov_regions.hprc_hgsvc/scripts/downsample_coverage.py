WINDOW = 5000  # bp, aggregate 100bp bins into 5kb windows via max (preserve spikes)

DIR = "/Users/xzhao/Downloads/analyses/gnomAD_LR_vcf/final_vcfs/low_cov_benchmark/low_cov_regions.hprc_hgsvc"
IN = f"{DIR}/hgsvc_hprc.cohort.coverage_counts.tsv"
OUT = "/private/tmp/claude-502/-Users-xzhao-Downloads-github/ce274515-e939-4f8a-9d3d-67c167dd0e2b/scratchpad/coverage_downsampled.tsv"

cur_chrom = None
cur_bin = None
cur_max_low = 0
cur_eligible = 0

with open(IN) as f, open(OUT, "w") as out:
    header = f.readline()
    for line in f:
        chrom, start, end, low, elig = line.rstrip("\n").split("\t")
        start = int(start)
        low = int(low)
        elig = int(elig)
        binid = start // WINDOW

        if chrom != cur_chrom or binid != cur_bin:
            if cur_chrom is not None:
                out.write(f"{cur_chrom}\t{cur_bin * WINDOW}\t{cur_max_low}\t{cur_eligible}\n")
            cur_chrom = chrom
            cur_bin = binid
            cur_max_low = low
            cur_eligible = elig
        else:
            if low > cur_max_low:
                cur_max_low = low
            cur_eligible = max(cur_eligible, elig)

    if cur_chrom is not None:
        out.write(f"{cur_chrom}\t{cur_bin * WINDOW}\t{cur_max_low}\t{cur_eligible}\n")

print("done ->", OUT)
