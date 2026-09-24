version 1.0

import "../wdl/utils/Structs.wdl"

task IndexVcf {
    input {
        File vcf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        ln -s ~{vcf} ~{prefix}.vcf.gz
        bcftools index -c ~{prefix}.vcf.gz
    >>>

    output {
        File vcf_csi = "~{prefix}.vcf.gz.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

# Rewrites a fully-phased diploid VCF into a pseudo-haploid VCF with 2x the
# samples: each original sample "S" becomes "S_hap1" and "S_hap2", each
# genotype encoded as a homozygous pseudo-diploid call (e.g. allele "1" on a
# haplotype becomes "1/1"; a missing allele becomes "./."). This lets plink
# and SAIGE - both diploid-oriented - operate on haplotypes unchanged.
# Column names match the *_hap1 / *_hap2 convention used in the per-haplotype
# methylation bed files (e.g. hprc_methylated.chr22.haplotype.bed.gz), so
# downstream matching is by sample-name string, not column order.
#
# Assumes genotypes are phased ("|"-separated). An unphased genotype ("/")
# is still split left/right allele but the hap1/hap2 assignment in that case
# is arbitrary, not a true haplotype - only meaningful for phased input.
task SplitPhasedVcfToHaplotypes {
    input {
        File vcf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} | awk -v OFS='\t' '
            /^##/ { print; next }
            /^#CHROM/ {
                printf "%s", $1
                for (i = 2; i <= 9; i++) { printf "\t%s", $i }
                for (i = 10; i <= NF; i++) { printf "\t%s_hap1\t%s_hap2", $i, $i }
                printf "\n"
                next
            }
            {
                nf = split($9, fmt, ":")
                gtidx = 1
                for (k = 1; k <= nf; k++) { if (fmt[k] == "GT") { gtidx = k; break } }

                printf "%s", $1
                for (i = 2; i <= 8; i++) { printf "\t%s", $i }
                printf "\tGT"
                for (i = 10; i <= NF; i++) {
                    split($i, subfields, ":")
                    gt = subfields[gtidx]
                    n = split(gt, alleles, /[|\/]/)
                    a = (n >= 1 ? alleles[1] : ".")
                    b = (n >= 2 ? alleles[2] : a)
                    h1 = (a == "." ? "./." : a "/" a)
                    h2 = (b == "." ? "./." : b "/" b)
                    printf "\t%s\t%s", h1, h2
                }
                printf "\n"
            }
        ' | bgzip -c > ~{prefix}.haplotypes.vcf.gz

        bcftools index -c ~{prefix}.haplotypes.vcf.gz
    >>>

    output {
        File haplotype_vcf = "~{prefix}.haplotypes.vcf.gz"
        File haplotype_vcf_csi = "~{prefix}.haplotypes.vcf.gz.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 4 * ceil(size(vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

# LD-prune the contig VCF (plink --indep-pairwise) and materialize the
# pruned marker set as a plink bfile. This bfile is reused, unchanged, as
# the --plinkFile for every locus's SAIGE step1 variance-ratio estimation
# on this contig, and as the input to CreateSparseGRM.
task LdPruneAndExtract {
    input {
        File vcf
        String prefix
        Int window_size_kb
        Int step_size
        Float r2_threshold
        String vcf_half_call
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        plink --vcf ~{vcf} --vcf-half-call ~{vcf_half_call} \
            --indep-pairwise ~{window_size_kb} ~{step_size} ~{r2_threshold} \
            --out ~{prefix}.prune

        plink --vcf ~{vcf} --vcf-half-call ~{vcf_half_call} \
            --extract ~{prefix}.prune.prune.in --keep-allele-order --make-bed \
            --out ~{prefix}.pruned
    >>>

    output {
        File bed = "~{prefix}.pruned.bed"
        File bim = "~{prefix}.pruned.bim"
        File fam = "~{prefix}.pruned.fam"
        File prune_in = "~{prefix}.prune.prune.in"
        File prune_out = "~{prefix}.prune.prune.out"
        File prune_log = "~{prefix}.prune.log"
        File extract_log = "~{prefix}.pruned.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 25,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CreateSparseGRM {
    input {
        File bed
        File bim
        File fam
        String prefix
        Int num_random_markers
        Float relatedness_cutoff
        Float min_maf_for_grm
        Float max_missing_rate_for_grm
        Int n_threads
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        ln -s ~{bed} ~{prefix}.bed
        ln -s ~{bim} ~{prefix}.bim
        ln -s ~{fam} ~{prefix}.fam

        createSparseGRM.R \
            --plinkFile=~{prefix} \
            --outputPrefix=~{prefix}.sparseGRM \
            --numRandomMarkerforSparseKin=~{num_random_markers} \
            --relatednessCutoff=~{relatedness_cutoff} \
            --minMAFforGRM=~{min_maf_for_grm} \
            --maxMissingRateforGRM=~{max_missing_rate_for_grm} \
            --nThreads=~{n_threads}
    >>>

    output {
        File sparse_grm = glob("~{prefix}.sparseGRM*.sparseGRM.mtx")[0]
        File sparse_grm_samples = glob("~{prefix}.sparseGRM*.sparseGRM.mtx.sampleIDs.txt")[0]
    }

    RuntimeAttr default_attr = object {
        cpu_cores: n_threads,
        mem_gb: 16,
        disk_gb: 5 * ceil(size(bed, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

# Filters a wide-format methylation bed file (chrom, start, end, then one
# value column per sample, "." for missing) down to rows whose per-row call
# rate meets the threshold. Emits a gzipped TSV with site_id/call_rate
# columns prepended, plus a plain-text count of qualifying rows so the
# workflow can compute how many chunks to scatter.
task FilterMethylationSites {
    input {
        File methylation_bed
        Float call_rate_threshold
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat <<'PYEOF' > filter_sites.py
import csv
import gzip
import sys

in_path = sys.argv[1]
threshold = float(sys.argv[2])
out_path = sys.argv[3]
count_path = sys.argv[4]

MISSING = {".", "NA", ""}

n_qualifying = 0
with gzip.open(in_path, "rt") as fin, gzip.open(out_path, "wt") as fout:
    reader = csv.reader(fin, delimiter="\t")
    header = next(reader)
    header[0] = header[0].lstrip("#")
    sample_ids = header[3:]
    writer = csv.writer(fout, delimiter="\t")
    writer.writerow(["chrom", "start", "end", "site_id", "call_rate"] + sample_ids)

    for row in reader:
        chrom, start, end = row[0], row[1], row[2]
        values = row[3:]
        total = len(values)
        if total == 0:
            continue
        missing = sum(1 for v in values if v in MISSING)
        call_rate = (total - missing) / total
        if call_rate >= threshold:
            site_id = f"{chrom}_{start}_{end}"
            writer.writerow([chrom, start, end, site_id, f"{call_rate:.6f}"] + values)
            n_qualifying += 1

with open(count_path, "w") as f:
    f.write(str(n_qualifying) + "\n")
PYEOF
        python3 filter_sites.py ~{methylation_bed} ~{call_rate_threshold} ~{prefix}.filtered_sites.tsv.gz ~{prefix}.n_sites.txt
    >>>

    output {
        File filtered_sites = "~{prefix}.filtered_sites.tsv.gz"
        Int n_sites = read_int("~{prefix}.n_sites.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(methylation_bed, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

# Runs every qualifying methylation site on the contig: for each, SAIGE
# step1_fitNULLGLMM.R (sparse-GRM null model for that site's phenotype)
# followed by step2_SPAtests.R restricted to a +/- cis_window region around
# the site via --rangestoIncludeFile. Sites are parallelized across CPU
# cores within this one task via Python multiprocessing, rather than as
# separate Cromwell-scheduled shards - a nested WDL scatter (one level for
# contigs, one for site chunks) hit a Cromwell input-resolution bug on
# optional inputs (covariates_file, runtime_attr_override) referenced two
# scatter levels deep ("Failed to lookup input value for required input"
# even though both are declared optional), so parallelism across sites is
# kept inside this task instead of a second WDL scatter level.
task RunCisMeQTLContig {
    input {
        File filtered_sites
        String contig
        File vcf
        File vcf_csi
        File pruned_bed
        File pruned_bim
        File pruned_fam
        File sparse_grm
        File sparse_grm_samples
        Float relatedness_cutoff
        Int cis_window
        File? covariates_file
        String covar_col_list
        String qcovar_col_list
        Int min_samples_per_site
        String vcf_field
        Boolean inv_normalize
        Int n_parallel_workers
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        ln -s ~{vcf} genotypes.vcf.gz
        ln -s ~{vcf_csi} genotypes.vcf.gz.csi
        ln -s ~{pruned_bed} variance_ratio_markers.bed
        ln -s ~{pruned_bim} variance_ratio_markers.bim
        ln -s ~{pruned_fam} variance_ratio_markers.fam
        ln -s ~{sparse_grm} sparseGRM.mtx
        ln -s ~{sparse_grm_samples} sparseGRM.sampleIDs.txt
        mkdir -p sites results

        cat <<'PYEOF' > run_contig.py
import argparse
import csv
import gzip
import multiprocessing
import os
import subprocess

p = argparse.ArgumentParser()
p.add_argument("--sites", required=True)
p.add_argument("--contig", required=True)
p.add_argument("--cis-window", type=int, required=True)
p.add_argument("--relatedness-cutoff", required=True)
p.add_argument("--covar-col-list", default="")
p.add_argument("--qcovar-col-list", default="")
p.add_argument("--min-samples", type=int, required=True)
p.add_argument("--vcf-field", default="GT")
p.add_argument("--inv-normalize", choices=["TRUE", "FALSE"], required=True)
p.add_argument("--covariates", default=None)
p.add_argument("--n-workers", type=int, required=True)
p.add_argument("--out", required=True)
p.add_argument("--skipped-log", required=True)
ARGS = p.parse_args()

MISSING = {".", "NA", ""}
SINGLE_THREAD_ENV = dict(
    os.environ,
    OMP_NUM_THREADS="1",
    OPENBLAS_NUM_THREADS="1",
    MKL_NUM_THREADS="1",
    VECLIB_MAXIMUM_THREADS="1",
)

COVARIATES = {}
COVAR_COLS = []
if ARGS.covariates:
    with open(ARGS.covariates) as cf:
        reader = csv.reader(cf, delimiter="\t")
        COVAR_COLS = next(reader)[1:]
        for row in reader:
            COVARIATES[row[0]] = row[1:]


def process_site(row):
    chrom, start, end, site_id, call_rate, sample_ids, values = row
    pos = int(start)

    rows_out = []
    for sid, val in zip(sample_ids, values):
        if val in MISSING:
            continue
        if ARGS.covariates and sid not in COVARIATES:
            continue
        rows_out.append((sid, val, COVARIATES.get(sid, [])))

    if len(rows_out) < ARGS.min_samples:
        return (site_id, "insufficient_samples", str(len(rows_out)), None)

    pheno_path = f"sites/{site_id}.pheno.txt"
    with open(pheno_path, "w") as pf:
        pf.write("\t".join(["person_id", "trait"] + COVAR_COLS) + "\n")
        for sid, val, cov_vals in rows_out:
            pf.write("\t".join([sid, val] + cov_vals) + "\n")

    window_start = max(0, pos - ARGS.cis_window)
    window_end = pos + ARGS.cis_window
    range_path = f"sites/{site_id}.range.bed"
    with open(range_path, "w") as rf:
        rf.write(f"{chrom}\t{window_start}\t{window_end}\n")

    step1_prefix = f"results/{site_id}.step1"
    step1_cmd = [
        "step1_fitNULLGLMM.R",
        "--plinkFile=variance_ratio_markers",
        "--useSparseGRMtoFitNULL=TRUE",
        "--sparseGRMFile=sparseGRM.mtx",
        "--sparseGRMSampleIDFile=sparseGRM.sampleIDs.txt",
        f"--phenoFile={pheno_path}",
        "--phenoCol=trait",
        f"--covarColList={ARGS.covar_col_list}",
        f"--qCovarColList={ARGS.qcovar_col_list}",
        f"--invNormalize={ARGS.inv_normalize}",
        "--sampleIDColinphenoFile=person_id",
        "--traitType=quantitative",
        "--IsOverwriteVarianceRatioFile=TRUE",
        "--isCateVarianceRatio=FALSE",
        f"--outputPrefix={step1_prefix}",
        "--maxiter=5000",
    ]
    r1 = subprocess.run(step1_cmd, capture_output=True, text=True, env=SINGLE_THREAD_ENV)
    if r1.returncode != 0 or not os.path.exists(step1_prefix + ".rda"):
        os.remove(pheno_path)
        os.remove(range_path)
        return (site_id, "step1_failed", r1.stderr[-500:].strip(), None)

    step2_out = f"results/{site_id}.step2.assoc.txt"
    step2_cmd = [
        "step2_SPAtests.R",
        "--vcfFile=genotypes.vcf.gz",
        "--vcfFileIndex=genotypes.vcf.gz.csi",
        f"--vcfField={ARGS.vcf_field}",
        f"--chrom={ARGS.contig}",
        f"--rangestoIncludeFile={range_path}",
        f"--GMMATmodelFile={step1_prefix}.rda",
        f"--varianceRatioFile={step1_prefix}.varianceRatio.txt",
        "--sparseGRMFile=sparseGRM.mtx",
        "--sparseGRMSampleIDFile=sparseGRM.sampleIDs.txt",
        f"--relatednessCutoff={ARGS.relatedness_cutoff}",
        "--LOCO=FALSE",
        f"--SAIGEOutputFile={step2_out}",
        "--is_output_moreDetails=FALSE",
    ]
    r2 = subprocess.run(step2_cmd, capture_output=True, text=True, env=SINGLE_THREAD_ENV)

    result = None
    status = None
    if r2.returncode != 0 or not os.path.exists(step2_out):
        status = ("step2_failed", r2.stderr[-500:].strip())
    else:
        with open(step2_out) as sf:
            out_header = sf.readline().rstrip("\n")
            lines = [f"{site_id}\t{chrom}\t{pos}\t{len(rows_out)}\t" + line for line in sf]
        result = (out_header, lines)

    os.remove(pheno_path)
    os.remove(range_path)
    for suffix in (".rda", ".varianceRatio.txt"):
        fp = step1_prefix + suffix
        if os.path.exists(fp):
            os.remove(fp)
    if os.path.exists(step2_out):
        os.remove(step2_out)

    if result is not None:
        return (site_id, "ok", None, result)
    return (site_id, status[0], status[1], None)


def main():
    rows = []
    with gzip.open(ARGS.sites, "rt") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        sample_ids = header[5:]
        for row in reader:
            chrom, start, end, site_id, call_rate = row[:5]
            rows.append((chrom, start, end, site_id, call_rate, sample_ids, row[5:]))

    n_tested = 0
    n_skipped = 0
    wrote_header = False
    with open(ARGS.out, "w") as combined_out, open(ARGS.skipped_log, "w") as skipped_log:
        with multiprocessing.Pool(ARGS.n_workers) as pool:
            for site_id, status, detail, result in pool.imap(process_site, rows, chunksize=1):
                if status != "ok":
                    n_skipped += 1
                    skipped_log.write(f"{site_id}\t{status}\t{detail}\n")
                    continue
                out_header, lines = result
                if not wrote_header:
                    combined_out.write(
                        "pheno_site_id\tpheno_chrom\tpheno_pos\tn_samples\t" + out_header + "\n"
                    )
                    wrote_header = True
                combined_out.writelines(lines)
                n_tested += 1

        if not wrote_header:
            combined_out.write("pheno_site_id\tpheno_chrom\tpheno_pos\tn_samples\n")

    print(f"tested={n_tested} skipped={n_skipped}")


if __name__ == "__main__":
    main()
PYEOF

        python3 run_contig.py \
            --sites ~{filtered_sites} \
            --contig ~{contig} \
            --cis-window ~{cis_window} \
            --relatedness-cutoff ~{relatedness_cutoff} \
            --covar-col-list "~{covar_col_list}" \
            --qcovar-col-list "~{qcovar_col_list}" \
            --min-samples ~{min_samples_per_site} \
            --vcf-field ~{vcf_field} \
            --inv-normalize ~{if inv_normalize then "TRUE" else "FALSE"} \
            --n-workers ~{n_parallel_workers} \
            --out ~{prefix}.assoc.txt \
            --skipped-log ~{prefix}.skipped_sites.log \
            ~{if defined(covariates_file) then "--covariates=" + covariates_file else ""}
    >>>

    output {
        File contig_assoc = "~{prefix}.assoc.txt"
        File skipped_sites_log = "~{prefix}.skipped_sites.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: n_parallel_workers,
        mem_gb: 4 * n_parallel_workers,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

# Concatenates an array of TSVs that all share the same header, keeping the
# header from the first file only. Used both to gather chunks into a
# per-contig result and to gather per-contig results into the final output.
task ConcatenateTsvs {
    input {
        Array[File] tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ ~{length(tsvs)} -eq 0 ]; then
            touch ~{prefix}.tsv
        else
            awk 'FNR==1 && NR!=1 { next } { print }' ~{sep=" " tsvs} > ~{prefix}.tsv
        fi
    >>>

    output {
        File merged = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 5 * ceil(size(tsvs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
