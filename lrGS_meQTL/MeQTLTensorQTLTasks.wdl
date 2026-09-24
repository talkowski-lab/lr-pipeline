version 1.0

import "../wdl/utils/Structs.wdl"

# Splits multiallelic records into biallelic ones (bcftools norm -m -any).
# tensorQTL's genotype reader, like most cis-QTL tools, expects biallelic
# variants; this also sidesteps plink2's ~254-ALT-allele import limit on
# the rare complex multiallelic site, without silently dropping it the way
# SAIGE's VCF reader does (see MeQTLTasks.wdl).
task NormalizeVcf {
    input {
        File vcf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools norm -m -any ~{vcf} -Oz -o ~{prefix}.norm.vcf.gz
    >>>

    output {
        File normalized_vcf = "~{prefix}.norm.vcf.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 20,
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

# Converts a (biallelic) VCF to plink2's pgen/pvar/psam format for tensorQTL.
# --output-chr chrM keeps the "chr"-prefixed contig naming (plink2's default
# strips it to bare "22"), which must match the phenotype bed's #chr column
# for tensorQTL to correctly identify cis variants.
task ConvertVcfToPgen {
    input {
        File vcf
        String vcf_half_call
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        plink2 --vcf ~{vcf} --vcf-half-call ~{vcf_half_call} --output-chr chrM \
            --make-pgen --out ~{prefix}
    >>>

    output {
        File pgen = "~{prefix}.pgen"
        File pvar = "~{prefix}.pvar"
        File psam = "~{prefix}.psam"
        File log = "~{prefix}.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 20,
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
# rate meets the threshold, mean-imputes any remaining missing values (a
# rectangular, complete matrix is required by tensorQTL, unlike SAIGE's
# per-site sample subsetting), and reshapes it into tensorQTL's expected
# phenotype-bed layout (#chr, start, end, phenotype_id, then sample columns).
task BuildPhenotypeBed {
    input {
        File methylation_bed
        Float call_rate_threshold
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat <<'PYEOF' > build_phenotype_bed.py
import csv
import gzip
import sys

in_path, threshold, out_path, count_path = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4]
MISSING = {".", "NA", ""}

n_qualifying = 0
with gzip.open(in_path, "rt") as fin, open(out_path, "w") as fout:
    reader = csv.reader(fin, delimiter="\t")
    header = next(reader)
    header[0] = header[0].lstrip("#")
    sample_ids = header[3:]
    writer = csv.writer(fout, delimiter="\t")
    writer.writerow(["#chr", "start", "end", "phenotype_id"] + sample_ids)

    for row in reader:
        chrom, start, end = row[0], row[1], row[2]
        values = row[3:]
        total = len(values)
        if total == 0:
            continue
        missing_idx = {i for i, v in enumerate(values) if v in MISSING}
        call_rate = (total - len(missing_idx)) / total
        if call_rate < threshold:
            continue
        present = [float(v) for i, v in enumerate(values) if i not in missing_idx]
        mean_val = sum(present) / len(present) if present else 0.0
        filled = [mean_val if i in missing_idx else float(v) for i, v in enumerate(values)]
        site_id = f"{chrom}_{start}_{end}"
        writer.writerow([chrom, start, end, site_id] + [f"{v:.6f}" for v in filled])
        n_qualifying += 1

with open(count_path, "w") as f:
    f.write(str(n_qualifying) + "\n")
PYEOF
        python3 build_phenotype_bed.py ~{methylation_bed} ~{call_rate_threshold} ~{prefix}.phenotype.bed ~{prefix}.n_sites.txt
    >>>

    output {
        File phenotype_bed_plain = "~{prefix}.phenotype.bed"
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

task BgzipTabixBed {
    input {
        File bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bgzip -c ~{bed} > ~{prefix}.bed.gz
        tabix -p bed ~{prefix}.bed.gz
    >>>

    output {
        File bed_gz = "~{prefix}.bed.gz"
        File bed_gz_index = "~{prefix}.bed.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(bed, "GB")) + 10,
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

# Builds tensorQTL's covariates file (first column covariate name, remaining
# columns one per sample - transposed at load time by tensorQTL itself) from
# our own long-format covariates_file (person_id + one column per covariate),
# reordering/subsetting samples to match the plink2 .psam sample order. If no
# covariates_file is supplied, emits a header-only (zero-covariate) file so
# the pipeline still runs with an intercept-only model.
task BuildCovariates {
    input {
        File psam
        File? covariates_file
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat <<'PYEOF' > build_covariates.py
import argparse
import csv

p = argparse.ArgumentParser()
p.add_argument("--psam", required=True)
p.add_argument("--out", required=True)
p.add_argument("--covariates", default=None)
args = p.parse_args()

with open(args.psam) as f:
    f.readline()
    sample_ids = [line.rstrip("\n").split("\t")[0] for line in f]

covar_rows = []
if args.covariates:
    with open(args.covariates) as cf:
        reader = csv.reader(cf, delimiter="\t")
        covar_names = next(reader)[1:]
        by_sample = {row[0]: row[1:] for row in reader}
    for i, name in enumerate(covar_names):
        covar_rows.append(
            [name] + [by_sample.get(sid, ["NA"] * len(covar_names))[i] for sid in sample_ids]
        )

with open(args.out, "w") as out:
    out.write("\t".join(["ID"] + sample_ids) + "\n")
    for row in covar_rows:
        out.write("\t".join(row) + "\n")
PYEOF
        python3 build_covariates.py \
            --psam ~{psam} \
            --out ~{prefix}.covariates.txt \
            ~{if defined(covariates_file) then "--covariates=" + covariates_file else ""}
    >>>

    output {
        File covariates = "~{prefix}.covariates.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(psam, "GB")) + 10,
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

# Modernized (WDL 1.0, this repo's RuntimeAttr/prefix conventions) port of
# AoU-Multiomics-Analysis/tensorQTL_cis_permutations's tensorqtl_cis_permutations
# task. Unlike SAIGE, tensorQTL processes every qualifying phenotype on the
# contig in a single vectorized/GPU call - no per-site chunking needed.
task TensorQTLCisPermutations {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        File phenotype_bed
        File phenotype_bed_index
        File covariates
        Int cis_window
        File? phenotype_groups
        Float? fdr
        Float? qvalue_lambda
        Float? pval_threshold
        Int? seed
        String? flags
        Int num_gpus
        String gpu_type
        Array[String] gpu_zones
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        ln -s ~{plink_pgen} ~{prefix}.pgen
        ln -s ~{plink_pvar} ~{prefix}.pvar
        ln -s ~{plink_psam} ~{prefix}.psam
        ln -s ~{phenotype_bed} ~{prefix}.phenotype.bed.gz
        ln -s ~{phenotype_bed_index} ~{prefix}.phenotype.bed.gz.tbi

        python3 -m tensorqtl \
            ~{prefix} ~{prefix}.phenotype.bed.gz ~{prefix} \
            --mode cis \
            --covariates ~{covariates} \
            --window ~{cis_window} \
            ~{"--phenotype_groups " + phenotype_groups} \
            ~{"--fdr " + fdr} \
            ~{"--pval_threshold " + pval_threshold} \
            ~{"--qvalue_lambda " + qvalue_lambda} \
            ~{"--seed " + seed} \
            ~{flags}
    >>>

    output {
        File cis_qtl = "~{prefix}.cis_qtl.txt.gz"
        File log = "~{prefix}.tensorQTL.cis.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 32,
        disk_gb: 3 * ceil(size(plink_pgen, "GB") + size(phenotype_bed, "GB")) + 20,
        boot_disk_gb: 25,
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
        gpuType: gpu_type
        gpuCount: num_gpus
        zones: gpu_zones
    }
}

# Concatenates an array of gzipped TSVs that all share the same header,
# keeping the header from the first file only. Used to gather per-contig
# cis_qtl.txt.gz outputs into one genome-wide file.
task ConcatenateGzippedTsvs {
    input {
        Array[File] gz_tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ ~{length(gz_tsvs)} -eq 0 ]; then
            touch ~{prefix}.tsv.gz
        else
            mkdir -p decompressed
            i=0
            for f in ~{sep=" " gz_tsvs}; do
                gunzip -c "$f" > decompressed/"$i".tsv
                i=$((i + 1))
            done
            awk 'FNR==1 && NR!=1 { next } { print }' decompressed/*.tsv | bgzip -c > ~{prefix}.tsv.gz
        fi
    >>>

    output {
        File merged = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 5 * ceil(size(gz_tsvs, "GB")) + 10,
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
