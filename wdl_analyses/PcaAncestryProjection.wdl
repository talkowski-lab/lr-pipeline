version 1.0

import "Structs.wdl"

# Infer the ancestry of study (case) samples by projecting them onto a PCA of labelled reference samples (e.g. HGDP+1kGP).
#
# 1. PrepPlink (case and reference): VCF -> PLINK, autosomal biallelic SNVs only, variant IDs set to CHROM:POS:REF:ALT
#    (CHROM without "chr"), genotype missingness <= max_missing and MAF >= min_maf. Reference samples can be restricted
#    with ref_samples (one sample ID per line, e.g. unrelated samples). With variants_tsv (columns chr, pos, end, ref, alt;
#    pos is the 1-based VCF position; optionally gzipped; a header line is skipped), only the listed SNVs are used in both
#    datasets.
# 2. MergeForPca: SNVs passing QC in both datasets, LD pruned on the reference (plink --indep-pairwise) when ld_prune is
#    true, then the two datasets are merged.
# 3. SmartPca: convertf to EIGENSTRAT, reference samples labelled from ref_ancestry_tsv (column ref_label_column, sample IDs
#    in the first column), study samples labelled STUDY. PCs are fit on the reference populations (poplistname) and
#    study samples are projected (lsqproject), with nsnpldregress neighbouring SNVs regressed out for LD (0 disables).
# 4. AssignAncestry: each study sample gets the majority label of its knn_k nearest reference samples in PC space; the
#    label is kept when that fraction is >= min_prob, otherwise "unassigned".

workflow PcaAncestryProjection {
    input {
        File case_vcf
        File case_vcf_idx
        File ref_vcf
        File ref_vcf_idx
        String prefix

        File ref_ancestry_tsv
        String ref_label_column = "genetic_region"
        File? ref_samples
        File? variants_tsv
        Float case_min_maf = 0.01
        Float ref_min_maf = 0.01
        Float max_missing = 0.05
        Boolean ld_prune = true
        Int ld_window = 200
        Int ld_step = 50
        Float ld_r2 = 0.2
        Int n_pcs = 10
        Int nsnpldregress = 2
        Int knn_k = 10
        Float min_prob = 0.75

        String plink2_docker = "quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0"
        String plink_docker = "quay.io/biocontainers/plink:1.90b7.7--h18e278d_1"
        String eigensoft_docker = "quay.io/biocontainers/eigensoft:8.0.0--h75d7a4a_6"
        String python_docker = "python:3.11-slim"

        RuntimeAttr? runtime_attr_prep_plink
        RuntimeAttr? runtime_attr_merge_for_pca
        RuntimeAttr? runtime_attr_smartpca
        RuntimeAttr? runtime_attr_assign_ancestry
    }

    call PrepPlink as PrepCase {
        input:
            vcf = case_vcf,
            vcf_idx = case_vcf_idx,
            variants_tsv = variants_tsv,
            min_maf = case_min_maf,
            max_missing = max_missing,
            prefix = "~{prefix}.case",
            docker = plink2_docker,
            runtime_attr_override = runtime_attr_prep_plink
    }

    call PrepPlink as PrepRef {
        input:
            vcf = ref_vcf,
            vcf_idx = ref_vcf_idx,
            keep_samples = ref_samples,
            variants_tsv = variants_tsv,
            min_maf = ref_min_maf,
            max_missing = max_missing,
            prefix = "~{prefix}.ref",
            docker = plink2_docker,
            runtime_attr_override = runtime_attr_prep_plink
    }

    call MergeForPca {
        input:
            case_bed = PrepCase.bed,
            case_bim = PrepCase.bim,
            case_fam = PrepCase.fam,
            ref_bed = PrepRef.bed,
            ref_bim = PrepRef.bim,
            ref_fam = PrepRef.fam,
            ld_prune = ld_prune,
            ld_window = ld_window,
            ld_step = ld_step,
            ld_r2 = ld_r2,
            prefix = prefix,
            docker = plink_docker,
            runtime_attr_override = runtime_attr_merge_for_pca
    }

    call SmartPca {
        input:
            bed = MergeForPca.bed,
            bim = MergeForPca.bim,
            fam = MergeForPca.fam,
            case_fam = PrepCase.fam,
            ref_ancestry_tsv = ref_ancestry_tsv,
            ref_label_column = ref_label_column,
            n_pcs = n_pcs,
            nsnpldregress = nsnpldregress,
            prefix = prefix,
            docker = eigensoft_docker,
            runtime_attr_override = runtime_attr_smartpca
    }

    call AssignAncestry {
        input:
            pcs_tsv = SmartPca.pcs_tsv,
            knn_k = knn_k,
            min_prob = min_prob,
            prefix = prefix,
            docker = python_docker,
            runtime_attr_override = runtime_attr_assign_ancestry
    }

    output {
        File case_ancestry = AssignAncestry.case_ancestry
        File pcs_tsv = SmartPca.pcs_tsv
        File eigenvalues = SmartPca.eval
        File snp_weights = SmartPca.snp_weights
        File smartpca_log = SmartPca.log
        File pca_variant_ids = MergeForPca.pca_variant_ids
        File variant_counts = MergeForPca.variant_counts
        File case_prep_log = PrepCase.log
        File ref_prep_log = PrepRef.log
    }
}

task PrepPlink {
    input {
        File vcf
        File vcf_idx
        File? keep_samples
        File? variants_tsv
        Float min_maf
        Float max_missing
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ -n "~{variants_tsv}" ]; then
            gzip -cdf ~{variants_tsv} | awk '$2 ~ /^[0-9]+$/ {sub(/^chr/, "", $1); print $1":"$2":"$4":"$5}' > variant_ids.txt
            echo "variant list: $(wc -l < variant_ids.txt) SNVs"
        fi

        plink2 \
            --vcf ~{vcf} \
            --vcf-half-call m \
            --autosome \
            --snps-only just-acgt \
            --max-alleles 2 \
            --set-all-var-ids '@:#:$r:$a' \
            --new-id-max-allele-len 1000 \
            ~{"--keep " + keep_samples} \
            ~{if defined(variants_tsv) then "--extract variant_ids.txt" else ""} \
            --geno ~{max_missing} \
            --maf ~{min_maf} \
            --make-bed \
            --out ~{prefix}

        echo "~{prefix}: $(wc -l < ~{prefix}.bim) SNVs, $(wc -l < ~{prefix}.fam) samples"
    >>>

    output {
        File bed = "~{prefix}.bed"
        File bim = "~{prefix}.bim"
        File fam = "~{prefix}.fam"
        File log = "~{prefix}.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
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

task MergeForPca {
    input {
        File case_bed
        File case_bim
        File case_fam
        File ref_bed
        File ref_bim
        File ref_fam
        Boolean ld_prune
        Int ld_window
        Int ld_step
        Float ld_r2
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cut -f2 ~{ref_bim} | sort > ref.ids
        cut -f2 ~{case_bim} | sort > case.ids
        comm -12 ref.ids case.ids > shared.ids

        if [ "~{ld_prune}" == "true" ]; then
            plink \
                --bed ~{ref_bed} --bim ~{ref_bim} --fam ~{ref_fam} \
                --extract shared.ids \
                --indep-pairwise ~{ld_window} ~{ld_step} ~{ld_r2} \
                --out ld
            sort ld.prune.in > ~{prefix}.pca_variant_ids.txt
        else
            cp shared.ids ~{prefix}.pca_variant_ids.txt
        fi

        plink --bed ~{ref_bed} --bim ~{ref_bim} --fam ~{ref_fam} --extract ~{prefix}.pca_variant_ids.txt --make-bed --out ref_pca
        plink --bed ~{case_bed} --bim ~{case_bim} --fam ~{case_fam} --extract ~{prefix}.pca_variant_ids.txt --make-bed --out case_pca
        plink --bfile ref_pca --bmerge case_pca --make-bed --out ~{prefix}.merged

        {
            echo -e "step\tn_snvs"
            echo -e "ref_qc\t$(wc -l < ref.ids)"
            echo -e "case_qc\t$(wc -l < case.ids)"
            echo -e "shared\t$(wc -l < shared.ids)"
            echo -e "pca\t$(wc -l < ~{prefix}.pca_variant_ids.txt)"
        } > ~{prefix}.variant_counts.tsv
        cat ~{prefix}.variant_counts.tsv
    >>>

    output {
        File bed = "~{prefix}.merged.bed"
        File bim = "~{prefix}.merged.bim"
        File fam = "~{prefix}.merged.fam"
        File pca_variant_ids = "~{prefix}.pca_variant_ids.txt"
        File variant_counts = "~{prefix}.variant_counts.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(case_bed, "GB") + size(ref_bed, "GB")) + 10,
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

task SmartPca {
    input {
        File bed
        File bim
        File fam
        File case_fam
        File ref_ancestry_tsv
        String ref_label_column
        Int n_pcs
        Int nsnpldregress
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat > convertf.par <<EOF
        genotypename: ~{bed}
        snpname: ~{bim}
        indivname: ~{fam}
        outputformat: EIGENSTRAT
        genotypeoutname: merged.geno
        snpoutname: merged.snp
        indivoutname: merged.raw.ind
        familynames: NO
        EOF
        convertf -p convertf.par

        awk -F'\t' -v col="~{ref_label_column}" \
            'NR == 1 {for (i = 1; i <= NF; i++) if ($i == col) c = i; if (!c) exit 1; next} {print $1"\t"$c}' \
            ~{ref_ancestry_tsv} > labels.tsv
        awk '{print $2"\tSTUDY"}' ~{case_fam} >> labels.tsv
        awk 'FNR == NR {lab[$1] = $2; next} {l = ($1 in lab) ? lab[$1] : "Ignore"; print $1, $2, l}' \
            labels.tsv merged.raw.ind > merged.ind
        awk '$3 != "STUDY" && $3 != "Ignore" {print $3}' merged.ind | sort -u > poplist.txt
        echo "reference populations: $(tr '\n' ' ' < poplist.txt)"
        awk '{print $3}' merged.ind | sort | uniq -c

        cat > smartpca.par <<EOF
        genotypename: merged.geno
        snpname: merged.snp
        indivname: merged.ind
        evecoutname: ~{prefix}.evec
        evaloutname: ~{prefix}.eval
        snpweightoutname: ~{prefix}.snp_weights.txt
        poplistname: poplist.txt
        lsqproject: YES
        numoutevec: ~{n_pcs}
        nsnpldregress: ~{nsnpldregress}
        numoutlieriter: 0
        altnormstyle: NO
        numthreads: 4
        EOF
        smartpca -p smartpca.par > ~{prefix}.smartpca.log

        awk -v n=~{n_pcs} 'BEGIN {printf "sample_id"; for (i = 1; i <= n; i++) printf "\tPC%d", i; print "\tlabel"}
            NR > 1 {printf "%s", $1; for (i = 2; i <= NF; i++) printf "\t%s", $i; print ""}' \
            ~{prefix}.evec > ~{prefix}.pcs.tsv
    >>>

    output {
        File pcs_tsv = "~{prefix}.pcs.tsv"
        File eval = "~{prefix}.eval"
        File snp_weights = "~{prefix}.snp_weights.txt"
        File log = "~{prefix}.smartpca.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: 5 * ceil(size(bed, "GB")) + 10,
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

task AssignAncestry {
    input {
        File pcs_tsv
        Int knn_k
        Float min_prob
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
        import collections
        import csv

        with open("~{pcs_tsv}") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        pc_cols = [c for c in rows[0] if c.startswith("PC")]
        ref = [(r["label"], [float(r[c]) for c in pc_cols]) for r in rows if r["label"] not in ("STUDY", "Ignore")]
        study = [r for r in rows if r["label"] == "STUDY"]

        with open("~{prefix}.case_ancestry.tsv", "w") as o:
            o.write("sample_id\tassigned_ancestry\tknn_label\tknn_fraction\t" + "\t".join(pc_cols) + "\n")
            for r in study:
                x = [float(r[c]) for c in pc_cols]
                nearest = sorted(ref, key=lambda p: sum((a - b) ** 2 for a, b in zip(p[1], x)))[:~{knn_k}]
                label, n = collections.Counter(p[0] for p in nearest).most_common(1)[0]
                frac = n / len(nearest)
                assigned = label if frac >= ~{min_prob} else "unassigned"
                o.write(f"{r['sample_id']}\t{assigned}\t{label}\t{frac:.2f}\t" + "\t".join(r[c] for c in pc_cols) + "\n")

        counts = collections.Counter(line.split("\t")[1] for line in open("~{prefix}.case_ancestry.tsv").readlines()[1:])
        print(f"{len(ref)} reference and {len(study)} study samples; assigned: {dict(counts)}")
        CODE
    >>>

    output {
        File case_ancestry = "~{prefix}.case_ancestry.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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
