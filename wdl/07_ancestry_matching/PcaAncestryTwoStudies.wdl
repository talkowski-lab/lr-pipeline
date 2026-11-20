version 1.0

import "../utils/Structs.wdl"
import "PcaAncestryProjection.wdl"

# Project two study cohorts (e.g. CDP WES cases and AoU long-read WGS controls) onto a PCA trained on a labelled reference
# (e.g. HGDP+1kGP) restricted to a shared SNV list, and assign each study sample an ancestry.
# Inputs are the snv_vcf outputs of SelectSnvsInBed (one VCF per call set, already restricted to the target regions and
# PASS/biallelic as chosen there).
# 1. FilterToSites (study1, study2, reference): keeps only SNVs whose CHROM:POS:REF:ALT is in sites_tsv (chr, pos, end,
#    ref, alt; optionally gzipped; header skipped). Local files, so this takes seconds to minutes.
# 2. MergeStudies: INFO and all FORMAT fields but GT are dropped, then the two studies are merged (bcftools merge -m none;
#    sample IDs must be unique across the studies). A site missing from one study is missing (./.) in its samples and is
#    removed by the PCA missingness filter.
# 3. PcaAncestryProjection (the PCA_for_DP102.sh recipe): plink QC on the merged studies and on the reference, SNVs
#    shared by both, LD pruned on the reference, merged, convertf, smartpca with PCs fit on the reference populations
#    only (poplistname) and study samples projected (lsqproject, nsnpldregress), then kNN ancestry labels.
# 4. AddCohort: adds each study sample's cohort (study1_name / study2_name) to the ancestry table.

workflow PcaAncestryTwoStudies {
    input {
        File study1_vcf
        File study1_vcf_idx
        File study2_vcf
        File study2_vcf_idx
        File ref_vcf
        File ref_vcf_idx
        File sites_tsv
        String prefix

        String study1_name = "CDP"
        String study2_name = "AoU"
        File ref_ancestry_tsv
        String ref_label_column = "genetic_region"
        File? ref_samples
        Float study_min_maf = 0.01
        Float ref_min_maf = 0.01
        Float max_missing = 0.05
        Boolean ld_prune = true
        Int n_pcs = 10
        Int nsnpldregress = 2
        Int knn_k = 10
        Float min_prob = 0.75

        String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
        String plink2_docker = "quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0"
        String plink_docker = "quay.io/biocontainers/plink:1.90b7.7--h18e278d_1"
        String eigensoft_docker = "quay.io/biocontainers/eigensoft:8.0.0--h75d7a4a_6"
        String python_docker = "python:3.11-slim"

        RuntimeAttr? runtime_attr_filter_to_sites
        RuntimeAttr? runtime_attr_merge_studies
        RuntimeAttr? runtime_attr_prep_plink
        RuntimeAttr? runtime_attr_merge_for_pca
        RuntimeAttr? runtime_attr_smartpca
        RuntimeAttr? runtime_attr_assign_ancestry
        RuntimeAttr? runtime_attr_add_cohort
    }

    call FilterToSites as FilterStudy1 {
        input:
            vcf = study1_vcf,
            vcf_idx = study1_vcf_idx,
            sites_tsv = sites_tsv,
            prefix = "~{prefix}.~{study1_name}",
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_filter_to_sites
    }

    call FilterToSites as FilterStudy2 {
        input:
            vcf = study2_vcf,
            vcf_idx = study2_vcf_idx,
            sites_tsv = sites_tsv,
            prefix = "~{prefix}.~{study2_name}",
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_filter_to_sites
    }

    call FilterToSites as FilterRef {
        input:
            vcf = ref_vcf,
            vcf_idx = ref_vcf_idx,
            sites_tsv = sites_tsv,
            prefix = "~{prefix}.ref",
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_filter_to_sites
    }

    call MergeStudies {
        input:
            study1_vcf = FilterStudy1.sites_vcf,
            study1_vcf_idx = FilterStudy1.sites_vcf_idx,
            study2_vcf = FilterStudy2.sites_vcf,
            study2_vcf_idx = FilterStudy2.sites_vcf_idx,
            study1_name = study1_name,
            study2_name = study2_name,
            prefix = prefix,
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_merge_studies
    }

    call PcaAncestryProjection.PcaAncestryProjection as Pca {
        input:
            case_vcf = MergeStudies.merged_vcf,
            case_vcf_idx = MergeStudies.merged_vcf_idx,
            ref_vcf = FilterRef.sites_vcf,
            ref_vcf_idx = FilterRef.sites_vcf_idx,
            prefix = prefix,
            ref_ancestry_tsv = ref_ancestry_tsv,
            ref_label_column = ref_label_column,
            ref_samples = ref_samples,
            variants_tsv = sites_tsv,
            case_min_maf = study_min_maf,
            ref_min_maf = ref_min_maf,
            max_missing = max_missing,
            ld_prune = ld_prune,
            n_pcs = n_pcs,
            nsnpldregress = nsnpldregress,
            knn_k = knn_k,
            min_prob = min_prob,
            plink2_docker = plink2_docker,
            plink_docker = plink_docker,
            eigensoft_docker = eigensoft_docker,
            python_docker = python_docker,
            runtime_attr_prep_plink = runtime_attr_prep_plink,
            runtime_attr_merge_for_pca = runtime_attr_merge_for_pca,
            runtime_attr_smartpca = runtime_attr_smartpca,
            runtime_attr_assign_ancestry = runtime_attr_assign_ancestry
    }

    call AddCohort {
        input:
            case_ancestry = Pca.case_ancestry,
            sample_cohorts = MergeStudies.sample_cohorts,
            prefix = prefix,
            docker = python_docker,
            runtime_attr_override = runtime_attr_add_cohort
    }

    output {
        File study_ancestry = AddCohort.study_ancestry
        File pcs_tsv = Pca.pcs_tsv
        File eigenvalues = Pca.eigenvalues
        File snp_weights = Pca.snp_weights
        File smartpca_log = Pca.smartpca_log
        File pca_variant_ids = Pca.pca_variant_ids
        File variant_counts = Pca.variant_counts
        File study_prep_log = Pca.case_prep_log
        File ref_prep_log = Pca.ref_prep_log
        File study1_site_count = FilterStudy1.site_count
        File study2_site_count = FilterStudy2.site_count
        File ref_site_count = FilterRef.site_count
        File merged_study_vcf = MergeStudies.merged_vcf
        File merged_study_vcf_idx = MergeStudies.merged_vcf_idx
        File ref_sites_vcf = FilterRef.sites_vcf
        File ref_sites_vcf_idx = FilterRef.sites_vcf_idx
    }
}

task FilterToSites {
    input {
        File vcf
        File vcf_idx
        File sites_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        gzip -cdf ~{sites_tsv} \
            | awk -F'\t' '$2 ~ /^[0-9]+$/ {print $1 "\t" $2 > "sites.pos"; print $1 ":" $2 ":" $4 ":" $5 > "sites.ids"}'
        sort -k1,1 -k2,2n -u sites.pos > sites.sorted.pos

        # -T streams the file and keeps listed positions; the ID match then requires the same REF and ALT
        bcftools view -T sites.sorted.pos ~{vcf} -Ou \
            | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Ou \
            | bcftools view -i 'ID=@sites.ids' -Oz -o ~{prefix}.sites.vcf.gz
        bcftools index -t ~{prefix}.sites.vcf.gz

        n_in=$(bcftools view -H ~{vcf} | wc -l)
        n_out=$(bcftools index -n ~{prefix}.sites.vcf.gz)
        printf 'listed_sites\tinput_snvs\tkept_snvs\n%s\t%s\t%s\n' "$(wc -l < sites.ids)" "$n_in" "$n_out" \
            | tee ~{prefix}.site_count.tsv
    >>>

    output {
        File sites_vcf = "~{prefix}.sites.vcf.gz"
        File sites_vcf_idx = "~{prefix}.sites.vcf.gz.tbi"
        File site_count = "~{prefix}.site_count.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
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

task MergeStudies {
    input {
        File study1_vcf
        File study1_vcf_idx
        File study2_vcf
        File study2_vcf_idx
        String study1_name
        String study2_name
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate -x INFO,^FORMAT/GT ~{study1_vcf} -Oz -o study1.gt.vcf.gz
        bcftools annotate -x INFO,^FORMAT/GT ~{study2_vcf} -Oz -o study2.gt.vcf.gz
        bcftools index -t study1.gt.vcf.gz
        bcftools index -t study2.gt.vcf.gz
        bcftools merge -m none study1.gt.vcf.gz study2.gt.vcf.gz -Oz -o ~{prefix}.studies_merged.vcf.gz
        bcftools index -t ~{prefix}.studies_merged.vcf.gz

        {
            echo -e "sample_id\tcohort"
            bcftools query -l study1.gt.vcf.gz | awk -v c="~{study1_name}" '{print $1"\t"c}'
            bcftools query -l study2.gt.vcf.gz | awk -v c="~{study2_name}" '{print $1"\t"c}'
        } > ~{prefix}.sample_cohorts.tsv
        echo "~{study1_name}: $(bcftools index -n study1.gt.vcf.gz) SNVs, $(bcftools query -l study1.gt.vcf.gz | wc -l) samples"
        echo "~{study2_name}: $(bcftools index -n study2.gt.vcf.gz) SNVs, $(bcftools query -l study2.gt.vcf.gz | wc -l) samples"
        echo "merged: $(bcftools index -n ~{prefix}.studies_merged.vcf.gz) SNVs"
    >>>

    output {
        File merged_vcf = "~{prefix}.studies_merged.vcf.gz"
        File merged_vcf_idx = "~{prefix}.studies_merged.vcf.gz.tbi"
        File sample_cohorts = "~{prefix}.sample_cohorts.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size(study1_vcf, "GB") + size(study2_vcf, "GB")) + 10,
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

task AddCohort {
    input {
        File case_ancestry
        File sample_cohorts
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
        import collections
        import csv

        cohort = {r["sample_id"]: r["cohort"] for r in csv.DictReader(open("~{sample_cohorts}"), delimiter="\t")}
        rows = list(csv.DictReader(open("~{case_ancestry}"), delimiter="\t"))
        cols = ["sample_id", "cohort"] + [c for c in rows[0] if c != "sample_id"]
        with open("~{prefix}.study_ancestry.tsv", "w") as o:
            w = csv.DictWriter(o, fieldnames=cols, delimiter="\t", lineterminator="\n")
            w.writeheader()
            for r in rows:
                r["cohort"] = cohort.get(r["sample_id"], "NA")
                w.writerow(r)
        counts = collections.Counter((r["cohort"], r["assigned_ancestry"]) for r in rows)
        for (c, a), n in sorted(counts.items()):
            print(f"{c}\t{a}\t{n}")
        CODE
    >>>

    output {
        File study_ancestry = "~{prefix}.study_ancestry.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
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
