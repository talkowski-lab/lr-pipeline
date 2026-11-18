version 1.0

import "Structs.wdl"
import "SelectSnvsInBed.wdl"
import "PcaAncestryTwoStudies.wdl"

# Pick ancestry-matched controls for a case cohort and annotate the case VCF with control allele counts, from raw call
# sets in one run (e.g. CDP WES cases, AoU long-read WGS controls, HGDP+1kGP reference).
# 1. SelectSnvsInBed, scattered over the three call sets (reference, case, control; each an array of gs:// VCF/BCF files
#    read in place with remote tabix): biallelic autosomal SNVs inside regions_bed with AF > min_af (AF recomputed from
#    the genotypes of that call set), PASS only where set.
# 2. IntersectSites: SNVs (same CHROM:POS:REF:ALT) common in all three call sets -> shared sites_tsv.
# 3-4. PcaAncestryTwoStudies: each step-1 VCF filtered to the shared sites (FilterToSites), case + control merged (GT
#    only), PCA fit on the reference and case/control projected (PcaAncestryProjection), kNN ancestry label per sample
#    with the fraction of neighbours from each reference population (frac_<label>), and a cohort column.
# 5. MatchControls: for each case, the n_controls control samples nearest in PC space (first match_n_pcs PCs, each
#    weighted by sqrt(eigenvalue) so distances follow the variance each PC explains). Matching is with replacement: a
#    control can be matched to several cases. Outputs case -> comma-separated controls, the long table with distances,
#    and how many times each control was matched.
# 6. Case annotation: every case variant (multiallelic records split into biallelic ones, left-aligned when ref_fa is
#    given) is looked up in the control call set by CHROM:POS:REF:ALT (control records split the same way, read with
#    remote tabix at the case positions). For a match, INFO gets AF, AC, AN, N_HET, N_HOMREF and N_HOMALT computed from
#    the control genotypes for three sets: all control samples (<control_name>_*), the matched controls counted once
#    (<control_name>_CTRL_*), and the matched controls counted once per case they were matched to (<control_name>_CTRL_DUP_*).
#    Case variants with no control record are left unannotated.

workflow MatchControlsByAncestry {
    input {
        Array[String] ref_vcfs
        Array[String] ref_vcf_idxs
        Array[String] case_vcfs
        Array[String] case_vcf_idxs
        Array[String] control_vcfs
        Array[String] control_vcf_idxs
        File regions_bed
        String prefix

        String case_name = "CDP"
        String control_name = "AoU"
        Float min_af = 0.01
        Boolean ref_pass_only = false
        Boolean case_pass_only = true
        Boolean control_pass_only = false
        Boolean ref_gvcf = false
        Boolean case_gvcf = false
        Boolean control_gvcf = false
        String? requester_pays_project
        Int regions_per_batch = 500
        Int n_parallel = 16
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
        Int n_controls = 5
        Int match_n_pcs = 10
        Int site_merge_gap = 1000

        File? ref_fa
        File? ref_fai

        String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
        String plink2_docker = "quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0"
        String plink_docker = "quay.io/biocontainers/plink:1.90b7.7--h18e278d_1"
        String eigensoft_docker = "quay.io/biocontainers/eigensoft:8.0.0--h75d7a4a_6"
        String python_docker = "python:3.11-slim"

        RuntimeAttr? runtime_attr_select_snvs_in_bed
        RuntimeAttr? runtime_attr_concat_snvs
        RuntimeAttr? runtime_attr_intersect_sites
        RuntimeAttr? runtime_attr_filter_to_sites
        RuntimeAttr? runtime_attr_merge_studies
        RuntimeAttr? runtime_attr_prep_plink
        RuntimeAttr? runtime_attr_merge_for_pca
        RuntimeAttr? runtime_attr_smartpca
        RuntimeAttr? runtime_attr_assign_ancestry
        RuntimeAttr? runtime_attr_add_cohort
        RuntimeAttr? runtime_attr_match_controls
        RuntimeAttr? runtime_attr_case_sites
        RuntimeAttr? runtime_attr_control_stats
        RuntimeAttr? runtime_attr_concat_stats
        RuntimeAttr? runtime_attr_annotate_case
        RuntimeAttr? runtime_attr_concat_case
    }

    Array[Array[String]] set_vcfs = [ref_vcfs, case_vcfs, control_vcfs]
    Array[Array[String]] set_vcf_idxs = [ref_vcf_idxs, case_vcf_idxs, control_vcf_idxs]
    Array[Boolean] set_pass_only = [ref_pass_only, case_pass_only, control_pass_only]
    Array[Boolean] set_gvcf = [ref_gvcf, case_gvcf, control_gvcf]
    Array[String] set_names = ["ref", case_name, control_name]

    scatter (i in range(3)) {
        call SelectSnvsInBed.SelectSnvsInBed as SelectCommonSnvs {
            input:
                vcfs = set_vcfs[i],
                vcf_idxs = set_vcf_idxs[i],
                regions_bed = regions_bed,
                prefix = "~{prefix}.~{set_names[i]}",
                pass_only = set_pass_only[i],
                gvcf = set_gvcf[i],
                min_af = min_af,
                requester_pays_project = requester_pays_project,
                regions_per_batch = regions_per_batch,
                n_parallel = n_parallel,
                bcftools_docker = bcftools_docker,
                runtime_attr_select_snvs_in_bed = runtime_attr_select_snvs_in_bed,
                runtime_attr_concat_snvs = runtime_attr_concat_snvs
        }
    }

    call IntersectSites {
        input:
            ref_sites_tsv = SelectCommonSnvs.snv_sites_tsv[0],
            case_sites_tsv = SelectCommonSnvs.snv_sites_tsv[1],
            control_sites_tsv = SelectCommonSnvs.snv_sites_tsv[2],
            case_name = case_name,
            control_name = control_name,
            prefix = prefix,
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_intersect_sites
    }

    call PcaAncestryTwoStudies.PcaAncestryTwoStudies as Pca {
        input:
            study1_vcf = SelectCommonSnvs.snv_vcf[1],
            study1_vcf_idx = SelectCommonSnvs.snv_vcf_idx[1],
            study2_vcf = SelectCommonSnvs.snv_vcf[2],
            study2_vcf_idx = SelectCommonSnvs.snv_vcf_idx[2],
            ref_vcf = SelectCommonSnvs.snv_vcf[0],
            ref_vcf_idx = SelectCommonSnvs.snv_vcf_idx[0],
            sites_tsv = IntersectSites.shared_sites_tsv,
            prefix = prefix,
            study1_name = case_name,
            study2_name = control_name,
            ref_ancestry_tsv = ref_ancestry_tsv,
            ref_label_column = ref_label_column,
            ref_samples = ref_samples,
            study_min_maf = study_min_maf,
            ref_min_maf = ref_min_maf,
            max_missing = max_missing,
            ld_prune = ld_prune,
            n_pcs = n_pcs,
            nsnpldregress = nsnpldregress,
            knn_k = knn_k,
            min_prob = min_prob,
            bcftools_docker = bcftools_docker,
            plink2_docker = plink2_docker,
            plink_docker = plink_docker,
            eigensoft_docker = eigensoft_docker,
            python_docker = python_docker,
            runtime_attr_filter_to_sites = runtime_attr_filter_to_sites,
            runtime_attr_merge_studies = runtime_attr_merge_studies,
            runtime_attr_prep_plink = runtime_attr_prep_plink,
            runtime_attr_merge_for_pca = runtime_attr_merge_for_pca,
            runtime_attr_smartpca = runtime_attr_smartpca,
            runtime_attr_assign_ancestry = runtime_attr_assign_ancestry,
            runtime_attr_add_cohort = runtime_attr_add_cohort
    }

    call MatchControls {
        input:
            study_ancestry = Pca.study_ancestry,
            eigenvalues = Pca.eigenvalues,
            case_name = case_name,
            control_name = control_name,
            n_controls = n_controls,
            match_n_pcs = match_n_pcs,
            prefix = prefix,
            docker = python_docker,
            runtime_attr_override = runtime_attr_match_controls
    }

    call CaseSites {
        input:
            case_vcfs = case_vcfs,
            case_vcf_idxs = case_vcf_idxs,
            site_merge_gap = site_merge_gap,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            prefix = prefix,
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_case_sites
    }

    scatter (i in range(length(control_vcfs))) {
        call ControlStats {
            input:
                vcf = control_vcfs[i],
                vcf_idx = control_vcf_idxs[i],
                case_regions_bed = CaseSites.regions_bed,
                case_sites_tsv = CaseSites.sites_tsv,
                control_weights = MatchControls.control_weights,
                pass_only = control_pass_only,
                gvcf = control_gvcf,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                requester_pays_project = requester_pays_project,
                regions_per_batch = regions_per_batch,
                n_parallel = n_parallel,
                prefix = "~{prefix}.~{control_name}.shard_~{i}",
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_control_stats
        }
    }

    call ConcatStats {
        input:
            stats_tsvs = ControlStats.stats_tsv,
            control_name = control_name,
            prefix = prefix,
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat_stats
    }

    scatter (i in range(length(case_vcfs))) {
        call AnnotateCase {
            input:
                vcf = case_vcfs[i],
                vcf_idx = case_vcf_idxs[i],
                stats_tsv = ConcatStats.stats_tsv,
                stats_tsv_idx = ConcatStats.stats_tsv_idx,
                control_name = control_name,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.~{case_name}.shard_~{i}",
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_annotate_case
        }
    }

    call ConcatCase {
        input:
            vcfs = AnnotateCase.annotated_vcf,
            vcf_idxs = AnnotateCase.annotated_vcf_idx,
            prefix = "~{prefix}.~{case_name}.~{control_name}_annotated",
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat_case
    }

    output {
        File annotated_case_vcf = ConcatCase.vcf
        File annotated_case_vcf_idx = ConcatCase.vcf_idx
        File matched_controls = MatchControls.matched_controls
        File matched_pairs = MatchControls.matched_pairs
        File control_weights = MatchControls.control_weights
        File study_ancestry = Pca.study_ancestry
        File pcs_tsv = Pca.pcs_tsv
        File eigenvalues = Pca.eigenvalues
        File snp_weights = Pca.snp_weights
        File smartpca_log = Pca.smartpca_log
        File pca_variant_ids = Pca.pca_variant_ids
        File variant_counts = Pca.variant_counts
        File shared_sites_tsv = IntersectSites.shared_sites_tsv
        File site_counts = IntersectSites.site_counts
        File control_stats_tsv = ConcatStats.stats_tsv
        Array[File] common_snv_vcfs = SelectCommonSnvs.snv_vcf
        Array[File] common_snv_vcf_idxs = SelectCommonSnvs.snv_vcf_idx
    }
}

task IntersectSites {
    input {
        File ref_sites_tsv
        File case_sites_tsv
        File control_sites_tsv
        String case_name
        String control_name
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        gzip -cdf ~{ref_sites_tsv} > ref.tsv
        gzip -cdf ~{case_sites_tsv} > case.tsv
        gzip -cdf ~{control_sites_tsv} > control.tsv

        # sites_tsv columns: chr, pos, end, ref, alt, then the AF in each call set
        (
            printf 'chr\tpos\tend\tref\talt\tAF_ref\tAF_~{case_name}\tAF_~{control_name}\n'
            awk -F'\t' -v OFS='\t' '
                FNR == 1 {f++}
                $2 !~ /^[0-9]+$/ {next}
                {k = $1 ":" $2 ":" $4 ":" $5}
                f == 1 {af1[k] = $6; next}
                f == 2 {if (k in af1) af2[k] = $6; next}
                (k in af2) {print $1, $2, $3, $4, $5, af1[k], af2[k], $6}' ref.tsv case.tsv control.tsv
        ) | gzip -c > ~{prefix}.shared_sites.tsv.gz

        n() { awk -F'\t' '$2 ~ /^[0-9]+$/' "$1" | wc -l; }
        printf 'set\tn_sites\nref\t%s\n~{case_name}\t%s\n~{control_name}\t%s\nshared\t%s\n' \
            "$(n ref.tsv)" "$(n case.tsv)" "$(n control.tsv)" "$(gzip -cd ~{prefix}.shared_sites.tsv.gz | sed 1d | wc -l)" \
            | tee ~{prefix}.site_counts.tsv
    >>>

    output {
        File shared_sites_tsv = "~{prefix}.shared_sites.tsv.gz"
        File site_counts = "~{prefix}.site_counts.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 20,
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

task MatchControls {
    input {
        File study_ancestry
        File eigenvalues
        String case_name
        String control_name
        Int n_controls
        Int match_n_pcs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
        import collections
        import csv
        import heapq
        import math

        rows = list(csv.DictReader(open("~{study_ancestry}"), delimiter="\t"))
        evals = [float(x) for x in open("~{eigenvalues}").read().split()]
        pcs = [c for c in rows[0] if c.startswith("PC")][:~{match_n_pcs}]
        w = [math.sqrt(max(evals[i], 0.0)) for i in range(len(pcs))]

        def coords(r):
            return [float(r[c]) * wi for c, wi in zip(pcs, w)]

        cases = [r for r in rows if r["cohort"] == "~{case_name}"]
        controls = [(r["sample_id"], r["assigned_ancestry"], coords(r)) for r in rows if r["cohort"] == "~{control_name}"]
        if len(controls) < ~{n_controls}:
            raise SystemExit(f"only {len(controls)} ~{control_name} samples, fewer than n_controls = ~{n_controls}")

        times = collections.Counter()
        with open("~{prefix}.matched_controls.tsv", "w") as m, open("~{prefix}.matched_pairs.tsv", "w") as p:
            m.write("~{case_name}_sample\tmatched_~{control_name}_samples\n")
            p.write("~{case_name}_sample\t~{case_name}_ancestry\trank\t~{control_name}_sample\t~{control_name}_ancestry\tpc_distance\n")
            for r in cases:
                x = coords(r)
                nearest = heapq.nsmallest(~{n_controls}, ((math.dist(x, c[2]), c) for c in controls), key=lambda t: t[0])
                m.write(r["sample_id"] + "\t" + ",".join(c[0] for _, c in nearest) + "\n")
                for rank, (d, c) in enumerate(nearest, 1):
                    p.write(f"{r['sample_id']}\t{r['assigned_ancestry']}\t{rank}\t{c[0]}\t{c[1]}\t{d:.6g}\n")
                    times[c[0]] += 1

        with open("~{prefix}.control_weights.tsv", "w") as o:
            o.write("~{control_name}_sample\ttimes_matched\n")
            for s, n in sorted(times.items(), key=lambda t: (-t[1], t[0])):
                o.write(f"{s}\t{n}\n")
        print(f"{len(cases)} cases, {len(controls)} controls, {len(times)} distinct controls matched, "
              f"max times one control is matched: {max(times.values()) if times else 0}; PCs used: {len(pcs)}")
        CODE
    >>>

    output {
        File matched_controls = "~{prefix}.matched_controls.tsv"
        File matched_pairs = "~{prefix}.matched_pairs.tsv"
        File control_weights = "~{prefix}.control_weights.tsv"
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

task CaseSites {
    input {
        Array[File] case_vcfs
        Array[File] case_vcf_idxs
        Int site_merge_gap
        File? ref_fa
        File? ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String norm_ref_arg = if defined(ref_fa) then "-f ~{ref_fa}" else ""

    command <<<
        set -euo pipefail

        # every case variant as biallelic CHROM, POS, REF, ALT (same splitting as in AnnotateCase)
        while read -r vcf; do
            bcftools norm ~{norm_ref_arg} -m-any -Ou "$vcf" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n'
        done < ~{write_lines(case_vcfs)} | sort -k1,1 -k2,2n -k3,3 -k4,4 -u | gzip -c > ~{prefix}.case_sites.tsv.gz

        # positions within site_merge_gap bp joined into one fetch region (0-based BED)
        gzip -cd ~{prefix}.case_sites.tsv.gz | cut -f1,2 | uniq \
            | awk -F'\t' -v OFS='\t' -v gap=~{site_merge_gap} '
                $1 != c || $2 - e > gap {if (c != "") print c, s, e; c = $1; s = $2 - 1}
                {e = $2}
                END {if (c != "") print c, s, e}' > ~{prefix}.case_regions.bed
        echo "$(gzip -cd ~{prefix}.case_sites.tsv.gz | wc -l) case variants in $(wc -l < ~{prefix}.case_regions.bed) regions"
    >>>

    output {
        File sites_tsv = "~{prefix}.case_sites.tsv.gz"
        File regions_bed = "~{prefix}.case_regions.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(case_vcfs, "GB")) + ceil(size(ref_fa, "GB")) + 20,
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

task ControlStats {
    input {
        String vcf
        String vcf_idx
        File case_regions_bed
        File case_sites_tsv
        File control_weights
        Boolean pass_only
        Boolean gvcf
        File? ref_fa
        File? ref_fai
        String? requester_pays_project
        Int regions_per_batch
        Int n_parallel
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String pass_arg = if pass_only then "-f PASS" else ""
    String gvcf_cmd = if gvcf then "bcftools view -a -Ou | " else ""
    String norm_ref_arg = if defined(ref_fa) then "-f ~{ref_fa}" else ""

    command <<<
        set -euo pipefail

        cat > get_token.sh <<'EOF'
        wget -q -O - --header "Metadata-Flavor: Google" \
            "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token" \
            | sed -E 's/.*"access_token" *: *"([^"]*)".*/\1/'
        EOF
        if [ -n "~{requester_pays_project}" ]; then
            export GCS_REQUESTER_PAYS_PROJECT="~{requester_pays_project}"
        fi
        export INPUT="~{vcf}##idx##~{vcf_idx}"
        export GCS_OAUTH_TOKEN=$(bash get_token.sh)

        # sample groups for bcftools +fill-tags: ALL = every control sample; CTRL = matched controls; W<m> = controls
        # matched to m cases (summed with weight m for the counts that count each match)
        bcftools query -l "$INPUT" > samples.txt
        awk -F'\t' 'FNR == NR {if (FNR > 1) w[$1] = $2; next}
            {g = "ALL"; if ($1 in w) g = g ",CTRL,W" w[$1]; print $1 "\t" g}' ~{control_weights} samples.txt > groups.txt
        echo "$(wc -l < samples.txt) samples, $(grep -c CTRL groups.txt || true) matched controls in this file"

        gzip -cd ~{case_sites_tsv} | awk -F'\t' '{print $1 ":" $2 ":" $3 ":" $4 > "case.ids"; print $1 "\t" $2 > "case.pos.tmp"}'
        uniq case.pos.tmp > case.pos
        bcftools index -s "$INPUT" | cut -f1 > file_contigs.txt
        awk 'FNR == NR {keep[$1] = 1; next} ($1 in keep)' file_contigs.txt ~{case_regions_bed} > regions.bed
        mkdir batches
        awk -v n=~{regions_per_batch} '
            $1 != prev || c == n {if (f) close(f); b++; c = 0; prev = $1; f = sprintf("batches/batch.%06d.bed", b)}
            {c++; print > f}' regions.bed
        echo "$(wc -l < regions.bed) case regions on contigs in this file, $(ls batches | wc -l) batches"

        # per variant: AF, AC, AN, N_HET, N_HOMREF, N_HOMALT for ALL, CTRL (each matched control once) and the weighted
        # W<m> groups (each matched control counted once per case); N_HOMREF = NS - N_HET - N_HOMALT - hemizygous ALT
        cat > summarize.awk <<'EOF'
        function stats(ac, an, ns, het, hom, hemi,    homalt, af) {
            ac += 0; an += 0; ns += 0; het += 0; hemi += 0; homalt = (hom + 0) / 2
            af = an > 0 ? sprintf("%.6g", ac / an) : "."
            return af OFS ac OFS an OFS het OFS (ns - het - homalt - hemi) OFS homalt
        }
        BEGIN {FS = OFS = "\t"}
        {
            delete v
            n = split($5, kv, ";")
            for (i = 1; i <= n; i++) {split(kv[i], p, "="); v[p[1]] = p[2]}
            ac = an = ns = het = hom = hemi = 0
            for (k in v) if (k ~ /^AN_W[0-9]+$/) {
                m = substr(k, 5) + 0; g = "_W" m
                an += m * v[k]; ac += m * v["AC" g]; ns += m * v["NS" g]
                het += m * v["AC_Het" g]; hom += m * v["AC_Hom" g]; hemi += m * v["AC_Hemi" g]
            }
            print $1, $2, $3, $4, \
                stats(v["AC_ALL"], v["AN_ALL"], v["NS_ALL"], v["AC_Het_ALL"], v["AC_Hom_ALL"], v["AC_Hemi_ALL"]), \
                stats(v["AC_CTRL"], v["AN_CTRL"], v["NS_CTRL"], v["AC_Het_CTRL"], v["AC_Hom_CTRL"], v["AC_Hemi_CTRL"]), \
                stats(ac, an, ns, het, hom, hemi)
        }
        EOF

        cat > stats_batch.sh <<'EOF'
        set -euo pipefail
        bed="$1"
        out="${bed%.bed}.tsv"
        for attempt in 1 2 3; do
            export GCS_OAUTH_TOKEN=$(bash get_token.sh)
            if bcftools view -R "$bed" -T case.pos ~{pass_arg} -Ou "$INPUT" \
                | ~{gvcf_cmd}bcftools annotate -x ID,INFO,^FORMAT/GT -Ou \
                | bcftools norm ~{norm_ref_arg} -m-any -Ou \
                | bcftools view -e 'ALT="*" || ALT="<*>" || ALT="<NON_REF>"' -Ou \
                | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Ou \
                | bcftools view -i 'ID=@case.ids' -Ou \
                | bcftools +fill-tags -Ou -- -S groups.txt -t AN,AC,NS,AC_Het,AC_Hom,AC_Hemi \
                | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\n' \
                | awk -f summarize.awk > "$out"; then
                exit 0
            fi
            echo "retrying $bed (attempt $attempt failed)" >&2
            sleep $((attempt * 30))
        done
        exit 1
        EOF
        if [ -s regions.bed ]; then
            ls batches/*.bed | xargs -P ~{n_parallel} -n 1 bash stats_batch.sh
        fi
        touch batches/none.tsv
        cat batches/*.tsv | gzip -c > ~{prefix}.control_stats.tsv.gz
        echo "$(gzip -cd ~{prefix}.control_stats.tsv.gz | wc -l) case variants found in this control file"
    >>>

    output {
        File stats_tsv = "~{prefix}.control_stats.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 50,
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

task ConcatStats {
    input {
        Array[File] stats_tsvs
        String control_name
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cols=""
        for set in "" _CTRL _CTRL_DUP; do
            for f in AF AC AN N_HET N_HOMREF N_HOMALT; do
                cols="$cols\t~{control_name}${set}_${f}"
            done
        done
        (
            printf "#CHROM\tPOS\tREF\tALT${cols}\n"
            while read -r f; do gzip -cd "$f"; done < ~{write_lines(stats_tsvs)} | sort -k1,1 -k2,2n -k3,3 -k4,4 -u
        ) | bgzip -c > ~{prefix}.~{control_name}_stats.tsv.gz
        tabix -s1 -b2 -e2 ~{prefix}.~{control_name}_stats.tsv.gz
        echo "$(bgzip -cd ~{prefix}.~{control_name}_stats.tsv.gz | grep -vc '^#') case variants with a ~{control_name} record"
    >>>

    output {
        File stats_tsv = "~{prefix}.~{control_name}_stats.tsv.gz"
        File stats_tsv_idx = "~{prefix}.~{control_name}_stats.tsv.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(stats_tsvs, "GB")) + 10,
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

task AnnotateCase {
    input {
        File vcf
        File vcf_idx
        File stats_tsv
        File stats_tsv_idx
        String control_name
        File? ref_fa
        File? ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String norm_ref_arg = if defined(ref_fa) then "-f ~{ref_fa}" else ""

    command <<<
        set -euo pipefail

        P="~{control_name}"
        cols="CHROM,POS,REF,ALT"
        : > header.txt
        for set in "" _CTRL _CTRL_DUP; do
            case "$set" in
                "") what="all $P samples" ;;
                _CTRL) what="$P samples matched as controls, each counted once" ;;
                _CTRL_DUP) what="$P samples matched as controls, each counted once per case it is matched to" ;;
            esac
            for f in AF AC AN N_HET N_HOMREF N_HOMALT; do
                case "$f" in
                    AF) type=Float; desc="Allele frequency" ;;
                    AC) type=Integer; desc="Alternate allele count" ;;
                    AN) type=Integer; desc="Number of called alleles" ;;
                    N_HET) type=Integer; desc="Number of heterozygous genotypes" ;;
                    N_HOMREF) type=Integer; desc="Number of homozygous reference genotypes" ;;
                    N_HOMALT) type=Integer; desc="Number of homozygous alternate genotypes" ;;
                esac
                echo "##INFO=<ID=${P}${set}_${f},Number=A,Type=${type},Description=\"${desc} in ${what}\">" >> header.txt
                cols="${cols},INFO/${P}${set}_${f}"
            done
        done

        bcftools norm ~{norm_ref_arg} -m-any -Ou ~{vcf} \
            | bcftools annotate -a ~{stats_tsv} -h header.txt -c "$cols" -Oz -o ~{prefix}.annotated.vcf.gz
        bcftools index -t ~{prefix}.annotated.vcf.gz
        n_all=$(bcftools index -n ~{prefix}.annotated.vcf.gz)
        n_hit=$(bcftools view -H -i "INFO/${P}_AN>=0" ~{prefix}.annotated.vcf.gz | wc -l)
        echo "${n_all} case variants (after splitting), ${n_hit} with a ${P} record"
    >>>

    output {
        File annotated_vcf = "~{prefix}.annotated.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.annotated.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + ceil(size(stats_tsv, "GB")) + ceil(size(ref_fa, "GB")) + 10,
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

task ConcatCase {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools concat -a --file-list ~{write_lines(vcfs)} --threads 2 -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
        echo "$(bcftools index -n ~{prefix}.vcf.gz) records"
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcfs, "GB")) + 10,
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
