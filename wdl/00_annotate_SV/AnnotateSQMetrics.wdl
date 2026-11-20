version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateSQMetrics {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        String utils_docker

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig as SubsetVcf {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset
            }
        }

        File contig_vcf = select_first([SubsetVcf.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcf.subset_vcf_idx, vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.site_quality",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [contig_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [contig_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call CalculateSiteMetrics {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    prefix = "~{prefix}.~{contig}.site_quality.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_annotate
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = CalculateSiteMetrics.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.site_quality",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, CalculateSiteMetrics.annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs as MergeAnnotations {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.site_quality",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotations_tsv_sq = select_first([MergeAnnotations.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task CalculateSiteMetrics {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            -x INFO/AS_VarDP,INFO/AS_QUALapprox,INFO/AS_QD,INFO/inbreeding_coeff,INFO/HWE,INFO/AS_pab_max \
            -Oz -o stripped.vcf.gz \
            "~{vcf}"
        tabix -p vcf stripped.vcf.gz

        python3 <<CODE
import pysam
import scipy.stats as stats

vcf_in = pysam.VariantFile("stripped.vcf.gz")

tsv_out = open("~{prefix}.annotations.tsv", "w")


def format_allele_values(values):
    formatted = ['.' if value is None else str(value) for value in values]
    return '.' if all(value == '.' for value in formatted) else ','.join(formatted)

for rec in vcf_in:
    if not rec.alts:
        continue

    n_alt_alleles = len(rec.alts)
    has_dp = 'DP' in rec.format
    has_pl = 'PL' in rec.format
    has_ad = 'AD' in rec.format

    as_vardp = [None] * n_alt_alleles
    as_qualapprox = [None] * n_alt_alleles
    as_qd_sum = [0.0] * n_alt_alleles
    as_qd_n = [0] * n_alt_alleles
    n_ref_per_allele = [0] * n_alt_alleles
    n_het_per_allele = [0] * n_alt_alleles
    n_alt_per_allele = [0] * n_alt_alleles
    n_ref_site = n_het_site = n_alt_site = 0
    as_pab_max = [None] * n_alt_alleles

    for sample in rec.samples.values():
        gt = sample.get('GT')
        if not gt or None in gt:
            continue

        n_non_ref = sum(1 for allele in gt if allele is not None and allele > 0)
        if n_non_ref == 0:
            n_ref_site += 1
        elif n_non_ref == 1:
            n_het_site += 1
        else:
            n_alt_site += 1

        pl = sample.get('PL') if has_pl else None
        hom_ref_pl = None
        if pl and len(pl) > 0 and pl[0] is not None:
            hom_ref_pl = int(pl[0])

        dp = sample.get('DP') if has_dp else None
        ad = sample.get('AD') if has_ad else None
        ref_copies = sum(1 for allele in gt if allele == 0)

        for allele_idx in range(1, n_alt_alleles + 1):
            allele_i = allele_idx - 1
            n_copies = sum(1 for allele in gt if allele == allele_idx)

            if n_copies == 0:
                n_ref_per_allele[allele_i] += 1
            elif n_copies == 1:
                n_het_per_allele[allele_i] += 1
            else:
                n_alt_per_allele[allele_i] += 1

            if n_copies > 0:
                if dp is not None:
                    if as_vardp[allele_i] is None:
                        as_vardp[allele_i] = 0
                    as_vardp[allele_i] += int(dp)
                if hom_ref_pl is not None:
                    if as_qualapprox[allele_i] is None:
                        as_qualapprox[allele_i] = 0
                    as_qualapprox[allele_i] += hom_ref_pl
                    if ad and len(ad) > allele_idx and ad[allele_idx] is not None and ad[allele_idx] > 0:
                        as_qd_sum[allele_i] += hom_ref_pl / float(ad[allele_idx])
                        as_qd_n[allele_i] += 1

            if n_copies == 1 and ref_copies == 1 and ad is not None and len(ad) > allele_idx:
                ref_dp = ad[0]
                alt_dp = ad[allele_idx]
                if ref_dp is not None and alt_dp is not None:
                    total_dp = ref_dp + alt_dp
                    if total_dp > 0:
                        pval = round(float(stats.binomtest(ref_dp, total_dp, 0.5).pvalue), 6)
                        if as_pab_max[allele_i] is None or pval > as_pab_max[allele_i]:
                            as_pab_max[allele_i] = pval

    inbreeding_coeff = []
    for i in range(n_alt_alleles):
        n_ref_allele = n_ref_per_allele[i]
        n_het_allele = n_het_per_allele[i]
        n_alt_allele = n_alt_per_allele[i]
        n_tot_allele = n_ref_allele + n_het_allele + n_alt_allele

        if n_tot_allele == 0:
            inbreeding_coeff.append(0.0)
            continue

        p_allele = (2.0 * n_ref_allele + n_het_allele) / (2.0 * n_tot_allele)
        q_allele = 1.0 - p_allele
        e_het_allele = 2.0 * p_allele * q_allele * n_tot_allele

        if e_het_allele > 0:
            inbreeding_coeff.append(round(1.0 - (n_het_allele / e_het_allele), 6))
        else:
            inbreeding_coeff.append(0.0)

    pab_max_vals = [v if v is not None else '.' for v in as_pab_max]

    hwe_val = '.'
    n_tot = n_ref_site + n_het_site + n_alt_site
    if n_tot > 0:
        p = (2.0 * n_ref_site + n_het_site) / (2.0 * n_tot)
        q = 1.0 - p
        e_het = 2.0 * p * q * n_tot
        e_ref = (p ** 2) * n_tot
        e_alt = (q ** 2) * n_tot
        chisq = sum(((o - e) ** 2) / e for o, e in zip([n_ref_site, n_het_site, n_alt_site], [e_ref, e_het, e_alt]) if e > 0)
        hwe_val = round(float(stats.chi2.sf(chisq, 1)), 6)

    as_qd = [round(as_qd_sum[i] / as_qd_n[i], 6) if as_qd_n[i] > 0 else None for i in range(n_alt_alleles)]

    row = [
        rec.chrom,
        str(rec.pos),
        rec.ref,
        ','.join(rec.alts),
        rec.id if rec.id else '.',
        format_allele_values(inbreeding_coeff),
        format_allele_values(pab_max_vals),
        format_allele_values(as_qualapprox),
        format_allele_values(as_qd),
        format_allele_values(as_vardp),
        str(hwe_val),
    ]
    tsv_out.write('\t'.join(row) + '\n')

tsv_out.close()
CODE
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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
