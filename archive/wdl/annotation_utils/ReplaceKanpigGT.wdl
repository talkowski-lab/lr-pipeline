version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ReplaceKanpigGT {
    input {
        File vcf
        File vcf_idx
        Array[String] sample_ids
        Array[File] sample_vcfs
        Array[File] sample_vcf_idxs
        String contig
        String prefix

        Int min_sv_length

        String utils_docker

        RuntimeAttr? runtime_attr_subset_cohort_to_contig
        RuntimeAttr? runtime_attr_subset_sample_to_contig
        RuntimeAttr? runtime_attr_replace_calls
    }

    scatter (i in range(length(sample_ids))) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = sample_vcfs[i],
                vcf_idx = sample_vcf_idxs[i],
                contig = contig,
                prefix = "~{prefix}.~{sample_ids[i]}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_sample_to_contig
        }
    }

    call ReplaceCalls {
        input:
            cohort_vcf = vcf,
            cohort_vcf_idx = vcf_idx,
            sample_ids = sample_ids,
            sample_vcfs = SubsetVcfToContig.subset_vcf,
            sample_vcf_idxs = SubsetVcfToContig.subset_vcf_idx,
            min_sv_length = min_sv_length,
            prefix = "~{prefix}.kanpig_replaced",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_replace_calls
    }

    output {
        File replaced_vcf = ReplaceCalls.replaced_vcf
        File replaced_vcf_idx = ReplaceCalls.replaced_vcf_idx
        File match_counts_tsv = ReplaceCalls.match_counts_tsv
    }
}

task ReplaceCalls {
    input {
        File cohort_vcf
        File cohort_vcf_idx
        Array[String] sample_ids
        Array[File] sample_vcfs
        Array[File] sample_vcf_idxs
        Int min_sv_length
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam
import math

def is_non_ref(gt):
    return any(a is not None and a > 0 for a in gt)

def count_non_ref(gt):
    return sum(1 for a in gt if a is not None and a > 0)

def calculate_pl(ref_reads, alt_reads):
    n = ref_reads + alt_reads
    if n == 0:
        return (0, 0, 0)
    means = [0.03, 0.50, 0.97]
    priors = [0.33, 0.34, 0.33]
    ll_raw = []
    for i in range(3):
        ll = math.log(priors[i]) + alt_reads * math.log(means[i]) + ref_reads * math.log(1.0 - means[i])
        ll_10 = ll / math.log(10)
        ll_raw.append(ll_10)
    max_ll = max(ll_raw)
    return tuple(int(round(-10 * (ll - max_ll))) for ll in ll_raw)

def calculate_gq(pls):
    return min(sorted(pls)[1], 99)

min_sv_length = ~{min_sv_length}

with open("~{write_lines(sample_ids)}") as fh:
    sample_ids = [line.strip() for line in fh]
with open("~{write_lines(sample_vcfs)}") as fh:
    sample_vcfs = [line.strip() for line in fh]

all_sample_data = {}
for sid, sf in zip(sample_ids, sample_vcfs):
    sample_in = pysam.VariantFile(sf)
    sample_name = list(sample_in.header.samples)[0]
    site_to_data = {}
    for r in sample_in:
        s = r.samples[sample_name]
        ref = r.ref.upper() if r.ref else r.ref
        alts = tuple(a.upper() for a in r.alts) if r.alts else r.alts
        site_to_data[(r.chrom, r.pos, ref, alts)] = {'GT': s['GT'], 'AD': s['AD']}
    sample_in.close()
    all_sample_data[sid] = site_to_data

cohort_in = pysam.VariantFile("~{cohort_vcf}")
vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=cohort_in.header)
match_counts = {sid: [0, 0, 0] for sid in sample_ids}

for rec in cohort_in:
    if abs(rec.info.get('allele_length', 0)) < min_sv_length:
        vcf_out.write(rec)
        continue

    ref_u = rec.ref.upper() if rec.ref else rec.ref
    alts_u = tuple(a.upper() for a in rec.alts) if rec.alts else rec.alts
    key = (rec.chrom, rec.pos, ref_u, alts_u)

    for sid in sample_ids:
        cs = rec.samples[sid]
        gt = cs['GT']
        if not is_non_ref(gt):
            continue
        match_counts[sid][0] += 1

        sd = all_sample_data[sid].get(key)
        if sd is None:
            continue
        match_counts[sid][1] += 1

        if count_non_ref(gt) != count_non_ref(sd['GT']):
            continue
        match_counts[sid][2] += 1

        ad = sd['AD']
        pls = calculate_pl(ad[0], ad[1])
        cs['AD'] = ad
        cs['PL'] = pls
        cs['GQ'] = calculate_gq(pls)

    vcf_out.write(rec)

cohort_in.close()
vcf_out.close()

with open("~{prefix}.match_counts.tsv", "w") as fh:
    fh.write("sample_id\tvariant_count\tvariant_match_count\tvariant_call_match_count\n")
    for sid in sample_ids:
        vc, vmc, vcmc = match_counts[sid]
        fh.write(f"{sid}\t{vc}\t{vmc}\t{vcmc}\n")
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File replaced_vcf = "~{prefix}.vcf.gz"
        File replaced_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File match_counts_tsv = "~{prefix}.match_counts.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(cohort_vcf, "GB") + size(sample_vcfs, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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
