version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow Kanpig {
    meta {
        description: [
            "This tool regenotypes a cohort SV VCF against each sample's aligned reads using Kanpig (https://github.com/ACEnglish/kanpig). It subsets the cohort to the target samples, runs Kanpig per sample with sex-aware ploidy beds, and merges the per-sample genotypes back into both a raw and a processed cohort VCF. It outputs the regenotyped (processed) and raw Kanpig VCFs."
        ]
    }

    parameter_meta {
        cohort_vcf: "Cohort SV VCF to regenotype."
        cohort_vcf_idx: "Index for the cohort VCF."
        bams: "Aligned reads, one per sample."
        bais: "Indexes for `bams`."
        ref_fa: "From references."
        ref_fai: "From references."
        ploidy_bed_male: "From references."
        ploidy_bed_female: "From references."
        sample_ids: "Samples to regenotype."
        sexes: "Sex of each sample in `sample_ids`."
        swap_samples: "Sample-ID swap map applied to the cohort VCF."
        merge_args: "Arguments controlling the per-sample genotype merge."
        kanpig_params: "Parameters passed to Kanpig."
        sv_kanpig_vcf: "Regenotyped (processed) cohort VCF."
        sv_kanpig_vcf_idx: "Index for the processed VCF."
        sv_kanpig_raw_vcf: "Raw Kanpig cohort VCF."
        sv_kanpig_raw_vcf_idx: "Index for the raw VCF."
    }

    input {
        File cohort_vcf
        File cohort_vcf_idx
        Array[File] bams
        Array[File] bais
        File ref_fa
        File ref_fai
        File ploidy_bed_male
        File ploidy_bed_female
        Array[String] sample_ids
        Array[String] sexes
        String prefix

        File? swap_samples
        String merge_args = "--merge id"
        String kanpig_params = "--neighdist 500 --gpenalty 0.04 --hapsim 0.97"

        String kanpig_docker
        String utils_docker

        RuntimeAttr? runtime_attr_swap_samples
        RuntimeAttr? runtime_attr_subset_cohort_to_samples
        RuntimeAttr? runtime_attr_subset_to_sample
        RuntimeAttr? runtime_attr_run_kanpig
        RuntimeAttr? runtime_attr_merge_genotypes
        RuntimeAttr? runtime_attr_merge_raw_vcfs
        RuntimeAttr? runtime_attr_merge_processed_vcfs
    }

    if (defined(swap_samples)) {
        call Helpers.SwapSampleIds {
            input:
                vcf = cohort_vcf,
                vcf_idx = cohort_vcf_idx,
                sample_swap_list = select_first([swap_samples]),
                prefix = "~{prefix}.swapped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_swap_samples
        }
    }

    File final_cohort_vcf = select_first([SwapSampleIds.swapped_vcf, cohort_vcf])
    File final_cohort_vcf_idx = select_first([SwapSampleIds.swapped_vcf_idx, cohort_vcf_idx])

    call Helpers.SubsetVcfToSamples as SubsetCohortToSamples {
        input:
            vcf = final_cohort_vcf,
            vcf_idx = final_cohort_vcf_idx,
            samples = sample_ids,
            prefix = "~{prefix}.subset",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_cohort_to_samples
    }

    scatter (i in range(length(sample_ids))) {
        File ploidy_bed = if sexes[i] == "M" then ploidy_bed_male else ploidy_bed_female

        call Helpers.SubsetVcfToSamples {
            input:
                vcf = SubsetCohortToSamples.subset_vcf,
                vcf_idx = SubsetCohortToSamples.subset_vcf_idx,
                samples = [sample_ids[i]],
                filter_to_sample = false,
                prefix = "~{prefix}.~{sample_ids[i]}.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_to_sample
        }

        call RunKanpig {
            input:
                input_vcf = SubsetVcfToSamples.subset_vcf,
                input_vcf_idx = SubsetVcfToSamples.subset_vcf_idx,
                bam = bams[i],
                bai = bais[i],
                sample_id = sample_ids[i],
                ploidy_bed = ploidy_bed,
                kanpig_params = kanpig_params,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.~{sample_ids[i]}.kanpig",
                docker = kanpig_docker,
                runtime_attr_override = runtime_attr_run_kanpig
        }

        call MergeGenotypes {
            input:
                base_vcf = SubsetVcfToSamples.subset_vcf,
                base_vcf_idx = SubsetVcfToSamples.subset_vcf_idx,
                kanpig_vcf = RunKanpig.regenotyped_vcf,
                kanpig_vcf_idx = RunKanpig.regenotyped_vcf_idx,
                sex = sexes[i],
                prefix = "~{prefix}.~{sample_ids[i]}.kanpig_merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_genotypes
        }
    }

    call Helpers.MergeVcfs as MergeRaw {
        input:
            vcfs = RunKanpig.regenotyped_vcf,
            vcf_idxs = RunKanpig.regenotyped_vcf_idx,
            prefix = "~{prefix}.kanpig_raw",
            extra_args = merge_args,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_raw_vcfs
    }

    call Helpers.MergeVcfs as MergeProcessed {
        input:
            vcfs = MergeGenotypes.merged_vcf,
            vcf_idxs = MergeGenotypes.merged_vcf_idx,
            prefix = "~{prefix}.kanpig_merged",
            extra_args = merge_args,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_processed_vcfs
    }

    output {
        File sv_kanpig_vcf = MergeProcessed.merged_vcf
        File sv_kanpig_vcf_idx = MergeProcessed.merged_vcf_idx
        File sv_kanpig_raw_vcf = MergeRaw.merged_vcf
        File sv_kanpig_raw_vcf_idx = MergeRaw.merged_vcf_idx
    }
}

task RunKanpig {
    input {
        File input_vcf
        File input_vcf_idx
        File bam
        File bai
        String sample_id
        File ploidy_bed
        String kanpig_params
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        kanpig gt \
            --input ~{input_vcf} \
            --out ~{prefix}.kanpig.vcf \
            --reads ~{bam} \
            --reference ~{ref_fa} \
            --ploidy-bed ~{ploidy_bed} \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --sample ~{sample_id} \
            ~{kanpig_params}

        bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -Oz -o ~{prefix}.vcf.gz \
            ~{prefix}.kanpig.vcf

        rm ~{prefix}.kanpig.vcf

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File regenotyped_vcf = "~{prefix}.vcf.gz"
        File regenotyped_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 4,
        disk_gb: ceil(size(bam, "GB")) + ceil(size(input_vcf, "GB")) + 20,
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

task MergeGenotypes {
    input {
        File base_vcf
        File base_vcf_idx
        File kanpig_vcf
        File kanpig_vcf_idx
        String sex
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam
import math
from math import comb

def is_called(gt):
    return any(a is not None and a > 0 for a in gt)

def is_missing(gt):
    return all(a is None for a in gt)

def calculate_pl(ref_reads, alt_reads):
    if ref_reads + alt_reads == 0:
        return (0, 0, 0)
    means = [0.05, 0.50, 0.95]
    log10 = math.log(10)
    ll = [(alt_reads * math.log(means[i]) + ref_reads * math.log(1.0 - means[i])) / log10 for i in range(3)]
    max_ll = max(ll)
    return tuple(int(round(-10 * (x - max_ll))) for x in ll)

def calculate_gq(pls):
    return min(sorted(pls)[1], 99)

def clear_format_fields(rec, sample, n_alleles):
    rec.samples[sample]['DP'] = None
    rec.samples[sample]['GQ'] = None
    rec.samples[sample]['AD'] = tuple(None for _ in range(n_alleles))
    rec.samples[sample]['PL'] = tuple(None for _ in range(comb(n_alleles + 1, 2)))
    rec.samples[sample]['GT'] = tuple(None for _ in range(n_alleles))

sex = "~{sex}"
is_male = (sex == "M")

kp_vcf = pysam.VariantFile("~{kanpig_vcf}")
base_vcf = pysam.VariantFile("~{base_vcf}")

header = base_vcf.header.copy()
if 'DP' not in header.formats:
    header.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">')
if 'AD' not in header.formats:
    header.add_line('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">')
if 'GQ' not in header.formats:
    header.add_line('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
if 'PL' not in header.formats:
    header.add_line('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">')

out = pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=header)
sample = list(base_vcf.header.samples)[0]

for rec in base_vcf:
    rec.translate(out.header)
    chrom = rec.chrom
    n_alleles = len(rec.alts) + 1
    base_gt = rec.samples[sample]['GT']
    kp_rec = next(
        r for r in kp_vcf.fetch(rec.chrom, rec.pos - 1, rec.pos)
        if r.id == rec.id and r.ref == rec.ref and r.alts == rec.alts
    )
    kp_gt = kp_rec.samples[sample]['GT']

    is_hemi = is_male and chrom in {"chrX", "chrY"}
    is_female_y = (not is_male) and chrom == "chrY"

    # Case 1: Ref/missing in base and ref in Kanpig
    if not is_called(base_gt) and not is_called(kp_gt) and not is_missing(kp_gt):
        if is_female_y:
            # Case 1a: Female on chrY, so clear everything
            clear_format_fields(rec, sample, n_alleles)
        else:
            # Case 1b: Autosome, male on chrX/chrY and female on chrX, so set AD/DP/PL/GQ from Kanpig
            rec.samples[sample]['DP'] = kp_rec.samples[sample]['DP']

            if kp_rec.samples[sample]['AD'] is not None and len(kp_rec.samples[sample]['AD']) == n_alleles:
                ad = kp_rec.samples[sample]['AD']
                pls = calculate_pl(ad[0], ad[1])
                rec.samples[sample]['AD'] = ad
                rec.samples[sample]['PL'] = pls
                rec.samples[sample]['GQ'] = calculate_gq(pls)

            if is_hemi:
                # Case 1b_i: Male on chrX/chrY, so set GT to 0/.
                rec.samples[sample]['GT'] = (0, None)
            else:
                # Case 1b_ii: Autosome or female on chrX, so set GT to 0/0
                rec.samples[sample]['GT'] = tuple(0 for _ in range(n_alleles))

    # Case 2: Ref/missing in base and missing/alt in Kanpig, so clear FORMAT fields
    if (not is_called(base_gt) and is_missing(kp_gt)) or (not is_called(base_gt) and is_called(kp_gt)):
        clear_format_fields(rec, sample, n_alleles)

    gt_current = rec.samples[sample]['GT']
    rec.samples[sample]['GT'] = tuple(sorted(gt_current, key=lambda a: (a is None, a if a is not None else 0)))
    rec.samples[sample].phased = False

    out.write(rec)

base_vcf.close()
kp_vcf.close()
out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(base_vcf, "GB") + size(kanpig_vcf, "GB")) + 10,
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
