version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ReplaceSampleCalls {
    input {
        Array[File] sample_vcfs
        Array[File] sample_vcf_idxs
        File cohort_vcf
        File cohort_vcf_idx
        String prefix

        String docker

        RuntimeAttr? runtime_attr_replace
    }

    call ReplaceCalls {
        input:
            sample_vcfs = sample_vcfs,
            sample_vcf_idxs = sample_vcf_idxs,
            cohort_vcf = cohort_vcf,
            cohort_vcf_idx = cohort_vcf_idx,
            prefix = prefix,
            docker = docker,
            runtime_attr_override = runtime_attr_replace
    }

    output {
        File replaced_vcf = ReplaceCalls.replaced_vcf
        File replaced_vcf_idx = ReplaceCalls.replaced_vcf_idx
    }
}

task ReplaceCalls {
    input {
        Array[File] sample_vcfs
        Array[File] sample_vcf_idxs
        File cohort_vcf
        File cohort_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

# Build per-sample lookup: {sample_name: {variant_id: data}}
all_sample_data = {}
for sf in [line.strip() for line in open("~{write_lines(sample_vcfs)}")]:
    sample_in = pysam.VariantFile(sf)
    sample_name = list(sample_in.header.samples)[0]
    sample_variant_ids = {}
    for record in sample_in:
        s = record.samples[sample_name]
        data = { 'GT': s['GT'], 'phased': s.phased, 'alts': record.alts }
        if 'PS' in record.format:
            data['PS'] = s['PS']
        if 'PF' in record.format:
            data['PF'] = s['PF']
        sample_variant_ids[record.id] = data
    sample_in.close()
    all_sample_data[sample_name] = sample_variant_ids

cohort_in = pysam.VariantFile("~{cohort_vcf}")
vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", 'w', header=cohort_in.header)

for record in cohort_in:
    for sample_name, sample_variant_ids in all_sample_data.items():
        match_data = sample_variant_ids.get(record.id)

        if match_data is not None:
            cohort_sample = record.samples[sample_name]
            sample_gt = match_data['GT']
            sample_phased = match_data['phased']

            # Only modify het calls
            real_alleles = [a for a in cohort_sample['GT'] if a is not None]
            if len(real_alleles) >= 2 and len(set(real_alleles)) > 1:
                if sample_gt is not None and len(sample_gt) == 2:
                    sample_alts = tuple(a.upper() for a in (match_data['alts'] or ()))
                    cohort_alt_idx = {seq.upper(): i + 1 for i, seq in enumerate(record.alts or ())}
                    def translate(a):
                        if a is None or a == 0:
                            return a
                        return cohort_alt_idx.get(sample_alts[a - 1])
                    cohort_sample['GT']  = (translate(sample_gt[0]), translate(sample_gt[1]))
                    cohort_sample.phased = sample_phased

            # Clear then copy PS/PF from sample_vcf
            if 'PS' in cohort_sample:
                cohort_sample['PS'] = None
            if 'PF' in cohort_sample:
                cohort_sample['PF'] = None
            if 'PS' in match_data:
                cohort_sample['PS'] = match_data['PS']
            if 'PF' in match_data:
                cohort_sample['PF'] = match_data['PF']

        elif record.info.get('allele_type') == 'trv':
            # Clear GT/PS/PF for unmatched trv variants where sample is called
            cohort_sample = record.samples[sample_name]
            gt = cohort_sample['GT']
            if not (gt is not None and len(gt) > 0 and all(a == 0 for a in gt)):
                cohort_sample['GT'] = tuple(None for _ in gt) if gt else (None, None)
                if 'PS' in record.format:
                    cohort_sample['PS'] = None
                if 'PF' in record.format:
                    cohort_sample['PF'] = None

    vcf_out.write(record)

cohort_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File replaced_vcf = "~{prefix}.vcf.gz"
        File replaced_vcf_idx = "~{prefix}.vcf.gz.tbi"
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
