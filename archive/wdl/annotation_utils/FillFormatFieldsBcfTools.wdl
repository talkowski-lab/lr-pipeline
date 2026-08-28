version 1.0

import "../utils/Structs.wdl"

workflow FillFormatFieldsBcfTools {
    input {
        File unfilled_vcf
        File unfilled_vcf_idx
        File filled_vcf
        File filled_vcf_idx
        String prefix

        Array[String] format_fields
        String? include_field
        String? include_value
        Boolean modify_ev_number = false
        Boolean unphase_gts = false
        Boolean add_pl = false

        String utils_docker

        RuntimeAttr? runtime_attr_fill
    }

    call FillVcfFormatFields {
        input:
            unfilled_vcf = unfilled_vcf,
            unfilled_vcf_idx = unfilled_vcf_idx,
            filled_vcf = filled_vcf,
            filled_vcf_idx = filled_vcf_idx,
            format_fields = format_fields,
            include_field = include_field,
            include_value = include_value,
            modify_ev_number = modify_ev_number,
            unphase_gts = unphase_gts,
            add_pl = add_pl,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_fill
    }

    output {
        File refilled_vcf = FillVcfFormatFields.output_vcf
        File refilled_vcf_idx = FillVcfFormatFields.output_vcf_idx
    }
}

task FillVcfFormatFields {
    input {
        File unfilled_vcf
        File unfilled_vcf_idx
        File filled_vcf
        File filled_vcf_idx
        Array[String] format_fields
        String? include_field
        String? include_value
        Boolean modify_ev_number = false
        Boolean unphase_gts = false
        Boolean add_pl = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        format_fields_file="~{write_lines(format_fields)}"
        filled_vcf_for_fill="~{filled_vcf}"
        threads=$(nproc)

        if [[ "~{modify_ev_number}" == "true" ]] && grep -Fxq "EV" "$format_fields_file"; then
            bcftools view -h "~{filled_vcf}" > filled.header.txt

            python3 <<'CODE'
import re

with open("filled.header.txt") as src:
    lines = src.readlines()

with open("filled.header.txt", "w") as dst:
    for line in lines:
        if line.startswith("##FORMAT=<ID=EV,"):
            line = re.sub(r"Number=[^,>]+", "Number=.", line, count=1)
        dst.write(line)
CODE

            bcftools reheader \
                -h filled.header.txt \
                -o filled.ev_number_fixed.vcf.gz \
                "~{filled_vcf}"
            tabix -p vcf filled.ev_number_fixed.vcf.gz
            filled_vcf_for_fill="filled.ev_number_fixed.vcf.gz"
        else
            if [[ "~{filled_vcf_idx}" != "~{filled_vcf}.tbi" ]]; then
                ln -sf "~{filled_vcf_idx}" "~{filled_vcf}.tbi"
            fi
        fi

        bcftools query -l "~{unfilled_vcf}" | sort > unfilled.samples.txt
        bcftools query -l "$filled_vcf_for_fill" | sort > filled.samples.txt
        comm -12 unfilled.samples.txt filled.samples.txt > common.samples.txt

        python3 <<CODE
with open("$format_fields_file") as fh:
    format_fields = [line.strip() for line in fh if line.strip()]

columns = [f"FORMAT/{field}" for field in format_fields]

with open("annotate.columns.txt", "w") as out:
    out.write(",".join(columns))
CODE

        annotate_columns=$(cat annotate.columns.txt)
        annotate_include_args=()
        if [[ -n "~{default="" include_field}" && -n "~{default="" include_value}" ]]; then
            annotate_include_args=(-i 'INFO/~{default="" include_field}="~{default="" include_value}"' -k)
        fi

        if [[ -n "$annotate_columns" ]]; then
            bcftools annotate \
                -a "$filled_vcf_for_fill" \
                -S common.samples.txt \
                -c "$annotate_columns" \
                --collapse none \
                --threads "$threads" \
                "${annotate_include_args[@]}" \
                -Oz -o annotated.vcf.gz \
                "~{unfilled_vcf}"
            tabix -p vcf -f annotated.vcf.gz
        else
            ln -sf "~{unfilled_vcf}" annotated.vcf.gz
            ln -sf "~{unfilled_vcf_idx}" annotated.vcf.gz.tbi
        fi

        python3 <<CODE
import math
import pysam

with open("$format_fields_file") as fh:
    format_fields = [l.strip() for l in fh if l.strip()]

with open("common.samples.txt") as fh:
    common_samples = [l.strip() for l in fh if l.strip()]

assert "GT" not in format_fields, "GT is not supported in FillFormatFields_V3"

include_field = "~{default="" include_field}" or None
include_value = "~{default="" include_value}" or None
unphase_gts = ~{true="True" false="False" unphase_gts}
add_pl = ~{true="True" false="False" add_pl}

annotated_in = pysam.VariantFile("annotated.vcf.gz")

out_header = annotated_in.header.copy()
if add_pl and "PL" not in out_header.formats:
    out_header.add_line('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">')

out = pysam.VariantFile("~{prefix}.vcf.gz", "w", header=out_header)
all_samples = list(out_header.samples)

def variant_passes_include(rec):
    if include_field is None:
        return True
    val = rec.info.get(include_field)
    if isinstance(val, tuple):
        val = val[0] if val else None
    return str(val) == include_value

def get_sample_ploidy(sample_data):
    gt = sample_data.get("GT")
    if gt is None:
        return None
    return len(gt)

def calculate_pl(ref_reads, alt_reads, ploidy):
    n = ref_reads + alt_reads
    if n == 0:
        return (0, 0) if ploidy == 1 else (0, 0, 0)

    if ploidy == 1:
        means = [0.03, 0.97]
        priors = [0.50, 0.50]
    else:
        means = [0.03, 0.50, 0.97]
        priors = [0.33, 0.34, 0.33]

    ll_raw = []
    for i in range(len(means)):
        ll = math.log(priors[i]) + alt_reads * math.log(means[i]) + ref_reads * math.log(1.0 - means[i])
        ll_10 = ll / math.log(10)
        ll_raw.append(ll_10)

    max_ll = max(ll_raw)
    return tuple(int(round(-10 * (ll - max_ll))) for ll in ll_raw)

def ad_is_populated(ad):
    return ad is not None and len(ad) == 2 and all(value is not None for value in ad)

def unphase_gt(gt):
    return tuple(sorted(gt, key=lambda allele: (allele is None, allele if allele is not None else 0)))

for annotated_rec in annotated_in:
    annotated_rec.translate(out_header)
    passes_include = variant_passes_include(annotated_rec)

    if passes_include:
        # Unphase genotypes
        if unphase_gts:
            for sample in all_samples:
                current_gt = annotated_rec.samples[sample].get("GT")
                if current_gt is not None:
                    annotated_rec.samples[sample]["GT"] = unphase_gt(current_gt)
                    annotated_rec.samples[sample].phased = False

        # Set PL
        if add_pl and "PL" not in annotated_rec.format and "AD" in annotated_rec.format:
            for sample in all_samples:
                ad = annotated_rec.samples[sample].get("AD")
                if ad_is_populated(ad):
                    ploidy = get_sample_ploidy(annotated_rec.samples[sample])
                    if ploidy is not None:
                        annotated_rec.samples[sample]["PL"] = calculate_pl(ad[0], ad[1], ploidy)

    out.write(annotated_rec)

out.close()
CODE

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File output_vcf = "~{prefix}.vcf.gz"
        File output_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(unfilled_vcf, "GB") + size(filled_vcf, "GB")) + 10,
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
