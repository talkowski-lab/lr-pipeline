version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CombineVcfsAcrossContigs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        String prefix

        Boolean drop_genotypes = false

        String utils_docker

        RuntimeAttr? runtime_attr_combine_vcfs
    }

    call CombineVcfsAcrossContigsTask {
        input:
            vcfs = vcfs,
            vcf_idxs = vcf_idxs,
            contigs = contigs,
            drop_genotypes = drop_genotypes,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_combine_vcfs
    }

    output {
        File concat_vcf = CombineVcfsAcrossContigsTask.concat_vcf
        File concat_vcf_idx = CombineVcfsAcrossContigsTask.concat_vcf_idx
    }
}

task CombineVcfsAcrossContigsTask {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        Boolean drop_genotypes = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import re
import sys
from collections import Counter


def natural_chunks(value):
    return tuple(
        (0, int(token)) if token.isdigit() else (1, token.lower())
        for token in re.findall(r"\d+|\D+", value)
    )


def contig_key(contig):
    normalized = re.sub(r"^chr", "", contig, flags=re.IGNORECASE)
    upper = normalized.upper()

    if upper.isdigit():
        return (0, int(upper))
    if upper == "X":
        return (1, 23)
    if upper == "Y":
        return (1, 24)
    if upper in {"M", "MT"}:
        return (1, 25)

    return (2, natural_chunks(normalized))


vcfs = [line.rstrip("\n") for line in open("~{write_lines(vcfs)}")]
vcf_idxs = [line.rstrip("\n") for line in open("~{write_lines(vcf_idxs)}")]
contigs = [line.rstrip("\n") for line in open("~{write_lines(contigs)}")]

if not vcfs:
    sys.stderr.write("ERROR: vcfs must contain at least one file.\n")
    sys.exit(1)

if len(vcfs) != len(vcf_idxs) or len(vcfs) != len(contigs):
    sys.stderr.write(
        f"ERROR: vcfs ({len(vcfs)}), vcf_idxs ({len(vcf_idxs)}), and contigs ({len(contigs)}) must have the same length.\n"
    )
    sys.exit(1)

duplicate_contigs = sorted(contig for contig, count in Counter(contigs).items() if count > 1)
if duplicate_contigs:
    sys.stderr.write(
        "ERROR: CombineVcfsAcrossContigs expects exactly one VCF per contig. "
        f"Duplicate contigs: {', '.join(duplicate_contigs)}\n"
    )
    sys.exit(1)

records = sorted(zip(contigs, vcfs, vcf_idxs), key=lambda record: contig_key(record[0]))

with open("sorted_inputs.tsv", "w") as out:
    for index, (contig, vcf, vcf_idx) in enumerate(records):
        out.write(f"{index:06d}\t{contig}\t{vcf}\t{vcf_idx}\n")
CODE

        mkdir -p prepared
        : > concat_vcfs.txt

        while IFS=$'\t' read -r sort_index contig vcf vcf_idx; do
            input_vcf="prepared/${sort_index}.input.vcf.gz"
            ln -sfn "$vcf" "$input_vcf"
            ln -sfn "$vcf_idx" "${input_vcf}.tbi"

            if [[ "~{drop_genotypes}" == "true" ]]; then
                stripped_vcf="prepared/${sort_index}.stripped.vcf.gz"
                bcftools view \
                    -G \
                    -Oz -o "$stripped_vcf" \
                    "$input_vcf"
                tabix -p vcf -f "$stripped_vcf"
                echo "$stripped_vcf" >> concat_vcfs.txt
            else
                echo "$input_vcf" >> concat_vcfs.txt
            fi
        done < sorted_inputs.tsv

        bcftools concat \
            --naive \
            --file-list concat_vcfs.txt \
            -Oz -o "~{prefix}.vcf.gz"

        tabix -p vcf -f "~{prefix}.vcf.gz"
    >>>

    output {
        File concat_vcf = "~{prefix}.vcf.gz"
        File concat_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcfs, "GB")) + 25,
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
