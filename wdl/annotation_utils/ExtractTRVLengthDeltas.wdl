version 1.0

import "../utils/Structs.wdl"

workflow ExtractTRVLengthDeltas {
    meta {
        description: [
            "This utility extracts ALT allele length deltas (`len(ALT) - len(REF)`) across the tandem repeat records (`INFO/allele_type=trv`) of one single-contig VCF whose shortest `INFO/MOTIFS` entry is 3 bp long. It outputs a per-allele table and a per-length-delta count table, both additive across contigs, so per-contig outputs aggregate by concatenation and by summing within each length delta respectively.",
            "Each count table row counts distinct ALT alleles (`n_alleles`) and their summed `INFO/AC` (`allele_count`) at one length delta. Genotypes are never decoded, so the task scales with the number of records rather than samples."
        ]
    }

    parameter_meta {
        vcf: "Single-contig VCF whose tandem repeat records are tabulated."
        vcf_idx: "Index for `vcf`."
        allele_deltas_tsv: "Gzipped per-allele table with columns `contig`, `pos`, `trid`, `motifs`, `ref_length`, `alt_length`, `length_delta` and `allele_count`."
        length_deltas_tsv: "Gzipped count table with columns `contig`, `length_delta`, `n_alleles` and `allele_count`."
    }

    input {
        File vcf
        File vcf_idx
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_extract
    }

    call ExtractLengthDeltas {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_extract
    }

    output {
        File allele_deltas_tsv = ExtractLengthDeltas.alleles_tsv
        File length_deltas_tsv = ExtractLengthDeltas.counts_tsv
    }
}

task ExtractLengthDeltas {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query \
            -i 'INFO/allele_type="trv"' \
            -f '%CHROM\t%POS\t%INFO/TRID\t%REF\t%ALT\t%INFO/MOTIFS\t%INFO/AC\n' \
            ~{vcf} > trv.tsv

        python3 <<'PYCODE'
import csv
import gzip
from collections import defaultdict

INPUT = "trv.tsv"
ALLELES_OUTPUT = "~{prefix}.trv_allele_deltas.tsv.gz"
COUNTS_OUTPUT = "~{prefix}.trv_length_deltas.tsv.gz"
MOTIF_LENGTH = 3

n_alleles = defaultdict(int)
allele_count = defaultdict(int)
with open(INPUT) as handle, gzip.open(ALLELES_OUTPUT, "wt", newline="") as alleles_handle:
    alleles_writer = csv.writer(alleles_handle, delimiter="\t")
    alleles_writer.writerow(["contig", "pos", "trid", "motifs", "ref_length", "alt_length", "length_delta", "allele_count"])
    for line in handle:
        contig, pos, trid, ref, alts, motifs, acs = line.rstrip("\n").split("\t")
        if min(len(motif) for motif in motifs.split(",")) != MOTIF_LENGTH:
            continue
        for alt, ac in zip(alts.split(","), acs.split(",")):
            length_delta = len(alt) - len(ref)
            alleles_writer.writerow([contig, pos, trid, motifs, len(ref), len(alt), length_delta, ac])
            key = (contig, length_delta)
            n_alleles[key] += 1
            allele_count[key] += int(ac)

with gzip.open(COUNTS_OUTPUT, "wt", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["contig", "length_delta", "n_alleles", "allele_count"])
    for key in sorted(n_alleles):
        writer.writerow([*key, n_alleles[key], allele_count[key]])
PYCODE
    >>>

    output {
        File alleles_tsv = "~{prefix}.trv_allele_deltas.tsv.gz"
        File counts_tsv = "~{prefix}.trv_length_deltas.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
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
