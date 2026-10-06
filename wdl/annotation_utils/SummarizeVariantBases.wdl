version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SummarizeVariantBases {
    meta {
        description: [
            "This utility measures how many bases are altered by variation in each of several callsets for one contig, producing a site-level table of the reference bases altered across the callset and a sample-level table of the bases altered per genome, each with one row per variant class and size bin and, for each entry of `vcf_names`, a column of bases and a column of those bases as a proportion of the contig length, all rounded to six decimal places.",
            "Records are classed by `INFO/allele_type` as SNV, DEL or INS, with any type containing 'dup' counted as INS, and DEL and INS are split by the absolute `INFO/allele_length` into 1-49, 50-499 and 500+ bp bins. Other types, such as tandem repeats, and records with an AC of zero are not counted.",
            "The site-level table measures how much of the reference is altered. Only records whose `INFO/allele_type` is snv or del alter reference bases, an SNV its own position and a deletion the `INFO/allele_length` bases after its anchor base, so insertion rows are always zero and a replaced anchor base is never counted. The value is the number of distinct reference bases altered by any record with a nonzero AC in the bin, so different ALT alleles or overlapping records at one base count it once.",
            "The sample-level table measures how much each genome differs from the reference. Every altered allele counts, with an SNV or insertion contributing 1 or its inserted length and a deletion its deleted length, and the value is the sum over records of AC times those bases divided by the number of samples, computed from `INFO/AC` alone rather than from the genotypes. Two different ALT alleles at one base in a sample, or a homozygous-alternate genotype, therefore count that base twice.",
            "The contig is cut into regions of `shard_bin_size` and each callset is streamed region by region straight from the bucket. Each reference base is counted only in the region containing it and each record's AC only in the region containing its position, so no base or record is counted twice across regions."
        ]
    }

    parameter_meta {
        vcfs: "Callsets to summarize. Read region by region, so each may hold any set of contigs."
        vcf_idxs: "Index for vcfs, each stored at its VCF's path with a '.tbi' suffix."
        vcf_names: "Name of each entry of `vcfs`, in the same order, used as its column header."
        contig: "Contig being summarized."
        subset_vcf_string: "`bcftools view` arguments applied to each callset before counting, such as an include expression or a sample list. Must not contain -r, -t, -G or -o."
        ref_fai: "From references."
        shard_bin_size: "Width in base pairs of the regions the contig is sharded into."
        site_bases_tsv: "TSV with one row per variant class and size bin and, for each entry of `vcf_names`, the number of distinct reference bases altered across the callset and its proportion of the contig length."
        sample_bases_tsv: "TSV laid out as `site_bases_tsv` holding the mean number of bases altered per sample, counting each altered allele, and its proportion of the contig length."
    }

    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] vcf_names
        String contig
        String prefix

        String subset_vcf_string = ""

        File ref_fai

        Int shard_bin_size = 5000000
        String utils_docker

        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_count_shard
        RuntimeAttr? runtime_attr_merge_counts
    }

    call Helpers.CreateContigShards {
        input:
            vcfs = vcfs,
            vcf_idxs = vcf_idxs,
            contig = contig,
            shard_bin_size = shard_bin_size,
            prefix = "~{prefix}.shards",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_create_shards
    }

    scatter (j in range(length(CreateContigShards.shard_regions))) {
        Boolean is_last_shard = j == length(CreateContigShards.shard_regions) - 1

        scatter (i in range(length(vcfs))) {
            call CountShardBases {
                input:
                    vcf = vcfs[i],
                    vcf_idx = vcf_idxs[i],
                    vcf_name = vcf_names[i],
                    region = CreateContigShards.shard_regions[j],
                    is_last_shard = is_last_shard,
                    prefix = "~{prefix}.shard_~{j}.callset_~{i}",
                    subset_vcf_string = subset_vcf_string,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_count_shard
            }
        }
    }

    call MergeVariantBaseCounts {
        input:
            counts_tsvs = flatten(CountShardBases.counts_tsv),
            vcf_names = vcf_names,
            contig = contig,
            ref_fai = ref_fai,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_counts
    }

    output {
        File site_bases_tsv = MergeVariantBaseCounts.site_bases_tsv
        File sample_bases_tsv = MergeVariantBaseCounts.sample_bases_tsv
    }
}

task CountShardBases {
    input {
        File vcf
        File vcf_idx
        String vcf_name
        String region
        Boolean is_last_shard
        String prefix
        String subset_vcf_string
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        vcf: { localization_optional: true }
        vcf_idx: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        n_samples=$(bcftools view -h ~{subset_vcf_string} ~{vcf} | bcftools query -l - | wc -l)

        # Omit -t so deletions starting in an earlier region still reach the bases they cover in this one
        bcftools view -r ~{region} ~{subset_vcf_string} ~{vcf} -Ou \
            | bcftools query -f '%POS\t%INFO/allele_type\t%INFO/allele_length\t%INFO/AC\n' - \
            > records.tsv

        python3 <<CODE
REGION = "~{region}"
IS_LAST_SHARD = "~{is_last_shard}".lower() == "true"
VCF_NAME = "~{vcf_name}"
N_SAMPLES = int("$n_samples")

CATEGORIES = ["SNV"]
for _type in ["DEL", "INS"]:
    CATEGORIES += [f"{_type} 1-49", f"{_type} 50-499", f"{_type} 500+"]


def size_label(length):
    if length < 50:
        return "1-49"
    if length < 500:
        return "50-499"
    return "500+"


region_start, region_end = (int(value) for value in REGION.rsplit(":", 1)[1].split("-"))

site_bases = {category: 0 for category in CATEGORIES}
allele_bases = {category: 0 for category in CATEGORIES}
altered_spans = {category: [] for category in CATEGORIES}

with open("records.tsv") as handle:
    for line in handle:
        pos, allele_type, allele_length, ac = line.rstrip("\n").split("\t")
        pos = int(pos)
        allele_type = allele_type.lower()
        if allele_type == "snv":
            category, length = "SNV", 1
        elif allele_type == "del" or allele_type == "ins" or "dup" in allele_type:
            length = abs(int(allele_length))
            category = f"{'DEL' if allele_type == 'del' else 'INS'} {size_label(length)}"
        else:
            continue

        alt_count = sum(int(value) for value in ac.split(",") if value != ".")
        if alt_count == 0:
            continue

        # Count a record's AC only in the region holding its position
        if region_start <= pos <= region_end:
            allele_bases[category] += alt_count * length

        # Keep the reference bases an SNV or deletion alters, clipped to the region, as insertions alter none
        if allele_type not in ("snv", "del"):
            continue
        span_start = pos if allele_type == "snv" else pos + 1
        span_end = span_start + length - 1
        span_start = max(span_start, region_start)
        if not IS_LAST_SHARD:
            span_end = min(span_end, region_end)
        if span_end >= span_start:
            altered_spans[category].append((span_start, span_end))

# Count each altered reference base once per bin, however many records alter it
for category, spans in altered_spans.items():
    covered_end = 0
    for span_start, span_end in sorted(spans):
        first_new = max(span_start, covered_end + 1)
        if span_end >= first_new:
            site_bases[category] += span_end - first_new + 1
            covered_end = span_end

with open("~{prefix}.counts.tsv", "w") as out:
    for category in CATEGORIES:
        out.write(f"{VCF_NAME}\t{category}\t{site_bases[category]}\t{allele_bases[category]}\t{N_SAMPLES}\n")
CODE
    >>>

    output {
        File counts_tsv = "~{prefix}.counts.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 20,
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

task MergeVariantBaseCounts {
    input {
        Array[File] counts_tsvs
        Array[String] vcf_names
        String contig
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{write_lines(counts_tsvs)} | xargs cat > counts.tsv

        contig_length=$(awk -v contig="~{contig}" '$1 == contig { print $2 }' ~{ref_fai})

        python3 <<CODE
from collections import defaultdict

NAMES = "~{sep='\t' vcf_names}".split("\t")
CONTIG_LENGTH = int("$contig_length")

CATEGORIES = ["SNV"]
for _type in ["DEL", "INS"]:
    CATEGORIES += [f"{_type} 1-49", f"{_type} 50-499", f"{_type} 500+"]

site_bases = defaultdict(int)
allele_bases = defaultdict(int)
n_samples = {}

with open("counts.tsv") as handle:
    for line in handle:
        vcf_name, category, sites, alleles, samples = line.rstrip("\n").split("\t")
        site_bases[(vcf_name, category)] += int(sites)
        allele_bases[(vcf_name, category)] += int(alleles)
        n_samples[vcf_name] = int(samples)

with open("~{prefix}.site_bases.tsv", "w") as sites_out, open("~{prefix}.sample_bases.tsv", "w") as samples_out:
    header = ["category"]
    for name in NAMES:
        header += [f"{name}_bases", f"{name}_proportion"]
    sites_out.write("\t".join(header) + "\n")
    samples_out.write("\t".join(header) + "\n")
    for category in CATEGORIES:
        site_values = []
        sample_values = []
        for name in NAMES:
            sites = site_bases[(name, category)]
            samples = allele_bases[(name, category)] / n_samples[name]
            site_values += [f"{sites:.6f}", f"{sites / CONTIG_LENGTH:.6f}"]
            sample_values += [f"{samples:.6f}", f"{samples / CONTIG_LENGTH:.6f}"]
        sites_out.write("\t".join([category] + site_values) + "\n")
        samples_out.write("\t".join([category] + sample_values) + "\n")
CODE
    >>>

    output {
        File site_bases_tsv = "~{prefix}.site_bases.tsv"
        File sample_bases_tsv = "~{prefix}.sample_bases.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(counts_tsvs, "GB")) + 10,
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
