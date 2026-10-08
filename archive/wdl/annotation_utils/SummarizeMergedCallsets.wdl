version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SummarizeMergedCallsets {
    meta {
        description: [
            "This utility counts the sites of a merged VCF from `MergeVcfs` by whether each site is supported by more than one callset, at the site level and per sample. Each sample is assigned to a callset by `sample_sources_tsv`; a site is matched when at least two callsets have a carrier among their own samples and unique to a callset when only that callset does, so a merged record whose other callset carries only reference or missing genotypes counts as unique. Sites with no carrier are skipped.",
            "Counts are binned by allele class and size from INFO/allele_type and INFO/allele_length as in `SummarizeAnnotations`, split by INFO/REGION, and repeated for sites lacking dbSNP_ID and dbGaP_ID, lacking gnomAD_V4_match_ID, or lacking both. The site table reports how many sites are matched and how many are unique to each callset; the per-sample table reports how many matched and unique sites each sample carries. The task fails if a callset has a carrier at a site whose SOURCE_NAMES does not list it, which means the input was not produced by `MergeVcfs`."
        ]
    }

    parameter_meta {
        merged_vcf: "Merged VCF from `MergeVcfs`, such as one contig. Pass only the merged VCF, not the callset VCFs merged into it."
        merged_vcf_idx: "Index for `merged_vcf`."
        sample_sources_tsv: "Two-column TSV without a header giving each sample ID and the SOURCE_NAMES callset name it belongs to. Every sample in `merged_vcf` must be listed."
        subset_vcf_string: "`bcftools view` arguments used to pre-subset `merged_vcf`."
        length_bins: "Size-bin edges used for the DEL and INS columns."
        records_per_shard: "Number of variants to keep within a single shard."
        site_counts_tsv: "Number of matched sites and of sites unique to each callset, per concordance group and region."
        sample_counts_tsv: "Per-sample counts of matched and unique sites per concordance group and region."
    }

    input {
        File merged_vcf
        File merged_vcf_idx
        String prefix

        File sample_sources_tsv
        String subset_vcf_string = ""
        Array[Int] length_bins = [0, 1, 50, 500]

        Int? records_per_shard
        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_count
        RuntimeAttr? runtime_attr_merge
    }

    if (defined(records_per_shard)) {
        call Helpers.ShardVcfByRecords {
            input:
                vcf = merged_vcf,
                vcf_idx = merged_vcf_idx,
                records_per_shard = select_first([records_per_shard]),
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_shard
        }
    }

    Array[File] shard_vcfs = select_first([ShardVcfByRecords.shards, [merged_vcf]])
    Array[File] shard_vcf_idxs = select_first([ShardVcfByRecords.shard_idxs, [merged_vcf_idx]])

    scatter (j in range(length(shard_vcfs))) {
        call Helpers.SubsetVcfByArgs {
            input:
                vcf = shard_vcfs[j],
                vcf_idx = shard_vcf_idxs[j],
                extra_args = subset_vcf_string,
                prefix = "~{prefix}.shard_~{j}.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call CountMergedCallsetShard {
            input:
                vcf = SubsetVcfByArgs.subset_vcf,
                vcf_idx = SubsetVcfByArgs.subset_vcf_idx,
                sample_sources_tsv = sample_sources_tsv,
                length_bins = length_bins,
                prefix = "~{prefix}.shard_~{j}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_count
        }
    }

    call MergeCallsetCountTables {
        input:
            site_tsvs = CountMergedCallsetShard.site_counts_tsv,
            sample_tsvs = CountMergedCallsetShard.sample_counts_tsv,
            sample_sources_tsv = sample_sources_tsv,
            length_bins = length_bins,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File site_counts_tsv = MergeCallsetCountTables.site_counts_tsv
        File sample_counts_tsv = MergeCallsetCountTables.sample_counts_tsv
    }
}

task CountMergedCallsetShard {
    input {
        File vcf
        File vcf_idx
        File sample_sources_tsv
        Array[Int] length_bins
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import subprocess
import numpy as np
import pandas as pd
import pysam

VCF_PATH = "~{vcf}"
SOURCES_PATH = "~{sample_sources_tsv}"
SITE_OUTPUT = "~{prefix}.sites.raw.tsv"
SAMPLE_OUTPUT = "~{prefix}.samples.raw.tsv"
LENGTH_BINS = [~{sep=", " length_bins}]
SIZE_LABELS = [f"{s}-{e - 1}" for s, e in zip(LENGTH_BINS, LENGTH_BINS[1:])] + [f"{LENGTH_BINS[-1]}+"]
REGION_ORDER = ["US", "RM", "SD", "SR"]
INFO_TAGS = ["allele_type", "allele_length", "REGION", "dbSNP_ID", "dbGaP_ID", "gnomAD_V4_match_ID", "SOURCE_NAMES"]
KEY_COLUMNS = ["dbsnp_missing", "gnomad_missing", "region", "status", "bucket"]
CHUNK_SIZE = 100000


def get_size_bucket(allele_length):
    size = abs(allele_length)
    for index, start in enumerate(LENGTH_BINS):
        if index + 1 == len(LENGTH_BINS) or size < LENGTH_BINS[index + 1]:
            return SIZE_LABELS[index]


def determine_column(allele_type, allele_length):
    allele_type = allele_type.lower()
    if allele_type == "snv": return "SNV"
    if allele_type == "trv": return "TRV"
    if allele_length == ".": return "Other"
    size_label = get_size_bucket(int(float(allele_length)))
    if allele_type in ("ins", "dup"): return f"INS {size_label}"
    if allele_type == "del": return f"DEL {size_label}"
    return "Other"


# Map every VCF sample to its callset, failing on a sample the TSV does not list
sources = pd.read_csv(SOURCES_PATH, sep="\t", header=None, names=["sample", "callset"], dtype=str)
callset_of = dict(zip(sources["sample"], sources["callset"]))
callsets = list(dict.fromkeys(sources["callset"]))
with pysam.VariantFile(VCF_PATH) as vcf_in:
    samples = list(vcf_in.header.samples)
    header_tags = set(vcf_in.header.info.keys())
missing = [s for s in samples if s not in callset_of]
if missing:
    raise SystemExit(f"Samples absent from sample_sources_tsv: {missing[:5]}")
sample_callset = np.array([callsets.index(callset_of[s]) for s in samples])
status_labels = np.array(["matched"] + [f"unique_{c}" for c in callsets])

# Stream only the carrier samples of each site, so bcftools drops sites with no carrier
query_format = "\t".join(f"%INFO/{t}" if t in header_tags else "." for t in INFO_TAGS) + "\t[%SAMPLE,]\n"
proc = subprocess.Popen(["bcftools", "query", "-i", 'GT="alt"', "-f", query_format, VCF_PATH], stdout=subprocess.PIPE)
chunks = pd.read_csv(
    proc.stdout, sep="\t", header=None, names=INFO_TAGS + ["carriers"],
    dtype=str, na_filter=False, quoting=3, chunksize=CHUNK_SIZE,
)

site_tables = []
sample_counts = {}
for chunk in chunks:
    if chunk.empty:
        continue
    n_records = len(chunk)
    carriers = chunk["carriers"].str.rstrip(",").str.split(",")
    record_idx = np.repeat(np.arange(n_records), carriers.str.len().to_numpy())
    sample_idx = pd.Categorical(np.concatenate(carriers.to_numpy()), categories=samples).codes

    # Call a site matched when carriers come from at least two callsets
    in_callset = np.zeros((n_records, len(callsets)), dtype=bool)
    in_callset[record_idx, sample_callset[sample_idx]] = True
    n_supporting = in_callset.sum(axis=1)
    status = np.where(n_supporting >= 2, 0, 1 + in_callset.argmax(axis=1))
    source_names = chunk["SOURCE_NAMES"].str.split(",", expand=True).to_numpy()
    listed = np.stack([(source_names == c).any(axis=1) for c in callsets], axis=1)
    if (in_callset & ~listed).any():
        raise SystemExit("A callset carries a site that SOURCE_NAMES does not list it in; pass only MergeVcfs outputs")

    type_lengths = chunk["allele_type"] + "|" + chunk["allele_length"]
    type_codes, type_uniques = pd.factorize(type_lengths)
    buckets = np.array([determine_column(*u.split("|")) for u in type_uniques])[type_codes]
    keys = pd.DataFrame({
        "dbsnp_missing": ((chunk["dbSNP_ID"] == ".") & (chunk["dbGaP_ID"] == ".")).astype(int).to_numpy(),
        "gnomad_missing": (chunk["gnomAD_V4_match_ID"] == ".").astype(int).to_numpy(),
        "region": chunk["REGION"].where(chunk["REGION"].isin(REGION_ORDER), "").to_numpy(),
        "status": status_labels[status],
        "bucket": buckets,
    })

    site_tables.append(keys.groupby(KEY_COLUMNS, as_index=False).size().rename(columns={"size": "sites"}))

    # Count each carrier sample once per site under a matched or unique status
    sample_keys = keys.assign(status=np.where(status == 0, "matched", "unique"))
    key_codes, key_uniques = pd.MultiIndex.from_frame(sample_keys).factorize()
    sample_sums = np.bincount(
        key_codes[record_idx] * len(samples) + sample_idx, minlength=len(key_uniques) * len(samples)
    ).reshape(len(key_uniques), len(samples))
    for key, sums in zip(key_uniques, sample_sums):
        sample_counts[key] = sample_counts.get(key, 0) + sums

proc.stdout.close()
if proc.wait() != 0:
    raise SystemExit("bcftools query failed")

# Write long-format shard tables that the merge task sums and pivots
if site_tables:
    site_table = pd.concat(site_tables).groupby(KEY_COLUMNS, as_index=False).sum()
else:
    site_table = pd.DataFrame(columns=KEY_COLUMNS + ["sites"])
site_table.to_csv(SITE_OUTPUT, sep="\t", index=False)

with open(SAMPLE_OUTPUT, "w") as handle:
    handle.write("\t".join(["sample", "callset"] + KEY_COLUMNS + ["count"]) + "\n")
    for key, sums in sample_counts.items():
        for index in np.flatnonzero(sums):
            sample = samples[index]
            handle.write("\t".join([sample, callset_of[sample]] + [str(k) for k in key] + [str(sums[index])]) + "\n")
CODE
    >>>

    output {
        File site_counts_tsv = "~{prefix}.sites.raw.tsv"
        File sample_counts_tsv = "~{prefix}.samples.raw.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 2 * ceil(size([vcf, vcf_idx], "GB")) + 10,
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

task MergeCallsetCountTables {
    input {
        Array[File] site_tsvs
        Array[File] sample_tsvs
        File sample_sources_tsv
        Array[Int] length_bins
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import pandas as pd

SITE_TSVS = [p for p in "~{sep=',' site_tsvs}".split(",") if p]
SAMPLE_TSVS = [p for p in "~{sep=',' sample_tsvs}".split(",") if p]
SOURCES_PATH = "~{sample_sources_tsv}"
SITE_OUTPUT = "~{prefix}.site_counts.tsv"
SAMPLE_OUTPUT = "~{prefix}.sample_counts.tsv"
LENGTH_BINS = [~{sep=", " length_bins}]
SIZE_LABELS = [f"{s}-{e - 1}" for s, e in zip(LENGTH_BINS, LENGTH_BINS[1:])] + [f"{LENGTH_BINS[-1]}+"]
COLUMN_BUCKETS = ["SNV"] + [f"DEL {_l}" for _l in SIZE_LABELS] + [f"INS {_l}" for _l in SIZE_LABELS] + ["TRV", "Other"]
CONCORDANCE_GROUPS = [
    ("All", None),
    ("dbSNP Missing", "dbsnp_missing"),
    ("gnomAD Missing", "gnomad_missing"),
    ("dbSNP/gnomAD Missing", "both_missing"),
]
REGION_ORDER = ["All", "US", "RM", "SD", "SR"]
RAW_KEYS = ["dbsnp_missing", "gnomad_missing", "region", "status", "bucket"]


def read_tables(paths, keys):
    frames = [pd.read_csv(p, sep="\t", dtype={"sample": str, "region": str}, keep_default_na=False) for p in paths]
    table = pd.concat([f for f in frames if not f.empty], ignore_index=True)
    table = table.groupby(keys, as_index=False).sum(numeric_only=True)
    table["both_missing"] = ((table["dbsnp_missing"] == 1) & (table["gnomad_missing"] == 1)).astype(int)
    return table


def expand_groups(table, group_keys, value_columns):
    frames = []
    for concordance, flag in CONCORDANCE_GROUPS:
        concordance_table = table if flag is None else table[table[flag] == 1]
        for region in REGION_ORDER:
            region_table = concordance_table if region == "All" else concordance_table[concordance_table["region"] == region]
            summed = region_table.groupby(group_keys, as_index=False)[value_columns].sum()
            summed.insert(0, "region", region)
            summed.insert(0, "concordance", concordance)
            frames.append(summed)
    return pd.concat(frames, ignore_index=True)


def pivot_buckets(table, index_columns, value_column, orders):
    wide = table.pivot_table(index=index_columns, columns="bucket", values=value_column, aggfunc="sum", fill_value=0)
    wide = wide.reindex(columns=COLUMN_BUCKETS, fill_value=0).reset_index()
    for column, order in orders.items():
        wide[column] = pd.Categorical(wide[column], order)
    return wide.sort_values(list(orders)).rename_axis(columns=None)


callsets = list(dict.fromkeys(pd.read_csv(SOURCES_PATH, sep="\t", header=None, dtype=str)[1]))
concordance_order = [g for g, _ in CONCORDANCE_GROUPS]

# Sum the site shards in each concordance and region group
site = expand_groups(read_tables(SITE_TSVS, RAW_KEYS), ["status", "bucket"], ["sites"])
site = pivot_buckets(site, ["concordance", "region", "status"], "sites", {
    "concordance": concordance_order,
    "region": REGION_ORDER,
    "status": [f"unique_{c}" for c in callsets] + ["matched"],
})
site.to_csv(SITE_OUTPUT, sep="\t", index=False)

# Sum the per-sample shards in the same concordance and region groups
sample_raw = read_tables(SAMPLE_TSVS, ["sample", "callset"] + RAW_KEYS)
sample = expand_groups(sample_raw, ["sample", "callset", "status", "bucket"], ["count"])
sample = pivot_buckets(sample, ["sample", "callset", "concordance", "region", "status"], "count", {
    "callset": callsets,
    "sample": list(dict.fromkeys(sample_raw["sample"])),
    "concordance": concordance_order,
    "region": REGION_ORDER,
    "status": ["matched", "unique"],
})
sample.to_csv(SAMPLE_OUTPUT, sep="\t", index=False)
CODE
    >>>

    output {
        File site_counts_tsv = "~{prefix}.site_counts.tsv"
        File sample_counts_tsv = "~{prefix}.sample_counts.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(site_tsvs, "GB") + size(sample_tsvs, "GB")) + 10,
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
