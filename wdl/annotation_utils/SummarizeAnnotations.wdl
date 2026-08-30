version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SummarizeAnnotations {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix

        Array[Int] length_bins_summary = [0, 1, 50, 500]
        Array[Int] length_bins_plotting = [0, 1, 50, 100, 500, 5000, 50000]
        Array[Float] af_bins_plotting = [0.0, 0.01, 0.05, 0.1, 0.5]

        Boolean create_per_sample = false
        Boolean create_per_allele = false
        Boolean create_list = false
        Boolean create_functional = false
        Boolean create_plotting = false

        Boolean use_ssd = false
        Boolean split_by_region = false
        String subset_vcf_string = ""
        Int max_length = -1
        Int min_length = -1

        Int? records_per_shard

        File? ped

        String utils_docker

        RuntimeAttr? runtime_attr_strip
        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_count
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_convert
        RuntimeAttr? runtime_attr_find_trios
    }

    Boolean sites_only = !(create_per_sample || create_per_allele || create_plotting)
    Boolean run_denovo = create_plotting && defined(ped)

    if (run_denovo) {
        call Helpers.FindTrios {
            input:
                vcf = vcfs[0],
                vcf_idx = vcf_idxs[0],
                ped = select_first([ped]),
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_find_trios
        }
    }

    scatter (i in range(length(vcfs))) {
        if (sites_only) {
            call Helpers.StripGenotypes {
                input:
                    vcf = vcfs[i],
                    vcf_idx = vcf_idxs[i],
                    prefix = "~{prefix}.input_~{i}.stripped",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_strip
            }
        }

        File vcf_for_shard = select_first([StripGenotypes.stripped_vcf, vcfs[i]])
        File vcf_idx_for_shard = select_first([StripGenotypes.stripped_vcf_idx, vcf_idxs[i]])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = vcf_for_shard,
                    vcf_idx = vcf_idx_for_shard,
                    records_per_shard = select_first([records_per_shard]),
                    use_ssd = use_ssd,
                    prefix = "~{prefix}.input_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] shard_vcfs = select_first([ShardVcfByRecords.shards, [vcf_for_shard]])
        Array[File] shard_vcf_idxs = select_first([ShardVcfByRecords.shard_idxs, [vcf_idx_for_shard]])

        scatter (j in range(length(shard_vcfs))) {
            call Helpers.SubsetVcfByArgs {
                input:
                    vcf = shard_vcfs[j],
                    vcf_idx = shard_vcf_idxs[j],
                    extra_args = subset_vcf_string,
                    prefix = "~{prefix}.input_~{i}.shard_~{j}.subset",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset
            }

            call CountAnnotationShard {
                input:
                    vcf = SubsetVcfByArgs.subset_vcf,
                    vcf_idx = SubsetVcfByArgs.subset_vcf_idx,
                    create_per_sample = create_per_sample,
                    create_per_allele = create_per_allele,
                    create_list = create_list,
                    create_plotting = create_plotting,
                    split_by_region = split_by_region,
                    trio_definitions = FindTrios.trio_definitions,
                    length_bins_summary = length_bins_summary,
                    length_bins_plotting = length_bins_plotting,
                    af_bins_plotting = af_bins_plotting,
                    max_length = max_length,
                    min_length = min_length,
                    prefix = "~{prefix}.input_~{i}.shard_~{j}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_count
            }
        }
    }

    Array[File] site_count_tables = flatten(CountAnnotationShard.site_counts_tsv)
    Array[File] variant_list_tables = flatten(CountAnnotationShard.plotting_variant_list_tsv)
    Array[File] sample_count_tables = flatten(CountAnnotationShard.sample_counts_tsv)
    Array[File] allele_count_tables = flatten(CountAnnotationShard.allele_counts_tsv)
    Array[File] list_tables = flatten(CountAnnotationShard.list_tsv)
    Array[File] sample_count_files = flatten(CountAnnotationShard.sample_count_file)
    Array[File] functional_site_count_tables = flatten(CountAnnotationShard.gene_counts_tsv)
    Array[File] functional_sample_count_tables = flatten(CountAnnotationShard.gene_sample_counts_tsv)
    Array[File] functional_allele_count_tables = flatten(CountAnnotationShard.gene_allele_counts_tsv)

    call MergeAnnotationCountTables as MergeSiteCounts {
        input:
            count_tsvs = site_count_tables,
            sample_count_files = sample_count_files,
            normalization_mode = "sites",
            split_by_region = split_by_region,
            drop_empty_types = true,
            prefix = "~{prefix}.summary_sites",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge
    }

    if (create_per_sample) {
        call MergeAnnotationCountTables as MergeSampleCounts {
            input:
                count_tsvs = sample_count_tables,
                sample_count_files = sample_count_files,
                normalization_mode = "samples",
                split_by_region = split_by_region,
                drop_empty_types = true,
                prefix = "~{prefix}.summary_samples",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    if (create_per_allele) {
        call MergeAnnotationCountTables as MergeAlleleCounts {
            input:
                count_tsvs = allele_count_tables,
                sample_count_files = sample_count_files,
                normalization_mode = "alleles",
                split_by_region = split_by_region,
                drop_empty_types = true,
                prefix = "~{prefix}.summary_alleles",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    if (create_list) {
        call MergeAnnotationListTables {
            input:
                list_tsvs = list_tables,
                prefix = "~{prefix}.summary_list",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    if (create_functional) {
        call MergeGeneCountTables as MergeGeneCounts {
            input:
                count_tsvs = functional_site_count_tables,
                prefix = "~{prefix}.summary_functional",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    if (create_per_sample && create_functional) {
        call MergeGeneCountTables as MergeGeneSampleCounts {
            input:
                count_tsvs = functional_sample_count_tables,
                prefix = "~{prefix}.summary_functional_samples",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    if (create_per_allele && create_functional) {
        call MergeGeneCountTables as MergeGeneAlleleCounts {
            input:
                count_tsvs = functional_allele_count_tables,
                prefix = "~{prefix}.summary_functional_alleles",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    if (create_plotting) {
        call MergeAnnotationCountTables as MergePlottingSites {
            input:
                count_tsvs = flatten(CountAnnotationShard.plotting_site_tsv),
                sample_count_files = sample_count_files,
                normalization_mode = "sites",
                split_by_region = split_by_region,
                prefix = "~{prefix}.plotting_sites",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }

        call Helpers.ConvertTsvToParquet as ConvertPlottingSites {
            input:
                tsv = MergePlottingSites.merged_counts_tsv,
                prefix = "~{prefix}.plotting_sites",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_convert
        }

        call MergeSampleSpecificTables as MergePlottingSamples {
            input:
                count_tsvs = flatten(CountAnnotationShard.plotting_sample_tsv),
                prefix = "~{prefix}.plotting_samples",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }

        call Helpers.ConvertTsvToParquet as ConvertPlottingSamples {
            input:
                tsv = MergePlottingSamples.merged_tsv,
                prefix = "~{prefix}.plotting_samples",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_convert
        }

        call MergeSampleSpecificTables as MergePlottingAlleles {
            input:
                count_tsvs = flatten(CountAnnotationShard.plotting_allele_tsv),
                prefix = "~{prefix}.plotting_alleles",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }

        call Helpers.ConvertTsvToParquet as ConvertPlottingAlleles {
            input:
                tsv = MergePlottingAlleles.merged_tsv,
                prefix = "~{prefix}.plotting_alleles",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_convert
        }

        if (run_denovo) {
            call MergeSampleSpecificTables as MergePlottingDenovo {
                input:
                    count_tsvs = select_all(flatten(CountAnnotationShard.plotting_denovo_tsv)),
                    prefix = "~{prefix}.plotting_denovo",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge
            }

            call Helpers.ConvertTsvToParquet as ConvertPlottingDenovo {
                input:
                    tsv = MergePlottingDenovo.merged_tsv,
                    prefix = "~{prefix}.plotting_denovo",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_convert
            }
        }

        call MergeVariantListTables {
            input:
                variant_list_tsvs = variant_list_tables,
                prefix = "~{prefix}.plotting_variant_list",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }

        call Helpers.ConvertTsvToParquet as ConvertPlottingVariantList {
            input:
                tsv = MergeVariantListTables.merged_tsv,
                prefix = "~{prefix}.plotting_variant_list",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_convert
        }
    }

    output {
        File summary_sites_tsv = MergeSiteCounts.merged_counts_tsv
        File? summary_samples_tsv = MergeSampleCounts.merged_counts_tsv
        File? summary_alleles_tsv = MergeAlleleCounts.merged_counts_tsv
        File? summary_list_tsv = MergeAnnotationListTables.merged_list_tsv
        File? summary_functional_tsv = MergeGeneCounts.merged_counts_tsv
        File? summary_functional_samples_tsv = MergeGeneSampleCounts.merged_counts_tsv
        File? summary_functional_alleles_tsv = MergeGeneAlleleCounts.merged_counts_tsv

        File? plotting_sites_parquet = ConvertPlottingSites.parquet
        File? plotting_samples_parquet = ConvertPlottingSamples.parquet
        File? plotting_alleles_parquet = ConvertPlottingAlleles.parquet
        File? plotting_denovo_parquet = ConvertPlottingDenovo.parquet
        File? plotting_variant_list_parquet = ConvertPlottingVariantList.parquet
    }
}

task CountAnnotationShard {
    input {
        File vcf
        File vcf_idx
        Boolean create_per_sample
        Boolean create_per_allele
        Boolean create_list
        Boolean create_plotting
        Boolean split_by_region
        File? trio_definitions
        Array[Int] length_bins_summary
        Array[Int] length_bins_plotting
        Array[Float] af_bins_plotting
        Int max_length
        Int min_length
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
# Imports
import csv
import re
from collections import defaultdict
import pysam


# Constants
VCF_PATH = "~{vcf}"
SITE_OUTPUT = "~{prefix}.sites.raw.tsv"
SAMPLE_OUTPUT = "~{prefix}.samples.raw.tsv"
ALLELE_OUTPUT = "~{prefix}.alleles.raw.tsv"
LIST_OUTPUT = "~{prefix}.list.raw.tsv"
GENE_OUTPUT = "~{prefix}.genes.raw.tsv"
GENE_SAMPLE_OUTPUT = "~{prefix}.genes.samples.raw.tsv"
GENE_ALLELE_OUTPUT = "~{prefix}.genes.alleles.raw.tsv"
SAMPLE_COUNT_OUTPUT = "~{prefix}.sample_count.txt"
PLOT_SITE_OUTPUT = "~{prefix}.plotting_sites.raw.tsv"
PLOT_SAMPLE_OUTPUT = "~{prefix}.plotting_samples.raw.tsv"
PLOT_ALLELE_OUTPUT = "~{prefix}.plotting_alleles.raw.tsv"
DENOVO_OUTPUT = "~{prefix}.plotting_denovo.raw.tsv"

CREATE_PER_SAMPLE = "~{create_per_sample}".lower() == "true"
CREATE_PER_ALLELE = "~{create_per_allele}".lower() == "true"
CREATE_LIST = "~{create_list}".lower() == "true"
CREATE_PLOTTING = "~{create_plotting}".lower() == "true"
SPLIT_BY_REGION = "~{split_by_region}".lower() == "true"
VARIANT_LIST_OUTPUT = "~{prefix}.plotting_variant_list.raw.tsv"
TRIO_DEF_PATH = "~{trio_definitions}"
MAX_LENGTH = ~{max_length}
MIN_LENGTH = ~{min_length}
MAX_LENGTH = MAX_LENGTH if MAX_LENGTH > 0 else None
MIN_LENGTH = MIN_LENGTH if MIN_LENGTH > 0 else None

LENGTH_BINS_SUMMARY = [~{sep=", " length_bins_summary}]
SIZE_LABELS_SUMMARY = [f"{s}-{e - 1}" for s, e in zip(LENGTH_BINS_SUMMARY, LENGTH_BINS_SUMMARY[1:])] + [f"{LENGTH_BINS_SUMMARY[-1]}+"]

COLUMN_BUCKETS_SUMMARY = ["SNV"]
for _label in SIZE_LABELS_SUMMARY:
    COLUMN_BUCKETS_SUMMARY.append(f"DEL {_label}")
for _label in SIZE_LABELS_SUMMARY:
    COLUMN_BUCKETS_SUMMARY.append(f"INS {_label}")
COLUMN_BUCKETS_SUMMARY += ["TRV", "Other"]

LENGTH_BINS_PLOTTING = [~{sep=", " length_bins_plotting}]
SIZE_LABELS_PLOTTING = [f"{s}-{e - 1}" for s, e in zip(LENGTH_BINS_PLOTTING, LENGTH_BINS_PLOTTING[1:])] + [f"{LENGTH_BINS_PLOTTING[-1]}+"]

COLUMN_BUCKETS_PLOTTING = ["SNV"]
for _label in SIZE_LABELS_PLOTTING:
    COLUMN_BUCKETS_PLOTTING.append(f"DEL {_label}")
for _label in SIZE_LABELS_PLOTTING:
    COLUMN_BUCKETS_PLOTTING.append(f"INS {_label}")
COLUMN_BUCKETS_PLOTTING += ["TRV", "Other"]
COLUMN_BUCKETS_PLOTTING_MATCHED = [f"{_bucket} (Matched)" for _bucket in COLUMN_BUCKETS_PLOTTING]

AF_BINS_PLOTTING = [~{sep=", " af_bins_plotting}]
if len(AF_BINS_PLOTTING) > 1:
    AF_LABELS_PLOTTING = [f"{AF_BINS_PLOTTING[i]}-{AF_BINS_PLOTTING[i + 1]}" for i in range(len(AF_BINS_PLOTTING) - 1)]
    AF_LABELS_PLOTTING.append(f">{AF_BINS_PLOTTING[-1]}")
else:
    AF_LABELS_PLOTTING = []

COLUMN_BUCKETS_PLOTTING_AF = []
for _type in ["SNV", "DEL", "INS", "TRV", "Other"]:
    for _af_label in AF_LABELS_PLOTTING:
        COLUMN_BUCKETS_PLOTTING_AF.append(f"{_type} [{_af_label}]")

def get_size_bucket(allele_length, length_bins, size_labels):
    size = abs(allele_length)
    for index, start in enumerate(length_bins):
        if index + 1 == len(length_bins) or size < length_bins[index + 1]:
            return size_labels[index]

INTERNAL_TOTAL_LABEL = "Total"
DISPLAY_ALL_LABEL = "All"
REGION_ORDER = ["US", "RM", "SD", "SR"]

CONSEQUENCE_PRIORITY = [
    ("LoF", {
        "transcript_ablation", "stop_gained", "frameshift_variant",
        "splice_donor_variant", "splice_acceptor_variant", "stop_lost",
        "start_lost", "transcript_amplification", "feature_elongation",
        "feature_truncation",
    }),
    ("Missense", {
        "missense_variant", "inframe_insertion", "inframe_deletion",
        "protein_altering_variant",
    }),
    ("Synonymous", {
        "synonymous_variant", "stop_retained_variant", "start_retained_variant",
        "incomplete_terminal_codon_variant", "coding_sequence_variant",
    }),
    ("Intronic", {
        "intron_variant", "splice_region_variant", "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant", "splice_donor_5th_base_variant",
        "NMD_transcript_variant", "non_coding_transcript_exon_variant",
        "non_coding_transcript_variant",
    }),
    ("Intergenic", {
        "intergenic_variant", "upstream_gene_variant", "downstream_gene_variant",
        "regulatory_region_variant", "regulatory_region_ablation",
        "regulatory_region_amplification", "TF_binding_site_variant",
        "TFBS_ablation", "TFBS_amplification", "3_prime_UTR_variant",
        "5_prime_UTR_variant",
    }),
]

SVANNOTATE_PRIORITY = [
    ("LoF", "PREDICTED_LOF"),
    ("Copy Gain", "PREDICTED_COPY_GAIN"),
    ("IED", "PREDICTED_INTRAGENIC_EXON_DUP"),
    ("PED", "PREDICTED_PARTIAL_EXON_DUP"),
    ("TSSD", "PREDICTED_TSS_DUP"),
    ("DP", "PREDICTED_DUP_PARTIAL"),
    ("UTR", "PREDICTED_UTR"),
]

VEP_LABELS = [label for label, _ in CONSEQUENCE_PRIORITY] + ["Other"]
SVANNOTATE_LABELS = [label for label, _ in SVANNOTATE_PRIORITY]

ROW_ORDER = [
    ("", INTERNAL_TOTAL_LABEL),
    ("", "ME"),
    ("ME", "ALU"),
    ("ME", "LINE"),
    ("ME", "SVA"),
    ("", "Duplication"),
    ("Duplication", "Filtered Tandem"),
    ("Duplication", "SVAN Tandem"),
    ("Duplication", "Interspersed"),
    ("Duplication", "Complex"),
    ("Duplication", "Inverted"),
    ("", "TR Parsed"),
    ("", "NUMT"),
    ("", "VEP"),
]
ROW_ORDER += [("VEP", label) for label in VEP_LABELS]
ROW_ORDER.append(("", "SVAnnotate"))
ROW_ORDER += [("SVAnnotate", label) for label in SVANNOTATE_LABELS]
ROW_ORDER.append(("", "Concordance"))
for missing_label in ["dbSNP Missing", "gnomAD Missing", "dbSNP/gnomAD Missing"]:
    ROW_ORDER.append(("Concordance", missing_label))
    ROW_ORDER += [("Concordance", f"{missing_label} + {label} (VEP)") for label in VEP_LABELS]
    ROW_ORDER += [("Concordance", f"{missing_label} + {label} (SVAnnotate)") for label in SVANNOTATE_LABELS]

PARENT_ROWS = {"Concordance", "ME", "Duplication", "VEP", "SVAnnotate"}

gene_counts = defaultdict(lambda: defaultdict(int))
gene_sample_counts = defaultdict(lambda: defaultdict(int))
gene_allele_counts = defaultdict(lambda: defaultdict(int))
PREDICTED_FIELDS = [
    "PREDICTED_LOF",
    "PREDICTED_COPY_GAIN",
    "PREDICTED_INTRAGENIC_EXON_DUP",
    "PREDICTED_PARTIAL_EXON_DUP",
    "PREDICTED_TSS_DUP",
    "PREDICTED_DUP_PARTIAL",
    "PREDICTED_INV_SPAN",
    "PREDICTED_UTR",
    "PREDICTED_INTRONIC",
    "PREDICTED_PROMOTER"
]


# Helper Functions
def init_table(column_buckets):
    cols = list(column_buckets)
    return defaultdict(lambda: {col: 0.0 for col in cols})

def first_value(value):
    if isinstance(value, (list, tuple)):
        return value[0] if value else None
    return value

def get_info_value(record, key):
    try:
        return record.info.get(key)
    except ValueError:
        return None

def get_string_info(record, key):
    value = first_value(get_info_value(record, key))
    return "" if value is None else str(value)

def get_int_info(record, key):
    value = first_value(get_info_value(record, key))
    if value is None or value == ".":
        return None
    return abs(int(value))

def get_float_info(record, key):
    value = first_value(get_info_value(record, key))
    if value is None or value == ".":
        return None
    return float(value)

def has_info(record, key):
    return key in record.info

def get_info_gene_names(record, key):
    value = get_info_value(record, key)
    if isinstance(value, str):
        values = [value]
    elif isinstance(value, (list, tuple)):
        values = value
    elif value is None:
        values = []
    else:
        values = [str(value)]

    genes = set()
    for item in values:
        if item is None:
            continue
        for gene_name in str(item).split(","):
            gene_name = gene_name.strip()
            if gene_name and gene_name != ".":
                genes.add(gene_name)
    return genes

def get_vep_field_idx(header):
    if "vep" not in header.info:
        return {}
    description = header.info["vep"].description or ""
    match = re.search(r"Format:\s*([^\"]+)", description)
    if match is None:
        return {}
    fields = [field.strip() for field in match.group(1).strip().rstrip(".").split("|")]
    return {field: idx for idx, field in enumerate(fields)}

def extract_all_consequences(record, vep_field_idx):
    consequence_idx = vep_field_idx.get("Consequence")
    if consequence_idx is None or "vep" not in record.info:
        return set()
    annotations = get_info_value(record, "vep")
    if isinstance(annotations, str):
        annotations = [annotations]
    consequences = set()
    for annotation in annotations:
        fields = annotation.split("|")
        if consequence_idx >= len(fields):
            continue
        for consequence in fields[consequence_idx].split("&"):
            consequence = consequence.strip()
            if consequence:
                consequences.add(consequence)
    return consequences

def determine_vep_label(consequences):
    for label, consequence_terms in CONSEQUENCE_PRIORITY:
        if consequence_terms & consequences:
            return label
    return "Other"

def determine_svannotate_label(record):
    for label, field in SVANNOTATE_PRIORITY:
        if has_info(record, field):
            return label
    return None

def determine_column(allele_type, allele_length, length_bins, size_labels):
    allele_type = (allele_type or "").lower()
    if allele_type == "snv": return "SNV"
    if allele_type == "trv": return "TRV"
    if allele_length is None: return "Other"
    size_label = get_size_bucket(allele_length, length_bins, size_labels)
    if allele_type in ("ins", "dup"): return f"INS {size_label}"
    if allele_type == "del": return f"DEL {size_label}"
    return "Other"

def get_af_bucket(af_value, af_bins, af_labels):
    if not af_labels or af_value is None:
        return None
    for i, label in enumerate(af_labels):
        if i == 0:
            if af_value <= af_bins[1]:
                return label
        elif i < len(af_bins) - 1:
            if af_bins[i] < af_value <= af_bins[i + 1]:
                return label
        else:
            if af_value > af_bins[-1]:
                return label
    return None

def determine_af_column(allele_type, af_value, af_bins, af_labels):
    if not af_labels or af_value is None:
        return None
    allele_type = (allele_type or "").lower()
    af_label = get_af_bucket(af_value, af_bins, af_labels)
    if af_label is None:
        return None
    if allele_type == "snv": return f"SNV [{af_label}]"
    if allele_type == "trv": return f"TRV [{af_label}]"
    if allele_type in ("ins", "dup"): return f"INS [{af_label}]"
    if allele_type == "del": return f"DEL [{af_label}]"
    return f"Other [{af_label}]"

def determine_row_weights(record, allele_type_value, allele_subtype_value, vep_field_idx):
    allele_type = (allele_type_value or "").lower()
    subtype = (allele_subtype_value or "").lower()
    consequences = extract_all_consequences(record, vep_field_idx)
    vep_label = determine_vep_label(consequences)
    svannotate_label = determine_svannotate_label(record)

    row_weights = {row: 0 for row in ROW_ORDER}
    row_weights[("", INTERNAL_TOTAL_LABEL)] = 1

    is_alu = subtype in ("alu_ins", "alu_del")
    is_line = subtype in ("line_ins", "line_del")
    is_sva = subtype in ("sva_ins", "sva_del")
    if is_alu or is_line or is_sva: row_weights[("", "ME")] = 1
    if is_alu: row_weights[("ME", "ALU")] = 1
    if is_line: row_weights[("ME", "LINE")] = 1
    if is_sva: row_weights[("ME", "SVA")] = 1

    is_filtered_tandem = subtype == "tandem_dup" and allele_type == "dup"
    is_svan_tandem = subtype == "tandem_dup" and allele_type == "ins"
    is_dup_interspersed = subtype == "dup_interspersed"
    is_dup_complex = subtype == "complex_dup"
    is_dup_inverted = subtype == "inv_dup"
    is_numt = subtype == "numt"
    is_tr_parsed = has_info(record, "TR_PARSED")

    if is_filtered_tandem or is_svan_tandem or is_dup_interspersed or is_dup_complex or is_dup_inverted:
        row_weights[("", "Duplication")] = 1
    if is_filtered_tandem: row_weights[("Duplication", "Filtered Tandem")] = 1
    if is_svan_tandem: row_weights[("Duplication", "SVAN Tandem")] = 1
    if is_dup_interspersed: row_weights[("Duplication", "Interspersed")] = 1
    if is_dup_complex: row_weights[("Duplication", "Complex")] = 1
    if is_dup_inverted: row_weights[("Duplication", "Inverted")] = 1
    if is_numt: row_weights[("", "NUMT")] = 1
    if is_tr_parsed: row_weights[("", "TR Parsed")] = 1

    row_weights[("", "VEP")] = 1
    row_weights[("VEP", vep_label)] = 1

    if svannotate_label is not None:
        row_weights[("", "SVAnnotate")] = 1
        row_weights[("SVAnnotate", svannotate_label)] = 1

    is_dbsnp = has_info(record, "dbSNP_ID") or has_info(record, "dbGaP_ID")
    is_gnomad_matched = has_info(record, "gnomAD_V4_match_ID")

    row_weights[("", "Concordance")] = 1
    if not is_dbsnp:
        row_weights[("Concordance", "dbSNP Missing")] = 1
        row_weights[("Concordance", f"dbSNP Missing + {vep_label} (VEP)")] = 1
        if svannotate_label is not None:
            row_weights[("Concordance", f"dbSNP Missing + {svannotate_label} (SVAnnotate)")] = 1

    if not is_gnomad_matched:
        row_weights[("Concordance", "gnomAD Missing")] = 1
        row_weights[("Concordance", f"gnomAD Missing + {vep_label} (VEP)")] = 1
        if svannotate_label is not None:
            row_weights[("Concordance", f"gnomAD Missing + {svannotate_label} (SVAnnotate)")] = 1

    if not is_dbsnp and not is_gnomad_matched:
        row_weights[("Concordance", "dbSNP/gnomAD Missing")] = 1
        row_weights[("Concordance", f"dbSNP/gnomAD Missing + {vep_label} (VEP)")] = 1
        if svannotate_label is not None:
            row_weights[("Concordance", f"dbSNP/gnomAD Missing + {svannotate_label} (SVAnnotate)")] = 1

    return {row_key: weight for row_key, weight in row_weights.items() if weight > 0}

def get_genotype_weights(record):
    carrier_count = 0
    alt_allele_count = 0
    for sample in record.samples.values():
        genotype = sample.get("GT")
        if genotype is None: continue
        alt_alleles = sum(1 for allele in genotype if allele is not None and allele > 0)
        alt_allele_count += alt_alleles
        if alt_alleles > 0: carrier_count += 1
    return carrier_count, alt_allele_count

def trim_alleles(ref, alt):
    i = 0
    while i < min(len(ref), len(alt)) and ref[i] == alt[i]:
        i += 1
    ref, alt = ref[i:], alt[i:]
    j = 0
    while j < min(len(ref), len(alt)) and ref[-1-j] == alt[-1-j]:
        j += 1
    return (ref[:-j], alt[:-j]) if j else (ref, alt)

def count_genotypes_for_alt(record, alt_idx_1based):
    n_ref = n_het = n_hom = n_missing = 0
    for s in record.samples.values():
        gt = s.get("GT")
        if gt is None or all(a is None for a in gt):
            n_missing += 1
            continue
        n = sum(1 for a in gt if a == alt_idx_1based)
        if n == 0: n_ref += 1
        elif n == 1: n_het += 1
        else: n_hom += 1
    return n_ref, n_het, n_hom, n_missing

def format_list_column(row_key):
    annotation_group, annotation = row_key
    if annotation_group: return annotation
    if annotation in PARENT_ROWS: return f"{annotation} - {INTERNAL_TOTAL_LABEL}"
    return annotation

def get_list_columns():
    default_labels = [format_list_column(row_key) for row_key in ROW_ORDER]
    label_counts = {}
    for label in default_labels:
        label_counts[label] = label_counts.get(label, 0) + 1
    list_columns = []
    for row_key, default_label in zip(ROW_ORDER, default_labels):
        annotation_group, annotation = row_key
        if annotation_group and label_counts[default_label] > 1:
            list_columns.append(f"{annotation_group} - {annotation}")
        else:
            list_columns.append(default_label)
    return list_columns

def write_table(path, table_data, integer_output, split_by_region, column_buckets):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        header = ["category", "sub_category", "tr_status"]
        if split_by_region:
            header.append("region")
        writer.writerow(header + column_buckets)
        for key in sorted(table_data.keys()):
            if split_by_region:
                category, sub_category, tr_status, region = key
            else:
                category, sub_category, tr_status = key
            display_category = DISPLAY_ALL_LABEL if category == INTERNAL_TOTAL_LABEL else category
            display_sub_category = DISPLAY_ALL_LABEL if sub_category == INTERNAL_TOTAL_LABEL else sub_category
            display_tr_status = DISPLAY_ALL_LABEL if tr_status == INTERNAL_TOTAL_LABEL else tr_status
            row_prefix = [display_category, display_sub_category, display_tr_status]
            if split_by_region:
                display_region = DISPLAY_ALL_LABEL if region == INTERNAL_TOTAL_LABEL else region
                row_prefix.append(display_region)
            values = []
            for col in column_buckets:
                value = table_data[key].get(col, 0.0)
                if integer_output: values.append(str(int(value)))
                else: values.append(str(value))
            writer.writerow(row_prefix + values)

def write_gene_counts(path, gene_count_data):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_name", "category", "count"])
        for gene in gene_count_data:
            for category, count in gene_count_data[gene].items():
                writer.writerow([gene, category, str(count)])


# Initialize data
site_table = init_table(COLUMN_BUCKETS_SUMMARY)
sample_table = init_table(COLUMN_BUCKETS_SUMMARY)
allele_table = init_table(COLUMN_BUCKETS_SUMMARY)

vcf_in = pysam.VariantFile(VCF_PATH)
vep_field_idx = get_vep_field_idx(vcf_in.header)
sample_count = len(vcf_in.header.samples)

with open(SAMPLE_COUNT_OUTPUT, "w") as handle:
    handle.write(f"{sample_count}\n")

trios = []
all_probands = []
if CREATE_PLOTTING and TRIO_DEF_PATH:
    with open(TRIO_DEF_PATH) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) >= 3:
                trios.append((fields[0], fields[1], fields[2]))
    all_probands = sorted(set(child for child, _, _ in trios))

if CREATE_PLOTTING:
    plotting_sample_names = list(vcf_in.header.samples)
    plot_site_table = defaultdict(lambda: {col: 0 for col in COLUMN_BUCKETS_PLOTTING + COLUMN_BUCKETS_PLOTTING_AF + COLUMN_BUCKETS_PLOTTING_MATCHED})
    plot_sample_table = {s: defaultdict(int) for s in plotting_sample_names}
    plot_allele_table = {s: defaultdict(int) for s in plotting_sample_names}

    if trios:
        DENOVO_SITE_TYPES = ['Proband', 'Mendelian', 'Paternal', 'Maternal', 'Paternal Total', 'Maternal Total']
        DENOVO_ALLELE_TYPES = ['Proband Alleles', 'Mendelian Alleles', 'Paternal Alleles', 'Maternal Alleles',
                               'Paternal Total Alleles', 'Maternal Total Alleles']
        DENOVO_COUNT_TYPES = DENOVO_SITE_TYPES + DENOVO_ALLELE_TYPES
        plot_denovo_table = {ct: {s: defaultdict(float) for s in all_probands} for ct in DENOVO_COUNT_TYPES}

# Populate outputs
with open(VARIANT_LIST_OUTPUT, 'w', newline='') as vl_handle:
    vl_writer = csv.writer(vl_handle, delimiter='\t')
    if CREATE_PLOTTING:
        vl_writer.writerow(['allele_class', 'svlen', 'af', 'n_ref', 'n_het', 'n_hom', 'n_missing', 'region'])

    with open(LIST_OUTPUT, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["variant_id", "classification"] + get_list_columns())

        for record in vcf_in:
            a_type = get_string_info(record, "allele_type")
            a_sub = get_string_info(record, "allele_subtype")
            a_length = get_int_info(record, "allele_length")

            if MAX_LENGTH is not None and (a_length is None or a_length > MAX_LENGTH):
                continue
            if MIN_LENGTH is not None and (a_length is None or a_length < MIN_LENGTH):
                continue

            carrier_count, alt_allele_count = get_genotype_weights(record)
            if carrier_count == 0:
                continue

            column_summary = determine_column(a_type, a_length, LENGTH_BINS_SUMMARY, SIZE_LABELS_SUMMARY)
            column_plotting = determine_column(a_type, a_length, LENGTH_BINS_PLOTTING, SIZE_LABELS_PLOTTING) if CREATE_PLOTTING else None
            row_weights = determine_row_weights(record, a_type, a_sub, vep_field_idx)

            tr_status = "TR" if has_info(record, "TR_ENVELOPED") else "Not TR"
            is_db_or_gnomad_matched = has_info(record, "dbSNP_ID") or has_info(record, "gnomAD_V4_match_ID")
            regions = [INTERNAL_TOTAL_LABEL]
            if SPLIT_BY_REGION:
                region = get_string_info(record, "REGION")
                if region in REGION_ORDER:
                    regions.append(region)

            af = get_float_info(record, "AF")
            ac = get_int_info(record, "AC")

            is_large = a_length is not None and a_length > 50000

            buckets = ["all"]
            if is_large:
                buckets.append("all_large")
            if ac == 1:
                buckets.append("singleton")
                if is_large:
                    buckets.append("singleton_large")
            if af is not None and af < 0.001:
                buckets.append("ultra_rare")
                if is_large:
                    buckets.append("ultra_rare_large")
            if af is not None and af < 0.01:
                buckets.append("rare")
                if is_large:
                    buckets.append("rare_large")
            if af is not None and af > 0.05:
                buckets.append("common")
                if is_large:
                    buckets.append("common_large")

            for field in PREDICTED_FIELDS:
                if has_info(record, field):
                    for gene_name in get_info_gene_names(record, field):
                        for bucket in buckets:
                            category = f"{field}.{bucket}"
                            gene_counts[gene_name][category] += 1
                            if CREATE_PER_SAMPLE:
                                gene_sample_counts[gene_name][category] += carrier_count
                            if CREATE_PER_ALLELE:
                                gene_allele_counts[gene_name][category] += alt_allele_count

            for row_key, weight in row_weights.items():
                cat = row_key[0] if row_key[0] else row_key[1]
                ann = row_key[1] if row_key[0] else INTERNAL_TOTAL_LABEL
                for current_tr_status in [INTERNAL_TOTAL_LABEL, tr_status]:
                    if SPLIT_BY_REGION:
                        for current_region in regions:
                            full_key = (cat, ann, current_tr_status, current_region)
                            site_table[full_key][column_summary] += weight
                            if CREATE_PER_SAMPLE: sample_table[full_key][column_summary] += carrier_count * weight
                            if CREATE_PER_ALLELE: allele_table[full_key][column_summary] += alt_allele_count * weight
                    else:
                        full_key = (cat, ann, current_tr_status)
                        site_table[full_key][column_summary] += weight
                        if CREATE_PER_SAMPLE: sample_table[full_key][column_summary] += carrier_count * weight
                        if CREATE_PER_ALLELE: allele_table[full_key][column_summary] += alt_allele_count * weight

            if CREATE_PLOTTING:
                # Update plotting site table (Total annotation row only)
                for current_tr_status in [INTERNAL_TOTAL_LABEL, tr_status]:
                    for current_region in regions:
                        if SPLIT_BY_REGION:
                            plot_site_key = (INTERNAL_TOTAL_LABEL, INTERNAL_TOTAL_LABEL, current_tr_status, current_region)
                        else:
                            plot_site_key = (INTERNAL_TOTAL_LABEL, INTERNAL_TOTAL_LABEL, current_tr_status)
                        plot_site_table[plot_site_key][column_plotting] += 1
                        if is_db_or_gnomad_matched:
                            plot_site_table[plot_site_key][f'{column_plotting} (Matched)'] += 1
                        if COLUMN_BUCKETS_PLOTTING_AF:
                            af_col = determine_af_column(a_type, af, AF_BINS_PLOTTING, AF_LABELS_PLOTTING)
                            if af_col is not None:
                                plot_site_table[plot_site_key][af_col] += 1

                # Update per-sample plotting tables and collect genotype counts
                vl_n_ref = vl_n_het = vl_n_hom = vl_n_miss = 0
                for sample_name, sample_data in record.samples.items():
                    genotype = sample_data.get("GT")
                    if genotype is None or all(a is None for a in genotype):
                        vl_n_miss += 1
                        continue
                    n_alleles = sum(1 for a in genotype if a is not None and a > 0)
                    if n_alleles == 0:
                        vl_n_ref += 1
                    elif n_alleles == 1:
                        vl_n_het += 1
                    else:
                        vl_n_hom += 1
                    if n_alleles == 0:
                        continue
                    for current_tr_status in [INTERNAL_TOTAL_LABEL, tr_status]:
                        for current_region in regions:
                            skey = (column_plotting, current_tr_status, current_region)
                            plot_sample_table[sample_name][skey] += 1
                            plot_allele_table[sample_name][skey] += n_alleles
                            if COLUMN_BUCKETS_PLOTTING_AF:
                                af_col = determine_af_column(a_type, af, AF_BINS_PLOTTING, AF_LABELS_PLOTTING)
                                if af_col is not None:
                                    af_skey = (af_col, current_tr_status, current_region)
                                    plot_sample_table[sample_name][af_skey] += 1
                                    plot_allele_table[sample_name][af_skey] += n_alleles

                # Write variant list row(s)
                rec_region = get_string_info(record, "REGION") if SPLIT_BY_REGION else ""
                if a_type.lower() == 'trv':
                    for alt_idx, alt in enumerate(record.alts or []):
                        tref, talt = trim_alleles(record.ref, alt)
                        net = len(talt) - len(tref)
                        svlen = abs(net)
                        cls = "TR_SNV" if net == 0 else ("TR_INS" if net > 0 else "TR_DEL")
                        n_r, n_h, n_ho, n_m = count_genotypes_for_alt(record, alt_idx + 1)
                        n_geno = n_r + n_h + n_ho
                        af_info = get_info_value(record, "AF")
                        if isinstance(af_info, (list, tuple)) and alt_idx < len(af_info):
                            af_v = af_info[alt_idx]
                        elif n_geno > 0:
                            af_v = (n_h + 2 * n_ho) / (n_geno * 2)
                        else:
                            af_v = ""
                        vl_writer.writerow([cls, svlen, af_v, n_r, n_h, n_ho, n_m, rec_region])
                else:
                    af_v = get_float_info(record, "AF")
                    vl_writer.writerow([a_type, a_length or 0, af_v if af_v is not None else "", vl_n_ref, vl_n_het, vl_n_hom, vl_n_miss, rec_region])

                # Update plotting denovo table
                if trios:
                    for child, father, mother in trios:
                        child_sample = record.samples.get(child)
                        child_gt = child_sample.get("GT") if child_sample is not None else None
                        child_carries = child_gt is not None and any(a is not None and a > 0 for a in child_gt)

                        father_sample = record.samples.get(father)
                        father_gt = father_sample.get("GT") if father_sample is not None else None
                        father_has = father_gt is not None and any(a is not None and a > 0 for a in father_gt)

                        mother_sample = record.samples.get(mother)
                        mother_gt = mother_sample.get("GT") if mother_sample is not None else None
                        mother_has = mother_gt is not None and any(a is not None and a > 0 for a in mother_gt)

                        is_mendelian = father_has or mother_has
                        n_child = sum(1 for a in child_gt if a is not None and a > 0) if child_gt else 0
                        n_father = sum(1 for a in father_gt if a is not None and a > 0) if father_gt else 0
                        n_mother = sum(1 for a in mother_gt if a is not None and a > 0) if mother_gt else 0
                        n_child_denovo = max(0, n_child - n_father - n_mother)
                        n_child_inherited = n_child - n_child_denovo
                        parental_total = n_father + n_mother
                        # Proportional inheritance formula
                        fa_transmitted = (n_child_inherited * n_father / parental_total) if parental_total > 0 else 0.0
                        mo_transmitted = (n_child_inherited * n_mother / parental_total) if parental_total > 0 else 0.0

                        for current_tr_status in [INTERNAL_TOTAL_LABEL, tr_status]:
                            for current_region in regions:
                                full_key = (column_plotting, current_tr_status, current_region)
                                if father_has:
                                    plot_denovo_table['Paternal Total'][child][full_key] += 1
                                if mother_has:
                                    plot_denovo_table['Maternal Total'][child][full_key] += 1
                                if child_carries:
                                    plot_denovo_table['Proband'][child][full_key] += 1
                                    if is_mendelian:
                                        plot_denovo_table['Mendelian'][child][full_key] += 1
                                    if father_has:
                                        plot_denovo_table['Paternal'][child][full_key] += 1
                                    if mother_has:
                                        plot_denovo_table['Maternal'][child][full_key] += 1
                                plot_denovo_table['Paternal Total Alleles'][child][full_key] += n_father
                                plot_denovo_table['Maternal Total Alleles'][child][full_key] += n_mother
                                plot_denovo_table['Proband Alleles'][child][full_key] += n_child
                                plot_denovo_table['Mendelian Alleles'][child][full_key] += n_child_inherited
                                plot_denovo_table['Paternal Alleles'][child][full_key] += fa_transmitted
                                plot_denovo_table['Maternal Alleles'][child][full_key] += mo_transmitted

            if CREATE_LIST:
                writer.writerow([str(record.id or "."), column_summary] + ["1" if row_key in row_weights else "0" for row_key in ROW_ORDER])

# Write output tables
write_table(SITE_OUTPUT, site_table, integer_output=True, split_by_region=SPLIT_BY_REGION, column_buckets=COLUMN_BUCKETS_SUMMARY)
write_table(SAMPLE_OUTPUT, sample_table, integer_output=True, split_by_region=SPLIT_BY_REGION, column_buckets=COLUMN_BUCKETS_SUMMARY)
write_table(ALLELE_OUTPUT, allele_table, integer_output=True, split_by_region=SPLIT_BY_REGION, column_buckets=COLUMN_BUCKETS_SUMMARY)
write_gene_counts(GENE_OUTPUT, gene_counts)
write_gene_counts(GENE_SAMPLE_OUTPUT, gene_sample_counts)
write_gene_counts(GENE_ALLELE_OUTPUT, gene_allele_counts)

if CREATE_PLOTTING:
    write_table(PLOT_SITE_OUTPUT, plot_site_table, integer_output=True, split_by_region=SPLIT_BY_REGION, column_buckets=COLUMN_BUCKETS_PLOTTING + COLUMN_BUCKETS_PLOTTING_AF + COLUMN_BUCKETS_PLOTTING_MATCHED)
    plot_col_keys = []
    plot_col_names = []
    for col_bucket in COLUMN_BUCKETS_PLOTTING + COLUMN_BUCKETS_PLOTTING_AF:
        plot_col_keys.append((col_bucket, INTERNAL_TOTAL_LABEL, INTERNAL_TOTAL_LABEL))
        plot_col_names.append(f'All - {col_bucket}')
        for tr in ['TR', 'Not TR']:
            plot_col_keys.append((col_bucket, tr, INTERNAL_TOTAL_LABEL))
            plot_col_names.append(f'All - {col_bucket} - {tr}')
        if SPLIT_BY_REGION:
            for region in REGION_ORDER:
                plot_col_keys.append((col_bucket, INTERNAL_TOTAL_LABEL, region))
                plot_col_names.append(f'All - {col_bucket} - {region}')

    with open(PLOT_SAMPLE_OUTPUT, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample_id'] + plot_col_names)
        for sample_name in plotting_sample_names:
            writer.writerow([sample_name] + [str(int(plot_sample_table[sample_name].get(k, 0))) for k in plot_col_keys])

    with open(PLOT_ALLELE_OUTPUT, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample_id'] + plot_col_names)
        for sample_name in plotting_sample_names:
            writer.writerow([sample_name] + [str(int(plot_allele_table[sample_name].get(k, 0))) for k in plot_col_keys])

    if trios:
        denovo_col_keys = []
        denovo_col_names = []
        denovo_combos = [(INTERNAL_TOTAL_LABEL, INTERNAL_TOTAL_LABEL, '')]
        denovo_combos += [('TR', INTERNAL_TOTAL_LABEL, ' - TR'), ('Not TR', INTERNAL_TOTAL_LABEL, ' - Not TR')]
        if SPLIT_BY_REGION:
            denovo_combos += [(INTERNAL_TOTAL_LABEL, r, f' - {r}') for r in REGION_ORDER]
        for col_bucket in COLUMN_BUCKETS_PLOTTING:
            for tr_stat, reg, suffix in denovo_combos:
                for count_type in DENOVO_COUNT_TYPES:
                    denovo_col_keys.append((col_bucket, tr_stat, reg, count_type))
                    denovo_col_names.append(f'All - {col_bucket}{suffix} ({count_type})')

        def fmt_val(v):
            return str(int(v)) if v % 1 == 0 else str(round(v, 6))

        with open(DENOVO_OUTPUT, 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['sample_id'] + denovo_col_names)
            for proband in all_probands:
                row = [proband]
                for col_bucket, tr_stat, reg, count_type in denovo_col_keys:
                    full_key = (col_bucket, tr_stat, reg)
                    row.append(fmt_val(plot_denovo_table[count_type][proband].get(full_key, 0)))
                writer.writerow(row)
else:
    with open(PLOT_SITE_OUTPUT, "w") as handle:
        pass
    with open(PLOT_SAMPLE_OUTPUT, "w") as handle:
        pass
    with open(PLOT_ALLELE_OUTPUT, "w") as handle:
        pass

if not CREATE_LIST:
    with open(LIST_OUTPUT, "w", newline="") as handle:
        pass
PYCODE
    >>>

    output {
        File site_counts_tsv = "~{prefix}.sites.raw.tsv"
        File sample_counts_tsv = "~{prefix}.samples.raw.tsv"
        File allele_counts_tsv = "~{prefix}.alleles.raw.tsv"
        File list_tsv = "~{prefix}.list.raw.tsv"
        File gene_counts_tsv = "~{prefix}.genes.raw.tsv"
        File gene_sample_counts_tsv = "~{prefix}.genes.samples.raw.tsv"
        File gene_allele_counts_tsv = "~{prefix}.genes.alleles.raw.tsv"
        File sample_count_file = "~{prefix}.sample_count.txt"
        File plotting_site_tsv = "~{prefix}.plotting_sites.raw.tsv"
        File plotting_sample_tsv = "~{prefix}.plotting_samples.raw.tsv"
        File plotting_allele_tsv = "~{prefix}.plotting_alleles.raw.tsv"
        File? plotting_denovo_tsv = "~{prefix}.plotting_denovo.raw.tsv"
        File plotting_variant_list_tsv = "~{prefix}.plotting_variant_list.raw.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 4 * ceil(size([vcf, vcf_idx], "GB")) + 10,
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

task MergeAnnotationCountTables {
    input {
        Array[File] count_tsvs
        Array[File] sample_count_files
        String normalization_mode
        Boolean split_by_region = false
        Boolean drop_empty_types = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
import csv

COUNT_FILES = [path for path in "~{sep=',' count_tsvs}".split(",") if path]
SAMPLE_COUNT_FILES = [path for path in "~{sep=',' sample_count_files}".split(",") if path]
MODE = "~{normalization_mode}"
SPLIT_BY_REGION = "~{split_by_region}".lower() == "true"
DROP_EMPTY_TYPES = "~{drop_empty_types}".lower() == "true"
OUTPUT = "~{prefix}.tsv"

INTERNAL_TOTAL_LABEL = "Total"
DISPLAY_ALL_LABEL = "All"
TR_ORDER = [INTERNAL_TOTAL_LABEL, "TR", "Not TR"]
REGION_ORDER = [INTERNAL_TOTAL_LABEL, "US", "RM", "SD", "SR"]

VEP_LABELS = ["LoF", "Missense", "Synonymous", "Intronic", "Intergenic", "Other"]
SVANNOTATE_LABELS = ["LoF", "Copy Gain", "IED", "PED", "TSSD", "DP", "UTR"]

ROW_ORDER = [
    ("", INTERNAL_TOTAL_LABEL),
    ("", "ME"),
    ("ME", "ALU"),
    ("ME", "LINE"),
    ("ME", "SVA"),
    ("", "Duplication"),
    ("Duplication", "Filtered Tandem"),
    ("Duplication", "SVAN Tandem"),
    ("Duplication", "Interspersed"),
    ("Duplication", "Complex"),
    ("Duplication", "Inverted"),
    ("", "TR Parsed"),
    ("", "NUMT"),
    ("", "VEP"),
]
ROW_ORDER += [("VEP", label) for label in VEP_LABELS]
ROW_ORDER.append(("", "SVAnnotate"))
ROW_ORDER += [("SVAnnotate", label) for label in SVANNOTATE_LABELS]
ROW_ORDER.append(("", "Concordance"))
for missing_label in ["dbSNP Missing", "gnomAD Missing", "dbSNP/gnomAD Missing"]:
    ROW_ORDER.append(("Concordance", missing_label))
    ROW_ORDER += [("Concordance", f"{missing_label} + {label} (VEP)") for label in VEP_LABELS]
    ROW_ORDER += [("Concordance", f"{missing_label} + {label} (SVAnnotate)") for label in SVANNOTATE_LABELS]

def get_cat_ann(row_key):
    group, ann = row_key
    if not group:
        return ann, INTERNAL_TOTAL_LABEL
    return group, ann

def normalize_input_label(value):
    return INTERNAL_TOTAL_LABEL if value == DISPLAY_ALL_LABEL else value

header = None
counts = {}
group_column_count = 4 if SPLIT_BY_REGION else 3

for path in COUNT_FILES:
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        current_header = next(reader, None)
        if current_header is None:
            continue
        if header is None:
            header = current_header
        elif current_header != header:
            raise ValueError(f"Mismatched headers while merging count tables: {path}")

        for row in reader:
            key = tuple(normalize_input_label(value) for value in row[:group_column_count])
            if key not in counts:
                counts[key] = [0.0] * (len(header) - group_column_count)
            for i, val in enumerate(row[group_column_count:]):
                counts[key][i] += float(val)

if header is None:
    raise ValueError("No count tables were provided for merging")

denominator = 1.0
if MODE != "sites":
    sample_counts = set()
    for path in SAMPLE_COUNT_FILES:
        with open(path, "r") as handle:
            sample_counts.add(int(handle.read().strip()))

    if len(sample_counts) != 1:
        raise ValueError("All input VCFs must contain the same number of samples for normalized count outputs")

    sample_count = sample_counts.pop()
    if sample_count <= 0:
        raise ValueError("Sample-normalized outputs require at least one sample in every input VCF")

    denominator = float(sample_count)

value_columns = header[group_column_count:]
value_count = len(value_columns)

def col_type(col):
    if col.startswith("DEL "): return "DEL"
    if col.startswith("INS "): return "INS"
    return col

if DROP_EMPTY_TYPES:
    type_totals = {}
    for row_values in counts.values():
        for col, v in zip(value_columns, row_values):
            type_totals[col_type(col)] = type_totals.get(col_type(col), 0.0) + v
    keep_idx = [i for i, col in enumerate(value_columns) if type_totals.get(col_type(col), 0.0) > 0]
else:
    keep_idx = list(range(value_count))

header_out = header[:group_column_count] + [value_columns[i] for i in keep_idx]

with open(OUTPUT, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(header_out)

    def get_row_prefix(c, a, t, region):
        display_category = DISPLAY_ALL_LABEL if c == INTERNAL_TOTAL_LABEL else c
        display_sub_category = DISPLAY_ALL_LABEL if a == INTERNAL_TOTAL_LABEL else a
        display_tr_status = DISPLAY_ALL_LABEL if t == INTERNAL_TOTAL_LABEL else t
        row_prefix = [display_category, display_sub_category, display_tr_status]
        if SPLIT_BY_REGION:
            display_region = DISPLAY_ALL_LABEL if region == INTERNAL_TOTAL_LABEL else region
            row_prefix.append(display_region)
        return row_prefix

    def write_count_row(c, a, t, region, values):
        kept = [values[i] for i in keep_idx]
        if MODE == "sites":
            formatted = [str(int(round(v))) for v in kept]
        else:
            formatted = [f"{v / denominator:.2f}" for v in kept]

        writer.writerow(get_row_prefix(c, a, t, region) + formatted)

    for row_key in ROW_ORDER:
        cat, ann = get_cat_ann(row_key)
        for tr in TR_ORDER:
            if SPLIT_BY_REGION:
                for region in REGION_ORDER:
                    values = counts.get((cat, ann, tr, region), [0.0] * value_count)
                    write_count_row(cat, ann, tr, region, values)
            else:
                values = counts.get((cat, ann, tr), [0.0] * value_count)
                write_count_row(cat, ann, tr, INTERNAL_TOTAL_LABEL, values)
PYCODE
    >>>

    output {
        File merged_counts_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size(count_tsvs, "GB")) + 10,
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

task MergeAnnotationListTables {
    input {
        Array[File] list_tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
import csv


LIST_FILES = [path for path in "~{sep=',' list_tsvs}".split(",") if path]
OUTPUT = "~{prefix}.tsv"


header = None
variant_order = []
variants = {}

for path in LIST_FILES:
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        current_header = next(reader)
        if header is None:
            header = current_header
        elif current_header != header:
            raise ValueError(f"Mismatched headers while merging annotation list tables: {path}")

        for row in reader:
            variant_id = row[0]
            classification = row[1]
            memberships = [int(value) for value in row[2:]]

            if variant_id not in variants:
                variants[variant_id] = [classification, memberships]
                variant_order.append(variant_id)
                continue

            existing_classification, existing_memberships = variants[variant_id]
            if existing_classification != classification:
                raise ValueError(
                    f"Variant {variant_id} has conflicting classifications while merging annotation list tables"
                )

            variants[variant_id][1] = [max(left, right) for left, right in zip(existing_memberships, memberships)]

if header is None:
    raise ValueError("No annotation list tables were provided for merging")

with open(OUTPUT, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(header)
    for variant_id in variant_order:
        classification, memberships = variants[variant_id]
        writer.writerow([variant_id, classification] + [str(value) for value in memberships])
PYCODE
    >>>

    output {
        File merged_list_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size(list_tsvs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: 4 + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeGeneCountTables {
    input {
        Array[File] count_tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
import csv
from collections import defaultdict

COUNT_FILES = [path for path in "~{sep=',' count_tsvs}".split(",") if path]
OUTPUT = "~{prefix}.tsv"

PREDICTED_FIELDS = [
    "PREDICTED_LOF",
    "PREDICTED_COPY_GAIN",
    "PREDICTED_INTRAGENIC_EXON_DUP",
    "PREDICTED_PARTIAL_EXON_DUP",
    "PREDICTED_TSS_DUP",
    "PREDICTED_DUP_PARTIAL",
    "PREDICTED_INV_SPAN",
    "PREDICTED_UTR",
    "PREDICTED_INTRONIC",
    "PREDICTED_PROMOTER"
]

BUCKETS = [
    "all",
    "common",
    "rare",
    "ultra_rare",
    "singleton",
    "all_large",
    "common_large",
    "rare_large",
    "ultra_rare_large",
    "singleton_large",
]

new_columns = []
for f in PREDICTED_FIELDS:
    for b in BUCKETS:
        new_columns.append(f"{f}.{b}")

gene_counts = defaultdict(lambda: {col: 0 for col in new_columns})

for path in COUNT_FILES:
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            next(reader) # skip header
        except StopIteration:
            pass
        for row in reader:
            if len(row) == 3:
                gene, cat, count = row
                gene_counts[gene][cat] += int(count)

with open(OUTPUT, "w", newline="") as out_handle:
    writer = csv.writer(out_handle, delimiter="\t")
    writer.writerow(["gene_name"] + new_columns)
    for gene in sorted(gene_counts.keys()):
        counts = gene_counts[gene]
        row = [gene] + [str(counts.get(col, 0)) for col in new_columns]
        writer.writerow(row)

PYCODE
    >>>

    output {
        File merged_counts_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size(count_tsvs, "GB")) + 10,
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

task MergeVariantListTables {
    input {
        Array[File] variant_list_tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
import csv

INPUT_FILES = [path for path in "~{sep=',' variant_list_tsvs}".split(",") if path]
OUTPUT = "~{prefix}.tsv"

header = None
with open(OUTPUT, "w", newline="") as out_handle:
    writer = csv.writer(out_handle, delimiter="\t")
    for path in INPUT_FILES:
        with open(path, "r", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            current_header = next(reader, None)
            if current_header is None:
                continue
            if header is None:
                header = current_header
                writer.writerow(header)
            for row in reader:
                writer.writerow(row)
PYCODE
    >>>

    output {
        File merged_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size(variant_list_tsvs, "GB")) + 10,
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

task MergeSampleSpecificTables {
    input {
        Array[File] count_tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
import csv
from collections import defaultdict

COUNT_FILES = [path for path in "~{sep=',' count_tsvs}".split(",") if path]
OUTPUT = "~{prefix}.tsv"

header = None
sample_counts = defaultdict(lambda: defaultdict(int))
sample_order = []

for path in COUNT_FILES:
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        current_header = next(reader, None)
        if current_header is None:
            continue
        if header is None:
            header = current_header
        elif current_header != header:
            raise ValueError(f"Mismatched headers while merging sample-specific tables: {path}")

        for row in reader:
            sample_id = row[0]
            if sample_id not in sample_counts:
                sample_order.append(sample_id)
            for col_name, val in zip(header[1:], row[1:]):
                sample_counts[sample_id][col_name] += float(val)

if header is None:
    with open(OUTPUT, "w") as handle:
        pass
else:
    active_columns = [col for col in header[1:] if any(sample_counts[s][col] > 0 for s in sample_order)]

    def fmt_val(v):
        return str(int(v)) if v % 1 == 0 else str(round(v, 6))

    with open(OUTPUT, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id"] + active_columns)
        for sample_id in sample_order:
            writer.writerow([sample_id] + [fmt_val(sample_counts[sample_id][col]) for col in active_columns])
PYCODE
    >>>

    output {
        File merged_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size(count_tsvs, "GB")) + 10,
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
