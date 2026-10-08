version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SummarizeTRVAllelesPerGenome {
    meta {
        description: [
            "This utility counts the tandem repeat variant (INFO/allele_type='trv') ALT alleles carried per genome across one or more VCFs, how many of them fall within coding sequences and how many distinct genes those coding alleles hit, producing the sentence 'An average of X alternative alleles were observed per individual genome, including Y that reside within coding sequences of Z genes'.",
            "Each sample counts every distinct ALT allele index in its genotype at a TRV record once, so 2/3 counts two alleles while 2/2 and 0/2 count one; missing alleles are ignored. This matches the per-ALT carrier counts in the `SummarizeAnnotations` plotting variant list, and no length filter is applied. Averages divide by the number of samples in the VCF header, which must be identical across all inputs.",
            "An ALT allele is coding when its changed span, found by trimming the bases REF and ALT share at the start and then at the end, overlaps a CDS line of `coding_gtf`; a pure insertion is coding only when both of its flanking bases lie within the same CDS line. Every CDS line counts regardless of gene type. A secondary metric instead calls every ALT allele coding when the whole REF span of its record overlaps a CDS line.",
            "Genes are the distinct `gene_name` values of the CDS lines hit by the coding ALT alleles a sample carries, unioned across all inputs before averaging."
        ]
    }

    parameter_meta {
        vcfs: "VCFs whose TRV alleles are counted."
        vcf_idxs: "Indexes for `vcfs`."
        coding_gtf: "From references."
        records_per_shard: "Number of records to keep within a single shard."
        summary_tsv: "TSV of metric and value rows, holding the record, allele and gene totals and their per-genome averages."
        summary_sentence: "Per-genome averages formatted into the summary sentence, to one decimal place."
        coding_alleles_tsv: "TSV with one row per coding TRV ALT allele, giving its record, allele index, carrier count and the genes whose CDS it overlaps."
        coding_genes_tsv: "TSV with one row per gene hit by a coding TRV ALT allele, giving the number of such alleles and of samples carrying at least one of them."
    }

    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix

        File coding_gtf

        Int? records_per_shard
        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_count
        RuntimeAttr? runtime_attr_merge
    }

    scatter (i in range(length(vcfs))) {
        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = vcfs[i],
                    vcf_idx = vcf_idxs[i],
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.input_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] shard_vcfs = select_first([ShardVcfByRecords.shards, [vcfs[i]]])
        Array[File] shard_vcf_idxs = select_first([ShardVcfByRecords.shard_idxs, [vcf_idxs[i]]])

        scatter (j in range(length(shard_vcfs))) {
            call CountTRVAllelesShard {
                input:
                    vcf = shard_vcfs[j],
                    vcf_idx = shard_vcf_idxs[j],
                    coding_gtf = coding_gtf,
                    prefix = "~{prefix}.input_~{i}.shard_~{j}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_count
            }
        }
    }

    call MergeTRVAlleleCounts {
        input:
            counts_tsvs = flatten(CountTRVAllelesShard.counts_tsv),
            coding_alleles_tsvs = flatten(CountTRVAllelesShard.coding_alleles_tsv),
            sample_genes_tsvs = flatten(CountTRVAllelesShard.sample_genes_tsv),
            sample_count_files = flatten(CountTRVAllelesShard.sample_count_file),
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File summary_tsv = MergeTRVAlleleCounts.summary_tsv
        String summary_sentence = MergeTRVAlleleCounts.summary_sentence
        File coding_alleles_tsv = MergeTRVAlleleCounts.coding_alleles_tsv
        File coding_genes_tsv = MergeTRVAlleleCounts.coding_genes_tsv
    }
}

task CountTRVAllelesShard {
    input {
        File vcf
        File vcf_idx
        File coding_gtf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf} > samples.txt
        wc -l < samples.txt | tr -d ' ' > ~{prefix}.sample_count.txt

        python3 <<'CODE'
import re
import subprocess
from bisect import bisect_right
from collections import Counter, defaultdict

GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')


def load_cds(path):
    intervals = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            intervals[fields[0]].append((int(fields[3]), int(fields[4]), GENE_NAME_RE.search(fields[8]).group(1)))
    index = {}
    for contig, rows in intervals.items():
        rows.sort()
        max_ends = []
        running = 0
        for _, end, _ in rows:
            running = max(running, end)
            max_ends.append(running)
        index[contig] = ([row[0] for row in rows], max_ends, rows)
    return index


def overlapping_genes(contig, start, end, contained=False):
    if contig not in CDS:
        return frozenset()
    starts, max_ends, rows = CDS[contig]
    genes = set()
    idx = bisect_right(starts, end) - 1
    while idx >= 0 and max_ends[idx] >= start:
        cds_start, cds_end, gene = rows[idx]
        if (cds_start <= start and cds_end >= end) if contained else cds_end >= start:
            genes.add(gene)
        idx -= 1
    return frozenset(genes)


def trim_alleles(ref, alt):
    i = 0
    while i < min(len(ref), len(alt)) and ref[i] == alt[i]:
        i += 1
    ref, alt = ref[i:], alt[i:]
    j = 0
    while j < min(len(ref), len(alt)) and ref[-1 - j] == alt[-1 - j]:
        j += 1
    return i, len(ref) - j


def alt_genes(contig, pos, ref, alt):
    prefix_len, ref_len = trim_alleles(ref, alt)
    start = pos + prefix_len
    if ref_len > 0:
        return overlapping_genes(contig, start, start + ref_len - 1)
    return overlapping_genes(contig, start - 1, start, contained=True)


GT_CACHE = {}


def gt_alts(gt):
    alts = GT_CACHE.get(gt)
    if alts is None:
        alts = frozenset(int(a) for a in re.split(r"[/|]", gt) if a not in (".", "0"))
        GT_CACHE[gt] = alts
    return alts


CDS = load_cds("~{coding_gtf}")
with open("samples.txt") as handle:
    samples = [line.strip() for line in handle if line.strip()]

counts = Counter()
sample_genes = defaultdict(set)
query = subprocess.Popen(
    ["bcftools", "query", "-i", 'INFO/allele_type="trv"', "-f", r"%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n", "~{vcf}"],
    stdout=subprocess.PIPE,
    text=True,
)

# Count distinct ALT alleles per sample, grouping identical genotypes and resolving per-sample genes only at coding records
with open("~{prefix}.coding_alleles.tsv", "w") as coding_out:
    coding_out.write("contig\tpos\tid\talt_index\tn_carriers\tgenes\n")
    for line in query.stdout:
        contig, pos, variant_id, ref, alts, *gts = line.rstrip("\n").split("\t")
        pos = int(pos)
        alts = alts.split(",")
        counts["n_trv_records"] += 1
        counts["n_trv_alt_alleles"] += len(alts)

        gt_counts = Counter(gts)
        n_carried = sum(len(gt_alts(gt)) * n for gt, n in gt_counts.items())
        counts["alt_allele_sum"] += n_carried
        if overlapping_genes(contig, pos, pos + len(ref) - 1):
            counts["coding_locus_alt_allele_sum"] += n_carried

        coding = {}
        for k, alt in enumerate(alts, start=1):
            genes = alt_genes(contig, pos, ref, alt)
            if genes:
                coding[k] = genes
        if not coding:
            continue

        carriers = Counter()
        for gt, n in gt_counts.items():
            for k in gt_alts(gt):
                if k in coding:
                    carriers[k] += n
        counts["n_coding_trv_alt_alleles"] += len(coding)
        counts["coding_alt_allele_sum"] += sum(carriers.values())
        for k, genes in coding.items():
            coding_out.write(f"{contig}\t{pos}\t{variant_id}\t{k}\t{carriers[k]}\t{','.join(sorted(genes))}\n")
        for sample, gt in zip(samples, gts):
            for k in gt_alts(gt):
                if k in coding:
                    sample_genes[sample] |= coding[k]

if query.wait() != 0:
    raise RuntimeError("bcftools query failed")

with open("~{prefix}.counts.tsv", "w") as handle:
    for key in ["n_trv_records", "n_trv_alt_alleles", "n_coding_trv_alt_alleles", "alt_allele_sum",
                "coding_alt_allele_sum", "coding_locus_alt_allele_sum"]:
        handle.write(f"{key}\t{counts[key]}\n")

with open("~{prefix}.sample_genes.tsv", "w") as handle:
    for sample in samples:
        for gene in sorted(sample_genes.get(sample, ())):
            handle.write(f"{sample}\t{gene}\n")
CODE
    >>>

    output {
        File counts_tsv = "~{prefix}.counts.tsv"
        File coding_alleles_tsv = "~{prefix}.coding_alleles.tsv"
        File sample_genes_tsv = "~{prefix}.sample_genes.tsv"
        File sample_count_file = "~{prefix}.sample_count.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, vcf_idx, coding_gtf], "GB")) + 10,
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

task MergeTRVAlleleCounts {
    input {
        Array[File] counts_tsvs
        Array[File] coding_alleles_tsvs
        Array[File] sample_genes_tsvs
        Array[File] sample_count_files
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
from collections import Counter, defaultdict

COUNT_FILES = "~{sep=',' counts_tsvs}".split(",")
CODING_ALLELE_FILES = "~{sep=',' coding_alleles_tsvs}".split(",")
SAMPLE_GENE_FILES = "~{sep=',' sample_genes_tsvs}".split(",")
SAMPLE_COUNT_FILES = "~{sep=',' sample_count_files}".split(",")

sample_counts = set()
for path in SAMPLE_COUNT_FILES:
    with open(path) as handle:
        sample_counts.add(int(handle.read().strip()))
if len(sample_counts) != 1:
    raise ValueError("All input VCFs must contain the same number of samples for per-genome averages")
n_samples = sample_counts.pop()
if n_samples <= 0:
    raise ValueError("Per-genome averages require at least one sample in every input VCF")

counts = Counter()
for path in COUNT_FILES:
    with open(path) as handle:
        for line in handle:
            key, value = line.rstrip("\n").split("\t")
            counts[key] += int(value)

# Concatenate the coding allele tables under one header and tally coding alleles per gene
gene_alleles = Counter()
with open("~{prefix}.coding_alleles.tsv", "w") as out:
    out.write("contig\tpos\tid\talt_index\tn_carriers\tgenes\n")
    for path in CODING_ALLELE_FILES:
        with open(path) as handle:
            next(handle)
            for line in handle:
                out.write(line)
                for gene in line.rstrip("\n").split("\t")[5].split(","):
                    gene_alleles[gene] += 1

# Union each sample's genes across shards before counting carriers per gene
sample_genes = defaultdict(set)
for path in SAMPLE_GENE_FILES:
    with open(path) as handle:
        for line in handle:
            sample, gene = line.rstrip("\n").split("\t")
            sample_genes[sample].add(gene)
gene_carriers = Counter(gene for genes in sample_genes.values() for gene in genes)

with open("~{prefix}.coding_genes.tsv", "w") as out:
    out.write("gene_name\tn_coding_alt_alleles\tn_carrier_samples\n")
    for gene in sorted(gene_alleles):
        out.write(f"{gene}\t{gene_alleles[gene]}\t{gene_carriers[gene]}\n")

alt_per_genome = counts["alt_allele_sum"] / n_samples
coding_per_genome = counts["coding_alt_allele_sum"] / n_samples
locus_per_genome = counts["coding_locus_alt_allele_sum"] / n_samples
genes_per_genome = sum(gene_carriers.values()) / n_samples

with open("~{prefix}.summary.tsv", "w") as out:
    out.write("metric\tvalue\n")
    out.write(f"n_samples\t{n_samples}\n")
    out.write(f"n_trv_records\t{counts['n_trv_records']}\n")
    out.write(f"n_trv_alt_alleles\t{counts['n_trv_alt_alleles']}\n")
    out.write(f"n_coding_trv_alt_alleles\t{counts['n_coding_trv_alt_alleles']}\n")
    out.write(f"n_coding_genes\t{len(gene_alleles)}\n")
    out.write(f"alt_alleles_per_genome\t{alt_per_genome:.4f}\n")
    out.write(f"coding_alt_alleles_per_genome\t{coding_per_genome:.4f}\n")
    out.write(f"coding_locus_alt_alleles_per_genome\t{locus_per_genome:.4f}\n")
    out.write(f"coding_genes_per_genome\t{genes_per_genome:.4f}\n")

with open("~{prefix}.summary_sentence.txt", "w") as out:
    out.write(
        f"An average of {alt_per_genome:,.1f} alternative alleles were observed per individual genome, "
        f"including {coding_per_genome:,.1f} that reside within coding sequences of {genes_per_genome:,.1f} genes\n"
    )
CODE
    >>>

    output {
        File summary_tsv = "~{prefix}.summary.tsv"
        String summary_sentence = read_string("~{prefix}.summary_sentence.txt")
        File coding_alleles_tsv = "~{prefix}.coding_alleles.tsv"
        File coding_genes_tsv = "~{prefix}.coding_genes.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(flatten([counts_tsvs, coding_alleles_tsvs, sample_genes_tsvs]), "GB")) + 10,
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
