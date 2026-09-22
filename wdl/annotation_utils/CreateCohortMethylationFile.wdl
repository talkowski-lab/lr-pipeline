version 1.0

import "../utils/Structs.wdl"

workflow CreateCohortMethylationFile {
    meta {
        description: [
            "This utility builds cohort-level CpG methylation matrices from per-sample `MethylationProfiling` BED outputs. For each contig, it merges every sample's combined and per-haplotype modification-score BEDs into a wide site-by-sample(/haplotype) matrix, filling `.` for sites missing in a given sample or haplotype. Samples can optionally be processed in shards (merged independently, then joined column-wise) to bound how many sample files are localized onto a single task at once."
        ]
    }

    parameter_meta {
        combined_beds: "Per-sample combined CpG pileup BEDs (`cpg_combined_bed` from MethylationProfiling)."
        combined_bed_idxs: "Indexes for `combined_beds`."
        hap1_beds: "Per-sample haplotype 1 CpG pileup BEDs (`cpg_hap1_bed` from MethylationProfiling)."
        hap1_bed_idxs: "Indexes for `hap1_beds`."
        hap2_beds: "Per-sample haplotype 2 CpG pileup BEDs (`cpg_hap2_bed` from MethylationProfiling)."
        hap2_bed_idxs: "Indexes for `hap2_beds`."
        sample_ids: "Sample IDs, parallel to `combined_beds`/`hap1_beds`/`hap2_beds`."
        contigs: "Contigs to process."
        samples_per_shard: "Maximum number of samples merged per shard before shards are joined column-wise into the final matrix."
        combined_methylation_beds: "Per-contig site-by-sample modification-score matrix BEDs."
        haplotype_methylation_beds: "Per-contig site-by-haplotype modification-score matrix BEDs."
    }

    input {
        Array[File] combined_beds
        Array[File] combined_bed_idxs
        Array[File] hap1_beds
        Array[File] hap1_bed_idxs
        Array[File] hap2_beds
        Array[File] hap2_bed_idxs
        Array[String] sample_ids
        Array[String] contigs
        String prefix

        Int? samples_per_shard

        String utils_docker

        RuntimeAttr? runtime_attr_merge_combined
        RuntimeAttr? runtime_attr_merge_haplotype
        RuntimeAttr? runtime_attr_join_combined
        RuntimeAttr? runtime_attr_join_haplotype
    }

    Int num_samples = length(sample_ids)

    scatter (sample_id in sample_ids) {
        String hap1_column_id = "~{sample_id}_hap1"
        String hap2_column_id = "~{sample_id}_hap2"
    }

    if (defined(samples_per_shard)) {
        Int samples_per_shard_value = select_first([samples_per_shard])
        Int num_shards = ceil((num_samples * 1.0) / samples_per_shard_value)

        scatter (shard_index in range(num_shards)) {
            Int shard_start = shard_index * samples_per_shard_value
            Int shard_end = if shard_start + samples_per_shard_value < num_samples then shard_start + samples_per_shard_value else num_samples

            scatter (offset in range(shard_end - shard_start)) {
                Int j = shard_start + offset
                File selected_combined_bed = combined_beds[j]
                File selected_combined_bed_idx = combined_bed_idxs[j]
                File selected_hap1_bed = hap1_beds[j]
                File selected_hap1_bed_idx = hap1_bed_idxs[j]
                File selected_hap2_bed = hap2_beds[j]
                File selected_hap2_bed_idx = hap2_bed_idxs[j]
                String selected_sample_id = sample_ids[j]
                String selected_hap1_column_id = hap1_column_id[j]
                String selected_hap2_column_id = hap2_column_id[j]
            }

            Array[String] shard_sample_ids = selected_sample_id
            Array[File] shard_combined_beds = selected_combined_bed
            Array[File] shard_combined_bed_idxs = selected_combined_bed_idx
            Array[File] shard_hap_beds = flatten([selected_hap1_bed, selected_hap2_bed])
            Array[File] shard_hap_bed_idxs = flatten([selected_hap1_bed_idx, selected_hap2_bed_idx])
            Array[String] shard_hap_column_ids = flatten([selected_hap1_column_id, selected_hap2_column_id])
        }
    }

    scatter (contig in contigs) {
        if (defined(samples_per_shard)) {
            Array[Array[File]] contig_shard_combined_beds = select_first([shard_combined_beds])
            Array[Array[File]] contig_shard_combined_bed_idxs = select_first([shard_combined_bed_idxs])
            Array[Array[String]] contig_shard_sample_ids = select_first([shard_sample_ids])
            Array[Array[File]] contig_shard_hap_beds = select_first([shard_hap_beds])
            Array[Array[File]] contig_shard_hap_bed_idxs = select_first([shard_hap_bed_idxs])
            Array[Array[String]] contig_shard_hap_column_ids = select_first([shard_hap_column_ids])

            scatter (i in range(length(contig_shard_combined_beds))) {
                call MergeMethylationShard as MergeCombinedShard {
                    input:
                        beds = contig_shard_combined_beds[i],
                        bed_idxs = contig_shard_combined_bed_idxs[i],
                        column_ids = contig_shard_sample_ids[i],
                        contig = contig,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.combined",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_merge_combined
                }

                call MergeMethylationShard as MergeHaplotypeShard {
                    input:
                        beds = contig_shard_hap_beds[i],
                        bed_idxs = contig_shard_hap_bed_idxs[i],
                        column_ids = contig_shard_hap_column_ids[i],
                        contig = contig,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.haplotype",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_merge_haplotype
                }
            }

            call JoinMethylationShards as JoinCombinedShards {
                input:
                    shard_beds = MergeCombinedShard.matrix_bed,
                    prefix = "~{prefix}.~{contig}.combined",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_join_combined
            }

            call JoinMethylationShards as JoinHaplotypeShards {
                input:
                    shard_beds = MergeHaplotypeShard.matrix_bed,
                    prefix = "~{prefix}.~{contig}.haplotype",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_join_haplotype
            }
        }

        if (!defined(samples_per_shard)) {
            call MergeMethylationShard as MergeCombinedContig {
                input:
                    beds = combined_beds,
                    bed_idxs = combined_bed_idxs,
                    column_ids = sample_ids,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.combined",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_combined
            }

            call MergeMethylationShard as MergeHaplotypeContig {
                input:
                    beds = flatten([hap1_beds, hap2_beds]),
                    bed_idxs = flatten([hap1_bed_idxs, hap2_bed_idxs]),
                    column_ids = flatten([hap1_column_id, hap2_column_id]),
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.haplotype",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_haplotype
            }
        }

        File contig_combined_bed = select_first([JoinCombinedShards.matrix_bed, MergeCombinedContig.matrix_bed])
        File contig_haplotype_bed = select_first([JoinHaplotypeShards.matrix_bed, MergeHaplotypeContig.matrix_bed])
    }

    output {
        Array[File] combined_methylation_beds = contig_combined_bed
        Array[File] haplotype_methylation_beds = contig_haplotype_bed
    }
}

task MergeMethylationShard {
    input {
        Array[File] beds
        Array[File] bed_idxs
        Array[String] column_ids
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import heapq
import subprocess

beds = "~{sep=',' beds}".split(',')
column_ids = "~{sep=',' column_ids}".split(',')

procs = {
    c: subprocess.Popen(['tabix', b, "~{contig}"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    for c, b in zip(column_ids, beds)
}

def next_site(proc):
    line = proc.stdout.readline()
    if not line:
        return None
    fields = line.rstrip('\n').split('\t')
    return int(fields[1]), int(fields[2]), fields[3]

heap = []
pending = {}
for c in column_ids:
    site = next_site(procs[c])
    if site is not None:
        start, end, score = site
        pending[c] = score
        heapq.heappush(heap, (start, end, c))

with open("~{prefix}.bed", "w") as out:
    out.write("\t".join(["#chrom", "start", "end"] + column_ids) + "\n")
    while heap:
        start, end, _ = heap[0]
        row_scores = {}
        while heap and heap[0][0] == start and heap[0][1] == end:
            _, _, c = heapq.heappop(heap)
            row_scores[c] = pending.pop(c)
            site = next_site(procs[c])
            if site is not None:
                next_start, next_end, next_score = site
                pending[c] = next_score
                heapq.heappush(heap, (next_start, next_end, c))
        row = ["~{contig}", str(start), str(end)] + [row_scores.get(c, ".") for c in column_ids]
        out.write("\t".join(row) + "\n")

for proc in procs.values():
    proc.stdout.close()
    proc.wait()
CODE

        bgzip ~{prefix}.bed
    >>>

    output {
        File matrix_bed = "~{prefix}.bed.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(beds, "GB")) + 20,
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

task JoinMethylationShards {
    input {
        Array[File] shard_beds
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import gzip
import heapq

shard_files = "~{sep=',' shard_beds}".split(',')
handles = [gzip.open(f, 'rt') for f in shard_files]
column_groups = [h.readline().rstrip('\n').split('\t')[3:] for h in handles]

def next_row(f):
    line = f.readline()
    if not line:
        return None
    fields = line.rstrip('\n').split('\t')
    return fields[0], int(fields[1]), int(fields[2]), fields[3:]

heap = []
pending = {}
for i, f in enumerate(handles):
    row = next_row(f)
    if row is not None:
        chrom, start, end, vals = row
        pending[i] = (chrom, vals)
        heapq.heappush(heap, (start, end, i))

all_columns = [c for group in column_groups for c in group]

with open("~{prefix}.bed", "w") as out:
    out.write("\t".join(["#chrom", "start", "end"] + all_columns) + "\n")
    while heap:
        start, end, _ = heap[0]
        row_values = {}
        chrom = None
        while heap and heap[0][0] == start and heap[0][1] == end:
            _, _, i = heapq.heappop(heap)
            chrom, vals = pending.pop(i)
            row_values[i] = vals
            next_row_val = next_row(handles[i])
            if next_row_val is not None:
                next_chrom, next_start, next_end, next_vals = next_row_val
                pending[i] = (next_chrom, next_vals)
                heapq.heappush(heap, (next_start, next_end, i))
        row = [chrom, str(start), str(end)]
        for i, group in enumerate(column_groups):
            row.extend(row_values.get(i, ["."] * len(group)))
        out.write("\t".join(row) + "\n")

for f in handles:
    f.close()
CODE

        bgzip ~{prefix}.bed
    >>>

    output {
        File matrix_bed = "~{prefix}.bed.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(shard_beds, "GB")) + 10,
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
