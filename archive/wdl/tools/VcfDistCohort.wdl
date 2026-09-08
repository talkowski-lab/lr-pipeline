version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "BackbonePhase.wdl"
import "VcfDist.wdl"

workflow VcfDistCohort {
    input {
        Array[File] eval_vcfs
        Array[File] eval_vcf_idxs
        Array[File] truth_vcfs
        Array[File] truth_vcf_idxs
        File ref_fa
        Array[String] contigs
        String prefix

        Array[String]? subset_samples
        String? vcfdist_args

        String utils_docker
        String vcfdist_docker

        RuntimeAttr? runtime_attr_subset_assignment_samples
        RuntimeAttr? runtime_attr_assign_samples
        RuntimeAttr? runtime_attr_extract_assigned_samples
        RuntimeAttr? runtime_attr_subset_truth_contig
        RuntimeAttr? runtime_attr_subset_truth_samples
        RuntimeAttr? runtime_attr_subset_eval_samples
        RuntimeAttr? runtime_attr_extract_eval_sample
        RuntimeAttr? runtime_attr_extract_truth_sample
        RuntimeAttr? runtime_attr_vcfdist
        RuntimeAttr? runtime_attr_annotate_vcfdist_results
        RuntimeAttr? runtime_attr_aggregate_results
    }

    if (defined(subset_samples)) {
        call Helpers.SubsetVcfToSamples as SubsetAssignmentVcf {
            input:
                vcf = eval_vcfs[0],
                vcf_idx = eval_vcf_idxs[0],
                samples = select_first([subset_samples]),
                prefix = "~{prefix}.assignment_subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_assignment_samples
        }
    }

    File assignment_vcf = select_first([SubsetAssignmentVcf.subset_vcf, eval_vcfs[0]])
    File assignment_vcf_idx = select_first([SubsetAssignmentVcf.subset_vcf_idx, eval_vcf_idxs[0]])

    call BackbonePhase.AssignSamplesToBaseVcfs {
        input:
            vcf = assignment_vcf,
            vcf_idx = assignment_vcf_idx,
            base_vcfs = truth_vcfs,
            prefix = "~{prefix}.assignments",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_assign_samples
    }

    scatter (j in range(length(truth_vcfs))) {
        call ExtractAssignedSamples {
            input:
                assignment_tsv = AssignSamplesToBaseVcfs.assignment_tsv,
                truth_vcf_index = j,
                prefix = "~{prefix}.truth_~{j}.assigned_samples",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_extract_assigned_samples
        }
    }

    scatter (i in range(length(eval_vcfs))) {
        scatter (j in range(length(truth_vcfs))) {
            call Helpers.SubsetVcfToSamples as SubsetTruthSamples {
                input:
                    vcf = truth_vcfs[j],
                    vcf_idx = truth_vcf_idxs[j],
                    samples = ExtractAssignedSamples.assigned_samples[j],
                    prefix = "~{prefix}.truth_~{j}.~{contigs[i]}.samples",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_truth_samples
            }

            call Helpers.SubsetVcfToSamples as SubsetEvalSamples {
                input:
                    vcf = eval_vcfs[i],
                    vcf_idx = eval_vcf_idxs[i],
                    samples = ExtractAssignedSamples.assigned_samples[j],
                    prefix = "~{prefix}.eval_~{contigs[i]}.truth_~{j}.samples",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_eval_samples
            }

            scatter (sample in ExtractAssignedSamples.assigned_samples[j]) {
                call Helpers.ExtractSample as ExtractEvalSample {
                    input:
                        vcf = SubsetEvalSamples.subset_vcf,
                        vcf_idx = SubsetEvalSamples.subset_vcf_idx,
                        sample = sample,
                        prefix = "~{prefix}.eval_~{contigs[i]}.truth_~{j}.~{sample}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_extract_eval_sample
                }

                call Helpers.ExtractSample as ExtractTruthSample {
                    input:
                        vcf = SubsetTruthSamples.subset_vcf,
                        vcf_idx = SubsetTruthSamples.subset_vcf_idx,
                        sample = sample,
                        prefix = "~{prefix}.truth_~{j}.~{contigs[i]}.~{sample}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_extract_truth_sample
                }

                call VcfDist.RunVcfDist {
                    input:
                        vcf_eval = ExtractEvalSample.subset_vcf,
                        vcf_eval_idx = ExtractEvalSample.subset_vcf_idx,
                        vcf_truth = ExtractTruthSample.subset_vcf,
                        vcf_truth_idx = ExtractTruthSample.subset_vcf_idx,
                        ref_fa = ref_fa,
                        vcfdist_args = select_first([vcfdist_args, ""]),
                        prefix = "~{prefix}.~{contigs[i]}.truth_~{j}.~{sample}",
                        docker = vcfdist_docker,
                        runtime_attr_override = runtime_attr_vcfdist
                }

                call AnnotateVcfDistResults {
                    input:
                        phasing_summary_tsv = RunVcfDist.phasing_summary_tsv,
                        precision_recall_summary_tsv = RunVcfDist.precision_recall_summary_tsv,
                        precision_recall_tsv = RunVcfDist.precision_recall_tsv,
                        switchflips_tsv = RunVcfDist.switchflips_tsv,
                        phase_blocks_tsv = RunVcfDist.phase_blocks_tsv,
                        contig = contigs[i],
                        sample = sample,
                        prefix = "~{prefix}.~{contigs[i]}.truth_~{j}.~{sample}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_annotate_vcfdist_results
                }
            }

            Array[File] shard_phasing_summary_tsvs = AnnotateVcfDistResults.annotated_phasing_summary_tsv
            Array[File] shard_precision_recall_summary_tsvs = AnnotateVcfDistResults.annotated_precision_recall_summary_tsv
            Array[File] shard_precision_recall_tsvs = AnnotateVcfDistResults.annotated_precision_recall_tsv
            Array[File] shard_switchflips_tsvs = AnnotateVcfDistResults.annotated_switchflips_tsv
            Array[File] shard_phase_blocks_tsvs = AnnotateVcfDistResults.annotated_phase_blocks_tsv
        }

        Array[File] contig_phasing_summary_tsvs = flatten(shard_phasing_summary_tsvs)
        Array[File] contig_precision_recall_summary_tsvs = flatten(shard_precision_recall_summary_tsvs)
        Array[File] contig_precision_recall_tsvs = flatten(shard_precision_recall_tsvs)
        Array[File] contig_switchflips_tsvs = flatten(shard_switchflips_tsvs)
        Array[File] contig_phase_blocks_tsvs = flatten(shard_phase_blocks_tsvs)
    }

    Array[File] all_phasing_summary_tsvs = flatten(contig_phasing_summary_tsvs)
    Array[File] all_precision_recall_summary_tsvs = flatten(contig_precision_recall_summary_tsvs)
    Array[File] all_precision_recall_tsvs = flatten(contig_precision_recall_tsvs)
    Array[File] all_switchflips_tsvs = flatten(contig_switchflips_tsvs)
    Array[File] all_phase_blocks_tsvs = flatten(contig_phase_blocks_tsvs)

    call AggregateVcfDistOutputs {
        input:
            phasing_summary_tsvs = all_phasing_summary_tsvs,
            precision_recall_summary_tsvs = all_precision_recall_summary_tsvs,
            precision_recall_tsvs = all_precision_recall_tsvs,
            switchflips_tsvs = all_switchflips_tsvs,
            phase_blocks_tsvs = all_phase_blocks_tsvs,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_aggregate_results
    }

    output {
        File vcfdist_phasing_summary_tsv = AggregateVcfDistOutputs.phasing_summary_tsv
        File vcfdist_precision_recall_summary_tsv = AggregateVcfDistOutputs.precision_recall_summary_tsv
        File vcfdist_precision_recall_tsv = AggregateVcfDistOutputs.precision_recall_tsv
        File vcfdist_switchflips_tsv = AggregateVcfDistOutputs.switchflips_tsv
        File vcfdist_phase_blocks_tsv = AggregateVcfDistOutputs.phase_blocks_tsv
        File vcfdist_missing_samples = AssignSamplesToBaseVcfs.missing_samples
    }
}

task ExtractAssignedSamples {
    input {
        File assignment_tsv
        Int truth_vcf_index
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
samples = []
with open("~{assignment_tsv}") as handle:
    for line in handle:
        sample, base_idx = line.rstrip("\n").split("\t")
        if int(base_idx) == ~{truth_vcf_index}:
            samples.append(sample)

samples = sorted(samples)
with open("~{prefix}.samples.txt", "w") as out:
    for sample in samples:
        out.write(sample + "\n")
CODE
    >>>

    output {
        File samples_file = "~{prefix}.samples.txt"
        Array[String] assigned_samples = read_lines("~{prefix}.samples.txt")
        Int sample_count = length(read_lines("~{prefix}.samples.txt"))
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10,
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

task AnnotateVcfDistResults {
    input {
        File phasing_summary_tsv
        File precision_recall_summary_tsv
        File precision_recall_tsv
        File switchflips_tsv
        File phase_blocks_tsv
        String contig
        String sample
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import csv


def annotate_tsv(input_path, output_path, contig, sample):
    with open(input_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        orig_fields = list(reader.fieldnames or [])
    has_contig_col = "CONTIG" in orig_fields
    if has_contig_col:
        cols = ["SAMPLE"] + orig_fields
    else:
        cols = ["CONTIG", "SAMPLE"] + orig_fields
    with open(output_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for row in rows:
            row["SAMPLE"] = sample
            if not has_contig_col:
                row["CONTIG"] = contig
            out.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")


contig = "~{contig}"
sample = "~{sample}"
prefix = "~{prefix}"

annotate_tsv("~{phasing_summary_tsv}", f"{prefix}.phasing_summary.tsv", contig, sample)
annotate_tsv("~{precision_recall_summary_tsv}", f"{prefix}.precision_recall_summary.tsv", contig, sample)
annotate_tsv("~{precision_recall_tsv}", f"{prefix}.precision_recall.tsv", contig, sample)
annotate_tsv("~{switchflips_tsv}", f"{prefix}.switchflips.tsv", contig, sample)
annotate_tsv("~{phase_blocks_tsv}", f"{prefix}.phase_blocks.tsv", contig, sample)
CODE
    >>>

    output {
        File annotated_phasing_summary_tsv = "~{prefix}.phasing_summary.tsv"
        File annotated_precision_recall_summary_tsv = "~{prefix}.precision_recall_summary.tsv"
        File annotated_precision_recall_tsv = "~{prefix}.precision_recall.tsv"
        File annotated_switchflips_tsv = "~{prefix}.switchflips.tsv"
        File annotated_phase_blocks_tsv = "~{prefix}.phase_blocks.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10,
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

task AggregateVcfDistOutputs {
    input {
        Array[File] phasing_summary_tsvs
        Array[File] precision_recall_summary_tsvs
        Array[File] precision_recall_tsvs
        Array[File] switchflips_tsvs
        Array[File] phase_blocks_tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import csv


def contig_sort_key(contig):
    value = contig[3:] if contig.lower().startswith("chr") else contig
    lower_value = value.lower()
    if lower_value.isdigit():
        return (0, int(lower_value), contig)
    if lower_value == "x":
        return (1, 23, contig)
    if lower_value == "y":
        return (1, 24, contig)
    if lower_value in ("m", "mt"):
        return (1, 25, contig)
    return (2, lower_value, contig)


def aggregate_tsv(manifest_file, output_file):
    rows = []
    header = None
    with open(manifest_file) as m:
        for path in m:
            path = path.strip()
            if not path:
                continue
            with open(path) as f:
                reader = csv.DictReader(f, delimiter="\t")
                if header is None and reader.fieldnames:
                    header = list(reader.fieldnames)
                rows.extend(reader)
    if not header:
        header = ["CONTIG", "SAMPLE"]
    rows.sort(key=lambda r: (contig_sort_key(r.get("CONTIG", "")), r.get("SAMPLE", "")))
    with open(output_file, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(str(row.get(c, "")) for c in header) + "\n")


prefix = "~{prefix}"

aggregate_tsv("~{write_lines(phasing_summary_tsvs)}", f"{prefix}.phasing_summary.tsv")
aggregate_tsv("~{write_lines(precision_recall_summary_tsvs)}", f"{prefix}.precision_recall_summary.tsv")
aggregate_tsv("~{write_lines(precision_recall_tsvs)}", f"{prefix}.precision_recall.tsv")
aggregate_tsv("~{write_lines(switchflips_tsvs)}", f"{prefix}.switchflips.tsv")
aggregate_tsv("~{write_lines(phase_blocks_tsvs)}", f"{prefix}.phase_blocks.tsv")
CODE
    >>>

    output {
        File phasing_summary_tsv = "~{prefix}.phasing_summary.tsv"
        File precision_recall_summary_tsv = "~{prefix}.precision_recall_summary.tsv"
        File precision_recall_tsv = "~{prefix}.precision_recall.tsv"
        File switchflips_tsv = "~{prefix}.switchflips.tsv"
        File phase_blocks_tsv = "~{prefix}.phase_blocks.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(phasing_summary_tsvs, "GiB")) + ceil(size(precision_recall_summary_tsvs, "GiB")) + ceil(size(precision_recall_tsvs, "GiB")) + ceil(size(switchflips_tsvs, "GiB")) + ceil(size(phase_blocks_tsvs, "GiB")) + 15,
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

