version 1.0

# Depth CNV calls from two callers (gCNV raw and the edge scan) -> breakpoint check (v8) -> hierarchical merge of overlapping
# calls of the two callers -> breakpoint check of the merged calls -> sites -> depth genotyping in the chosen samples ->
# breakpoint check of every genotyped carrier. The cohort RD matrix is read in place over gs:// (pysam, OAuth token from the
# metadata server). Scripts are passed as files (score_cnv_breakpoints.py, unify_caller_pairs.py, cluster_depth_sites.py,
# genotype_depth_sites.py).

import "Structs.wdl"

workflow DepthCnvMergeGenotype {
    input {
        File gcnv_calls
        File edge_calls
        File samples
        Array[String] regions
        String prefix

        Int min_size = 10000
        Float close_frac = 0.1
        Float site_ro = 0.8
        Float del_max = 1.5
        Float hom_del_max = 0.5
        Float dup_min = 2.5
        String score_args = ""

        String rd_matrix
        File median_coverage
        File segdup_bed
        File sharpness_model
        File score_script
        File unify_script
        File cluster_script
        File genotype_script

        String docker

        RuntimeAttr? runtime_attr_prepare_calls
        RuntimeAttr? runtime_attr_score_calls
        RuntimeAttr? runtime_attr_merge_callers
        RuntimeAttr? runtime_attr_combine_filtered
        RuntimeAttr? runtime_attr_cluster_sites
        RuntimeAttr? runtime_attr_genotype_sites
        RuntimeAttr? runtime_attr_finalize_genotypes
    }

    parameter_meta {
        gcnv_calls: "gCNV raw per-sample call table (sample chrom start end SVTYPE ...)"
        edge_calls: "Edge-scan per-sample call table (sample chrom start end SVTYPE ...)"
        samples: "Samples to call and genotype, one per line"
        regions: "Regions (chrom or chrom:start-end), one scoring / genotyping task each; their contigs are the contigs analysed"
        close_frac: "Merge: two callers' ends closer than this fraction of the group size are 'close'"
        site_ro: "Reciprocal overlap linking calls into one site"
        del_max: "DEL carrier when site CN < this"
        hom_del_max: "DEL 1/1 when site CN < this"
        dup_min: "DUP carrier when site CN > this"
        score_args: "Extra arguments for score_script (default v8 rule)"
        rd_matrix: "Cohort RD matrix ('#Chr Start End <samples>', bgzipped + tabix-indexed), gs:// path, read in place"
        median_coverage: "Sample median coverage, 2 rows (sample IDs / medians)"
        sharpness_model: "Sharpness model JSON"
    }

    call PrepareCalls {
        input:
            gcnv_calls = gcnv_calls,
            edge_calls = edge_calls,
            samples = samples,
            regions = regions,
            min_size = min_size,
            prefix = prefix,
            docker = docker,
            runtime_attr_override = runtime_attr_prepare_calls
    }

    scatter (region in regions) {
        call ScoreCalls as ScoreRawCalls {
            input:
                calls = PrepareCalls.union_calls,
                rd_matrix = rd_matrix,
                region = region,
                score_args = score_args,
                median_coverage = median_coverage,
                segdup_bed = segdup_bed,
                sharpness_model = sharpness_model,
                score_script = score_script,
                prefix = prefix + ".raw." + sub(region, "[:-]", "_"),
                docker = docker,
                runtime_attr_override = runtime_attr_score_calls
        }
    }

    call MergeCallers {
        input:
            gcnv_calls = PrepareCalls.gcnv_subset,
            edge_calls = PrepareCalls.edge_subset,
            decisions = ScoreRawCalls.scores,
            samples = samples,
            contigs = PrepareCalls.contigs,
            min_size = min_size,
            close_frac = close_frac,
            sharpness_model = sharpness_model,
            unify_script = unify_script,
            prefix = prefix,
            docker = docker,
            runtime_attr_override = runtime_attr_merge_callers
    }

    scatter (region in regions) {
        call ScoreCalls as ScoreMergedCalls {
            input:
                calls = MergeCallers.merged_calls,
                rd_matrix = rd_matrix,
                region = region,
                score_args = score_args,
                median_coverage = median_coverage,
                segdup_bed = segdup_bed,
                sharpness_model = sharpness_model,
                score_script = score_script,
                prefix = prefix + ".merged." + sub(region, "[:-]", "_"),
                docker = docker,
                runtime_attr_override = runtime_attr_score_calls
        }
    }

    call CombineFiltered {
        input:
            union_calls = PrepareCalls.union_calls,
            raw_scores = ScoreRawCalls.scores,
            members = MergeCallers.members,
            merged_calls = MergeCallers.merged_calls,
            merged_scores = ScoreMergedCalls.scores,
            prefix = prefix,
            docker = docker,
            runtime_attr_override = runtime_attr_combine_filtered
    }

    call ClusterSites {
        input:
            calls = CombineFiltered.filtered_calls,
            site_ro = site_ro,
            cluster_script = cluster_script,
            prefix = prefix,
            docker = docker,
            runtime_attr_override = runtime_attr_cluster_sites
    }

    scatter (i in range(length(regions))) {
        call GenotypeSites {
            input:
                sites = ClusterSites.sites,
                samples = samples,
                rd_matrix = rd_matrix,
                region = regions[i],
                del_max = del_max,
                hom_del_max = hom_del_max,
                dup_min = dup_min,
                median_coverage = median_coverage,
                genotype_script = genotype_script,
                prefix = prefix + ".genotyped." + sub(regions[i], "[:-]", "_"),
                docker = docker,
                runtime_attr_override = runtime_attr_genotype_sites
        }

        call ScoreCalls as ScoreCarriers {
            input:
                calls = GenotypeSites.carriers,
                rd_matrix = rd_matrix,
                region = regions[i],
                score_args = score_args,
                median_coverage = median_coverage,
                segdup_bed = segdup_bed,
                sharpness_model = sharpness_model,
                score_script = score_script,
                prefix = prefix + ".carriers." + sub(regions[i], "[:-]", "_"),
                docker = docker,
                runtime_attr_override = runtime_attr_score_calls
        }
    }

    call FinalizeGenotypes {
        input:
            genotypes = GenotypeSites.genotypes,
            carrier_scores = ScoreCarriers.scores,
            prefix = prefix,
            docker = docker,
            runtime_attr_override = runtime_attr_finalize_genotypes
    }

    output {
        File union_calls = PrepareCalls.union_calls
        Array[File] raw_scores = ScoreRawCalls.scores
        File merge_groups = MergeCallers.groups
        File merged_calls = MergeCallers.merged_calls
        Array[File] merged_scores = ScoreMergedCalls.scores
        File filtered_calls = CombineFiltered.filtered_calls
        File sites = ClusterSites.sites
        Array[File] carrier_scores = ScoreCarriers.scores
        File final_genotypes = FinalizeGenotypes.final_genotypes
        File final_calls = FinalizeGenotypes.final_calls
        File summary = FinalizeGenotypes.summary
    }
}

task PrepareCalls {
    input {
        File gcnv_calls
        File edge_calls
        File samples
        Array[String] regions
        Int min_size
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import csv
import gzip

samples = set(open("~{samples}").read().split())
contigs = sorted({r.split(":")[0] for r in "~{sep=',' regions}".split(",")})


def opn(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


out = {}
union = []
for label, path, name in (("gCNV", "~{gcnv_calls}", "gcnv"), ("edge", "~{edge_calls}", "edge")):
    rows = []
    with opn(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            s, e = int(r["start"]), int(r["end"])
            if r["sample"] in samples and r["chrom"] in contigs and e - s >= ~{min_size}:
                rows.append((r["sample"], r["chrom"], s, e, r["SVTYPE"]))
    with gzip.open("~{prefix}." + name + ".calls.tsv.gz", "wt") as o:
        o.write("sample\tchrom\tstart\tend\tSVTYPE\tSVLEN\tFILTER\tCN\n")
        for x in rows:
            o.write(f"{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]}\t{x[3] - x[2]}\tPASS\t\n")
    union += [x + (label,) for x in rows]
    print(label, len(rows), "calls")
with gzip.open("~{prefix}.union.calls.tsv.gz", "wt") as o:
    o.write("sample\tchrom\tstart\tend\tSVTYPE\tSVLEN\tFILTER\tCN\tcaller\n")
    for x in sorted(union, key=lambda y: (y[1], y[2], y[3], y[0])):
        o.write(f"{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]}\t{x[3] - x[2]}\tPASS\t\t{x[5]}\n")
open("contigs.txt", "w").write(",".join(contigs))
CODE
    >>>

    output {
        File gcnv_subset = "~{prefix}.gcnv.calls.tsv.gz"
        File edge_subset = "~{prefix}.edge.calls.tsv.gz"
        File union_calls = "~{prefix}.union.calls.tsv.gz"
        String contigs = read_string("contigs.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(gcnv_calls, "GB") + size(edge_calls, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task ScoreCalls {
    input {
        File calls
        String rd_matrix
        String region
        String score_args
        File median_coverage
        File segdup_bed
        File sharpness_model
        File score_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(curl -sS -H "Metadata-Flavor: Google" \
            http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
            | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')

        python3 ~{score_script} \
            --calls ~{calls} \
            --rd ~{rd_matrix} \
            --medians ~{median_coverage} \
            --segdup ~{segdup_bed} \
            --model ~{sharpness_model} \
            --region ~{region} \
            ~{score_args} \
            --out ~{prefix}.decisions.tsv.gz
    >>>

    output {
        File scores = "~{prefix}.decisions.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(size(calls, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task MergeCallers {
    input {
        File gcnv_calls
        File edge_calls
        Array[File] decisions
        File samples
        String contigs
        Int min_size
        Float close_frac
        File sharpness_model
        File unify_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 ~{unify_script} \
            --calls-a gCNV=~{gcnv_calls} \
            --calls-b edge=~{edge_calls} \
            --decisions ~{sep=' ' decisions} \
            --model ~{sharpness_model} \
            --samples ~{samples} \
            --chroms ~{contigs} \
            --min-size ~{min_size} \
            --min-diff -1 \
            --close-frac ~{close_frac} \
            --out-groups ~{prefix}.merge_groups.tsv \
            --out-calls ~{prefix}.merged.calls.tsv.gz \
            --out-members ~{prefix}.merge_members.tsv.gz
    >>>

    output {
        File groups = "~{prefix}.merge_groups.tsv"
        File merged_calls = "~{prefix}.merged.calls.tsv.gz"
        File members = "~{prefix}.merge_members.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(decisions, "GB") + size(edge_calls, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task CombineFiltered {
    input {
        File union_calls
        Array[File] raw_scores
        File members
        File merged_calls
        Array[File] merged_scores
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import csv
import gzip


def rows(path):
    with gzip.open(path, "rt") as f:
        yield from csv.DictReader(f, delimiter="\t")


def key(r):
    return (r["sample"], r["chrom"], int(r["start"]), int(r["end"]), r["SVTYPE"])


caller = {}
for r in rows("~{union_calls}"):
    caller.setdefault(key(r), set()).add(r["caller"])
grouped = {key(r) for r in rows("~{members}")}
conf = {key(r): r["confidence"] for r in rows("~{merged_calls}")}
out = {}
n_raw = n_merged = 0
for path in "~{sep=',' raw_scores}".split(","):
    for r in rows(path):
        k = key(r)
        if k in grouped or r["kept"] != "True":
            continue
        n_raw += 1
        nk = (r["sample"], r["chrom"], int(r["new_start"]), int(r["new_end"]), r["SVTYPE"])
        out[nk] = ("standard", "+".join(sorted(caller.get(k, {"?"}))))
for path in "~{sep=',' merged_scores}".split(","):
    for r in rows(path):
        if r["kept"] != "True":
            continue
        n_merged += 1
        nk = (r["sample"], r["chrom"], int(r["new_start"]), int(r["new_end"]), r["SVTYPE"])
        out[nk] = (conf.get(key(r), "standard"), "merged")
with gzip.open("~{prefix}.filtered.calls.tsv.gz", "wt") as o:
    o.write("sample\tchrom\tstart\tend\tSVTYPE\tSVLEN\tFILTER\tCN\tconfidence\tsource\n")
    for k in sorted(out, key=lambda x: (x[1], x[2], x[3], x[0])):
        o.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{k[4]}\t{k[3] - k[2]}\tPASS\t\t{out[k][0]}\t{out[k][1]}\n")
print(f"kept unmerged calls {n_raw}, kept merged calls {n_merged}, filtered call set {len(out)}")
CODE
    >>>

    output {
        File filtered_calls = "~{prefix}.filtered.calls.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(raw_scores, "GB") + size(union_calls, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task ClusterSites {
    input {
        File calls
        Float site_ro
        File cluster_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 ~{cluster_script} \
            --calls ~{calls} \
            --ro ~{site_ro} \
            --prefix ~{prefix} \
            --out ~{prefix}.sites.tsv
    >>>

    output {
        File sites = "~{prefix}.sites.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(calls, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task GenotypeSites {
    input {
        File sites
        File samples
        String rd_matrix
        String region
        Float del_max
        Float hom_del_max
        Float dup_min
        File median_coverage
        File genotype_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(curl -sS -H "Metadata-Flavor: Google" \
            http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
            | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')

        python3 ~{genotype_script} \
            --sites ~{sites} \
            --rd ~{rd_matrix} \
            --medians ~{median_coverage} \
            --samples ~{samples} \
            --region ~{region} \
            --del-max ~{del_max} \
            --hom-del-max ~{hom_del_max} \
            --dup-min ~{dup_min} \
            --out-genotypes ~{prefix}.genotypes.tsv.gz \
            --out-carriers ~{prefix}.carriers.tsv.gz
    >>>

    output {
        File genotypes = "~{prefix}.genotypes.tsv.gz"
        File carriers = "~{prefix}.carriers.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(sites, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task FinalizeGenotypes {
    input {
        Array[File] genotypes
        Array[File] carrier_scores
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import csv
import gzip


def rows(path):
    with gzip.open(path, "rt") as f:
        yield from csv.DictReader(f, delimiter="\t")


kept = {}
for path in "~{sep=',' carrier_scores}".split(","):
    for r in rows(path):
        kept[(r["sample"], r["chrom"], r["start"], r["end"], r["SVTYPE"])] = (r["kept"], r["fail_reason"], r["new_start"], r["new_end"])
n = {"sites": set(), "carriers": 0, "kept": 0}
with gzip.open("~{prefix}.final_genotypes.tsv.gz", "wt") as og, gzip.open("~{prefix}.final.calls.tsv.gz", "wt") as oc:
    oc.write("sample\tchrom\tstart\tend\tSVTYPE\tSVLEN\tFILTER\tCN\tsite_id\tconfidence\tdiscovered\tGT\n")
    w = None
    for path in "~{sep=',' genotypes}".split(","):
        for r in rows(path):
            if w is None:
                w = csv.DictWriter(og, fieldnames=list(r) + ["carrier_kept", "fail_reason", "new_start", "new_end"], delimiter="\t",
                                   lineterminator="\n")
                w.writeheader()
            n["sites"].add(r["site_id"])
            k = kept.get((r["sample"], r["chrom"], r["start"], r["end"], r["SVTYPE"]))
            if r["GT"] == "0/0":
                r.update(carrier_kept=".", fail_reason=".", new_start=".", new_end=".")
            else:
                n["carriers"] += 1
                r.update(carrier_kept=k[0] if k else "NA", fail_reason=k[1] if k else "NA", new_start=k[2] if k else ".",
                         new_end=k[3] if k else ".")
                if k and k[0] == "True":
                    n["kept"] += 1
                    oc.write(f"{r['sample']}\t{r['chrom']}\t{k[2]}\t{k[3]}\t{r['SVTYPE']}\t{int(k[3]) - int(k[2])}\tPASS\t{r['CN']}\t"
                             f"{r['site_id']}\t{r['confidence']}\t{r['discovered']}\t{r['GT']}\n")
            w.writerow(r)
with open("~{prefix}.summary.tsv", "w") as o:
    o.write(f"sites\t{len(n['sites'])}\ncarriers\t{n['carriers']}\ncarriers_kept\t{n['kept']}\n")
CODE
    >>>

    output {
        File final_genotypes = "~{prefix}.final_genotypes.tsv.gz"
        File final_calls = "~{prefix}.final.calls.tsv.gz"
        File summary = "~{prefix}.summary.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(genotypes, "GB") + size(carrier_scores, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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
