version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FillBackbonePhasedGenotypes {
    input {
        File backbone_phased_vcf
        File backbone_phased_vcf_idx
        File backbone_phased_notrgt_vcf
        File backbone_phased_notrgt_vcf_idx
        String contig
        String prefix

        Int? shard_bin_size

        String utils_docker

        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_subset_backbone
        RuntimeAttr? runtime_attr_subset_notrgt
        RuntimeAttr? runtime_attr_pull_phase
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_aggregate_stats
    }

    if (defined(shard_bin_size)) {
        call Helpers.CreateContigShards {
            input:
                vcfs = [backbone_phased_vcf, backbone_phased_notrgt_vcf],
                vcf_idxs = [backbone_phased_vcf_idx, backbone_phased_notrgt_vcf_idx],
                contig = contig,
                shard_bin_size = select_first([shard_bin_size]),
                prefix = "~{prefix}.shards",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_shards
        }

        scatter (i in range(length(CreateContigShards.shard_regions))) {
            String shard_region = CreateContigShards.shard_regions[i]

            call Helpers.SubsetVcfToRegion as SubsetBackbone {
                input:
                    vcf = backbone_phased_vcf,
                    vcf_idx = backbone_phased_vcf_idx,
                    region = shard_region,
                    prefix = "~{prefix}.shard_~{i}.backbone",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_backbone
            }

            call Helpers.SubsetVcfToRegion as SubsetNotrgt {
                input:
                    vcf = backbone_phased_notrgt_vcf,
                    vcf_idx = backbone_phased_notrgt_vcf_idx,
                    region = shard_region,
                    prefix = "~{prefix}.shard_~{i}.notrgt",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_notrgt
            }

            call PullPhasedGenotypes {
                input:
                    backbone_phased_vcf = SubsetBackbone.subset_vcf,
                    backbone_phased_vcf_idx = SubsetBackbone.subset_vcf_idx,
                    backbone_phased_notrgt_vcf = SubsetNotrgt.subset_vcf,
                    backbone_phased_notrgt_vcf_idx = SubsetNotrgt.subset_vcf_idx,
                    prefix = "~{prefix}.shard_~{i}.pulled",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_pull_phase
            }
        }

        call Helpers.ConcatVcfs {
            input:
                vcfs = PullPhasedGenotypes.pulled_vcf,
                vcf_idxs = PullPhasedGenotypes.pulled_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.pulled",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    if (!defined(shard_bin_size)) {
        call PullPhasedGenotypes as PullPhasedGenotypesNoSharding {
            input:
                backbone_phased_vcf = backbone_phased_vcf,
                backbone_phased_vcf_idx = backbone_phased_vcf_idx,
                backbone_phased_notrgt_vcf = backbone_phased_notrgt_vcf,
                backbone_phased_notrgt_vcf_idx = backbone_phased_notrgt_vcf_idx,
                prefix = "~{prefix}.pulled",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_pull_phase
        }
    }

    Array[File] stats_shards = if defined(shard_bin_size)
        then select_first([PullPhasedGenotypes.counts_tsv])
        else [select_first([PullPhasedGenotypesNoSharding.counts_tsv])]

    call AggregatePhasePullStats {
        input:
            counts_tsvs = stats_shards,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_aggregate_stats
    }

    output {
        File backbone_merged_vcf = select_first([ConcatVcfs.concat_vcf, PullPhasedGenotypesNoSharding.pulled_vcf])
        File backbone_merged_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, PullPhasedGenotypesNoSharding.pulled_vcf_idx])
        File backbone_merged_tsv = AggregatePhasePullStats.stats_tsv
    }
}

task PullPhasedGenotypes {
    input {
        File backbone_phased_vcf
        File backbone_phased_vcf_idx
        File backbone_phased_notrgt_vcf
        File backbone_phased_notrgt_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [[ "~{backbone_phased_notrgt_vcf_idx}" != "~{backbone_phased_notrgt_vcf}.tbi" ]]; then
            ln -sf "~{backbone_phased_notrgt_vcf_idx}" "~{backbone_phased_notrgt_vcf}.tbi"
        fi

        python3 <<CODE
import pysam

backbone_in = pysam.VariantFile("~{backbone_phased_vcf}")
notrgt_in = pysam.VariantFile("~{backbone_phased_notrgt_vcf}")

notrgt_samples = set(notrgt_in.header.samples)
samples = list(backbone_in.header.samples)

n_het = {s: 0 for s in samples}
n_unphased = {s: 0 for s in samples}
n_unphased_post_pull = {s: 0 for s in samples}

out_header = backbone_in.header.copy()
out = pysam.VariantFile("~{prefix}.vcf.gz", "w", header=out_header)

def alleles_key(rec):
    ref = rec.ref.upper() if rec.ref else rec.ref
    alts = tuple(a.upper() for a in rec.alts) if rec.alts else ()
    return (rec.chrom, rec.pos, ref, alts)

for rec in backbone_in:
    rec.translate(out_header)
    match = None
    rec_key = alleles_key(rec)
    for cand in notrgt_in.fetch(rec.chrom, rec.start, rec.stop):
        if alleles_key(cand) == rec_key:
            match = cand
            break

    for sample in samples:
        gt_data = rec.samples[sample]
        gt = gt_data.get("GT")
        if gt is None or None in gt or len(gt) != 2:
            continue
        if gt[0] == gt[1]:
            continue
        n_het[sample] += 1
        if gt_data.phased:
            continue
        n_unphased[sample] += 1

        pulled = False
        if match is not None and sample in notrgt_samples:
            notrgt_gt_data = match.samples[sample]
            notrgt_gt = notrgt_gt_data.get("GT")
            if notrgt_gt is not None and notrgt_gt_data.phased:
                rec.samples[sample]["GT"] = notrgt_gt
                rec.samples[sample].phased = True
                notrgt_ps = notrgt_gt_data.get("PS")
                if notrgt_ps is not None:
                    rec.samples[sample]["PS"] = notrgt_ps
                pulled = True
        if not pulled:
            n_unphased_post_pull[sample] += 1

    out.write(rec)

out.close()

with open("~{prefix}.tsv", "w") as out_tsv:
    out_tsv.write("SAMPLE\tN_HET\tN_HET_UNPHASED\tN_HET_UNPHASED_POST_PULL\n")
    for sample in samples:
        out_tsv.write(f"{sample}\t{n_het[sample]}\t{n_unphased[sample]}\t{n_unphased_post_pull[sample]}\n")
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File pulled_vcf = "~{prefix}.vcf.gz"
        File pulled_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File counts_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(backbone_phased_vcf, "GB") + size(backbone_phased_notrgt_vcf, "GB")) + 10,
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

task AggregatePhasePullStats {
    input {
        Array[File] counts_tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
n_het = {}
n_unphased = {}
n_unphased_post_pull = {}
order = []

counts_tsv_paths = []
with open("~{write_lines(counts_tsvs)}") as f:
    for line in f:
        path = line.strip()
        if path:
            counts_tsv_paths.append(path)

for path in counts_tsv_paths:
    with open(path) as fh:
        fh.readline()
        for line in fh:
            sample, het, unphased, post_pull = line.rstrip("\n").split("\t")
            if sample not in n_het:
                order.append(sample)
                n_het[sample] = 0
                n_unphased[sample] = 0
                n_unphased_post_pull[sample] = 0
            n_het[sample] += int(het)
            n_unphased[sample] += int(unphased)
            n_unphased_post_pull[sample] += int(post_pull)

def proportion(numerator, denominator):
    return f"{(numerator / denominator):.6f}" if denominator else "0.000000"

with open("~{prefix}.tsv", "w") as out:
    out.write(
        "SAMPLE\tN_HET\tN_HET_UNPHASED\tPROP_HET_UNPHASED\t"
        "N_HET_UNPHASED_POST_PULL\tPROP_HET_UNPHASED_POST_PULL\n"
    )
    for sample in order:
        het = n_het[sample]
        unphased = n_unphased[sample]
        post_pull = n_unphased_post_pull[sample]
        out.write(
            f"{sample}\t{het}\t{unphased}\t{proportion(unphased, het)}\t"
            f"{post_pull}\t{proportion(post_pull, het)}\n"
        )
CODE
    >>>

    output {
        File stats_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
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
