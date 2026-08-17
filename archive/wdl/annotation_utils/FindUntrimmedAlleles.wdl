version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FindUntrimmedAlleles {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        String utils_docker

        RuntimeAttr? runtime_attr_subset_contig
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_find_untrimmed
        RuntimeAttr? runtime_attr_concat_shard_vcfs
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_contig
            }
        }

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] shard_vcfs = select_first([ShardVcfByRecords.shards, [contig_vcf]])
        Array[File] shard_vcf_idxs = select_first([ShardVcfByRecords.shard_idxs, [contig_vcf_idx]])

        scatter (i in range(length(shard_vcfs))) {
            call FindUntrimmedInShard {
                input:
                    vcf = shard_vcfs[i],
                    vcf_idx = shard_vcf_idxs[i],
                    prefix = "~{prefix}.~{contig}.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_find_untrimmed
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatVcfs as ConcatShardVcfs {
                input:
                    vcfs = FindUntrimmedInShard.subset_vcf,
                    vcf_idxs = FindUntrimmedInShard.subset_vcf_idx,
                    allow_overlaps = true,
                    naive = false,
                    prefix = "~{prefix}.~{contig}.untrimmed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shard_vcfs
            }
        }

        File contig_subset_vcf = select_first([ConcatShardVcfs.concat_vcf, FindUntrimmedInShard.subset_vcf[0]])
        File contig_subset_vcf_idx = select_first([ConcatShardVcfs.concat_vcf_idx, FindUntrimmedInShard.subset_vcf_idx[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = contig_subset_vcf,
                vcf_idxs = contig_subset_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.untrimmed",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcfs
        }
    }

    output {
        File untrimmed_vcf = select_first([ConcatVcfs.concat_vcf, contig_subset_vcf[0]])
        File untrimmed_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, contig_subset_vcf_idx[0]])
    }
}

task FindUntrimmedInShard {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<EOF
import pysam
import sys

vcf_in = pysam.VariantFile("~{vcf}")
if "original_ID" not in vcf_in.header.info:
    vcf_in.header.add_line(
        '##INFO=<ID=original_ID,Number=1,Type=String,'
        'Description="Variant ID prior to REF/ALT trimming">'
    )

vcf_out = pysam.VariantFile("unsorted.vcf.gz", "wz", header=vcf_in.header)

for rec in vcf_in:
    if rec.info.get("allele_type") == "trv":
        continue
    if rec.alts is None or rec.id is None:
        continue
    if len(rec.alts) != 1:
        sys.stderr.write(f"Skipping multi-allelic record {rec.id} at {rec.chrom}:{rec.pos}\n")
        continue

    alt = rec.alts[0]
    if alt is None or alt.startswith("<") or "*" in alt:
        continue

    ref, a = rec.ref, alt
    while len(ref) > 1 and len(a) > 1 and ref[-1] == a[-1]:
        ref = ref[:-1]
        a = a[:-1]
    pre = 0
    while len(ref) > 1 and len(a) > 1 and ref[0] == a[0]:
        ref = ref[1:]
        a = a[1:]
        pre += 1

    if ref == rec.ref and a == alt:
        continue

    original_id = rec.id
    new_rec = rec.copy()
    new_rec.pos = rec.pos + pre
    new_rec.ref = ref
    new_rec.alts = (a,)
    new_rec.info["original_ID"] = original_id

    a_type_raw = new_rec.info.get("allele_type")
    a_len_raw = new_rec.info.get("allele_length")
    a_type = str(a_type_raw).upper()
    a_len = abs(int(a_len_raw))
    if a_type == "SNV":
        new_rec.id = f"{new_rec.chrom}-{new_rec.pos}-{new_rec.ref}-{new_rec.alts[0]}"
    else:
        new_rec.id = f"{new_rec.chrom}-{new_rec.pos}-{a_type}-{a_len}"

    vcf_out.write(new_rec)

vcf_in.close()
vcf_out.close()
EOF

        bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o "~{prefix}.vcf.gz" \
            unsorted.vcf.gz

        tabix -p vcf "~{prefix}.vcf.gz"
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, vcf_idx], "GB")) + 5,
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
