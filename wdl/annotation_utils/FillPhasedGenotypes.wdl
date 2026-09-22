version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FillPhasedGenotypes {
    meta {
        description: [
            "This utility transfers phasing information from a phased VCF onto the genotypes of an unphased VCF over matching sites, optionally sharding each contig by region. It outputs the phased VCF."
        ]
    }

    parameter_meta {
        phased_vcf: "VCF providing the phasing information."
        phased_vcf_idx: "Index for `phased_vcf`."
        unphased_vcf: "VCF whose genotypes are phased."
        unphased_vcf_idx: "Index for `unphased_vcf`."
        contigs: "Contigs to process."
        shard_bin_size: "Region-bin size, in bp, used when sharding each contig."
        hiphase_phased_vcf: "Phased VCF."
        hiphase_phased_vcf_idx: "Index for the phased VCF."
    }

    input {
        File phased_vcf
        File phased_vcf_idx
        File unphased_vcf
        File unphased_vcf_idx
        Array[String] contigs
        String prefix

        Int? shard_bin_size

        String utils_docker

        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_fill
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetPhased {
            input:
                vcf = phased_vcf,
                vcf_idx = phased_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.phased",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        call Helpers.SubsetVcfToContig as SubsetUnphased {
            input:
                vcf = unphased_vcf,
                vcf_idx = unphased_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.unphased",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        if (defined(shard_bin_size)) {
            call Helpers.CreateContigShards {
                input:
                    vcfs = [SubsetPhased.subset_vcf, SubsetUnphased.subset_vcf],
                    vcf_idxs = [SubsetPhased.subset_vcf_idx, SubsetUnphased.subset_vcf_idx],
                    contig = contig,
                    shard_bin_size = select_first([shard_bin_size]),
                    prefix = "~{prefix}.~{contig}.shards",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_create_shards
            }

            scatter (i in range(length(CreateContigShards.shard_regions))) {
                String shard_region = CreateContigShards.shard_regions[i]

                call Helpers.SubsetVcfToRegion as SubsetPhasedShard {
                    input:
                        vcf = SubsetPhased.subset_vcf,
                        vcf_idx = SubsetPhased.subset_vcf_idx,
                        region = shard_region,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.phased",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset
                }

                call Helpers.SubsetVcfToRegion as SubsetUnphasedShard {
                    input:
                        vcf = SubsetUnphased.subset_vcf,
                        vcf_idx = SubsetUnphased.subset_vcf_idx,
                        region = shard_region,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.unphased",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset
                }

                call FillGenotypes {
                    input:
                        phased_vcf = SubsetPhasedShard.subset_vcf,
                        phased_vcf_idx = SubsetPhasedShard.subset_vcf_idx,
                        unphased_vcf = SubsetUnphasedShard.subset_vcf,
                        unphased_vcf_idx = SubsetUnphasedShard.subset_vcf_idx,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.filled",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_fill
                }
            }

            call Helpers.ConcatVcfs as ConcatContig {
                input:
                    vcfs = FillGenotypes.filled_vcf,
                    vcf_idxs = FillGenotypes.filled_vcf_idx,
                    allow_overlaps = false,
                    naive = true,
                    prefix = "~{prefix}.~{contig}.filled",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat
            }
        }

        if (!defined(shard_bin_size)) {
            call FillGenotypes as FillGenotypesNoSharding {
                input:
                    phased_vcf = SubsetPhased.subset_vcf,
                    phased_vcf_idx = SubsetPhased.subset_vcf_idx,
                    unphased_vcf = SubsetUnphased.subset_vcf,
                    unphased_vcf_idx = SubsetUnphased.subset_vcf_idx,
                    prefix = "~{prefix}.~{contig}.filled",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_fill
            }
        }

        File contig_filled_vcf = select_first([ConcatContig.concat_vcf, FillGenotypesNoSharding.filled_vcf])
        File contig_filled_vcf_idx = select_first([ConcatContig.concat_vcf_idx, FillGenotypesNoSharding.filled_vcf_idx])
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = contig_filled_vcf,
            vcf_idxs = contig_filled_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.filled",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File hiphase_phased_vcf = ConcatVcfs.concat_vcf
        File hiphase_phased_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task FillGenotypes {
    input {
        File phased_vcf
        File phased_vcf_idx
        File unphased_vcf
        File unphased_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

phased_in = pysam.VariantFile("~{phased_vcf}")
unphased_in = pysam.VariantFile("~{unphased_vcf}")

NULL_GT = [(None, None), (None, 0), (0, None), (0, ), (None, ), None]

extra_fmt_keys = []
for key in unphased_in.header.formats:
    if key not in phased_in.header.formats:
        phased_in.header.add_record(unphased_in.header.formats[key].record)
        extra_fmt_keys.append(key)

out = pysam.VariantFile("~{prefix}.vcf.gz", "w", header=phased_in.header)

def alleles_key(rec):
    ref = rec.ref.upper() if rec.ref else rec.ref
    alts = tuple(a.upper() for a in rec.alts) if rec.alts else ()
    return (rec.chrom, rec.pos, ref, alts)

for record in phased_in:
    match = None
    rec_key = alleles_key(record)
    for cand in unphased_in.fetch(record.chrom, record.start, record.stop):
        if alleles_key(cand) == rec_key:
            match = cand
            break

    if match:
        for sample in record.samples:
            for fmt_key in record.samples[sample].keys():
                if fmt_key == "GT":
                    if record.samples[sample]['GT'] in NULL_GT:
                        try:
                            record.samples[sample]['GT'] = match.samples[sample]['GT']
                        except Exception:
                            pass
                else:
                    if fmt_key in match.samples[sample]:
                        try:
                            record.samples[sample][fmt_key] = match.samples[sample][fmt_key]
                        except Exception:
                            pass

            for fmt_key in extra_fmt_keys:
                if fmt_key in match.samples[sample]:
                    try:
                        record.samples[sample][fmt_key] = match.samples[sample][fmt_key]
                    except Exception:
                        pass
    out.write(record)

out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File filled_vcf = "~{prefix}.vcf.gz"
        File filled_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(phased_vcf, "GB") + size(unphased_vcf, "GB")) + 25,
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
