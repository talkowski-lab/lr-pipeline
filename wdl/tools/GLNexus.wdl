# Derived from broadinstitute/long-read-pipelines:
# https://github.com/broadinstitute/long-read-pipelines/blob/main/wdl/pipelines/TechAgnostic/VariantCalling/LRJointCallGVCFs.wdl

version 1.0

import "../utils/Helpers.wdl"

workflow GLNexus {
    meta {
        description: "Joint-call gVCFs with GLNexus and convert resulting callset to a Hail MatrixTable."
    }

    parameter_meta {
        gvcfs: "GVCF files to joint-call."
        gvcf_idxs: "GVCF index files, ordered to correspond to gvcfs."
        ref_map_file: "Table mapping reference sequence names to auxiliary files; must include a dict entry."
        prefix: "Prefix for joint-called VCF, its index, and MatrixTable archive."
        background_sample_gvcfs: "Nested arrays of background GVCFs for joint calling."
        background_sample_gvcf_idxs: "Nested arrays of indexes corresponding to background_sample_gvcfs."
        force_add_missing_dp: "Add missing DP fields to gVCFs before joint calling."
        remove_duplicate_zero_depth_reference_blocks: "Remove exact duplicate input gVCF records with a non-alt GT and MIN_DP=0 from each shard."
        bed: "Intervals to which joint calling should be restricted."
        config: "GLNexus configuration preset or .yml filename."
        config_file: "Custom GLNexus configuration file; overrides config."
        more_PL: "Include PL from reference bands and other cases omitted by default."
        squeeze: "Reduce pVCF size by suppressing detail derived from reference bands."
        trim_uncalled_alleles: "Remove alleles with no output GT calls in postprocessing."
        num_cpus: "Number of CPUs to use."
        max_cpus: "Maximum number of CPUs to allow."
        reference: "Reference assembly label for MatrixTable conversion."
        ref_fa: "Reference FASTA file for MatrixTable conversion."
        ref_fai: "Reference FASTA index file for MatrixTable conversion."
        glnexus_docker: "Docker image for GLNexus tasks."
        hail_docker: "Docker image for Hail MatrixTable conversion."
        runtime_attr_get_ranges: "Override runtime attributes for reference range discovery."
        runtime_attr_shard_vcf_by_ranges: "Override runtime attributes for gVCF sharding."
        runtime_attr_call: "Override runtime attributes for GLNexus calling."
        runtime_attr_concat_variants: "Override runtime attributes for VCF concatenation."
        runtime_attr_convert_to_hail_mt: "Override runtime attributes for Hail MatrixTable conversion."
    }

    input {
        Array[File] gvcfs
        Array[File] gvcf_idxs
        File ref_map_file
        String prefix

        Array[Array[File]]? background_sample_gvcfs
        Array[Array[File]]? background_sample_gvcf_idxs
        Boolean force_add_missing_dp = false
        Boolean remove_duplicate_zero_depth_reference_blocks = false
        File? bed
        String config = "DeepVariantWGS"
        File? config_file
        Boolean more_PL = false
        Boolean squeeze = false
        Boolean trim_uncalled_alleles = false
        Int? num_cpus
        Int max_cpus = 64
        String reference = "GRCh38"
        String? ref_fa
        String? ref_fai

        String glnexus_docker
        String hail_docker

        RuntimeAttr? runtime_attr_get_ranges
        RuntimeAttr? runtime_attr_shard_vcf_by_ranges
        RuntimeAttr? runtime_attr_call
        RuntimeAttr? runtime_attr_concat_variants
        RuntimeAttr? runtime_attr_convert_to_hail_mt
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    Int cpus_exp = if defined(num_cpus) then select_first([num_cpus]) else 2 * length(gvcfs)
    Int cpus_act = if cpus_exp < max_cpus then cpus_exp else max_cpus

    call GetRanges {
        input:
            dict = ref_map["dict"],
            bed = bed,
            docker = glnexus_docker,
            runtime_attr_override = runtime_attr_get_ranges
    }

    if (defined(background_sample_gvcfs)) {
        Array[File] flattened_background_sample_gvcfs = flatten(select_first([background_sample_gvcfs]))
        Array[File] flattened_background_sample_gvcf_idxs = flatten(select_first([background_sample_gvcf_idxs]))
        Array[File] gvcfs_and_background_samples = flatten([gvcfs, flattened_background_sample_gvcfs])
        Array[File] gvcf_idxs_and_background_samples = flatten([gvcf_idxs, flattened_background_sample_gvcf_idxs])
    }
    Array[File] final_gvcfs = select_first([gvcfs_and_background_samples, gvcfs])
    Array[File] final_gvcf_idxs = select_first([gvcf_idxs_and_background_samples, gvcf_idxs])

    scatter (p in zip(final_gvcfs, final_gvcf_idxs)) {
        call ShardVCFByRanges {
            input:
                gvcf = p.left,
                tbi = p.right,
                ranges = GetRanges.ranges,
                remove_duplicate_zero_depth_reference_blocks = remove_duplicate_zero_depth_reference_blocks,
                docker = glnexus_docker,
                runtime_attr_override = runtime_attr_shard_vcf_by_ranges
        }
    }

    scatter (i in range(length(ShardVCFByRanges.sharded_gvcfs[0]))) {
        Array[File] per_contig_gvcfs = transpose(ShardVCFByRanges.sharded_gvcfs)[i]

        call GLNexusJointCall {
            input:
                gvcfs = per_contig_gvcfs,
                config = config,
                config_file = config_file,
                more_PL = more_PL,
                squeeze = squeeze,
                trim_uncalled_alleles = trim_uncalled_alleles,
                force_add_missing_dp = force_add_missing_dp,
                num_cpus = cpus_act,
                prefix = prefix,
                docker = glnexus_docker,
                runtime_attr_override = runtime_attr_call
        }
    }

    call ConcatVariants {
        input:
            variant_files = GLNexusJointCall.joint_bcf,
            is_gvcf = false,
            num_cpus = 4,
            prefix = prefix,
            docker = glnexus_docker,
            runtime_attr_override = runtime_attr_concat_variants
    }

    call Helpers.ConvertToHailMT {
        input:
            gvcf = ConcatVariants.combined_vcf,
            tbi = ConcatVariants.combined_vcf_tbi,
            prefix = prefix,
            docker = hail_docker,
            reference = reference,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            runtime_attr_override = runtime_attr_convert_to_hail_mt
    }

    output {
        File joint_vcf = ConcatVariants.combined_vcf
        File joint_vcf_idx = ConcatVariants.combined_vcf_tbi
        File joint_mt = ConvertToHailMT.mt_tar
    }
}

task GetRanges {
    meta {
        description: "Select loci over which to parallelize downstream operations."
    }

    input {
        File dict
        File? bed
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(size(dict, "GB"))

    command <<<
        set -euxo pipefail

        if [[ "~{defined(bed)}" == "true" ]]; then
            cat ~{bed} | awk '{ print $1 ":" $2 "-" $3 }' > ranges.txt
        else
            grep '^@SQ' ~{dict} | \
                awk '{ print $2, $3 }' | \
                sed 's/[SL]N://g' | \
                awk '{ print $1 ":0-" $2 }' \
                > ranges.txt
        fi
    >>>

    output {
        Array[String] ranges = read_lines("ranges.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: disk_size,
        boot_disk_gb: 25,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task ShardVCFByRanges {
    meta {
        description: "Split VCF into smaller ranges for parallelization."
    }

    input {
        File gvcf
        File tbi
        Array[String] ranges
        Boolean remove_duplicate_zero_depth_reference_blocks = false
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2 * ceil(size(gvcf, "GB"))

    command <<<
        set -euxo pipefail

        mkdir per_contig

        INDEX=0
        for RANGE in ~{sep=' ' ranges}
        do
            PINDEX=$(printf "%06d" $INDEX)
            FRANGE=$(echo $RANGE | sed 's/[:-]/___/g')
            OUTFILE="per_contig/$PINDEX.~{basename(gvcf, ".g.vcf.gz")}.locus_$FRANGE.g.vcf.gz"

            if [[ "~{remove_duplicate_zero_depth_reference_blocks}" == "true" ]]; then
                bcftools view ~{gvcf} $RANGE | awk -F$'\t' '
                    BEGIN { removed_count = 0 }
                    function clear_group(    i) {
                        for (i = 1; i <= row_count; i++) delete rows[i]
                        for (i in duplicate_count) delete duplicate_count[i]
                        for (i in removable) delete removable[i]
                        row_count = 0
                    }
                    function flush_group(    i, line) {
                        for (i = 1; i <= row_count; i++) {
                            line = rows[i]
                            if (removable[line] && duplicate_count[line] > 1) {
                                removed_count++
                            } else {
                                print line
                            }
                        }
                    }
                    function is_removable_record(    format_fields, sample_fields, field_count, sample_count, i, gt_index, min_dp_index, gt) {
                        field_count = split($9, format_fields, ":")
                        gt_index = 0
                        min_dp_index = 0
                        for (i = 1; i <= field_count; i++) {
                            if (format_fields[i] == "GT") gt_index = i
                            if (format_fields[i] == "MIN_DP") min_dp_index = i
                        }
                        if (gt_index == 0 || min_dp_index == 0 || NF < 10) return 0
                        sample_count = split($10, sample_fields, ":")
                        if (sample_count < gt_index || sample_count < min_dp_index) return 0
                        gt = sample_fields[gt_index]
                        return sample_fields[min_dp_index] == "0" && gt ~ /^(0|\.)([\/|](0|\.))*$/
                    }
                    /^#/ { print; next }
                    {
                        coordinate = $1 SUBSEP $2
                        if (row_count > 0 && coordinate != current_coordinate) {
                            flush_group()
                            clear_group()
                        }
                        current_coordinate = coordinate
                        rows[++row_count] = $0
                        duplicate_count[$0]++
                        if (is_removable_record()) removable[$0] = 1
                    }
                    END {
                        if (row_count > 0) flush_group()
                        print "Removed " removed_count " duplicate zero-depth non-alt gVCF records" > "/dev/stderr"
                    }
                ' | bgzip > $OUTFILE
            else
                bcftools view ~{gvcf} $RANGE | bgzip > $OUTFILE
            fi

            INDEX=$(($INDEX+1))
        done
    >>>

    output {
        Array[File] sharded_gvcfs = glob("per_contig/*")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: disk_size,
        boot_disk_gb: 25,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task GLNexusJointCall {
    meta {
        description: "Joint-call gVCFs with GLNexus."
        note: "This task includes code to force add missing DP field to gVCFs. This is included here to avoid localization overhead for moving the files around 2x."
    }

    parameter_meta {
        gvcfs: "gVCF files to perform joint calling upon"
        config: "GLNexus configuration preset. One of: gatk, gatk_unfiltered, xAtlas, xAtlas_unfiltered, weCall, weCall_unfiltered, DeepVariant, DeepVariantWGS, DeepVariantWES, DeepVariantWES_MED_DP, DeepVariant_unfiltered, Strelka2, GxS."
        config_file: "Custom configuration file override for GLNexus. If provided, this will override the config parameter."
        more_PL: "Include PL from reference bands and other cases omitted by default"
        squeeze: "Reduce pVCF size by suppressing detail in cells derived from reference bands"
        trim_uncalled_alleles: "Remove alleles with no output GT calls in postprocessing"
        force_add_missing_dp: "Adds DP field from INFO to sample-level data in gVCFs. This is required to enable GLNexus calling on GATK-called GVCFs. This is included as part of this task so that the gVCF files do not need to be localized twice (otherwise joint calling won't scale well to thousands of samples)."
        num_cpus: "Number of CPUs to use"
        prefix: "Output prefix for joined-called BCF and GVCF files"
        runtime_attr_override: "Runtime attributes override struct"
    }

    input {
        Array[File] gvcfs
        String config = "DeepVariantWGS"
        File? config_file
        Boolean more_PL = false
        Boolean squeeze = false
        Boolean trim_uncalled_alleles = false
        Boolean force_add_missing_dp = false
        Int num_cpus = 96
        String prefix = "out"
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 5 * ceil(size(gvcfs, "GB")) * if (force_add_missing_dp) then 2 else 1
    Int mem = 4 * num_cpus

    command <<<
        set -x

        # For guidance on performance settings, see https://github.com/dnanexus-rnd/GLnexus/wiki/Performance
        ulimit -Sn 65536

        echo ~{gvcfs[0]} | sed 's/.*locus_//' | sed 's/.g.vcf.bgz//' | sed 's/___/\t/g' > range.bed

        gvcf_file_list=~{write_lines(gvcfs)}

        if [[ "~{force_add_missing_dp}" == "true" ]]; then
            echo "Force adding missing DP field to gVCFs"
            fixed_gvcf_file_list="fixed_gvcf_file_list.txt"
            while read gvcf ; do
                echo "Processing ${gvcf}"
                bn=$(basename ${gvcf} | sed -e 's/\.gz$//' -e 's/\.bgz$//' -e 's/\.vcf$//' -e 's@\.g$@@')
                nn=${bn}.missing_dp_added.g.vcf.gz
                SM=$(bcftools query -l ${gvcf})
                bcftools view ${gvcf} | grep -v ':DP:' | grep -v '^#' | awk -F$'\t' 'BEGIN{OFS="\t"}{print $1,$2,"0"}' | bgzip -c > annot.txt.gz
                tabix -s1 -b2 -e2 annot.txt.gz
                bcftools annotate -Oz2 -o ${nn} -s "${SM}" -a annot.txt.gz -c CHROM,POS,FORMAT/DP ${gvcf}
                bcftools index -t ${nn}
                echo ${nn} >> ${fixed_gvcf_file_list}
            done < ${gvcf_file_list}

            gvcf_file_list=${fixed_gvcf_file_list}
        fi

        glnexus_cli \
            --config ~{if (defined(config_file)) then "~{config_file}" else "~{config}"} \
            --bed range.bed \
            ~{if more_PL then "--more-PL" else ""} \
            ~{if squeeze then "--squeeze" else ""} \
            ~{if trim_uncalled_alleles then "--trim-uncalled-alleles" else ""} \
            --list ${gvcf_file_list} \
            > ~{prefix}.bcf
    >>>

    output {
        File joint_bcf = "~{prefix}.bcf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: num_cpus,
        mem_gb: mem,
        disk_gb: disk_size,
        boot_disk_gb: 25,
        preemptible_tries: 0,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

task ConcatVariants {
    meta {
        description: "Concatenate VCFs/BCFs into a single .vcf.bgz file and index it."
    }

    input {
        Array[File] variant_files
        Boolean is_gvcf = false
        Int num_cpus = 4
        String prefix = "out"
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2 * ceil(size(variant_files, "GB"))
    String file_suffix = if is_gvcf then "g.vcf.bgz" else "vcf.bgz"

    command <<<
        set -euxo pipefail

        bcftools concat -n ~{sep=' ' variant_files} | bcftools view | bgzip -@ ~{num_cpus} -c > ~{prefix}.~{file_suffix}
        tabix -p vcf ~{prefix}.~{file_suffix}
    >>>

    output {
        File combined_vcf = "~{prefix}.~{file_suffix}"
        File combined_vcf_tbi = "~{prefix}.~{file_suffix}.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: num_cpus,
        mem_gb: 8,
        disk_gb: disk_size,
        boot_disk_gb: 25,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}
