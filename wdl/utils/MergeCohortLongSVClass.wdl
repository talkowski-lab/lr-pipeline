# Derived from yuliamostovoy/lr_callset_integration, a fork of fabio-cunial/callset_integration_phase2:
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_WorkflowD_PerSuffix.wdl
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_Workpackage12.wdl
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_Workpackage13.wdl
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_Workpackage14.wdl
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_Workpackage15.wdl

version 1.0

import "Structs.wdl"

workflow MergeCohortLongSVClass {
    meta {
        description: [
            "This sub-workflow merges one class of per-sample calls, either the SVs longer than the MergeSampleSVCallsets length range or the breakends, into a cohort callset. The per-sample calls are merged with bcftools merge, each chromosome is resharded into truvari-collapse shards, matching sites within each shard are collapsed with Truvari (https://github.com/ACEnglish/truvari) collapse, and the collapsed shards are concatenated per chromosome and then genome-wide."
        ]
    }

    parameter_meta {
        suffix: "Class to merge: 'ultralong' or 'bnd'."
        remote_indir: "GCS directory holding the per-sample calls of the class."
        remote_outdir_suffix: "GCS directory the class's merge, shard, collapse and concatenation stages are written under."
        chromosomes: "Chromosomes to process, in output order."
        n_expected_samples: "Number of samples to merge. Derived from the per-sample files of the class when omitted."
        truvari_chunk_min_records: "Minimum number of records in each truvari-collapse shard."
        truvari_collapse_refdist: "Distance, in bp, that shard boundaries keep from any record, so that records Truvari could collapse together fall in one shard."
        consistency_checks: "Whether to verify that sharding kept every record: 1 for yes, 0 for no."
        truvari_matching_parameters: "Truvari collapse matching arguments."
        max_resolve: "Maximum length of a symbolic SV whose sequence Truvari resolves from the reference before collapsing."
        use_bed: "Whether Truvari collapse is restricted to each shard's intervals with a BED."
        chunk_ids_per_file: "Number of truvari-collapse shards processed on each VM."
        concat_all_naive: "Whether the genome-wide concatenation uses bcftools concat --naive: 1 for yes, 0 for no."
        ref_fa: "From references."
        ref_fai: "From references."
        done: "Completion signal of the genome-wide concatenation."
        cohort_dir: "GCS directory holding the class's cohort callset, 'truvari_collapsed.bcf'."
    }

    input {
        String suffix
        String remote_indir
        String remote_outdir_suffix
        Array[String] chromosomes
        Int? n_expected_samples

        Int truvari_chunk_min_records = 2000
        Int truvari_collapse_refdist = 1000
        Int consistency_checks = 1
        String truvari_matching_parameters = "--refdist 500 --pctseq 0.95 --pctsize 0.95 --pctovl 0.0"
        Int max_resolve = 100000
        Boolean use_bed = false
        Int chunk_ids_per_file = 100
        Int concat_all_naive = 1

        File ref_fa
        File ref_fai

        String sv_integration_docker

        RuntimeAttr? runtime_attr_write_sample_list
        RuntimeAttr? runtime_attr_merge_calls
        RuntimeAttr? runtime_attr_shard_chromosome
        RuntimeAttr? runtime_attr_derive_chunk_ids
        RuntimeAttr? runtime_attr_collapse_shards
        RuntimeAttr? runtime_attr_concat_chromosome_shards
        RuntimeAttr? runtime_attr_concat_chromosomes
    }

    String merge_dir = remote_outdir_suffix + "/12_merge"
    String shard_dir = remote_outdir_suffix + "/13_shard"
    String collapse_dir = remote_outdir_suffix + "/14_collapse"
    String concat_dir = remote_outdir_suffix + "/15_concat"

    call WriteSampleList {
        input:
            remote_indir = remote_indir,
            suffix = suffix,
            docker = sv_integration_docker,
            runtime_attr_override = runtime_attr_write_sample_list
    }
    Int n_samples = select_first([n_expected_samples, WriteSampleList.n_samples])

    call MergeCalls {
        input:
            sample_ids = WriteSampleList.sample_ids_file,
            suffix = suffix,
            chromosomes = chromosomes,
            remote_indir = remote_indir,
            n_expected_samples = n_samples,
            remote_outdir = merge_dir,
            docker = sv_integration_docker,
            runtime_attr_override = runtime_attr_merge_calls
    }

    scatter (chromosome in chromosomes) {
        call ShardChromosome {
            input:
                chromosome_id = chromosome,
                suffix = suffix,
                truvari_chunk_min_records = truvari_chunk_min_records,
                truvari_collapse_refdist = truvari_collapse_refdist,
                consistency_checks = consistency_checks,
                remote_indir = merge_dir,
                remote_outdir = shard_dir,
                upstream_signal = [MergeCalls.done],
                docker = sv_integration_docker,
                runtime_attr_override = runtime_attr_shard_chromosome
        }

        call DeriveChunkIds {
            input:
                regions_txt = ShardChromosome.regions_txt,
                chunk_ids_per_file = chunk_ids_per_file,
                docker = sv_integration_docker,
                runtime_attr_override = runtime_attr_derive_chunk_ids
        }

        scatter (chunk_ids_file in DeriveChunkIds.chunk_id_files) {
            call CollapseShards {
                input:
                    remote_indir = shard_dir,
                    chromosome_id = chromosome,
                    chunks_ids = chunk_ids_file,
                    remote_outdir = collapse_dir,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    truvari_matching_parameters = truvari_matching_parameters,
                    max_resolve = max_resolve,
                    use_bed = use_bed,
                    docker = sv_integration_docker,
                    runtime_attr_override = runtime_attr_collapse_shards
            }
        }

        call ConcatChromosomeShards {
            input:
                chromosome = chromosome,
                remote_indir = collapse_dir,
                remote_outdir = concat_dir,
                upstream_signal = CollapseShards.done,
                docker = sv_integration_docker,
                runtime_attr_override = runtime_attr_concat_chromosome_shards
        }
    }

    call ConcatChromosomes {
        input:
            chromosomes = chromosomes,
            out_txt = ConcatChromosomeShards.out_txt,
            remote_outdir = concat_dir,
            naive = concat_all_naive,
            docker = sv_integration_docker,
            runtime_attr_override = runtime_attr_concat_chromosomes
    }

    output {
        String done = ConcatChromosomes.done
        String cohort_dir = concat_dir
    }
}

task WriteSampleList {
    input {
        String remote_indir
        String suffix
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        gcloud storage ls ~{remote_indir}/'*_~{suffix}.bcf' | sed 's#.*/##; s#_~{suffix}\.bcf$##' | sort -u > sample_ids.txt
        if [ ! -s sample_ids.txt ]; then
            echo "ERROR: no <sample>_~{suffix}.bcf files found under ~{remote_indir}."
            exit 1
        fi
        wc -l < sample_ids.txt > n.txt
        cat sample_ids.txt 1>&2
    >>>

    output {
        File sample_ids_file = "sample_ids.txt"
        Int n_samples = read_int("n.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 16,
        boot_disk_gb: 10,
        preemptible_tries: 3,
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

task MergeCalls {
    input {
        File sample_ids
        String suffix
        Array[String] chromosomes

        String remote_indir
        Int n_expected_samples

        Int n_files_per_merge = 100
        String remote_outdir
        String docker

        RuntimeAttr? runtime_attr_override
    }

    String docker_dir = "/callset_integration"

    command <<<
        set -euo pipefail

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        function LocalizeFiles() {
            local N_SAMPLES=$(wc -l < ~{sample_ids})
            if [ ${N_SAMPLES} -ne ~{n_expected_samples} ]; then
                echo "ERROR: sample_ids has ${N_SAMPLES} samples != ~{n_expected_samples}"
                exit 1
            fi

            rm -f all_remote_files.txt uri_list.txt
            while read -u 5 SAMPLE_ID; do
                echo "~{remote_indir}/${SAMPLE_ID}_~{suffix}.bcf" >> uri_list.txt
                echo "~{remote_indir}/${SAMPLE_ID}_~{suffix}.bcf.csi" >> uri_list.txt
            done 5< ~{sample_ids}
            while read -u 6 URI; do
                gcloud storage ls -l "${URI}" | grep -E '^[[:space:]]*[0-9]+' | sed 's/^[ ]*//' >> all_remote_files.txt
            done 6< uri_list.txt

            # Failing immediately if the files are too large WRT the available disk. Otherwise the VM may get stuck forever, and this gets worse with preemption.
            local AVAILABLE_GB=$(df -h | grep "cromwell_root" | tr -s ' ' | cut -d ' ' -f 4)
            AVAILABLE_GB=${AVAILABLE_GB%G}
            AVAILABLE_GB=${AVAILABLE_GB%.*}
            local REMOTE_GB=$(java -cp ~{docker_dir} SumFileSizes all_remote_files.txt)
            local SLACK_GB="5"
            REMOTE_GB=$(( ${REMOTE_GB} + ${SLACK_GB} ))
            if [ ${REMOTE_GB} -gt ${AVAILABLE_GB} ]; then
                echo "ERROR: the remote files are larger than the available disk space. Remote files + slack: ${REMOTE_GB}GB. Available disk: ${AVAILABLE_GB}GB."
                exit 1
            fi

            # Localizing all the single-sample VCFs.
            date 1>&2
            mkdir ./input_files/
            ${TIME_COMMAND} gcloud storage cp -I ./input_files/ < uri_list.txt
            date 1>&2
            local N_DOWNLOADED_SAMPLES=$(ls ./input_files/*.bcf | wc -l)
            if [ ${N_DOWNLOADED_SAMPLES} -lt ${N_SAMPLES} ]; then
                echo "ERROR: The number of downloaded samples (${N_DOWNLOADED_SAMPLES}) is smaller than the number of samples specified (${N_SAMPLES})."
                exit 1
            elif [ ${N_DOWNLOADED_SAMPLES} -gt ${N_SAMPLES} ]; then
                echo "ERROR: The number of downloaded samples (${N_DOWNLOADED_SAMPLES}) is larger than the number of samples specified (${N_SAMPLES})."
                exit 1
            fi
            df -h 1>&2
        }

        cat << 'END' > chunk_by_chr.sh
#!/bin/bash

INPUT_BCF=$1
CHROMOSOME=$2
mkdir -p ./${CHROMOSOME}/
bcftools view --output-type b ${INPUT_BCF} ${CHROMOSOME} --output ./${CHROMOSOME}/${INPUT_BCF}
bcftools index -f ./${CHROMOSOME}/${INPUT_BCF}
END
        chmod +x chunk_by_chr.sh

        cat ~{write_lines(chromosomes)} > chromosomes.txt
        LocalizeFiles

        # Trivial "hierarchical" bcftools merge with just two steps. Step 1: merging a few samples at a time over the whole genome. Reheader each per-sample BCF to its sample_ids name before merging so the cohort uses the same sample names as the main branch (WP3).
        rm -f list.txt
        while read -u 3 SAMPLE_ID; do
            echo ${SAMPLE_ID} > ${SAMPLE_ID}.sample_name.txt
            ${TIME_COMMAND} bcftools reheader --samples ${SAMPLE_ID}.sample_name.txt --output ./input_files/${SAMPLE_ID}_~{suffix}.reheader.bcf ./input_files/${SAMPLE_ID}_~{suffix}.bcf
            mv ./input_files/${SAMPLE_ID}_~{suffix}.reheader.bcf ./input_files/${SAMPLE_ID}_~{suffix}.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ./input_files/${SAMPLE_ID}_~{suffix}.bcf
            echo ./input_files/${SAMPLE_ID}_~{suffix}.bcf >> list.txt
        done 3< ~{sample_ids}
        split -l ~{n_files_per_merge} -d -a 4 list.txt list_
        N_LIST_FILES=$(ls list_* | wc -l)
        for LIST_FILE in $(ls list_* | sort -V); do
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --force-single --merge none --file-list ${LIST_FILE} --output-type b --output ${LIST_FILE}_merged.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_merged.bcf
            xargs --arg-file=${LIST_FILE} --max-lines=1 --max-procs=${N_THREADS} rm -f
            rm -f ${LIST_FILE}
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ${LIST_FILE}_merged.bcf --output ${LIST_FILE}_normed.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ${LIST_FILE}_normed.bcf
            rm -f ${LIST_FILE}_merged.bcf*
            ${TIME_COMMAND} xargs --arg-file=chromosomes.txt --max-lines=1 --max-procs=${N_THREADS} ./chunk_by_chr.sh ./${LIST_FILE}_normed.bcf
            rm -f ${LIST_FILE}_normed.bcf*
        done
        rm -rf ./input_files/

        # Step 2: merging all samples over each chromosome.
        rm -f files_list.txt
        while read -u 4 CHROMOSOME; do
            ls ./${CHROMOSOME}/*.bcf | sort -V > list.txt
            # --force-single lets the merge proceed when there is only one file (small cohorts produce a single batch per chromosome).
            ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --force-single --merge none --file-list list.txt --output-type b --output ./${CHROMOSOME}/merged.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ./${CHROMOSOME}/merged.bcf
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type b ./${CHROMOSOME}/merged.bcf --output ./${CHROMOSOME}/normed.bcf
            ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f ./${CHROMOSOME}/normed.bcf
            mv ./${CHROMOSOME}/normed.bcf ${CHROMOSOME}.bcf
            mv ./${CHROMOSOME}/normed.bcf.csi ${CHROMOSOME}.bcf.csi
            echo "${CHROMOSOME}.bcf" >> files_list.txt
            echo "${CHROMOSOME}.bcf.csi" >> files_list.txt
            rm -rf ./${CHROMOSOME}/
        done 4< chromosomes.txt
        df -h 1>&2
        ls -laht 1>&2

        # Uploading
        date 1>&2
        cat files_list.txt | gcloud storage mv -I ~{remote_outdir}/
        date 1>&2

        # Completion signal for orchestrator ordering. Ignored standalone.
        echo "done" > wp12.signal
    >>>

    output {
        String done = read_string("wp12.signal")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 50,
        boot_disk_gb: 10,
        preemptible_tries: 4,
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

task ShardChromosome {
    input {
        String chromosome_id
        String suffix
        Int truvari_chunk_min_records
        Int truvari_collapse_refdist
        Int consistency_checks

        String remote_indir
        String remote_outdir

        Array[String]? upstream_signal
        String docker

        RuntimeAttr? runtime_attr_override
    }

    String docker_dir = "/callset_integration"

    command <<<
        set -euo pipefail

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]))} - 1 ))

        cat << 'END' > chunk_by_region.sh
#!/bin/bash

REGION=$1
CHUNK_ID=$2
bcftools view --regions ${REGION} --regions-overlap pos --output-type b ~{chromosome_id}.bcf --output chunk_${CHUNK_ID}.bcf
bcftools index -f chunk_${CHUNK_ID}.bcf
END
        chmod +x chunk_by_region.sh

        gcloud storage cp ~{remote_indir}/~{chromosome_id}.'bcf*' .

        # Splitting
        N_RECORDS=$(bcftools index --nrecords ~{chromosome_id}.bcf)
        if [ ~{consistency_checks} -eq 1 ]; then
            ${TIME_COMMAND} bcftools query --format '%ID\n' ~{chromosome_id}.bcf > ids_truth.txt
        fi
        ${TIME_COMMAND} bcftools query --format '%POS\t%REF\t%ALT\n' ~{chromosome_id}.bcf > pos_ref_alt.tsv
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB}G TruvariDivide2Ultralong pos_ref_alt.tsv ~{truvari_collapse_refdist} ~{truvari_chunk_min_records} ~{chromosome_id} ${N_RECORDS} ~{suffix} > regions.txt
        rm -f pos_ref_alt.tsv
        ${TIME_COMMAND} xargs --arg-file=regions.txt --max-lines=1 --max-procs=${N_THREADS} ./chunk_by_region.sh
        ls -laht 1>&2
        df -h  1>&2
        rm -f ~{chromosome_id}.bcf*

        # Simple consistency checks
        if [ ~{consistency_checks} -eq 1 ]; then
            N_RECORDS_CHUNKED="0"
            for FILE in $(ls chunk_*.bcf.csi | sort -V); do
                N=$( bcftools index --nrecords ${FILE} )
                N_RECORDS_CHUNKED=$(( ${N_RECORDS_CHUNKED} + ${N} ))
            done
            if [ ${N_RECORDS_CHUNKED} -ne ${N_RECORDS} ]; then
                echo "ERROR: The truvari collapse chunks contain ${N_RECORDS_CHUNKED} total records, but the chromosome VCF contains ${N_RECORDS} records."
                exit 1
            fi
            rm -f ids_test.txt
            for FILE in $(ls chunk_*.bcf | sort -V); do
                bcftools query --format '%ID\n' ${FILE} >> ids_test.txt
            done
            diff --brief ids_test.txt ids_truth.txt
        fi

        # Uploading
        ${TIME_COMMAND} gcloud storage mv 'chunk_*.bcf*' ~{remote_outdir}/~{chromosome_id}/
    >>>

    output {
        # col 1 = region, col 2 = chunk id. Delocalized so an orchestrator can derive the truvari-collapse chunk-id lists in-graph (this task does not upload regions.txt to GCS).
        File regions_txt = "regions.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 50,
        boot_disk_gb: 10,
        preemptible_tries: 4,
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

task DeriveChunkIds {
    input {
        File regions_txt
        Int chunk_ids_per_file
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        awk 'NF>=2 {print $2}' ~{regions_txt} > ids.txt
        if [ ! -s ids.txt ]; then
            echo "ERROR: no chunk ids found in regions.txt."
            exit 1
        fi
        split -l ~{chunk_ids_per_file} -d -a 4 ids.txt chunk_ids_
        ls chunk_ids_* 1>&2
    >>>

    output {
        Array[File] chunk_id_files = glob("chunk_ids_*")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 16,
        boot_disk_gb: 10,
        preemptible_tries: 3,
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

task CollapseShards {
    input {
        String remote_indir
        String chromosome_id
        File chunks_ids
        String remote_outdir
        File ref_fa
        File ref_fai

        String truvari_matching_parameters
        Int max_resolve
        Boolean use_bed
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]))} - 1 ))

        # Removes SVLEN from symbolic ALTs, in order not to interfere with `truvari collapse`.
        function ResetAlts() {
            local CHUNK_ID=$1

            date 1>&2
            ( bcftools view --header-only chunk_${CHUNK_ID}.bcf ; bcftools view --no-header chunk_${CHUNK_ID}.bcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                if (substr($0,1,1)!="#" && substr($5,1,1)=="<") $5 = substr($5,1,4) ">"; \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' ) | bcftools view --output-type b --output out.bcf
            date 1>&2
            rm -f chunk_${CHUNK_ID}.bcf* ; mv out.bcf chunk_${CHUNK_ID}.bcf ; bcftools index --threads ${N_THREADS} chunk_${CHUNK_ID}.bcf
        }

        # Sets QUAL to the number of samples where a record was discovered, to simulate `--keep common` (which is slow on 10k samples) with `--keep maxqual` in `truvari collapse`. See e.g.: https://github.com/ACEnglish/truvari/issues/220#issuecomment- 2830920205 Remark: we use the number of samples rather than AC, since we don't trust genotypes at this stage, and since a record being discovered independently in more samples is more informative than it being genotyped more times in fewer samples.
        function CopyNSamplesToQual() {
            local CHUNK_ID=$1

            mv chunk_${CHUNK_ID}.bcf chunk_${CHUNK_ID}_in.bcf
            mv chunk_${CHUNK_ID}.bcf.csi chunk_${CHUNK_ID}_in.bcf.csi

            # Remark: we cannot join annotations just by ID at this stage, since the IDs in the output of bcftools merge are not necessarily all distinct.
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%COUNT(GT="alt")\n' chunk_${CHUNK_ID}_in.bcf | bgzip -c > chunk_${CHUNK_ID}_annotations.tsv.gz
            tabix -@ ${N_THREADS} -s1 -b2 -e2 chunk_${CHUNK_ID}_annotations.tsv.gz
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --annotations chunk_${CHUNK_ID}_annotations.tsv.gz --columns CHROM,POS,~ID,REF,ALT,QUAL --output-type z chunk_${CHUNK_ID}_in.bcf --output chunk_${CHUNK_ID}_out.vcf.gz
            rm -f chunk_${CHUNK_ID}_in.bcf* ; mv chunk_${CHUNK_ID}_out.vcf.gz chunk_${CHUNK_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t chunk_${CHUNK_ID}_in.vcf.gz
            rm -f chunk_${CHUNK_ID}_annotations.tsv.gz

            mv chunk_${CHUNK_ID}_in.vcf.gz chunk_${CHUNK_ID}_annotated.vcf.gz
            mv chunk_${CHUNK_ID}_in.vcf.gz.tbi chunk_${CHUNK_ID}_annotated.vcf.gz.tbi
        }

        # Remark: in theory we should set `--gt all` to avoid collapsing records that are present in the same sample, since we assume that intra- sample merging has already done that upstream. In practice `--gt all` is too slow on 10k samples. Remark: to further improve speed we could think of dropping genotypes before running truvari collapse. See e.g.: https://github.com/ACEnglish/truvari/issues/220#issuecomment- 2830920205 However, this would also discard e.g. SUPP fields that were copied to FORMAT upstream, so it is not correct for our setup. It would also make it impossible e.g. to compare precision/recall after collapse to precision/recall after cohort re-genotyping.
        function Collapse() {
            local CHUNK_ID=$1

            mv chunk_${CHUNK_ID}_annotated.vcf.gz chunk_${CHUNK_ID}_in.vcf.gz
            mv chunk_${CHUNK_ID}_annotated.vcf.gz.tbi chunk_${CHUNK_ID}_in.vcf.gz.tbi

            # Remark: we do not store `removed.vcf` since it's not needed and it can be much bigger than the collapsed output.
            ${TIME_COMMAND} truvari collapse --sizemin 0 --sizemax ${INFINITY} --keep maxqual --gt off --reference ~{ref_fa} --max-resolve ~{max_resolve} --dup-to-ins ~{truvari_matching_parameters} ${BED_FLAGS} --input chunk_${CHUNK_ID}_in.vcf.gz --output chunk_${CHUNK_ID}_out.vcf --removed-output /dev/null
            df -h 1>&2
            ls -laht 1>&2
            rm -f chunk_${CHUNK_ID}_in.vcf.gz* ; mv chunk_${CHUNK_ID}_out.vcf chunk_${CHUNK_ID}_in.vcf

            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type b chunk_${CHUNK_ID}_in.vcf --output chunk_${CHUNK_ID}_out.bcf
            df -h 1>&2
            ls -laht 1>&2
            rm -f chunk_${CHUNK_ID}_in.vcf ; mv chunk_${CHUNK_ID}_out.bcf chunk_${CHUNK_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f chunk_${CHUNK_ID}_in.bcf

            # Dropping the IDs written by truvari collapse, since they can be very long on a large cohort and needlessly inflate output size. IDs are reassigned by record order (chunk id, running index) via a streaming rewrite of the ID column. This is done in a single pass instead of a position-based `bcftools annotate`, because the latter matches records on CHROM,POS,REF,ALT and therefore cannot tell apart same-position symbolic records (identical CHROM/POS/REF/ALT but different END/SVLEN), assigning them the same ID.
            ${TIME_COMMAND} bcftools view chunk_${CHUNK_ID}_in.bcf | awk -v id=${CHUNK_ID} 'BEGIN { FS="\t"; OFS="\t"; i=0; } /^#/ { print; next } { $3=sprintf("%s_%d",id,i++); print $0 }' | bcftools view --output-type b --output chunk_${CHUNK_ID}_out.bcf
            rm -f chunk_${CHUNK_ID}_in.bcf ; mv chunk_${CHUNK_ID}_out.bcf chunk_${CHUNK_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f chunk_${CHUNK_ID}_in.bcf

            mv chunk_${CHUNK_ID}_in.bcf chunk_${CHUNK_ID}_truvari.bcf
            mv chunk_${CHUNK_ID}_in.bcf.csi chunk_${CHUNK_ID}_truvari.bcf.csi
        }

        INFINITY="1000000000"
        truvari --help 1>&2
        ls ~{ref_fai} 1>&2
        df -h 1>&2

        if ~{use_bed} ; then
            gcloud storage cp ~{remote_indir}/~{chromosome_id}/included.bed .
            BED_FLAGS="--bed included.bed"
        else
            BED_FLAGS=" "
        fi
        while read -u 3 CHUNK_ID; do
            # Skipping the chunk if it has already been processed
            TEST=$( gsutil ls ~{remote_outdir}/~{chromosome_id}/chunk_${CHUNK_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi

            # Collapsing
            gcloud storage cp ~{remote_indir}/~{chromosome_id}/chunk_${CHUNK_ID}.'bcf*' .
            ResetAlts ${CHUNK_ID}
            CopyNSamplesToQual ${CHUNK_ID}
            Collapse ${CHUNK_ID}

            # Uploading
            gcloud storage mv chunk_${CHUNK_ID}_truvari.bcf'*' ~{remote_outdir}/~{chromosome_id}/
            touch chunk_${CHUNK_ID}.done
            gcloud storage mv chunk_${CHUNK_ID}.done ~{remote_outdir}/~{chromosome_id}/
            rm -rf chunk_${CHUNK_ID}*
            ls -laht 1>&2
        done 3< ~{chunks_ids}

        # Completion signal for orchestrator ordering. Ignored standalone.
        echo "done" > wp14.signal
    >>>

    output {
        String done = read_string("wp14.signal")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 4,
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

task ConcatChromosomeShards {
    input {
        String chromosome
        String remote_indir
        String remote_outdir

        Array[String]? upstream_signal
        String docker

        RuntimeAttr? runtime_attr_override
    }

    String docker_dir = "/callset_integration"

    command <<<
        set -euo pipefail

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        TEST=$( gcloud storage ls ~{remote_outdir}/~{chromosome}/~{chromosome}.done || echo "0" )
        if [ ${TEST} != "0" ]; then
            # Skipping the chromosome if it has already been processed
            :
        else
            # Localizing all chunks
            gcloud storage ls ~{remote_indir}/~{chromosome}/chunk_'*.bcf' > test.txt
            if grep -q '.bcf' test.txt ; then
                :
            else
                echo "ERROR: ~{chromosome} has no truvari collapse chunks."
                exit
            fi
            ${TIME_COMMAND} gcloud storage cp ~{remote_indir}/~{chromosome}/chunk_'*.bcf*' .
            ls chunk_*.bcf | sort -V > chunk_list.txt
            cat chunk_list.txt
            df -h 1>&2

            # Concatenating all chunks
            N_CHUNKS=$(wc -l < chunk_list.txt)
            if [ ${N_CHUNKS} -gt 1 ]; then
                ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list chunk_list.txt --output-type b --output out.bcf
                df -h 1>&2
                rm -rf chunk_* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf
            else
                CHUNK_FILE=$(head -n 1 chunk_list.txt)
                mv ${CHUNK_FILE} in.bcf
                mv ${CHUNK_FILE}.csi in.bcf.csi
            fi

            # Enforcing a distinct ID in every record, and annotating every record with the number of samples it occurs in. Note that the latter is not equal to the QUAL field in input to truvari collapse upstream, so we have to recompute this number.
            CHR=~{chromosome}
            CHR=${CHR#chr}
            # Enforcing a distinct ID in every record, and annotating every record with the number of samples it occurs in. Note that the latter is not equal to the QUAL field in input to truvari collapse upstream, so we have to recompute this number. Both operations are keyed on record ORDER, not on CHROM/POS/REF/ALT. A position-based `bcftools annotate` matches on CHROM,POS,REF,ALT and cannot distinguish same-position symbolic records (identical CHROM/POS/REF/ALT but different END/SVLEN); it would assign such records the same ID and copy the same N_DISCOVERY_SAMPLES to all of them. Matching on ID (`-c ~ID`) is not usable either (unimplemented in bcftools). `bcftools query` and `bcftools view` both emit records in file order, so counts.txt aligns 1:1 with the streamed records.
            ${TIME_COMMAND} bcftools query --format '%COUNT(GT="alt")\n' in.bcf > counts.txt
            ${TIME_COMMAND} bcftools view in.bcf | awk -v id=${CHR} 'BEGIN { FS="\t"; OFS="\t"; i=0; j=0; while ((getline c < "counts.txt") > 0) { cnt[j++]=c; } } /^#CHROM/ { print "##INFO=<ID=N_DISCOVERY_SAMPLES,Number=1,Type=Integer,Description=\"Number of samples where the record was discovered\">"; print $0; next } /^#/ { print $0; next } { $3=sprintf("%s_%d",id,i); if ($8==".") { $8="N_DISCOVERY_SAMPLES=" cnt[i] } else { $8=$8 ";N_DISCOVERY_SAMPLES=" cnt[i] } i++; print $0 }' | bcftools view --output-type b --output out.bcf
            df -h 1>&2
            rm -f in.bcf* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf
            gcloud storage cp in.bcf ~{remote_outdir}/~{chromosome}/truvari_collapsed.bcf
            gcloud storage cp in.bcf.csi ~{remote_outdir}/~{chromosome}/truvari_collapsed.bcf.csi
        fi
        echo "~{chromosome}" > out.txt
    >>>

    output {
        File out_txt = "out.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 4,
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

task ConcatChromosomes {
    input {
        Array[String] chromosomes
        Array[File] out_txt
        String remote_outdir
        Int naive
        String docker

        RuntimeAttr? runtime_attr_override
    }

    String docker_dir = "/callset_integration"

    command <<<
        set -euo pipefail

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        # Localizing
        CHROMOSOMES=~{sep=',' chromosomes}
        echo ${CHROMOSOMES} | tr ',' '\n' > chr_list.txt
        rm -f file_list.txt
        while read -u 3 CHROMOSOME; do
            TEST=$( gcloud storage ls ~{remote_outdir}/${CHROMOSOME}/truvari_collapsed.bcf || echo 1 )
            if [ ${TEST} -eq 1 ]; then
                echo "ERROR: ${CHROMOSOME} has not been truvari collapsed."
                exit
            fi
            gcloud storage cp ~{remote_outdir}/${CHROMOSOME}/'*.bcf*' .
            mv truvari_collapsed.bcf ${CHROMOSOME}_truvari_collapsed.bcf
            mv truvari_collapsed.bcf.csi ${CHROMOSOME}_truvari_collapsed.bcf.csi
            echo ${CHROMOSOME}_truvari_collapsed.bcf >> file_list.txt
        done 3< chr_list.txt

        # Concatenating
        if [ ~{naive} -eq 1 ]; then
            CONCAT_FLAGS="--naive"
        else
            CONCAT_FLAGS=" "
        fi
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} ${CONCAT_FLAGS} --file-list file_list.txt --output-type b --output truvari_collapsed.bcf
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f truvari_collapsed.bcf

        # Uploading
        gcloud storage mv truvari_collapsed.'bcf*' ~{remote_outdir}/

        # Completion signal for orchestrator ordering. Ignored standalone.
        echo "done" > allchr.signal
    >>>

    output {
        String done = read_string("allchr.signal")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 4,
        disk_gb: 200,
        boot_disk_gb: 10,
        preemptible_tries: 4,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
