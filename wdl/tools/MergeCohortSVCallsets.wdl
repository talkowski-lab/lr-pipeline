# Derived from yuliamostovoy/lr_callset_integration, a fork of fabio-cunial/callset_integration_phase2:
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_WorkflowB_Merge_Collapse.wdl
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_Workpackage4_Main_shard.wdl
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_Workpackage5_Main_truvari_collapse.wdl
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_Workpackage6_Main_concat_shards.wdl

version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergeCohortSVCallsets {
    meta {
        description: [
            "This tool merges the scored per-sample SV calls written by MergeSampleSVCallsets into one cohort SV callset. Each chunk of the interval CSV is merged across samples with bcftools merge, each chromosome is resharded into truvari-collapse shards, matching sites within each shard are collapsed with Truvari (https://github.com/ACEnglish/truvari) collapse, and the collapsed shards are concatenated per chromosome and then genome-wide.",
            "Nothing is returned as a workflow output. Each stage writes under 'remote_outdir', and the cohort callset is '06_concat/truvari_collapsed.bcf' there, which RegenotypeFamilySVCallsets reads."
        ]
    }

    parameter_meta {
        remote_indir: "The 'remote_outdir' of MergeSampleSVCallsets, whose '02_scoring' subdirectory is read."
        remote_outdir: "GCS directory the merge, shard, collapse and concatenation stages are written under, without a trailing slash."
        sample_ids_file: "Samples to merge, one per line, in the column order of the merged VCF. Derived from the per-sample marker files in 'remote_indir' when omitted."
        chromosomes: "Chromosomes to process, in output order."
        merge_mode: "How bcftools merge matches records: 1 by CHROM, POS, REF and ALT, or 2 by ID."
        truvari_chunk_min_records: "Minimum number of records in each truvari-collapse shard."
        truvari_collapse_refdist: "Distance, in bp, that shard boundaries keep from any record, so that records Truvari could collapse together fall in one shard."
        consistency_checks: "Whether to verify that sharding kept every record: 1 for yes, 0 for no."
        truvari_matching_parameters: "Truvari collapse matching arguments."
        use_bed: "Whether Truvari collapse is restricted to each shard's intervals with a BED."
        chunk_ids_per_file: "Number of truvari-collapse shards processed on each VM."
        concat_all_naive: "Whether the genome-wide concatenation uses bcftools concat --naive: 1 for yes, 0 for no."
        split_for_bcftools_merge_csv: "The interval CSV MergeSampleSVCallsets split the scored calls into."
    }

    input {
        String remote_indir
        String remote_outdir
        File? sample_ids_file
        Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

        Int merge_mode = 1
        Int truvari_chunk_min_records = 2000
        Int truvari_collapse_refdist = 1000
        Int consistency_checks = 1
        String truvari_matching_parameters = "--refdist 500 --pctseq 0.95 --pctsize 0.95 --pctovl 0.0"
        Boolean use_bed = false
        Int chunk_ids_per_file = 100
        Int concat_all_naive = 1

        File split_for_bcftools_merge_csv

        String sv_integration_docker

        RuntimeAttr? runtime_attr_write_sample_list
        RuntimeAttr? runtime_attr_merge_chunk
        RuntimeAttr? runtime_attr_chunks_for_chromosome
        RuntimeAttr? runtime_attr_shard_chromosome
        RuntimeAttr? runtime_attr_derive_chunk_ids
        RuntimeAttr? runtime_attr_collapse_shards
        RuntimeAttr? runtime_attr_concat_chromosome_shards
        RuntimeAttr? runtime_attr_concat_chromosomes
    }

    String indir = sub(remote_indir, "/+$", "") + "/02_scoring"
    String outdir = sub(remote_outdir, "/+$", "")
    String merge_dir = outdir + "/03_merge"
    String shard_dir = outdir + "/04_shard"
    String collapse_dir = outdir + "/05_collapse"
    String concat_dir = outdir + "/06_concat"

    Int n_chunks = length(read_lines(split_for_bcftools_merge_csv))

    if (!defined(sample_ids_file)) {
        call WriteSampleList {
            input:
                remote_indir = indir,
                docker = sv_integration_docker,
                runtime_attr_override = runtime_attr_write_sample_list
        }
    }
    File sample_ids = select_first([sample_ids_file, WriteSampleList.sample_ids_file])

    scatter (chunk_id in range(n_chunks)) {
        call Helpers.BcftoolsMergeChunk as MergeChunk {
            input:
                chunk_id = chunk_id,
                sample_ids = sample_ids,
                remote_indir = indir,
                merge_mode = merge_mode,
                remote_outdir = merge_dir,
                docker = sv_integration_docker,
                runtime_attr_override = runtime_attr_merge_chunk
        }
    }

    scatter (chromosome in chromosomes) {
        call ChunksForChromosome {
            input:
                split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
                chromosome_id = chromosome,
                docker = sv_integration_docker,
                runtime_attr_override = runtime_attr_chunks_for_chromosome
        }

        call ShardChromosome {
            input:
                chromosome_id = chromosome,
                bcftools_chunks = ChunksForChromosome.chunks,
                truvari_chunk_min_records = truvari_chunk_min_records,
                truvari_collapse_refdist = truvari_collapse_refdist,
                consistency_checks = consistency_checks,
                remote_indir = merge_dir,
                remote_outdir = shard_dir,
                upstream_signal = MergeChunk.done,
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
                    truvari_matching_parameters = truvari_matching_parameters,
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
    }
}

task WriteSampleList {
    input {
        String remote_indir
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        gcloud storage ls ~{remote_indir}/'*.done' | sed 's#.*/##; s#\.done$##' | sort -u > sample_ids.txt
        wc -l sample_ids.txt 1>&2
        if [ ! -s sample_ids.txt ]; then
            echo "ERROR: no <sample>.done markers found under ~{remote_indir}."
            exit 1
        fi
    >>>

    output {
        File sample_ids_file = "sample_ids.txt"
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

task ChunksForChromosome {
    input {
        File split_for_bcftools_merge_csv
        String chromosome_id
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        awk -F, -v chr=~{chromosome_id} '$1==chr {print NR-1}' ~{split_for_bcftools_merge_csv} | paste -sd, - > chunks.txt
        if [ ! -s chunks.txt ]; then
            echo "ERROR: chromosome ~{chromosome_id} not present in the interval CSV."
            exit 1
        fi
        cat chunks.txt 1>&2
    >>>

    output {
        String chunks = read_string("chunks.txt")
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

task ShardChromosome {
    input {
        String chromosome_id
        String bcftools_chunks

        Int truvari_chunk_min_records
        Int truvari_collapse_refdist

        String remote_indir
        String remote_outdir
        Int consistency_checks = 1

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

        # Localizing all the bcftools merge chunks of the chromosome
        rm -f uri_list.txt file_list.txt
        for CHUNK in $(echo ~{bcftools_chunks} | tr ',' ' '); do
            echo ~{remote_indir}/chunk_${CHUNK}.bcf >> uri_list.txt
            echo ~{remote_indir}/chunk_${CHUNK}.bcf.csi >> uri_list.txt
            echo chunk_${CHUNK}.bcf >> file_list.txt
        done
        date 1>&2
        cat uri_list.txt | gcloud storage cp -I .
        date 1>&2
        df -h 1>&2
        rm -f uri_list.txt

        # Concatenating all the bcftools merge chunks to build a whole- chromosome VCF. This is because a truvari collapse chunk may straddle multiple bcftools merge chunks.
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list file_list.txt --output-type b --output ~{chromosome_id}.bcf
        bcftools index --threads ${N_THREADS} -f ~{chromosome_id}.bcf
        df -h 1>&2
        rm -f chunk_*.bcf* file_list.txt
        N_RECORDS=$(bcftools index --nrecords ~{chromosome_id}.bcf.csi)
        if [ ~{consistency_checks} -eq 1 ]; then
            ${TIME_COMMAND} bcftools query --format '%ID\n' ~{chromosome_id}.bcf > ids_truth.txt
        fi

        # Chunking the chromosome for truvari collapse
        ${TIME_COMMAND} bcftools query --format '%POS\t%REF\t%ALT\n' ~{chromosome_id}.bcf > pos_ref_alt.tsv
        ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB}G TruvariDivide2 pos_ref_alt.tsv ~{truvari_collapse_refdist} ~{truvari_chunk_min_records} ~{chromosome_id} ${N_RECORDS} > regions.txt
        rm -f pos_ref_alt.tsv
        cat << 'END' > chunk_by_region.sh
#!/bin/bash

INPUT_BCF=$1
REGION=$2
CHUNK_ID=$3
bcftools view --regions ${REGION} --regions-overlap pos --output-type b ${INPUT_BCF} --output chunk_${CHUNK_ID}.bcf
bcftools index -f chunk_${CHUNK_ID}.bcf
df -h 1>&2
END
        chmod +x chunk_by_region.sh
        ${TIME_COMMAND} xargs --arg-file=regions.txt --max-lines=1 --max-procs=${N_THREADS} ./chunk_by_region.sh ~{chromosome_id}.bcf
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
        ls chunk_*.bcf* > file_list.txt
        cat file_list.txt | gcloud storage cp -I ~{remote_outdir}/~{chromosome_id}/
        gcloud storage cp regions.txt ~{remote_outdir}/~{chromosome_id}/
    >>>

    output {
        # regions.txt (col 1 = region, col 2 = chunk id) is delocalized so an orchestrator can derive the truvari-collapse chunk-id lists in-graph, replacing make_workpackage7_chunk_id_files.sh.
        File regions_txt = "regions.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 12,
        disk_gb: 50,
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

        String truvari_matching_parameters
        Boolean use_bed
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

            # Remark: we do not store `removed.vcf` since it's several GBs per chunk and much bigger than the collapsed output.
            ${TIME_COMMAND} truvari collapse --sizemin 0 --sizemax ${INFINITY} --keep maxqual --gt off ~{truvari_matching_parameters} ${BED_FLAGS} --input chunk_${CHUNK_ID}_in.vcf.gz --output chunk_${CHUNK_ID}_out.vcf --removed-output /dev/null
            df -h 1>&2
            ls -laht 1>&2
            rm -f chunk_${CHUNK_ID}_in.vcf.gz* ; mv chunk_${CHUNK_ID}_out.vcf chunk_${CHUNK_ID}_in.vcf

            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type b chunk_${CHUNK_ID}_in.vcf --output chunk_${CHUNK_ID}_out.bcf
            df -h 1>&2
            ls -laht 1>&2
            rm -f chunk_${CHUNK_ID}_in.vcf ; mv chunk_${CHUNK_ID}_out.bcf chunk_${CHUNK_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f chunk_${CHUNK_ID}_in.bcf

            # Setting to `.` every ID written by truvari collapse, since these can be very long in a large cohort and needlessly inflate output size. Unique IDs genome-wide will be assigned downstream.
            ${TIME_COMMAND} bcftools annotate --remove ID --output-type b chunk_${CHUNK_ID}_in.bcf --output chunk_${CHUNK_ID}_out.bcf
            rm -f chunk_${CHUNK_ID}_in.bcf* ; mv chunk_${CHUNK_ID}_out.bcf chunk_${CHUNK_ID}_in.bcf ; bcftools index --threads ${N_THREADS} -f chunk_${CHUNK_ID}_in.bcf

            mv chunk_${CHUNK_ID}_in.bcf chunk_${CHUNK_ID}_truvari.bcf
            mv chunk_${CHUNK_ID}_in.bcf.csi chunk_${CHUNK_ID}_truvari.bcf.csi
        }

        INFINITY="1000000000"
        truvari --help 1>&2
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
        echo "done" > wp5.signal
    >>>

    output {
        String done = read_string("wp5.signal")
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
            ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list chunk_list.txt --output-type b --output out.bcf
            df -h 1>&2
            rm -rf chunk_* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf

            # Enforcing a distinct ID in every record, and annotating every record with the number of samples it occurs in. Note that the latter is not equal to the QUAL field in input to truvari collapse upstream, so we have to recompute this number.
            CHR=~{chromosome}
            CHR=${CHR#chr}
            ${TIME_COMMAND} bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%COUNT(GT="alt")\n' in.bcf | awk -v id=${CHR} 'BEGIN { FS="\t"; OFS="\t"; i=0; } { $3=sprintf("%s_%d",id,i++); print $0 }' | bgzip -c > annotations.tsv.gz
            tabix -@ ${N_THREADS} -s1 -b2 -e2 annotations.tsv.gz
            echo '##INFO=<ID=N_DISCOVERY_SAMPLES,Number=1,Type=Integer,Description="Number of samples where the record was discovered">' > header.txt
            ${TIME_COMMAND} bcftools annotate --header-lines header.txt --annotations annotations.tsv.gz --columns CHROM,POS,ID,REF,ALT,N_DISCOVERY_SAMPLES --output-type b in.bcf --output out.bcf
            df -h 1>&2
            rm -f in.bcf* ; mv out.bcf in.bcf ; bcftools index --threads ${N_THREADS} -f in.bcf
            gcloud storage cp in.bcf ~{remote_outdir}/~{chromosome}/truvari_collapsed.bcf
            gcloud storage cp in.bcf.csi ~{remote_outdir}/~{chromosome}/truvari_collapsed.bcf.csi

            touch ~{chromosome}.done
            gcloud storage mv ~{chromosome}.done ~{remote_outdir}/~{chromosome}/
        fi
        echo "~{chromosome}" > out.txt
    >>>

    output {
        File out_txt = "out.txt"
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
        bcftools index --threads ${N_THREADS} -f truvari_collapsed.bcf

        # Uploading
        gcloud storage mv truvari_collapsed.'bcf*' ~{remote_outdir}/
    >>>

    output {
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
