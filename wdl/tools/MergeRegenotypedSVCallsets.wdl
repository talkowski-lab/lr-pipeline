# Derived from yuliamostovoy/lr_callset_integration, a fork of fabio-cunial/callset_integration_phase2:
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_WorkflowE_Regenotype_Merge.wdl
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_Workpackage8_Main_concat_regenotyped_shards.wdl

version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergeRegenotypedSVCallsets {
    meta {
        description: [
            "This tool merges the per-family regenotyped SV calls written by RegenotypeFamilySVCallsets into one regenotyped cohort SV callset. Each chunk of the interval CSV is merged across every sample of every family with bcftools merge, matching records by ID and dropping records that are REF in every sample, and the merged chunks are concatenated in genome order.",
            "Nothing is returned as a workflow output. The merged chunks are written under 'merge' in 'remote_outdir', and the regenotyped cohort callset is 'concat/merged.bcf' there."
        ]
    }

    parameter_meta {
        remote_indir: "The 'remote_outdir' shared by every RegenotypeFamilySVCallsets run."
        remote_outdir: "GCS directory the merged chunks and the concatenated callset are written under, without a trailing slash."
        sample_ids_file: "Samples to merge, one per line, in the column order of the merged VCF. Derived from the per-sample marker files in 'remote_indir' when omitted."
        merge_mode: "How bcftools merge matches records: 1 by CHROM, POS, REF and ALT, or 2 by ID. Regenotyped calls share the cohort IDs, so 2 is correct here."
        split_for_bcftools_merge_csv: "The interval CSV RegenotypeFamilySVCallsets split the regenotyped calls into."
    }

    input {
        String remote_indir
        String remote_outdir
        File? sample_ids_file

        Int merge_mode = 2

        File split_for_bcftools_merge_csv

        String sv_integration_docker

        RuntimeAttr? runtime_attr_write_sample_list
        RuntimeAttr? runtime_attr_make_chunk_ids_csv
        RuntimeAttr? runtime_attr_merge_chunk
        RuntimeAttr? runtime_attr_concat_regenotyped_chunks
    }

    String indir = sub(remote_indir, "/+$", "")
    String outdir = sub(remote_outdir, "/+$", "")
    String merge_dir = outdir + "/merge"
    String concat_dir = outdir + "/concat"

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

    call MakeChunkIdsCsv {
        input:
            n_chunks = n_chunks,
            docker = sv_integration_docker,
            runtime_attr_override = runtime_attr_make_chunk_ids_csv
    }

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

    call ConcatRegenotypedChunks {
        input:
            chunk_ids = MakeChunkIdsCsv.csv,
            remote_indir = merge_dir,
            remote_outdir = concat_dir,
            upstream_signal = MergeChunk.done,
            docker = sv_integration_docker,
            runtime_attr_override = runtime_attr_concat_regenotyped_chunks
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

task MakeChunkIdsCsv {
    input {
        Int n_chunks
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        seq 0 $(( ~{n_chunks} - 1 )) | paste -sd, - > csv.txt
    >>>

    output {
        String csv = read_string("csv.txt")
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

task ConcatRegenotypedChunks {
    input {
        String chunk_ids

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
        EFFECTIVE_RAM_GB=$(( ~{ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]))} - 2 ))

        # Localizing
        for CHUNK in $(echo ~{chunk_ids} | tr ',' ' '); do
            echo ~{remote_indir}/chunk_${CHUNK}.bcf >> uri_list.txt
            echo ~{remote_indir}/chunk_${CHUNK}.bcf.csi >> uri_list.txt
            echo chunk_${CHUNK}.bcf >> file_list.txt
        done
        date 1>&2
        cat uri_list.txt | gcloud storage cp -I .
        date 1>&2
        df -h 1>&2

        # Concatenating
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list file_list.txt --output-type b --output merged.bcf
        ${TIME_COMMAND} bcftools index --threads ${N_THREADS} -f merged.bcf
        df -h 1>&2

        # Uploading
        gcloud storage mv merged.'bcf*' ~{remote_outdir}/
    >>>

    output {
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
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
