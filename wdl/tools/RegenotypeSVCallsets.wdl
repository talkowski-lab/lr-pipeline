# Derived from yuliamostovoy/lr_callset_integration wdl/SV_Integration_WorkflowC_Regenotype.wdl
# Derived from yuliamostovoy/lr_callset_integration wdl/SV_Integration_Workpackage7_Main_joint_genotype_families.wdl
# Fork of fabio-cunial/callset_integration_phase2
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_WorkflowC_Regenotype.wdl
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_Workpackage7_Main_joint_genotype_families.wdl

version 1.0

import "../utils/Structs.wdl"

workflow RegenotypeSVCallsets {
    meta {
        description: [
            "This tool regenotypes every sample at every site of the cohort SV callset from MergeCohortSVCallsets with Kanpig (https://github.com/ACEnglish/kanpig), against the sample's reads, so that each sample carries a called genotype at each cohort site. Each sample's regenotyped calls are split into the chunks of the interval CSV.",
            "It is run once over the cohort, as one Terra sample set, with the samples scattered. A sample whose marker file already exists in 'remote_outdir' is skipped, so a resubmission only reruns the samples that failed. Nothing is returned as a workflow output; MergeRegenotypedSVCallsets merges what the samples wrote into the regenotyped cohort callset."
        ]
    }

    parameter_meta {
        sample_ids: "Samples to regenotype, aligned by index with every other per-sample array."
        sample_sexes: "Sex of each sample. 'M' selects the male ploidy BED for Kanpig; any other value selects the female one."
        aligned_bams: "GCS paths to each sample's aligned reads."
        aligned_bais: "GCS paths to the indexes for 'aligned_bams'."
        remote_indir: "The 'remote_outdir' of MergeCohortSVCallsets, whose '06_concat' subdirectory holds the cohort callset as '<prefix>.bcf'. The workflow 'prefix' must match the one given to MergeCohortSVCallsets."
        remote_outdir: "GCS directory the per-sample regenotyped chunks are written under, without a trailing slash."
        requester_pays_project: "Project billed for reads from requester-pays buckets. Leave empty when none are read."
        kanpig_params_cohort: "Kanpig arguments for regenotyping against the cohort callset."
        split_for_bcftools_merge_csv: "The interval CSV the rest of the SV integration used."
        ref_fa: "From references."
        ref_fai: "From references."
        ploidy_bed_female: "From references."
        ploidy_bed_male: "From references."
        autosomes_bed: "Autosomes, used to report each member's heterozygous-call rate."
    }

    input {
        Array[String] sample_ids
        Array[String] sample_sexes
        Array[String] aligned_bams
        Array[String] aligned_bais
        String remote_indir
        String remote_outdir
        String prefix

        String requester_pays_project = ""
        String kanpig_params_cohort = "--neighdist 500 --gpenalty 0.04 --hapsim 0.97"
        File split_for_bcftools_merge_csv

        File ref_fa
        File ref_fai
        File ploidy_bed_female
        File ploidy_bed_male
        File autosomes_bed

        String sv_integration_docker

        RuntimeAttr? runtime_attr_regenotype_sample
    }

    String cohort_indir = sub(remote_indir, "/+$", "") + "/06_concat"
    String regenotyped_dir = sub(remote_outdir, "/+$", "")

    scatter (i in range(length(sample_ids))) {
        call RegenotypeSample {
            input:
                sample_id = sample_ids[i],
                sample_sex = sample_sexes[i],
                aligned_bam = aligned_bams[i],
                aligned_bai = aligned_bais[i],
                split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
                remote_indir = cohort_indir,
                remote_outdir = regenotyped_dir,
                prefix = prefix,
                requester_pays_project = requester_pays_project,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                ploidy_bed_female = ploidy_bed_female,
                ploidy_bed_male = ploidy_bed_male,
                autosomes_bed = autosomes_bed,
                kanpig_params_cohort = kanpig_params_cohort,
                docker = sv_integration_docker,
                runtime_attr_override = runtime_attr_regenotype_sample
        }
    }

    output {
    }
}

task RegenotypeSample {
    input {
        String sample_id
        String sample_sex
        String aligned_bam
        String aligned_bai
        File split_for_bcftools_merge_csv

        String remote_indir
        String remote_outdir
        String prefix
        String requester_pays_project

        File ref_fa
        File ref_fai
        File ploidy_bed_female
        File ploidy_bed_male
        File autosomes_bed

        String kanpig_params_cohort

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
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        export RUST_BACKTRACE="full"
        GCLOUD_STORAGE_BILLING_FLAGS=""
        if [ -n "~{requester_pays_project}" ]; then
            GCLOUD_STORAGE_BILLING_FLAGS="--billing-project=~{requester_pays_project}"
        fi

        function LocalizeSample() {
            local SAMPLE_ID=$1

            local ALIGNED_BAI="~{aligned_bai}"
            local ALIGNED_BAM="~{aligned_bam}"

            local AVAILABLE_GB=$(df -h | grep "cromwell_root" | tr -s ' ' | cut -d ' ' -f 4)
            AVAILABLE_GB=${AVAILABLE_GB%G}
            AVAILABLE_GB=${AVAILABLE_GB%.*}
            local BAM_BYTES=$(gcloud storage ls -l ${GCLOUD_STORAGE_BILLING_FLAGS} "${ALIGNED_BAM}" | awk '$1 ~ /^[0-9]+$/ { print $1; exit }')
            if [ -z "${BAM_BYTES}" ]; then
                echo "ERROR: could not determine BAM size for ${ALIGNED_BAM}."
                exit 1
            fi
            local BAM_GB=${BAM_BYTES}
            BAM_GB=$(( (${BAM_GB} + 1073741823) / 1073741824 + 5 ))
            if [ ${BAM_GB} -gt ${AVAILABLE_GB} ]; then
                echo "ERROR: the BAM is larger than the available disk space. BAM size + slack: ${BAM_GB}GB. Available disk: ${AVAILABLE_GB}GB."
                exit 1
            fi

            date 1>&2
            gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} "${ALIGNED_BAM}" ./${SAMPLE_ID}_aligned.bam
            gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} "${ALIGNED_BAI}" ./${SAMPLE_ID}_aligned.bam.bai
            date 1>&2
            touch ${SAMPLE_ID}_aligned.bam.bai
        }

        function DelocalizeSample() {
            local SAMPLE_ID=$1

            rm -f ${SAMPLE_ID}_*
        }

        function Kanpig() {
            local SAMPLE_ID=$1
            local SEX=$2

            local PLOIDY_BED
            if [ ${SEX} == "M" ]; then
                PLOIDY_BED=$(echo ~{ploidy_bed_male})
            else
                PLOIDY_BED=$(echo ~{ploidy_bed_female})
            fi

            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ${SAMPLE_ID} --output-type z cohort.bcf --output ${SAMPLE_ID}_personalized.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_personalized.vcf.gz

            # Remark: kanpig needs --sizemin >= --kmer.
            ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --sizemin 10 --sizemax ${INFINITY} ~{kanpig_params_cohort} --reference ~{ref_fa} --ploidy-bed ${PLOIDY_BED} --input ${SAMPLE_ID}_personalized.vcf.gz --reads ${SAMPLE_ID}_aligned.bam --out ${SAMPLE_ID}_out.vcf --sample ${SAMPLE_ID}
            rm -f ${SAMPLE_ID}_personalized.vcf.gz*

            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_out.vcf --output ${SAMPLE_ID}_kanpig.vcf.gz
            rm -f ${SAMPLE_ID}_out.vcf
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_kanpig.vcf.gz

            local N_RECORDS=$(bcftools index --nrecords ${SAMPLE_ID}_kanpig.vcf.gz)
            local N_PRESENT_RECORDS=$(bcftools query --format '%ID\n' --include 'GT="alt"' ${SAMPLE_ID}_kanpig.vcf.gz | wc -l)
            local PERCENT="0"
            if [ ${N_RECORDS} -gt 0 ]; then
                PERCENT=$(echo "scale=2; 100 * ${N_PRESENT_RECORDS} / ${N_RECORDS}" | bc)
            fi
            echo "${N_PRESENT_RECORDS},${N_RECORDS},${PERCENT},Number of records that are marked as ALT by kanpig" >> ${SAMPLE_ID}_kanpig.csv

            local N_HETS_IN_AUTOSOMES=$(bcftools query --format '%ID\n' --include 'GT="het"' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_kanpig.vcf.gz | wc -l)
            local N_PRESENT_RECORDS_IN_AUTOSOMES=$(bcftools query --format '%ID\n' --include 'GT="alt"' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_kanpig.vcf.gz | wc -l)
            PERCENT="0"
            if [ ${N_PRESENT_RECORDS_IN_AUTOSOMES} -gt 0 ]; then
                PERCENT=$(echo "scale=2; 100 * ${N_HETS_IN_AUTOSOMES} / ${N_PRESENT_RECORDS_IN_AUTOSOMES}" | bc)
            fi
            echo "${N_HETS_IN_AUTOSOMES},${N_PRESENT_RECORDS_IN_AUTOSOMES},${PERCENT},Number of records in autosomes that are marked as HET by kanpig" >> ${SAMPLE_ID}_kanpig.csv
            ${TIME_COMMAND} java -cp ~{docker_dir} GetKanpigWindows ${SAMPLE_ID}_kanpig.vcf.gz | bgzip > ${SAMPLE_ID}_kanpig.bed.gz
        }

        # Kanpig writes a genotype-level FORMAT/FT with integer values (0, 1, ...). FT is a reserved per-sample genotype-filter key, so htsjdk/GATK (and other strict parsers used on the final callset) reject values like "0" as invalid filter names. Rename FT -> FTK (values preserved) so downstream tools parse cleanly. Same operation as utils/FixKanpigFT.wdl. No-op if FT is absent.
        function FixKanpigFT() {
            local SAMPLE_ID=$1

            bcftools view --header-only ${SAMPLE_ID}_kanpig.vcf.gz > ${SAMPLE_ID}_ft_header.txt
            if grep -q '^##FORMAT=<ID=FT,' ${SAMPLE_ID}_ft_header.txt; then
                printf 'FORMAT/FT\tFTK\n' > ${SAMPLE_ID}_ft_rename.txt
                ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --rename-annots ${SAMPLE_ID}_ft_rename.txt --output-type z ${SAMPLE_ID}_kanpig.vcf.gz --output ${SAMPLE_ID}_ftfix.vcf.gz
                rm -f ${SAMPLE_ID}_kanpig.vcf.gz* ; mv ${SAMPLE_ID}_ftfix.vcf.gz ${SAMPLE_ID}_kanpig.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_kanpig.vcf.gz
                rm -f ${SAMPLE_ID}_ft_rename.txt
            fi
            rm -f ${SAMPLE_ID}_ft_header.txt
        }

        function ChunkAndUpload() {
            local SAMPLE_ID=$1

            local i="0"
            local INTERVAL
            while read -u 5 INTERVAL; do
                echo ${INTERVAL} | tr ',' '\t' > ${SAMPLE_ID}.bed
                ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --regions-file ${SAMPLE_ID}.bed --regions-overlap pos --output-type b ${SAMPLE_ID}_kanpig.vcf.gz --output ${SAMPLE_ID}_chunk_${i}.bcf
                bcftools index --threads ${N_THREADS} -f ${SAMPLE_ID}_chunk_${i}.bcf
                gcloud storage cp ${SAMPLE_ID}_chunk_${i}.bcf ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf
                gcloud storage cp ${SAMPLE_ID}_chunk_${i}.bcf.csi ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf.csi
                i=$(( ${i} + 1 ))
            done 5< ~{split_for_bcftools_merge_csv}
            gcloud storage cp ${SAMPLE_ID}_kanpig.bed.gz ${SAMPLE_ID}_kanpig.csv ~{remote_outdir}/
        }

        INFINITY="1000000000"
        ~{docker_dir}/kanpig --version 1>&2

        # Localizing the WP6 cohort VCF.
        ${TIME_COMMAND} gcloud storage cp ~{remote_indir}/~{prefix}.bcf cohort.bcf
        gcloud storage cp ~{remote_indir}/~{prefix}.bcf.csi cohort.bcf.csi

        TEST=$(gcloud storage ls ~{remote_outdir}/~{sample_id}.done || echo "0")
        if [ ${TEST} == "0" ]; then
            LocalizeSample ~{sample_id}
            Kanpig ~{sample_id} ~{sample_sex}
            FixKanpigFT ~{sample_id}
            ChunkAndUpload ~{sample_id}

            touch ~{sample_id}.done
            gcloud storage mv ~{sample_id}.done ~{remote_outdir}/
            DelocalizeSample ~{sample_id}
            ls -laht 1>&2
        fi

        # Completion signal for orchestrator ordering. Ignored standalone.
        echo "done" > wp7.signal
    >>>

    output {
        String done = read_string("wp7.signal")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 6,
        mem_gb: 8,
        disk_gb: 100,
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
