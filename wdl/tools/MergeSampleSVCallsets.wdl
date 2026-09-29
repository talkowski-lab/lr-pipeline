# Derived from yuliamostovoy/lr_callset_integration, a fork of fabio-cunial/callset_integration_phase2:
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_WorkflowA_Intrasample_Scoring.wdl
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_Workpackage1_intrasample_merge.wdl
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_Workpackage2_Main_scoring.wdl

version 1.0

import "../utils/Structs.wdl"

workflow MergeSampleSVCallsets {
    meta {
        description: [
            "This tool merges each sample's PAV, pbsv and Sniffles calls into one scored VCF per sample. Each caller's VCF is normalized and split into SVs between the minimum and maximum length, longer SVs and breakends, and the callers are merged within each class with Truvari (https://github.com/ACEnglish/truvari) collapse. The SVs within the length range are then genotyped against the sample's reads with Kanpig (https://github.com/ACEnglish/kanpig), keeping only calls with ALT support, and records matching the training resource are marked for XGBoost training.",
            "Each sample is then scored with an XGBoost model trained on those records, through GATK ScoreVariantAnnotations, and split into the chunks of the interval CSV for the cohort merge. Samples are processed in batches, one batch per VM.",
            "Nothing is returned as a workflow output. The per-sample merged, longer-SV and breakend calls are written under '01_intrasample' in 'remote_outdir', which MergeCohortLongSVCallsets reads, and the scored chunks with one marker file per sample are written under '02_scoring', which MergeCohortSVCallsets reads. A sample whose outputs already exist is skipped, so a resubmission only reruns the samples that failed."
        ]
    }

    parameter_meta {
        sample_ids: "Samples to process, aligned by index with every other per-sample array."
        sample_sexes: "Sex of each sample. 'M' selects the male ploidy BED for Kanpig; any other value selects the female one."
        aligned_bams: "GCS paths to each sample's aligned reads, streamed rather than localized."
        aligned_bais: "GCS paths to the indexes for 'aligned_bams'."
        pbsv_vcfs: "GCS paths to each sample's pbsv calls."
        pbsv_vcf_idxs: "GCS paths to the indexes for 'pbsv_vcfs'."
        sniffles_vcfs: "GCS paths to each sample's Sniffles calls."
        sniffles_vcf_idxs: "GCS paths to the indexes for 'sniffles_vcfs'."
        pav_vcfs: "GCS paths to each sample's PAV calls. Required only when 'has_pav' is true."
        pav_vcf_idxs: "GCS paths to the indexes for 'pav_vcfs'. Required only when 'has_pav' is true."
        pav_beds: "GCS paths to each sample's PAV callable-region BEDs. Required only when 'has_pav' is true."
        has_pav: "Whether to merge PAV calls alongside pbsv and Sniffles. When false, only pbsv and Sniffles are merged and PAV support is recorded as zero."
        remote_outdir: "GCS directory the per-sample outputs are written under, without a trailing slash."
        region: "Region to restrict every caller's calls to, or 'all' to keep the whole genome."
        requester_pays_project: "Project billed for reads from requester-pays buckets. Leave empty when none are read."
        min_sv_length: "Minimum SV length kept."
        max_sv_length: "Maximum length of an SV genotyped with Kanpig and scored. Longer SVs are integrated separately by MergeCohortLongSVCallsets."
        kanpig_params_singlesample: "Kanpig arguments for genotyping a single sample."
        ultralong_collapse_mode: "Whether Truvari collapse uses sequence similarity when merging SVs longer than 'max_sv_length': 0 for no, 1 for yes."
        filter_string: "bcftools expression for records to keep after scoring, or 'none' to keep every record."
        annotations: "INFO fields the XGBoost model is trained and scored on."
        training_resource_vcf: "Truth SV calls whose matches in each sample are marked as XGBoost training records."
        training_resource_vcf_idx: "Index for 'training_resource_vcf'."
        training_resource_bed: "Regions in which 'training_resource_vcf' is considered complete."
        ref_fa: "From references."
        ref_fai: "From references."
        standard_chromosomes_bed: "Chromosomes calls are restricted to."
        autosomes_bed: "Autosomes, used to report each sample's heterozygous-call rate."
        ref_agp: "Assembly gap layout of the reference, whose gaps are removed from the calls and the training regions."
        ploidy_bed_female: "From references."
        ploidy_bed_male: "From references."
        split_for_bcftools_merge_csv: "Intervals the scored calls are split into for the cohort merge, one chunk per line."
        training_python_script: "Python script GATK TrainVariantAnnotationsModel runs to train the XGBoost model."
        scoring_python_script: "Python script GATK ScoreVariantAnnotations runs to score with the XGBoost model."
        hyperparameters_json: "XGBoost hyperparameters."
        batch_size: "Number of samples processed one after another on each VM."
    }

    input {
        Array[String] sample_ids
        Array[String] sample_sexes
        Array[String] aligned_bams
        Array[String] aligned_bais
        Array[String] pbsv_vcfs
        Array[String] pbsv_vcf_idxs
        Array[String] sniffles_vcfs
        Array[String] sniffles_vcf_idxs
        Array[String] pav_vcfs = []
        Array[String] pav_vcf_idxs = []
        Array[String] pav_beds = []
        Boolean has_pav = true
        String remote_outdir

        String region = "all"
        String requester_pays_project = ""
        Int min_sv_length = 20
        Int max_sv_length = 2000
        String kanpig_params_singlesample = "--neighdist 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"
        Int ultralong_collapse_mode = 0
        String filter_string = "none"
        Array[String] annotations = ["KS_1", "KS_2", "SQ", "GQ", "DP", "AD_NON_ALT", "AD_ALL", "GT_COUNT", "SUPP_PAV", "SUPP_SNIFFLES", "SUPP_PBSV", "SVLEN"]

        File training_resource_vcf
        File training_resource_vcf_idx
        File training_resource_bed
        File ref_fa
        File ref_fai
        File standard_chromosomes_bed
        File autosomes_bed
        File ref_agp
        File ploidy_bed_female
        File ploidy_bed_male
        File split_for_bcftools_merge_csv
        File training_python_script
        File scoring_python_script
        File hyperparameters_json

        Int batch_size = 20
        String sv_integration_docker
        String xgb_scoring_docker

        RuntimeAttr? runtime_attr_make_manifests
        RuntimeAttr? runtime_attr_merge_sample_calls
        RuntimeAttr? runtime_attr_score_sample_calls
    }

    String outdir = sub(remote_outdir, "/+$", "")
    String intrasample_dir = outdir + "/01_intrasample"
    String scoring_dir = outdir + "/02_scoring"

    call MakeManifests {
        input:
            sample_ids = sample_ids,
            sample_sexes = sample_sexes,
            aligned_bais = aligned_bais,
            aligned_bams = aligned_bams,
            pbsv_vcf_idxs = pbsv_vcf_idxs,
            pbsv_vcfs = pbsv_vcfs,
            sniffles_vcf_idxs = sniffles_vcf_idxs,
            sniffles_vcfs = sniffles_vcfs,
            pav_beds = pav_beds,
            pav_vcf_idxs = pav_vcf_idxs,
            pav_vcfs = pav_vcfs,
            has_pav = has_pav,
            batch_size = batch_size,
            docker = sv_integration_docker,
            runtime_attr_override = runtime_attr_make_manifests
    }

    scatter (manifest in MakeManifests.manifests) {
        call MergeSampleCallsBatch {
            input:
                sv_integration_chunk_tsv = manifest,
                has_pav = has_pav,
                region = region,
                remote_outdir = intrasample_dir,
                requester_pays_project = requester_pays_project,
                min_sv_length = min_sv_length,
                max_sv_length = max_sv_length,
                kanpig_params_singlesample = kanpig_params_singlesample,
                ultralong_collapse_mode = ultralong_collapse_mode,
                training_resource_vcf = training_resource_vcf,
                training_resource_vcf_idx = training_resource_vcf_idx,
                training_resource_bed = training_resource_bed,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                standard_chromosomes_bed = standard_chromosomes_bed,
                autosomes_bed = autosomes_bed,
                ref_agp = ref_agp,
                ploidy_bed_female = ploidy_bed_female,
                ploidy_bed_male = ploidy_bed_male,
                docker = sv_integration_docker,
                runtime_attr_override = runtime_attr_merge_sample_calls
        }

        call ScoreSampleCallsBatch {
            input:
                sv_integration_chunk_tsv = manifest,
                split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
                filter_string = filter_string,
                remote_indir = intrasample_dir,
                remote_outdir = scoring_dir,
                training_resource_bed = training_resource_bed,
                annotations = annotations,
                training_python_script = training_python_script,
                scoring_python_script = scoring_python_script,
                hyperparameters_json = hyperparameters_json,
                upstream_signal = MergeSampleCallsBatch.done,
                docker = xgb_scoring_docker,
                runtime_attr_override = runtime_attr_score_sample_calls
        }
    }

    output {
    }
}

# Manifest columns are sliced by position inside the containers, so their order must not change
task MakeManifests {
    input {
        Array[String] sample_ids
        Array[String] sample_sexes
        Array[String] aligned_bais
        Array[String] aligned_bams
        Array[String] pbsv_vcf_idxs
        Array[String] pbsv_vcfs
        Array[String] sniffles_vcf_idxs
        Array[String] sniffles_vcfs
        Array[String] pav_beds
        Array[String] pav_vcf_idxs
        Array[String] pav_vcfs
        Boolean has_pav
        Int batch_size
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        N_ID=$(wc -l < ~{write_lines(sample_ids)})
        N_SEX=$(wc -l < ~{write_lines(sample_sexes)})
        N_BAI=$(wc -l < ~{write_lines(aligned_bais)})
        N_BAM=$(wc -l < ~{write_lines(aligned_bams)})
        N_PBSV_TBI=$(wc -l < ~{write_lines(pbsv_vcf_idxs)})
        N_PBSV_VCF=$(wc -l < ~{write_lines(pbsv_vcfs)})
        N_SNIF_TBI=$(wc -l < ~{write_lines(sniffles_vcf_idxs)})
        N_SNIF_VCF=$(wc -l < ~{write_lines(sniffles_vcfs)})
        for V in ${N_SEX} ${N_BAI} ${N_BAM} ${N_PBSV_TBI} ${N_PBSV_VCF} ${N_SNIF_TBI} ${N_SNIF_VCF}; do
            if [ ${V} -ne ${N_ID} ]; then
                echo "ERROR: a per-sample column has ${V} rows != ${N_ID} sample_ids."
                exit 1
            fi
        done

        if [ ~{true="1" false="0" has_pav} -eq 1 ]; then
            N_PAV_BED=$(wc -l < ~{write_lines(pav_beds)})
            N_PAV_TBI=$(wc -l < ~{write_lines(pav_vcf_idxs)})
            N_PAV_VCF=$(wc -l < ~{write_lines(pav_vcfs)})
            for V in ${N_PAV_BED} ${N_PAV_TBI} ${N_PAV_VCF}; do
                if [ ${V} -ne ${N_ID} ]; then
                    echo "ERROR: has_pav=true but a PAV column has ${V} rows != ${N_ID} sample_ids."
                    exit 1
                fi
            done
            paste ~{write_lines(sample_ids)} ~{write_lines(sample_sexes)} \
                  ~{write_lines(aligned_bais)} ~{write_lines(aligned_bams)} \
                  ~{write_lines(pav_beds)} ~{write_lines(pav_vcf_idxs)} ~{write_lines(pav_vcfs)} \
                  ~{write_lines(pbsv_vcf_idxs)} ~{write_lines(pbsv_vcfs)} \
                  ~{write_lines(sniffles_vcf_idxs)} ~{write_lines(sniffles_vcfs)} > all.tsv
        else
            paste ~{write_lines(sample_ids)} ~{write_lines(sample_sexes)} \
                  ~{write_lines(aligned_bais)} ~{write_lines(aligned_bams)} \
                  ~{write_lines(pbsv_vcf_idxs)} ~{write_lines(pbsv_vcfs)} \
                  ~{write_lines(sniffles_vcf_idxs)} ~{write_lines(sniffles_vcfs)} > all.tsv
        fi

        split --lines=~{batch_size} --numeric-suffixes=0 --suffix-length=6 all.tsv batch_
        ls batch_* 1>&2
    >>>

    output {
        Array[File] manifests = glob("batch_*")
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

task MergeSampleCallsBatch {
    input {
        File sv_integration_chunk_tsv
        Boolean has_pav
        String region
        String remote_outdir
        String requester_pays_project

        Int min_sv_length
        Int max_sv_length
        String kanpig_params_singlesample
        Int ultralong_collapse_mode

        File training_resource_vcf
        File training_resource_vcf_idx
        File training_resource_bed

        File ref_fa
        File ref_fai
        File standard_chromosomes_bed
        File autosomes_bed
        File ref_agp
        File ploidy_bed_female
        File ploidy_bed_male
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
        HAS_PAV="~{has_pav}"
        GCLOUD_STORAGE_BILLING_FLAGS=""
        if [ -n "~{requester_pays_project}" ]; then
            GCLOUD_STORAGE_BILLING_FLAGS="--billing-project=~{requester_pays_project}"
        fi

        # Applies `bcftools annotate`, transferring tags matched on the record ID. `~ID` only takes effect when REF,ALT are ALSO present in the annotation file and in --columns; otherwise bcftools silently falls back to CHROM,POS-only matching and mis-assigns tags across records that share a start coordinate. This wrapper injects REF,ALT (columns 4,5, looked up by the unique record ID) into the pre-built annotation TSV and splices REF,ALT into the column list, so the ~ID match engages. @param 1 bgzip'd annotation TSV: CHROM POS ID <tags...> (position-sorted). @param 2 header-lines file.  @param 3 --columns string containing ",~ID,". @param 4 input VCF (must have UNIQUE ids).  @param 5 output. @param 6 --output-type value (v/z/b).
        function AnnotateById() {
            local TSV=$1
            local HDR=$2
            local COLS=$3
            local IN=$4
            local OUT=$5
            local OFMT=$6

            bcftools query --format '%ID\t%REF\t%ALT\n' ${IN} | sort -t $'\t' -k1,1 > refalt_by_id.tsv
            # Always emit REF and ALT as two separate fields (default '.' if an ID is somehow absent from the map). A collapsed/empty ra[$3] would shift the columns and can segfault `bcftools annotate` below.
            bgzip -dc ${TSV} | awk 'BEGIN { FS="\t"; OFS="\t"; } NR==FNR { ra[$1]=$2 FS $3; next } { ref="."; alt="."; if ($3 in ra) { split(ra[$3], a, "\t"); ref=a[1]; alt=a[2]; } line=$1 FS $2 FS $3 FS ref FS alt; for (k=4; k<=NF; k++) line=line FS $k; print line }' refalt_by_id.tsv - | bgzip > refalt.${TSV}
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 refalt.${TSV}
            # Single-threaded annotate: `bcftools annotate --threads` with ~ID matching has been observed to segfault on large inputs (ONT samples). The annotate core is single-threaded anyway, so --threads only risks the crash without a speedup.
            bcftools annotate --annotations refalt.${TSV} --header-lines ${HDR} --columns "${COLS/,~ID,/,~ID,REF,ALT,}" --output-type ${OFMT} ${IN} --output ${OUT}
            rm -f refalt_by_id.tsv refalt.${TSV} refalt.${TSV}.tbi
        }

        # Remark: localizing the BAM could be avoided by running kanpig on the remote BAM (likely only on the windows it considers). This does not work on Google Cloud yet. @param $2 1=Localizes everything except the BAM. 2=Localizes just the BAM. $3 A row of `sv_integration_chunk_tsv`.
        function LocalizeSample() {
            local SAMPLE_ID=$1
            local MODE=$2
            local LINE=$3

            local ALIGNED_BAI=$(echo ${LINE} | cut -d , -f 3)
            local ALIGNED_BAM=$(echo ${LINE} | cut -d , -f 4)
            local PAV_BED=$(echo ${LINE} | cut -d , -f 5)
            local PAV_TBI=$(echo ${LINE} | cut -d , -f 6)
            local PAV_VCF_GZ=$(echo ${LINE} | cut -d , -f 7)
            local PBSV_TBI=$(echo ${LINE} | cut -d , -f 8)
            local PBSV_VCF_GZ=$(echo ${LINE} | cut -d , -f 9)
            local SNIFFLES_TBI=$(echo ${LINE} | cut -d , -f 10)
            local SNIFFLES_VCF_GZ=$(echo ${LINE} | cut -d , -f 11)
            if [ ${HAS_PAV} = "false" ]; then
                local N_FIELDS=$(echo ${LINE} | awk -F ',' '{ print NF }')
                if [ ${N_FIELDS} -lt 11 ]; then
                    PBSV_TBI=$(echo ${LINE} | cut -d , -f 5)
                    PBSV_VCF_GZ=$(echo ${LINE} | cut -d , -f 6)
                    SNIFFLES_TBI=$(echo ${LINE} | cut -d , -f 7)
                    SNIFFLES_VCF_GZ=$(echo ${LINE} | cut -d , -f 8)
                fi
            fi

            if [ ${MODE} -eq 2 ]; then
                date 1>&2
                gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${ALIGNED_BAM} ./${SAMPLE_ID}_aligned.bam
                date 1>&2
                gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${ALIGNED_BAI} ./${SAMPLE_ID}_aligned.bam.bai
            else
                if [ ${HAS_PAV} = "true" ]; then
                    gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${PAV_VCF_GZ} ./${SAMPLE_ID}_pav.vcf.gz
                    gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${PAV_TBI} ./${SAMPLE_ID}_pav.vcf.gz.tbi
                fi
                gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${PBSV_VCF_GZ} ./${SAMPLE_ID}_pbsv.vcf.gz
                gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${PBSV_TBI} ./${SAMPLE_ID}_pbsv.vcf.gz.tbi
                gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${SNIFFLES_VCF_GZ} ./${SAMPLE_ID}_sniffles.vcf.gz
                gcloud storage cp ${GCLOUD_STORAGE_BILLING_FLAGS} ${SNIFFLES_TBI} ./${SAMPLE_ID}_sniffles.vcf.gz.tbi
            fi
        }

        # Deletes all and only the files downloaded by `LocalizeSample()`.
        function DelocalizeSample() {
            local SAMPLE_ID=$1

            rm -f ${SAMPLE_ID}_aligned.bam* ${SAMPLE_ID}_pav.vcf.gz* ${SAMPLE_ID}_pbsv.vcf.gz* ${SAMPLE_ID}_sniffles.vcf.gz*
        }

        # Builds a BED file that excludes every gap from the AGP file of the reference.
        function GetReferenceGaps() {
            # Computing non-gap regions
            awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    if ( ( $1=="chr1" || $1=="chr2" || $1=="chr3" || $1=="chr4" || $1=="chr5" || $1=="chr6" || $1=="chr7" || $1=="chr8" || $1=="chr9" || $1=="chr10" || \
                           $1=="chr11" || $1=="chr12" || $1=="chr13" || $1=="chr14" || $1=="chr15" || $1=="chr16" || $1=="chr17" || $1=="chr18" || $1=="chr19" || $1=="chr20" || \
                           $1=="chr21" || $1=="chr22" || $1=="chrX" || $1=="chrY" || $1=="chrM" \
                         ) && $5=="N" \
                       ) print $0 \
                 }' ~{ref_agp} > gaps_unsorted.bed
            bedtools sort -i gaps_unsorted.bed -faidx ~{ref_fai} > gaps.bed
            bedtools complement -L -i gaps.bed -g ~{ref_fai} > not_gaps.bed

            # Intersecting non-gap regions with the training BED
            bedtools sort -i ~{training_resource_bed} -faidx ~{ref_fai} > training_resource_sorted.bed
            rm -f training_not_gaps_beds.wsv
            local ID="0"
            local ROW
            while read -u 4 ROW; do
                ID=$(( ${ID} + 1 ))
                echo "${ROW}" > ${ID}.bed
                bedtools intersect -a ${ID}.bed -b training_resource_sorted.bed -sorted -g ~{ref_fai} > training_not_gaps_${ID}.bed
                if [ -s training_not_gaps_${ID}.bed ]; then
                    echo "${ID} training_not_gaps_${ID}.bed" >> training_not_gaps_beds.wsv
                else
                    rm -f training_not_gaps_${ID}.bed
                fi
                rm -f ${ID}.bed
            done 4< not_gaps.bed
            ls -lht *.bed 1>&2

            # Removing temporary files
            rm -f gaps_unsorted.bed training_resource_sorted.bed
        }

        # Puts in canonical form a raw VCF from an SV caller. The procedure creates sorted output files `SAMPLEID_CALLERID_X.vcf.gz`, where X is: sv: non-BND records with length in [MIN_SV_LENGTH..MAX_SV_LENGTH], in canonical form; sv_ultralong: non-BND records with length >MAX_SV_LENGTH, devoid of sequence where possible to save space; bnd: BND records, in their original form. Remark: the funtion outputs indexed `.vcf.gz` files, since they are needed by `bcftools merge`.
        function CanonizeVcf() {
            local INPUT_VCF_GZ=$1
            local INPUT_TBI=$2
            local SAMPLE_ID=$3
            local CALLER_ID=$4
            local MIN_SV_LENGTH=$5
            local MAX_SV_LENGTH=$6
            local STANDARD_CHROMOSOMES_BED=$7
            local NOT_GAPS_BED=$8

            # QUAL is used by truvari collapse to select a representation. We assign values based on which representations we observed to be more accurate in a few test examples.
            if [ ${CALLER_ID} = 'pav' ]; then
                local QUAL="4"
            elif [ ${CALLER_ID} = 'pbsv' ]; then
                local QUAL="3"
            elif [ ${CALLER_ID} = 'sniffles' ]; then
                local QUAL="2"
            fi

            mv ${INPUT_VCF_GZ} ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            mv ${INPUT_TBI} ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi

            # Subsetting to standard chromosomes, or to a given region if any.
            if [ ~{region} != "all" ]; then
                ${TIME_COMMAND} bcftools view --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ~{region} --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            else
                ${TIME_COMMAND} bcftools filter --regions-file ${STANDARD_CHROMOSOMES_BED} --regions-overlap pos --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            fi

            # Removing records in reference gaps
            ${TIME_COMMAND} bcftools filter --regions-file ${NOT_GAPS_BED} --regions-overlap pos --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # Removing uncalled ALT alleles before SVLEN is treated as Number=A. This prevents malformed multiallelic records with too few SVLEN values from crashing bcftools norm.
            ${TIME_COMMAND} bcftools view --output-type u --min-ac 1 --trim-alt-alleles ${SAMPLE_ID}_${CALLER_ID}_in.vcf | bcftools +fill-tags --output-type v --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf -- -t AC,AN
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # Ensuring that SVLEN has the correct type for bcftools norm
            bcftools view --header-only ${SAMPLE_ID}_${CALLER_ID}_in.vcf | sed 's/ID=SVLEN,Number=.,/ID=SVLEN,Number=A,/g' > ${SAMPLE_ID}_${CALLER_ID}_header.txt
            ${TIME_COMMAND} bcftools reheader --header ${SAMPLE_ID}_${CALLER_ID}_header.txt --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # Splitting multiallelic records into biallelic records
            ${TIME_COMMAND} bcftools norm --multiallelics -any --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # Removing any remaining uncalled ALT alleles after the multiallelic split
            ${TIME_COMMAND} bcftools view --output-type v -i 'GT=="alt"' ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # Removing SNVs, if any.
            if [ ${CALLER_ID} = 'pav' ]; then
                ${TIME_COMMAND} bcftools filter --exclude 'SVTYPE="SNV"' --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            fi

            # Making sure SVLEN and SVTYPE are consistently annotated
            ${TIME_COMMAND} java -cp ~{docker_dir} AddSvtypeSvlen ${SAMPLE_ID}_${CALLER_ID}_in.vcf > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # Isolating BNDs
            ${TIME_COMMAND} bcftools filter --include 'SVTYPE="BND"' --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf
            ${TIME_COMMAND} bcftools filter --exclude 'SVTYPE="BND"' --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # Isolating ultra-long records and discarding short records
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>'${MAX_SV_LENGTH} --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf
            ${TIME_COMMAND} bcftools filter --include 'ABS(SVLEN)>='${MIN_SV_LENGTH}' && ABS(SVLEN)<='${MAX_SV_LENGTH} --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # 1. Main VCF ------------------------------------------------------

            # 1.1 Sorting
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # 1.2 Fixing symbolic records
            ${TIME_COMMAND} java -cp ~{docker_dir} -Xmx${EFFECTIVE_RAM_GB}G FixSymbolicRecords ${SAMPLE_ID}_${CALLER_ID}_in.vcf ~{ref_fa} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # 1.3 Fixing REF
            ${TIME_COMMAND} bcftools norm --check-ref s --fasta-ref ~{ref_fa} --do-not-normalize --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # 1.4 Cleaning REF, ALT, QUAL, FILTER. - REF and ALT must be uppercase for XGBoost scoring downstream to work. - QUAL is used by truvari collapse to select a representation. Symbolic records are NOT given low quality (it was 1 in Phase 1) since e.g. all DEL records made by Sniffles are symbolic. - We force every record to PASS, to rule out any filter-dependent effect in downstream tools.
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanRefAltQual ${SAMPLE_ID}_${CALLER_ID}_in.vcf ${QUAL} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # 1.5 Removing END, since its values may be inconsistent and make GATK crash downstream.
            ${TIME_COMMAND} bcftools annotate --remove INFO/END --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # 1.6 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --remove-duplicates --output-type z ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz

            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_sv.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_sv.vcf.gz.tbi

            # 2. BND VCF -------------------------------------------------------

            # 2.1 Sorting
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # Remark: we do not run the following command, since it seems to destroy BNDs ALTs (example: N]chr5:181473415] -> GNcNNNNNNNNNNNNNN ): bcftools norm --check-ref s --fasta-ref ~{ref_fa} --do-not-normalize

            # 2.2 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --rm-dup exact --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # 2.3 Forcing every record to PASS and adding QUAL, since it is used by `truvari collapse` to select a representation.
            ${TIME_COMMAND} java -cp ~{docker_dir} CleanQual ${SAMPLE_ID}_${CALLER_ID}_in.vcf ${QUAL} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # 2.4 Setting to 0/1 every non-ALT record, otherwise the corresponding truvari collapse SUPP field becomes zero. REMOVING THIS - we don't want to keep these variants ${TIME_COMMAND} bcftools +setGT --output-type z --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf -- --target-gt q --include 'GT="ref" || GT="mis"' --new-gt c:0/1 rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz

            # Step 2.4 (setGT) above used to compress+index _in.vcf into _in.vcf.gz; with it removed, bgzip+index the plain VCF here so the rename below (and the indexed .vcf.gz that bcftools merge needs) still works.
            bgzip ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_bnd.vcf.gz.tbi

            # 3. Ultralong VCF -------------------------------------------------

            # 3.1 Sorting
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # 3.2 Removing duplicated records
            ${TIME_COMMAND} bcftools norm --remove-duplicates --output-type v ${SAMPLE_ID}_${CALLER_ID}_in.vcf --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf

            # 3.3 Removing sequence (lossless), forcing every record to PASS, and setting QUAL, since it is used by `truvari collapse` to select a representation.
            if [ ~{ultralong_collapse_mode} -eq 0 ]; then
                ${TIME_COMMAND} java -cp ~{docker_dir} RemoveRefAlt ${SAMPLE_ID}_${CALLER_ID}_in.vcf ${QUAL} ~{ref_fai} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            elif [ ~{ultralong_collapse_mode} -eq 1 ]; then
                ${TIME_COMMAND} java -cp ~{docker_dir} CleanQual ${SAMPLE_ID}_${CALLER_ID}_in.vcf ${QUAL} > ${SAMPLE_ID}_${CALLER_ID}_out.vcf
                rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf ${SAMPLE_ID}_${CALLER_ID}_in.vcf
            fi

            # 3.4 Setting to 0/1 every non-ALT record, otherwise the corresponding truvari collapse SUPP field becomes zero.
            ${TIME_COMMAND} bcftools +setGT --output-type z --output ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf -- --target-gt q --include 'GT="ref" || GT="mis"' --new-gt c:0/1
            rm -f ${SAMPLE_ID}_${CALLER_ID}_in.vcf ; mv ${SAMPLE_ID}_${CALLER_ID}_out.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz

            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz
            mv ${SAMPLE_ID}_${CALLER_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_${CALLER_ID}_ultralong.vcf.gz.tbi
        }

        # Collapses with truvari all files `SAMPLEID_CALLERID_sv.vcf.gz`, creating an output file `SAMPLEID_sv.vcf.gz`. Remark: the funtion's inputs are indexed `.vcf.gz`, since they are needed by `bcftools merge`. It outputs a `.vcf.gz` since it's needed downstream.
        function IntrasampleMerge_sv() {
            local SAMPLE_ID=$1

            # Remark: the order of the callers in `bcftools merge` affects the value of the SAMPLE column emitted by `truvari collapse --intra`.
            if [ ${HAS_PAV} = "true" ]; then
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type z ${SAMPLE_ID}_pav_sv.vcf.gz ${SAMPLE_ID}_pbsv_sv.vcf.gz ${SAMPLE_ID}_sniffles_sv.vcf.gz --output ${SAMPLE_ID}_out.vcf
            else
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type z ${SAMPLE_ID}_pbsv_sv.vcf.gz ${SAMPLE_ID}_sniffles_sv.vcf.gz --output ${SAMPLE_ID}_out.vcf
            fi
            rm -f ${SAMPLE_ID}_*_sv.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics -any --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz

            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_in.vcf.gz --intra --keep maxqual --refdist 500 --pctseq 0.90 --pctsize 0.90 --sizemin 0 --sizemax ${INFINITY} --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            # Ensuring that every record has a unique ID, to enable joining by CHROM,POS,ID in downstream calls to `bcftools annotate`. Using CHROM,POS,REF,ALT can make `bcftools annotate` segfault, and the speed of joining by CHROM,POS,ID is independent of SVLEN. Remark: we preserve the original ID just for debugging reasons.
            (bcftools view --header-only ${SAMPLE_ID}_in.vcf ; bcftools view --no-header ${SAMPLE_ID}_in.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { printf("%s\t%s\t%d-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,++i,$3,$4,$5,$6,$7,$8,$9,$10); }') | bgzip --compress-level 1 > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            (bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | head -n 1 || echo "0") 1>&2

            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_sv.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_sv.vcf.gz.tbi
        }

        # Collapses with truvari all files `SAMPLEID_CALLERID_ultralong.vcf.gz`, creating an output file `SAMPLEID_ultralong.vcf.gz`. Remark: if `truvari collapse` is run without taking sequence similarity into account, different INS/DUP/CNV sequences of similar length at similar POS may be wrongly collapsed. We tolerate this for speed reasons. Remark: the function's inputs are indexed `.vcf.gz`, since they are needed by `bcftools merge`. It outputs a `.vcf.gz` since it's needed downstream.
        function IntrasampleMerge_ultralong() {
            local SAMPLE_ID=$1

            # Remark: the order of the callers in `bcftools merge` affects the value of the SAMPLE column emitted by `truvari collapse --intra`.
            if [ ${HAS_PAV} = "true" ]; then
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type v ${SAMPLE_ID}_pav_ultralong.vcf.gz ${SAMPLE_ID}_pbsv_ultralong.vcf.gz ${SAMPLE_ID}_sniffles_ultralong.vcf.gz --output ${SAMPLE_ID}_out.vcf
            else
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type v ${SAMPLE_ID}_pbsv_ultralong.vcf.gz ${SAMPLE_ID}_sniffles_ultralong.vcf.gz --output ${SAMPLE_ID}_out.vcf
            fi
            rm -f ${SAMPLE_ID}_*_ultralong.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics -any --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            # Removing SVLEN from symbolic ALTs, in order not to interfere with `truvari collapse`.
            local PCTSEQ_VALUE
            if [ ~{ultralong_collapse_mode} -eq 0 ]; then
                bcftools view --header-only ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
                ${TIME_COMMAND} bcftools view --no-header ${SAMPLE_ID}_in.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    if (substr($0,1,1)!="#" && substr($5,1,1)=="<") $5 = substr($5,1,4) ">"; \
                    printf("%s",$1); \
                    for (i=2; i<=NF; i++) printf("\t%s",$i); \
                    printf("\n"); \
                }' >> ${SAMPLE_ID}_out.vcf
                rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
                PCTSEQ_VALUE="0"
            elif [ ~{ultralong_collapse_mode} -eq 1 ]; then
                PCTSEQ_VALUE="0.90"
            fi

            bgzip --compress-level 1 ${SAMPLE_ID}_in.vcf ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_in.vcf.gz --intra --keep maxqual --refdist 500 --pctseq ${PCTSEQ_VALUE} --pctsize 0.90 --sizemin 0 --sizemax ${INFINITY} --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            # Adding SVLEN back into symbolic ALTs, to avoid overcollapse in the cohort-level bcftools merge downstream.
            if [ ~{ultralong_collapse_mode} -eq 0 ]; then
                ${TIME_COMMAND} java -cp ~{docker_dir} AddSvlenToSymbolicAlt ${SAMPLE_ID}_in.vcf > ${SAMPLE_ID}_out.vcf
                rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf
            fi

            # Ensuring that every record has a unique ID. Remark: we preserve the original ID just for debugging reasons.
            (bcftools view --header-only ${SAMPLE_ID}_in.vcf ; bcftools view --no-header ${SAMPLE_ID}_in.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { printf("%s\t%s\t%d-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,++i,$3,$4,$5,$6,$7,$8,$9,$10); }') | bgzip --compress-level 1 > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            (bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | head -n 1 || echo "0") 1>&2

            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_ultralong.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_ultralong.vcf.gz.tbi
        }

        # Collapses with truvari all files `SAMPLEID_CALLERID_bnd.vcf.gz`, creating an output file `SAMPLEID_bnd.vcf.gz`. Remark: the funtion's inputs are indexed `.vcf.gz`, since they are needed by `bcftools merge`. It outputs a `.vcf.gz` since it's needed downstream.
        function IntrasampleMerge_bnd() {
            local SAMPLE_ID=$1

            # Remark: the order of the callers in `bcftools merge` affects the value of the SAMPLE column emitted by `truvari collapse --intra`.
            if [ ${HAS_PAV} = "true" ]; then
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type v ${SAMPLE_ID}_pav_bnd.vcf.gz ${SAMPLE_ID}_pbsv_bnd.vcf.gz ${SAMPLE_ID}_sniffles_bnd.vcf.gz --output ${SAMPLE_ID}_out.vcf
            else
                ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --output-type v ${SAMPLE_ID}_pbsv_bnd.vcf.gz ${SAMPLE_ID}_sniffles_bnd.vcf.gz --output ${SAMPLE_ID}_out.vcf
            fi
            rm -f ${SAMPLE_ID}_*_bnd.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics -any --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz

            ${TIME_COMMAND} truvari collapse --input ${SAMPLE_ID}_in.vcf.gz --intra --keep maxqual --refdist 500 --pctseq 0.90 --pctsize 0.90 --sizemin 0 --sizemax ${INFINITY} --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type v ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            # Ensuring that every record has a unique ID Remark: we preserve the original ID just for debugging reasons.
            (bcftools view --header-only ${SAMPLE_ID}_in.vcf ; bcftools view --no-header ${SAMPLE_ID}_in.vcf | awk 'BEGIN { FS="\t"; OFS="\t"; i=0; } { printf("%s\t%s\t%d-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,++i,$3,$4,$5,$6,$7,$8,$9,$10); }') | bgzip --compress-level 1 > ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            (bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | head -n 1 || echo "0") 1>&2

            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_bnd.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_bnd.vcf.gz.tbi
        }

        # Copies truvari's SUPP field from SAMPLE to three tags in INFO. This is necessary, since kanpig overwrites the SAMPLE column. Remark: the funtion requires an indexed `.vcf.gz` in input, and it outputs an indexed `.vcf.gz` of `bcf`, depending on `OUTPUT_FORMAT` (`z` or `b`).
        function CopySuppToInfo() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local OUTPUT_FORMAT=$3
            local OUTPUT_VCF_GZ=$4

            if [ ${HAS_PAV} = "true" ]; then
                bcftools query --format '%CHROM\t%POS\t%ID\t[%SUPP]\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    printf("%s",$1); \
                    for (i=2; i<=NF-1; i++) printf("\t%s",$i); \
                    if ($4=="0") printf("\t0\t0\t0");
                    else if ($4=="1") printf("\t0\t0\t1");
                    else if ($4=="2") printf("\t0\t1\t0");
                    else if ($4=="3") printf("\t0\t1\t1");
                    else if ($4=="4") printf("\t1\t0\t0");
                    else if ($4=="5") printf("\t1\t0\t1");
                    else if ($4=="6") printf("\t1\t1\t0");
                    else if ($4=="7") printf("\t1\t1\t1");
                    printf("\n"); \
                }' | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
            else
                bcftools query --format '%CHROM\t%POS\t%ID\t[%SUPP]\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                    printf("%s",$1); \
                    for (i=2; i<=NF-1; i++) printf("\t%s",$i); \
                    if ($4=="0") printf("\t0\t0\t0");
                    else if ($4=="1") printf("\t1\t0\t0");
                    else if ($4=="2") printf("\t0\t1\t0");
                    else if ($4=="3") printf("\t1\t1\t0");
                    printf("\n"); \
                }' | bgzip -c > ${SAMPLE_ID}_annotations.tsv.gz
            fi
            tabix -@ ${N_THREADS} -f -s1 -b2 -e2 ${SAMPLE_ID}_annotations.tsv.gz
            echo '##INFO=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by pav">' > ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by sniffles">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> ${SAMPLE_ID}_header.txt
            # Remark: the order of the callers is now the reverse of the one in which they were bcftools-merged.
            AnnotateById ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt "CHROM,POS,~ID,INFO/SUPP_SNIFFLES,INFO/SUPP_PBSV,INFO/SUPP_PAV" ${INPUT_VCF_GZ} ${OUTPUT_VCF_GZ} ${OUTPUT_FORMAT}
            if [ ${OUTPUT_FORMAT} = z ]; then
                bcftools index --threads ${N_THREADS} -f -t ${OUTPUT_VCF_GZ}
            elif [ ${OUTPUT_FORMAT} = b ]; then
                bcftools index --threads ${N_THREADS} -f -c ${OUTPUT_VCF_GZ}
            fi

            rm -f ${SAMPLE_ID}_annotations.tsv.gz ${SAMPLE_ID}_header.txt ${INPUT_VCF_GZ}*
        }

        # Remark: the function outputs a `.vcf.gz`, since it's needed by the following steps.
        function Kanpig() {
            local SAMPLE_ID=$1
            local SEX=$2
            local INPUT_VCF=$3
            local ALIGNMENTS_BAM=$4

            if [ ${SEX} == "M" ]; then
                PLOIDY_BED=$(echo ~{ploidy_bed_male})
            else
                PLOIDY_BED=$(echo ~{ploidy_bed_female})
            fi

            # Remark: kanpig needs --sizemin >= --kmer
            ${TIME_COMMAND} ~{docker_dir}/kanpig gt --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_singlesample} --sizemin 10 --sizemax ${INFINITY} --reference ~{ref_fa} --input ${INPUT_VCF} --reads ${ALIGNMENTS_BAM} --out ${SAMPLE_ID}_out.vcf
            rm -f ${INPUT_VCF} ; mv ${SAMPLE_ID}_out.vcf ${SAMPLE_ID}_in.vcf

            # Sorting
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${SAMPLE_ID}_in.vcf --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz

            # Discarding records that are not marked as present by kanpig
            local N_RECORDS_BEFORE_KANPIG=$( bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz.tbi )
            ${TIME_COMMAND} bcftools filter --include 'GT="alt"' --output-type z ${SAMPLE_ID}_in.vcf.gz --output ${SAMPLE_ID}_out.vcf.gz
            rm -f ${SAMPLE_ID}_in.vcf.gz* ; mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz ; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            local N_RECORDS_AFTER_KANPIG=$( bcftools index --nrecords ${SAMPLE_ID}_in.vcf.gz.tbi )

            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_kanpig.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_kanpig.vcf.gz.tbi

            # Printing debug information
            local PERCENT=$( echo "scale=2; 100 * ${N_RECORDS_AFTER_KANPIG} / ${N_RECORDS_BEFORE_KANPIG}" | bc )
            echo "${N_RECORDS_AFTER_KANPIG},${N_RECORDS_BEFORE_KANPIG},${PERCENT},Number of records that are marked as ALT by kanpig" > ${SAMPLE_ID}_kanpig.csv
            local N_HETS_IN_AUTOSOMES=$( bcftools query --format '%ID' --include 'GT="het"' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_kanpig.vcf.gz | wc -l )
            local N_RECORDS_IN_AUTOSOMES=$( bcftools query --format '%ID' --regions-file ~{autosomes_bed} --regions-overlap pos ${SAMPLE_ID}_kanpig.vcf.gz | wc -l )
            local PERCENT=$( echo "scale=2; 100 * ${N_HETS_IN_AUTOSOMES} / ${N_RECORDS_IN_AUTOSOMES}" | bc )
            echo "${N_HETS_IN_AUTOSOMES},${N_RECORDS_IN_AUTOSOMES},${PERCENT},Number of records in autosomes that are marked as HET by kanpig" >> ${SAMPLE_ID}_kanpig.csv
            ${TIME_COMMAND} java -cp ~{docker_dir} GetKanpigWindows ${SAMPLE_ID}_kanpig.vcf.gz | bgzip > ${SAMPLE_ID}_kanpig.bed.gz
        }

        # Copies the following kanpig fields from SAMPLE to INFO: KS_1, KS_2, SQ, GQ, DP, AD_NON_ALT, AD_ALL This is necessary, since XGBoost downstream uses only INFO fields. Remark: the funtion requires an indexed `.vcf.gz` in input. It outputs a `.vcf.gz` since it's needed by the following steps.
        function CopyKanpigFieldsToInfo() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2

            # Creating new header lines
            touch ${SAMPLE_ID}_header.txt
            for FIELD in SQ GQ DP
            do
                bcftools view --header-only ${INPUT_VCF_GZ} | grep ID="${FIELD}," | sed -e 's/FORMAT/INFO/g' >> ${SAMPLE_ID}_header.txt
            done
            echo '##INFO=<ID=AD_NON_ALT,Number=1,Type=Integer,Description="Coverage for non-alternate alleles">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=AD_ALL,Number=1,Type=Integer,Description="Coverage for all alleles">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=KS_1,Number=1,Type=Integer,Description="Kanpig score 1">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=KS_2,Number=1,Type=Integer,Description="Kanpig score 2">' >> ${SAMPLE_ID}_header.txt
            echo '##INFO=<ID=GT_COUNT,Number=1,Type=Integer,Description="GT converted to an integer in {0,1,2}.">' >> ${SAMPLE_ID}_header.txt

            # Copying fields from FORMAT to INFO. Every record is assumed to have a distinct ID, which is enforced by the steps of the pipeline upstream.
            bcftools query -f '%CHROM\t%POS\t%ID\t[%KS]\t[%SQ]\t[%GQ]\t[%DP]\t[%AD]\t[%GT]\t%INFO/SUPP_PBSV\t%INFO/SUPP_SNIFFLES\t%INFO/SUPP_PAV\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                KS_1=-1; KS_2=-1; \
                p=0; \
                for (i=1; i<=length($4); i++) { \
                    if (substr($4,i,1)==",") { p=i; break; } \
                } \
                if (p==0) { KS_1=$4; KS_2=$4; } \
                else { KS_1=substr($4,1,p-1); KS_2=substr($4,p+1); } \
                if (KS_1==".") KS_1=-1; \
                if (KS_2==".") KS_2=-1; \
                \
                SQ=$5; \
                if (SQ==".") SQ=-1; \
                \
                GQ=$6; \
                if (GQ==".") GQ=-1; \
                \
                DP=$7; \
                if (DP==".") DP=-1; \
                \
                AD_NON_ALT=-1; AD_ALL=1; \
                p=0; \
                for (i=1; i<=length($8); i++) { \
                    if (substr($8,i,1)==",") { p=i; break; } \
                } \
                if (p==0) { AD_NON_ALT=$8; AD_ALL=$8; } \
                else { AD_NON_ALT=substr($8,1,p-1); AD_ALL=substr($8,p+1); } \
                if (AD_NON_ALT==".") AD_NON_ALT=-1; \
                if (AD_ALL==".") AD_ALL=-1; \
                \
                GT_COUNT=-1; \
                if ($9=="0/0" || $9=="0|0" || $9=="./."  || $9==".|." || $9=="./0" || $9==".|0" || $9=="0/." || $9=="0|." || $9=="0" || $9==".") GT_COUNT=0; \
                else if ($9=="0/1" || $9=="0|1" || $9=="1/0" || $9=="1|0" || $9=="./1" || $9==".|1" || $9=="1/." || $9=="1|." || $9=="1") GT_COUNT=1; \
                else if ($9=="1/1" || $9=="1|1") GT_COUNT=2; \
                \
                printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$1,$2,$3,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT,$10,$11,$12); \
            }' | bgzip -c > ${SAMPLE_ID}_format.tsv.gz
            tabix -@ ${N_THREADS} -s1 -b2 -e2 ${SAMPLE_ID}_format.tsv.gz
            AnnotateById ${SAMPLE_ID}_format.tsv.gz ${SAMPLE_ID}_header.txt "CHROM,POS,~ID,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT,SUPP_PBSV,SUPP_SNIFFLES,SUPP_PAV" ${INPUT_VCF_GZ} ${SAMPLE_ID}_out.vcf.gz z
            mv ${SAMPLE_ID}_out.vcf.gz ${SAMPLE_ID}_in.vcf.gz; bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_in.vcf.gz
            (bcftools view --no-header ${SAMPLE_ID}_in.vcf.gz | head -n 1 || echo "0") 1>&2

            rm -f ${INPUT_VCF_GZ}*
            mv ${SAMPLE_ID}_in.vcf.gz ${SAMPLE_ID}_kanpig.vcf.gz
            mv ${SAMPLE_ID}_in.vcf.gz.tbi ${SAMPLE_ID}_kanpig.vcf.gz.tbi

            # Removing temporary files
            rm -f ${SAMPLE_ID}_header.txt ${SAMPLE_ID}_format.tsv.gz
        }

        # Kanpig writes a genotype-level FORMAT/FT with integer values (0, 1, ...). FT is a reserved per-sample genotype-filter key, so htsjdk/GATK (used by the XGBoost scoring in the next workpackage) reject values like "0" as invalid filter names. Rename FT -> FTK (values preserved, no longer treated as a genotype filter) so every downstream tool parses cleanly. Same operation as utils/FixKanpigFT.wdl. No-op if FT is absent.
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

        cat << 'END' > truvari_bench.sh
#!/bin/bash

SAMPLE_ID=$1
INPUT_VCF_GZ=$2
TRAINING_RESOURCE_VCF_GZ=$3
INFINITY=$4
CHUNK_ID=$5
INCLUDE_BED=$6
${TIME_COMMAND} truvari bench -b ${TRAINING_RESOURCE_VCF_GZ} -c ${INPUT_VCF_GZ} --includebed ${INCLUDE_BED} --sizemin 1 --sizemax ${INFINITY} --sizefilt 1 --pctsize 0.9 --pctseq 0.9 --pick single -o ${SAMPLE_ID}_truvari_${CHUNK_ID}/
END
        chmod +x truvari_bench.sh

        # Extracts every record that has a stringent `truvari bench` match with some records in the resource. Remark: we use `--pick single` to force every resource record to be matched with at most one sample record, which is hopefully the most similar to it. This is because we assume that using a contaminated training set in XGBoost downstream is worse than using a slightly smaller training set. With `--pick multi` e.g. two records in the sample VCF might be matched to the same record in the resource VCF (probably not good) and vice versa (good). Remark: multiple instances of `truvari bench` are run in parallel using `not_gaps.bed`. Remark: in few anecdotal tests, `--pick multi` seems a bit faster than `--pick single` (4m vs 5m with 6 hyperthreading cores). Remark: both the inputs and the output of the function are indexed `.vcf.gz`, since they are needed by `truvari bench`.
        function GetTrainingRecords() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2

            # Running in parallel
            ${TIME_COMMAND} xargs --arg-file=training_not_gaps_beds.wsv --max-lines=1 --max-procs=${N_THREADS} ./truvari_bench.sh ${SAMPLE_ID} ${INPUT_VCF_GZ} ~{training_resource_vcf} ${INFINITY}

            # Concatenating outputs
            rm -f ${SAMPLE_ID}_outputs.txt
            while read -u 4 ROW; do
                ID=$(echo ${ROW} | cut -d ' ' -f 1)
                echo ${SAMPLE_ID}_truvari_${ID}/tp-comp.vcf.gz >> ${SAMPLE_ID}_outputs.txt
            done 4< training_not_gaps_beds.wsv
            ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --naive --file-list ${SAMPLE_ID}_outputs.txt --output-type z --output ${SAMPLE_ID}_training.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${SAMPLE_ID}_training.vcf.gz

            # Removing temporary files
            rm -rf ${SAMPLE_ID}_script.sh ${SAMPLE_ID}_outputs.txt ./${SAMPLE_ID}_truvari_*/
        }

        INFINITY="1000000000"
        truvari --help 1>&2
        ~{docker_dir}/kanpig --version 1>&2

        GetReferenceGaps ~{ref_agp} not_gaps.bed
        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)
            SEX=$(echo ${LINE} | cut -d , -f 2)

            # Skipping the sample if it has already been processed
            TEST=$( gcloud storage ls ${GCLOUD_STORAGE_BILLING_FLAGS} ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi

            # Merging
            LocalizeSample ${SAMPLE_ID} 1 ${LINE}
            if [ ${HAS_PAV} = "true" ]; then
                CanonizeVcf ${SAMPLE_ID}_pav.vcf.gz ${SAMPLE_ID}_pav.vcf.gz.tbi ${SAMPLE_ID} pav ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
            fi
            CanonizeVcf ${SAMPLE_ID}_pbsv.vcf.gz ${SAMPLE_ID}_pbsv.vcf.gz.tbi ${SAMPLE_ID} pbsv ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
            CanonizeVcf ${SAMPLE_ID}_sniffles.vcf.gz ${SAMPLE_ID}_sniffles.vcf.gz.tbi ${SAMPLE_ID} sniffles ~{min_sv_length} ~{max_sv_length} ~{standard_chromosomes_bed} not_gaps.bed
            IntrasampleMerge_sv ${SAMPLE_ID}
            IntrasampleMerge_ultralong ${SAMPLE_ID}
            IntrasampleMerge_bnd ${SAMPLE_ID}

            # Genotyping and marking training records
            LocalizeSample ${SAMPLE_ID} 2 ${LINE}
            CopySuppToInfo ${SAMPLE_ID} ${SAMPLE_ID}_sv.vcf.gz z ${SAMPLE_ID}_sv_supp.vcf.gz
            Kanpig ${SAMPLE_ID} ${SEX} ${SAMPLE_ID}_sv_supp.vcf.gz ${SAMPLE_ID}_aligned.bam
            CopyKanpigFieldsToInfo ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz
            FixKanpigFT ${SAMPLE_ID}
            GetTrainingRecords ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz

            # Copying SUPP fields to INFO in the BND and ultralong VCFs as well, just for uniformity. The original SUPP in FORMAT remains there, and since these VCFs won't be re-genotyped with kanpig, it will be correctly preserved by cohort-level truvari collapse.
            CopySuppToInfo ${SAMPLE_ID} ${SAMPLE_ID}_bnd.vcf.gz b ${SAMPLE_ID}_bnd_supp.bcf
            mv ${SAMPLE_ID}_bnd_supp.bcf ${SAMPLE_ID}_bnd.bcf
            mv ${SAMPLE_ID}_bnd_supp.bcf.csi ${SAMPLE_ID}_bnd.bcf.csi
            CopySuppToInfo ${SAMPLE_ID} ${SAMPLE_ID}_ultralong.vcf.gz b ${SAMPLE_ID}_ultralong_supp.bcf
            mv ${SAMPLE_ID}_ultralong_supp.bcf ${SAMPLE_ID}_ultralong.bcf
            mv ${SAMPLE_ID}_ultralong_supp.bcf.csi ${SAMPLE_ID}_ultralong.bcf.csi

            # Uploading
            gcloud storage mv ${GCLOUD_STORAGE_BILLING_FLAGS} ${SAMPLE_ID}_kanpig.vcf.'gz*' ${SAMPLE_ID}_kanpig.bed.gz ${SAMPLE_ID}_kanpig.csv ${SAMPLE_ID}_training.vcf.'gz*' ${SAMPLE_ID}_ultralong.'bcf*' ${SAMPLE_ID}_bnd.'bcf*' ~{remote_outdir}/
            touch ${SAMPLE_ID}.done
            gcloud storage mv ${GCLOUD_STORAGE_BILLING_FLAGS} ${SAMPLE_ID}.done ~{remote_outdir}/
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done 3< chunk.csv

        # Batch-completion signal, consumed by orchestrator workflows to order a downstream step after this one. Ignored by standalone runs.
        echo "done" > wp1.signal
    >>>

    output {
        String done = read_string("wp1.signal")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 6,
        mem_gb: 8,
        disk_gb: 256,
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

task ScoreSampleCallsBatch {
    input {
        File sv_integration_chunk_tsv
        File split_for_bcftools_merge_csv
        String filter_string

        String remote_indir
        String remote_outdir

        File training_resource_bed

        Array[String] annotations
        File training_python_script
        File scoring_python_script
        File hyperparameters_json

        String upstream_signal = ""
        String docker

        RuntimeAttr? runtime_attr_override
    }

    String docker_dir = "/root"

    command <<<
        set -euo pipefail

        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]))} - 1 ))
        GSUTIL_DELAY_S="600"
        export GATK_LOCAL_JAR="/root/gatk.jar"

        function LocalizeSample() {
            local SAMPLE_ID=$1
            local REMOTE_DIR=$2

            gsutil cp ${REMOTE_DIR}/${SAMPLE_ID}_kanpig.vcf.'gz*' ${REMOTE_DIR}/${SAMPLE_ID}_training.vcf.'gz*' .
        }

        # Deletes all files and directories related to the sample.
        function DelocalizeSample() {
            local SAMPLE_ID=$1

            rm -rf ./${SAMPLE_ID}_*
        }

        # Remark: the procedure's input and output are indexed `.vcf.gz`.
        function JointVcfFiltering() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2
            local RESOURCE_VCF_GZ=$3

            gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" ExtractVariantAnnotations -V ${INPUT_VCF_GZ} -O ${SAMPLE_ID}_extract -A ~{sep=" -A " annotations} --resource:resource,training=true,calibration=true ${RESOURCE_VCF_GZ} --maximum-number-of-unlabeled-variants 10000000 --mode INDEL --mnp-type INDEL -L ~{training_resource_bed}
            ls -laht
            # Output: ${SAMPLE_ID}_extract.annot.hdf5 ${SAMPLE_ID}_extract.unlabeled.annot.hdf5 ${SAMPLE_ID}_extract.vcf.gz ${SAMPLE_ID}_extract.vcf.gz.tbi
            gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" TrainVariantAnnotationsModel --annotations-hdf5 ${SAMPLE_ID}_extract.annot.hdf5 --unlabeled-annotations-hdf5 ${SAMPLE_ID}_extract.unlabeled.annot.hdf5 --model-backend PYTHON_SCRIPT --python-script ~{training_python_script} --hyperparameters-json ~{hyperparameters_json} -O ${SAMPLE_ID}.train --mode INDEL --verbosity DEBUG
            ls -laht
            # Output: ${SAMPLE_ID}.train.*
            gatk --java-options "-Xmx${EFFECTIVE_RAM_GB}G" ScoreVariantAnnotations -V ${INPUT_VCF_GZ} -O ${SAMPLE_ID}_score -A ~{sep=" -A " annotations} --resource:resource,training=true,calibration=true ${RESOURCE_VCF_GZ} --resource:extracted,extracted=true ${SAMPLE_ID}_extract.vcf.gz --model-prefix ${SAMPLE_ID}.train --model-backend PYTHON_SCRIPT --python-script ~{scoring_python_script} --mode INDEL --mnp-type INDEL --ignore-all-filters --verbosity DEBUG
            ls -laht
            # Output: ${SAMPLE_ID}_score.vcf.gz ${SAMPLE_ID}_score.vcf.gz.tbi ${SAMPLE_ID}_score.annot.hdf5 ${SAMPLE_ID}_score.scores.hdf5

            # Removing temporary files
            rm -f ${SAMPLE_ID}_extract.annot.hdf5 ${SAMPLE_ID}_extract.unlabeled.annot.hdf5 ${SAMPLE_ID}_extract.vcf.gz* ${SAMPLE_ID}.train.* ${SAMPLE_ID}_score.annot.hdf5 ${SAMPLE_ID}_score.scores.hdf5
        }

        # Copies the following fields from INFO to FORMAT, so that they are preserved by the inter-sample merge downstream: SUPP_*, SCORE, CALIBRATION_SENSITIVITY Remark: the procedure outputs an indexed `.bcf`. @param 2 A VCF where all IDs are distinct. This is guaranteed by workpackages upstream.
        function CopyInfoToFormat() {
            local SAMPLE_ID=$1
            local INPUT_VCF_GZ=$2

            echo '##FORMAT=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> ${SAMPLE_ID}_header.txt
            echo '##FORMAT=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by sniffles">' >> ${SAMPLE_ID}_header.txt
            echo '##FORMAT=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by pav">' >> ${SAMPLE_ID}_header.txt
            echo '##FORMAT=<ID=SCORE,Number=1,Type=Float,Description="Score according to the XGBoost model">' >> ${SAMPLE_ID}_header.txt
            echo '##FORMAT=<ID=CALIBRATION_SENSITIVITY,Number=1,Type=Float,Description="Calibration sensitivity according to the model applied by ScoreVariantAnnotations">' >> ${SAMPLE_ID}_header.txt
            # REF,ALT are emitted (cols 4,5) and added to --columns so the `~ID` match actually engages. Without REF,ALT present bcftools ignores `~ID` and matches on CHROM,POS only, mis-assigning FORMAT values across records sharing a start coordinate. IDs are distinct (see @param note above).
            bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%SUPP_PBSV\t%SUPP_SNIFFLES\t%SUPP_PAV\t%SCORE\t%CALIBRATION_SENSITIVITY\n' ${INPUT_VCF_GZ} | bgzip -c > ${SAMPLE_ID}_format.tsv.gz
            tabix -f -s1 -b2 -e2 ${SAMPLE_ID}_format.tsv.gz
            bcftools annotate --threads ${N_THREADS} --header-lines ${SAMPLE_ID}_header.txt --annotations ${SAMPLE_ID}_format.tsv.gz --columns CHROM,POS,~ID,REF,ALT,FORMAT/SUPP_PBSV,FORMAT/SUPP_SNIFFLES,FORMAT/SUPP_PAV,FORMAT/SCORE,FORMAT/CALIBRATION_SENSITIVITY --output-type b ${INPUT_VCF_GZ} --output ${SAMPLE_ID}_scored.bcf
            bcftools index --threads ${N_THREADS} ${SAMPLE_ID}_scored.bcf
            (bcftools view --no-header ${SAMPLE_ID}_scored.bcf | head -n 1 || echo "0") 1>&2

            # Removing temporary files
            rm -f ${SAMPLE_ID}_header.txt ${SAMPLE_ID}_format.tsv.gz*
        }

        # Assumes that `CopyInfoToFormat()` has already been executed.
        function PrintDebugInformation() {
            local SAMPLE_ID=$1
            local INPUT_BCF=$2

            rm -rf ${SAMPLE_ID}_xgboost.csv
            local N_RECORDS_BEFORE_FILTERING=$(bcftools index --nrecords ${INPUT_BCF})
            for THRESHOLD in 0.7 0.8 0.9 0.95 ; do
                local N_RECORDS_AFTER_FILTERING=$( bcftools query --format '%ID' --include "FORMAT/CALIBRATION_SENSITIVITY<=${THRESHOLD}" ${INPUT_BCF} | wc -l )
                local PERCENT=$( echo "scale=2; 100 * ${N_RECORDS_AFTER_FILTERING} / ${N_RECORDS_BEFORE_FILTERING}" | bc )
                echo "${N_RECORDS_AFTER_FILTERING},${N_RECORDS_BEFORE_FILTERING},${PERCENT},Number of records with CALIBRATION_SENSITIVITY<=${THRESHOLD}" >> ${SAMPLE_ID}_xgboost.csv
            done
            if [ "~{filter_string}" != "none" ]; then
                local N_RECORDS_AFTER_FILTERING=$( bcftools query --format '%ID' --include "~{filter_string}" ${INPUT_BCF} | wc -l )
                local PERCENT=$( echo "scale=2; 100 * ${N_RECORDS_AFTER_FILTERING} / ${N_RECORDS_BEFORE_FILTERING}" | bc )
                echo "${N_RECORDS_AFTER_FILTERING},${N_RECORDS_BEFORE_FILTERING},${PERCENT},Number of records that pass the specified filter" >> ${SAMPLE_ID}_xgboost.csv
            fi
        }

        # Remark: the procedure's input and output are indexed `.bcf`.
        function FilterChunkUpload() {
            local SAMPLE_ID=$1
            local INPUT_BCF=$2

            i="0"
            local INTERVAL
            while read -u 4 INTERVAL; do
                echo ${INTERVAL} | tr ',' '\t' > ${SAMPLE_ID}.bed
                # Remark: we use `targets` rather than `regions` because the former considers just the POS coordinate for overlaps.
                bcftools view --threads ${N_THREADS} --targets-file ${SAMPLE_ID}.bed --output-type b ${INPUT_BCF} --output ${SAMPLE_ID}_chunk_${i}.bcf
                bcftools index --threads ${N_THREADS} ${SAMPLE_ID}_chunk_${i}.bcf
                gsutil mv ${SAMPLE_ID}_chunk_${i}.bcf ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf
                gsutil mv ${SAMPLE_ID}_chunk_${i}.bcf.csi ~{remote_outdir}/chunk_${i}/${SAMPLE_ID}.bcf.csi
                i=$(( ${i} + 1 ))
            done 4< ~{split_for_bcftools_merge_csv}
            touch ${SAMPLE_ID}.done
            gsutil mv ${SAMPLE_ID}.done ~{remote_outdir}/ && echo 0 || echo 1
        }

        cat ~{sv_integration_chunk_tsv} | tr '\t' ',' > chunk.csv
        N_OUTPUT_CHUNKS=$(wc -l < ~{split_for_bcftools_merge_csv})
        while read -u 3 LINE; do
            SAMPLE_ID=$(echo ${LINE} | cut -d , -f 1)

            # Skipping the sample if it has already been processed
            TEST=$( gsutil ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0" )
            if [ ${TEST} != "0" ]; then
                continue
            fi

            # Filtering
            LocalizeSample ${SAMPLE_ID} ~{remote_indir}
            JointVcfFiltering ${SAMPLE_ID} ${SAMPLE_ID}_kanpig.vcf.gz ${SAMPLE_ID}_training.vcf.gz
            CopyInfoToFormat ${SAMPLE_ID} ${SAMPLE_ID}_score.vcf.gz
            PrintDebugInformation ${SAMPLE_ID} ${SAMPLE_ID}_scored.bcf
            FilterChunkUpload ${SAMPLE_ID} ${SAMPLE_ID}_scored.bcf
            DelocalizeSample ${SAMPLE_ID}
            ls -laht
        done 3< chunk.csv

        # Batch-completion signal for orchestrator ordering. Ignored standalone.
        echo "done" > wp2.signal
    >>>

    output {
        String done = read_string("wp2.signal")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 3,
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
