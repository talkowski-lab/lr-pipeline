# Derived from yuliamostovoy/lr_callset_integration, a fork of fabio-cunial/callset_integration_phase2:
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_WorkflowC_Regenotype.wdl
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_Workpackage7_Main_joint_genotype_families.wdl

version 1.0

import "../utils/Structs.wdl"

workflow RegenotypeFamilySVCallsets {
    meta {
        description: [
            "This tool regenotypes the cohort SV callset from MergeCohortSVCallsets jointly within one family with Kanpig (https://github.com/ACEnglish/kanpig). The family's members are first checked against the PED, then a family-specific candidate VCF is built from the cohort callset and each member is regenotyped against its reads, and each member's regenotyped calls are split into the chunks of the interval CSV.",
            "It is run once per family, as one Terra sample set per family, and every family's run writes into the same 'remote_outdir'. Nothing is returned as a workflow output; MergeRegenotypedSVCallsets merges what the family runs wrote into the regenotyped cohort callset."
        ]
    }

    parameter_meta {
        family_id: "Family to regenotype, which must match a family ID in column 1 of the PED."
        sample_ids: "Members of the family, aligned by index with every other per-sample array. Must match the family's members in the PED exactly."
        sample_sexes: "Sex of each member. 'M' selects the male ploidy BED for Kanpig; any other value selects the female one."
        aligned_bams: "GCS paths to each member's aligned reads."
        aligned_bais: "GCS paths to the indexes for 'aligned_bams'."
        ped: "Six-column pedigree defining each family's members."
        remote_indir: "The 'remote_outdir' of MergeCohortSVCallsets, whose '06_concat' subdirectory holds the cohort callset."
        remote_outdir: "GCS directory shared by every family's run, which the per-sample regenotyped chunks are written under, without a trailing slash."
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
        String family_id
        Array[String] sample_ids
        Array[String] sample_sexes
        Array[String] aligned_bams
        Array[String] aligned_bais
        File ped
        String remote_indir
        String remote_outdir

        String requester_pays_project = ""
        String kanpig_params_cohort = "--neighdist 500 --gpenalty 0.04 --hapsim 0.97"
        File split_for_bcftools_merge_csv

        File ref_fa
        File ref_fai
        File ploidy_bed_female
        File ploidy_bed_male
        File autosomes_bed

        String sv_integration_docker

        RuntimeAttr? runtime_attr_make_family_inputs
        RuntimeAttr? runtime_attr_joint_genotype_family
    }

    String cohort_indir = sub(remote_indir, "/+$", "") + "/06_concat"
    String regenotyped_dir = sub(remote_outdir, "/+$", "")

    call MakeFamilyInputs {
        input:
            family_id = family_id,
            sample_ids = sample_ids,
            ped = ped,
            docker = sv_integration_docker,
            runtime_attr_override = runtime_attr_make_family_inputs
    }

    call JointGenotypeFamily {
        input:
            family_ids = [family_id],
            ped = ped,
            sample_ids = sample_ids,
            sample_sexes = sample_sexes,
            aligned_bais = aligned_bais,
            aligned_bams = aligned_bams,
            split_for_bcftools_merge_csv = split_for_bcftools_merge_csv,
            remote_indir = cohort_indir,
            remote_outdir = regenotyped_dir,
            requester_pays_project = requester_pays_project,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male,
            autosomes_bed = autosomes_bed,
            kanpig_params_cohort = kanpig_params_cohort,
            upstream_signal = [MakeFamilyInputs.done],
            docker = sv_integration_docker,
            runtime_attr_override = runtime_attr_joint_genotype_family
    }

    output {
    }
}

task MakeFamilyInputs {
    input {
        String family_id
        Array[String] sample_ids
        File ped
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Members from the Terra sample_set.
        cat > set_members.raw <<'EOF_SAMPLE_IDS'
~{sep="\n" sample_ids}
EOF_SAMPLE_IDS
        grep -v '^[[:space:]]*$' set_members.raw | sort -u > set_members.txt
        if [ ! -s set_members.txt ]; then
            echo "ERROR: sample set ~{family_id} has no members."
            exit 1
        fi

        # Members of this family per the PED.
        awk -v fam="~{family_id}" 'BEGIN { FS="[ \t]+" } $1==fam && $2!="0" && $2!="." { print $2 }' ~{ped} | sort -u > ped_members.txt
        if [ ! -s ped_members.txt ]; then
            echo "ERROR: family_id '~{family_id}' not found in PED column 1 (or it lists no members). The Terra sample_set id must match the PED family id."
            exit 1
        fi

        # The set must EXACTLY match the PED family -- guards against a mis-built set.
        if ! diff -q set_members.txt ped_members.txt >/dev/null; then
            echo "ERROR: sample set '~{family_id}' does not match its PED family." 1>&2
            echo "  In the set but NOT in the PED family:" 1>&2; comm -23 set_members.txt ped_members.txt 1>&2
            echo "  In the PED family but NOT in the set:" 1>&2; comm -13 set_members.txt ped_members.txt 1>&2
            exit 1
        fi

        echo "validated" > validated.txt
    >>>

    output {
        String done = read_string("validated.txt")
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

task JointGenotypeFamily {
    input {
        Array[String] family_ids
        File ped
        Array[String] sample_ids
        Array[String] sample_sexes
        Array[String] aligned_bais
        Array[String] aligned_bams
        File split_for_bcftools_merge_csv

        String remote_indir
        String remote_outdir
        String requester_pays_project

        File ref_fa
        File ref_fai
        File ploidy_bed_female
        File ploidy_bed_male
        File autosomes_bed

        String kanpig_params_cohort

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
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.22/plugins"
        export RUST_BACKTRACE="full"
        GCLOUD_STORAGE_BILLING_FLAGS=""
        if [ -n "~{requester_pays_project}" ]; then
            GCLOUD_STORAGE_BILLING_FLAGS="--billing-project=~{requester_pays_project}"
        fi

        function LocalizeSample() {
            local SAMPLE_ID=$1
            local LINE=$2

            local ALIGNED_BAI=$(echo "${LINE}" | cut -f 3)
            local ALIGNED_BAM=$(echo "${LINE}" | cut -f 4)

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

        function BuildFamilyCandidateVcf() {
            local FAMILY_ID=$1

            awk -v family="${FAMILY_ID}" 'BEGIN { FS="[ \t]+" } $1==family && $2!="0" && $2!="." { print $2 }' ~{ped} | sort -u > ${FAMILY_ID}.samples.txt
            local N_FAMILY_SAMPLES=$(wc -l < ${FAMILY_ID}.samples.txt)
            if [ ${N_FAMILY_SAMPLES} -eq 0 ]; then
                echo "ERROR: family ${FAMILY_ID} has no samples in the PED."
                exit 1
            fi

            local SAMPLE_ID
            while read -u 4 SAMPLE_ID; do
                if ! awk -v sample="${SAMPLE_ID}" 'BEGIN { FS="\t"; found=0 } $1==sample { found=1 } END { exit(found ? 0 : 1) }' sample_metadata.tsv; then
                    echo "ERROR: sample ${SAMPLE_ID} from family ${FAMILY_ID} is missing from sample_ids."
                    exit 1
                fi
            done 4< ${FAMILY_ID}.samples.txt

            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples-file ${FAMILY_ID}.samples.txt --output-type z cohort.bcf --output ${FAMILY_ID}_all.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${FAMILY_ID}_all.vcf.gz
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --include 'COUNT(GT="alt")>0' --output-type z ${FAMILY_ID}_all.vcf.gz --output ${FAMILY_ID}_present.vcf.gz
            bcftools index --threads ${N_THREADS} -f -t ${FAMILY_ID}_present.vcf.gz
            rm -f ${FAMILY_ID}_all.vcf.gz*

            local N_RECORDS=$(bcftools index --nrecords ${FAMILY_ID}_present.vcf.gz)
            echo "${FAMILY_ID},${N_FAMILY_SAMPLES},${N_RECORDS},Number of family samples and family-present records" > ${FAMILY_ID}_family.csv
        }

        function Kanpig() {
            local FAMILY_ID=$1
            local SAMPLE_ID=$2
            local SEX=$3

            local PLOIDY_BED
            if [ ${SEX} == "M" ]; then
                PLOIDY_BED=$(echo ~{ploidy_bed_male})
            else
                PLOIDY_BED=$(echo ~{ploidy_bed_female})
            fi

            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --samples ${SAMPLE_ID} --output-type z ${FAMILY_ID}_present.vcf.gz --output ${SAMPLE_ID}_personalized.vcf.gz
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

        cat > sample_ids.txt <<'EOF_SAMPLE_IDS'
~{sep="\n" sample_ids}
EOF_SAMPLE_IDS
        cat > sample_sexes.txt <<'EOF_SAMPLE_SEXES'
~{sep="\n" sample_sexes}
EOF_SAMPLE_SEXES
        cat > aligned_bais.txt <<'EOF_ALIGNED_BAIS'
~{sep="\n" aligned_bais}
EOF_ALIGNED_BAIS
        cat > aligned_bams.txt <<'EOF_ALIGNED_BAMS'
~{sep="\n" aligned_bams}
EOF_ALIGNED_BAMS
        cat > family_ids.txt <<'EOF_FAMILY_IDS'
~{sep="\n" family_ids}
EOF_FAMILY_IDS
        grep -v '^[[:space:]]*$' family_ids.txt | sort -u > family_ids.unique.txt
        mv family_ids.unique.txt family_ids.txt

        N_SAMPLE_IDS=$(wc -l < sample_ids.txt)
        N_SAMPLE_SEXES=$(wc -l < sample_sexes.txt)
        N_ALIGNED_BAIS=$(wc -l < aligned_bais.txt)
        N_ALIGNED_BAMS=$(wc -l < aligned_bams.txt)
        if [ ${N_SAMPLE_IDS} -ne ${N_SAMPLE_SEXES} ] || [ ${N_SAMPLE_IDS} -ne ${N_ALIGNED_BAIS} ] || [ ${N_SAMPLE_IDS} -ne ${N_ALIGNED_BAMS} ]; then
            echo "ERROR: sample_ids, sample_sexes, aligned_bais, and aligned_bams must have the same length."
            echo "sample_ids=${N_SAMPLE_IDS}, sample_sexes=${N_SAMPLE_SEXES}, aligned_bais=${N_ALIGNED_BAIS}, aligned_bams=${N_ALIGNED_BAMS}"
            exit 1
        fi
        paste sample_ids.txt sample_sexes.txt aligned_bais.txt aligned_bams.txt > sample_metadata.tsv

        # Localizing the WP6 cohort VCF.
        ${TIME_COMMAND} gcloud storage cp ~{remote_indir}/truvari_collapsed.'bcf*' .
        mv truvari_collapsed.bcf cohort.bcf
        mv truvari_collapsed.bcf.csi cohort.bcf.csi

        while read -u 3 FAMILY_ID; do
            if [ -z "${FAMILY_ID}" ] || [[ "${FAMILY_ID}" == \#* ]]; then
                continue
            fi

            BuildFamilyCandidateVcf ${FAMILY_ID}

            while read -u 4 SAMPLE_ID; do
                TEST=$(gcloud storage ls ~{remote_outdir}/${SAMPLE_ID}.done || echo "0")
                if [ ${TEST} != "0" ]; then
                    continue
                fi

                LINE=$(awk -v sample="${SAMPLE_ID}" 'BEGIN { FS="\t" } $1==sample { print; exit }' sample_metadata.tsv)
                SEX=$(echo "${LINE}" | cut -f 2)
                LocalizeSample ${SAMPLE_ID} "${LINE}"
                Kanpig ${FAMILY_ID} ${SAMPLE_ID} ${SEX}
                FixKanpigFT ${SAMPLE_ID}
                ChunkAndUpload ${SAMPLE_ID}

                touch ${SAMPLE_ID}.done
                gcloud storage mv ${SAMPLE_ID}.done ~{remote_outdir}/
                DelocalizeSample ${SAMPLE_ID}
                ls -laht 1>&2
            done 4< ${FAMILY_ID}.samples.txt

            gcloud storage cp ${FAMILY_ID}_family.csv ~{remote_outdir}/
            rm -f ${FAMILY_ID}.samples.txt ${FAMILY_ID}_present.vcf.gz* ${FAMILY_ID}_family.csv
        done 3< family_ids.txt

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
