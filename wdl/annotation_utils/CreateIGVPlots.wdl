# Derived from wdl/RunVisualizePlotsLongRead.wdl of yuliamostovoy/variant-interpretation
# https://github.com/yuliamostovoy/variant-interpretation/blob/main/wdl/RunVisualizePlotsLongRead.wdl
# Derived from wdl/IgvBamPlots.wdl of yuliamostovoy/variant-interpretation
# https://github.com/yuliamostovoy/variant-interpretation/blob/main/wdl/IgvBamPlots.wdl
# Derived from wdl/LongReadDepthPlot.wdl of yuliamostovoy/variant-interpretation
# https://github.com/yuliamostovoy/variant-interpretation/blob/main/wdl/LongReadDepthPlot.wdl

version 1.0

import "../utils/Structs.wdl"

workflow CreateIGVPlots {
    meta {
        description: [
            "This utility takes IGV snapshots of the long reads at a curated list of variants, so calls can be reviewed by eye, and packs them into one tarball with one PNG per variant per family.",
            "A family is plotted when any of its members carries at least one listed variant, and each of its plots shows every family member's reads, labelled and ordered by carrier and affected status from `ped`. Only the reads around the listed variants are streamed from each BAM. A variant whose window, its length plus 1.5 times its length, is shorter than `igv_max_window` gets one snapshot spanning it; a longer one gets a snapshot of each breakpoint, placed side by side. A DEL or DUP also gets a read-depth track below the snapshot, each family member's mosdepth depth normalized to its median over the flanks, with carriers in red. A header on each plot gives the variant ID, length, type and location.",
            "This is a port of the long-read visualization workflow of variant-interpretation (https://github.com/yuliamostovoy/variant-interpretation), by Yulia Mostovoy and the Talkowski lab."
        ]
    }

    parameter_meta {
        variant_list: "Tab-separated list of variants to plot, with columns `chrom`, 0-based `start`, `end`, `ID`, `svtype` and comma-separated carrier `samples`; a header line is optional. Generate it with `bcftools query -f '%CHROM\\t%POS0\\t%END\\t%ID\\t%INFO/SVTYPE\\t[%SAMPLE,]\\n'` over the carriers' records."
        ped: "Pedigree file for the cohort, without a header; families and carrier labels come from it."
        sample_ids: "Sample IDs, matching `ped` and the `samples` column of `variant_list`."
        bams: "Aligned BAMs, in `sample_ids` order; only the plotted regions are streamed."
        bais: "Indexes for `bams`."
        igv_genome: "IGV-hosted genome ID that IGV loads for its reference sequence and gene track."
        igv_max_window: "Longest window, in bp, shown in a single snapshot; longer variants get one snapshot per breakpoint."
        breakpoint_window: "Span, in bp, shown on each side of a breakpoint in a breakpoint snapshot."
        depth_flank: "Minimum flank, in bp, shown on each side of a DEL or DUP in its depth track."
        depth_flank_fraction: "Flank shown on each side of a DEL or DUP in its depth track, as a share of its length, when above `depth_flank`."
        depth_window: "Window size, in bp, of the depth track."
        annotation_beds: "BED files drawn as feature tracks under the reads and shaded on the depth tracks, such as segmental duplications."
        annotation_names: "Track names for `annotation_beds`, in the same order."
        igv_plots_tar: "Tarball of `igv_plots/<family>/<variant_id>_<family>.png`, one plot per family per variant it carries."
    }

    input {
        File variant_list
        File ped
        Array[String] sample_ids
        Array[File] bams
        Array[File] bais
        String prefix

        String igv_genome = "hg38"
        Int igv_max_window = 150000
        Int breakpoint_window = 500
        Int depth_flank = 5000
        Float depth_flank_fraction = 0.1
        Int depth_window = 250
        Array[File] annotation_beds = []
        Array[String] annotation_names = []

        String utils_docker
        String igv_docker

        RuntimeAttr? runtime_attr_prepare_igv_variants
        RuntimeAttr? runtime_attr_plot_family_igv
        RuntimeAttr? runtime_attr_merge_igv_plots
    }

    call PrepareIGVVariants {
        input:
            variant_list = variant_list,
            ped = ped,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_prepare_igv_variants
    }

    scatter (family in PrepareIGVVariants.families) {
        call PlotFamilyIGV {
            input:
                variants_bed = PrepareIGVVariants.variants_bed,
                ped = ped,
                family = family,
                sample_ids = sample_ids,
                bams = bams,
                bais = bais,
                igv_genome = igv_genome,
                igv_max_window = igv_max_window,
                breakpoint_window = breakpoint_window,
                depth_flank = depth_flank,
                depth_flank_fraction = depth_flank_fraction,
                depth_window = depth_window,
                annotation_beds = annotation_beds,
                annotation_names = annotation_names,
                docker = igv_docker,
                runtime_attr_override = runtime_attr_plot_family_igv
        }
    }

    call MergeIGVPlots {
        input:
            plot_tars = PlotFamilyIGV.plots_tar,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_igv_plots
    }

    output {
        File igv_plots_tar = MergeIGVPlots.merged_tar
    }
}

task PrepareIGVVariants {
    input {
        File variant_list
        File ped
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/scripts/qc/igv/reformat_variants_for_visualization.py \
            --input ~{variant_list} \
            --output ~{prefix}.igv_variants.bed

        # Families with a member carrying any listed variant
        tail -n +2 ~{prefix}.igv_variants.bed | cut -f6 | tr ',' '\n' | sort -u > carriers.txt
        awk 'NR == FNR { carrier[$1]; next } $2 in carrier { print $1 }' carriers.txt ~{ped} | sort -u > families.txt
    >>>

    output {
        File variants_bed = "~{prefix}.igv_variants.bed"
        Array[String] families = read_lines("families.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
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

task PlotFamilyIGV {
    input {
        File variants_bed
        File ped
        String family
        Array[String] sample_ids
        Array[File] bams
        Array[File] bais
        String igv_genome
        Int igv_max_window
        Int breakpoint_window
        Int depth_flank
        Float depth_flank_fraction
        Int depth_window
        Array[File] annotation_beds
        Array[String] annotation_names
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bams: { localization_optional: true }
        bais: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        awk -v family="~{family}" '$1 == family { print $2 }' ~{ped} > members.txt
        paste ~{write_lines(sample_ids)} ~{write_lines(bams)} ~{write_lines(bais)} \
            | awk 'NR == FNR { member[$1]; next } $1 in member' members.txt - > family_bams.tsv

        # Variants carried by a family member
        awk -F'\t' 'NR == FNR { member[$1]; next } FNR > 1 {
            n = split($6, carriers, ",")
            for (i = 1; i <= n; i++) if (carriers[i] in member) { print; next }
        }' members.txt ~{variants_bed} > family_variants.bed

        # Stream the reads covering every snapshot and depth track
        awk -v OFS='\t' -v flank=~{depth_flank} -v fraction=~{depth_flank_fraction} -v breakpoint=~{breakpoint_window} '{
            length_bp = $3 - $2
            pad = flank
            if (fraction * length_bp > pad) pad = int(fraction * length_bp)
            if (0.25 * length_bp > pad) pad = int(0.25 * length_bp)
            if (breakpoint > pad) pad = breakpoint
            start = $2 - pad
            print $1, (start < 0 ? 0 : start), $3 + pad
        }' family_variants.bed | sort -k1,1 -k2,2n | bedtools merge > regions.bed
        bedtools makewindows -b regions.bed -w ~{depth_window} > windows.bed

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        while read -r SAMPLE BAM BAI; do
            samtools view -h -b -M -X -L regions.bed -o "${SAMPLE}.bam" "${BAM}" "${BAI}"
            samtools index "${SAMPLE}.bam"
            mosdepth --no-per-base --by windows.bed "${SAMPLE}" "${SAMPLE}.bam"
        done < family_bams.tsv
        cut -f1 family_bams.tsv | sed 's/$/.bam/' > bams.txt

        mkdir igv_plots depth_plots
        python3 /opt/scripts/qc/igv/makeigvpesr.py \
            -v family_variants.bed \
            -fam_id ~{family} \
            -samples "$(cut -f1 family_bams.tsv | paste -sd, -)" \
            -crams bams.txt \
            -p ~{ped} \
            -o igv_plots \
            -b ~{breakpoint_window} \
            -i igv_batch.txt \
            -bam igv_setup.sh \
            -m ~{igv_max_window} \
            --genome ~{igv_genome} \
            --long_read \
            --status_labels \
            --annotation_beds ~{sep=' ' annotation_beds} \
            --annotation_names ~{sep=' ' annotation_names}

        # Size the IGV heap to the VM, leaving room for the virtual display
        MEM_KB=$(awk '/MemTotal/ { print $2 }' /proc/meminfo)
        sed -i -E "s/-Xmx[0-9]+[mMgG]/-Xmx$(( MEM_KB * 80 / 100 / 1024 ))m/" /opt/igv/igv.sh
        # Exit on running out of memory, which IGV otherwise catches and hangs on
        export _JAVA_OPTIONS="-XX:+ExitOnOutOfMemoryError"
        # IGV fetches its genome over HTTP, so retry transient failures
        for ATTEMPT in 1 2 3; do
            if xvfb-run --auto-servernum --server-args="-screen 0 1920x1080x24" /opt/igv/igv.sh -g ~{igv_genome} -b igv_batch.txt; then
                break
            fi
            if [ "${ATTEMPT}" -eq 3 ]; then
                echo "IGV failed after ${ATTEMPT} attempts" >&2
                exit 1
            fi
            sleep $(( ATTEMPT * 30 ))
        done

        python3 /opt/scripts/qc/igv/plot_longread_depth.py \
            --bed family_variants.bed \
            --ped ~{ped} \
            --family ~{family} \
            --flank ~{depth_flank} \
            --flank-frac ~{depth_flank_fraction} \
            --depth-dir . \
            --outdir depth_plots \
            --annotation-beds ~{sep=' ' annotation_beds} \
            --annotation-names ~{sep=' ' annotation_names}

        mkdir -p igv_plots_merged
        python3 /opt/scripts/qc/igv/concat_igv_depth.py \
            --igv-dir igv_plots \
            --depth-dir depth_plots \
            --varfile family_variants.bed \
            --outdir igv_plots_merged/~{family}
        tar -czf ~{family}.igv_plots.tar.gz -C igv_plots_merged ~{family}
    >>>

    output {
        File plots_tar = "~{family}.igv_plots.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 12,
        disk_gb: 50,
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

task MergeIGVPlots {
    input {
        Array[File] plot_tars
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir igv_plots
        for TAR in ~{sep=' ' plot_tars}; do
            tar -xzf "$TAR" -C igv_plots
        done

        tar -czf ~{prefix}.igv_plots.tar.gz igv_plots
    >>>

    output {
        File merged_tar = "~{prefix}.igv_plots.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(plot_tars, "GB")) * 3 + 10,
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
