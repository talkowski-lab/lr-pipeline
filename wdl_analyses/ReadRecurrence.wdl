version 1.0

# Recurrence (dot) plots of the reads or assembly contigs of two haplotype BAMs (aligned to hg38) across one region,
# against hg38 and against T2T-CHM13.
# 1. PrepareRegion (once): hg38 region +/- flank, and the region lifted to T2T through an hg38 -> CHM13 chain file
#    +/- flank, as two small FASTAs (one sequence named chrom:start-end each); the optional SD and gene tracks are cut
#    to the two segments (empty files when not given).
# 2. ListRegionReads (per haplotype): reads whose non-secondary alignments overlap the region by >= min_overlap bp
#    (spanning_only: also covering the whole region), at most max_reads per haplotype; one small BAM per read with all its
#    primary and supplementary alignments (SA tag followed).
# 3. Recurrence (one job per read and reference): the read is cut by its CIGAR to the part aligned to the region +/- flank,
#    split into window bp windows (step bp apart), and every window is matched against the reference segment with edit
#    distance <= max_ed_frac x window, on both strands. Reads matching mostly on the reverse strand (reverse-strand
#    primary or supplementary alignments) are reverse-complemented so that their diagonal runs forward. One PDF per read
#    and reference: reference on x (SD and gene tracks underneath when given), read on y.
# 4. BundleOutputs: all figures in one tar.gz (<prefix>_recurrence_figures/<haplotype>/<reference>/<read>.pdf), all tables
#    in another, and the per-haplotype read lists.

import "Structs.wdl"

workflow ReadRecurrence {
    input {
        File hap1_bam
        File hap1_bam_idx
        File hap2_bam
        File hap2_bam_idx
        String region
        String prefix

        Int flank = 5000
        Int window = 100
        Int step = 100
        Float max_ed_frac = 0.1
        Int min_overlap = 1000
        Boolean spanning_only = false
        Int max_reads = 100
        String? t2t_region

        File hg38_fa
        File hg38_fai
        File? hg38_gzi
        File t2t_fa
        File t2t_fai
        File? t2t_gzi
        File chain
        File? hg38_segdup_bed
        File? hg38_genes_gtf
        File? t2t_segdup_bed
        File? t2t_genes_gtf

        String recurrence_docker
        Int recurrence_threads = 4

        RuntimeAttr? runtime_attr_prepare_region
        RuntimeAttr? runtime_attr_list_region_reads
        RuntimeAttr? runtime_attr_recurrence
        RuntimeAttr? runtime_attr_bundle_outputs
    }

    parameter_meta {
        hap1_bam: "Haplotype 1 reads or contigs aligned to hg38 (BAM)"
        hap2_bam: "Haplotype 2 reads or contigs aligned to hg38 (BAM)"
        region: "hg38 region, chrom:start-end"
        prefix: "Prefix for output files"
        flank: "bp added on each side of the region, on the references and on the reads"
        window: "Window size (bp)"
        step: "bp between window starts (= window: windows tile the read)"
        max_ed_frac: "Maximum edit distance of a match, as a fraction of the window"
        min_overlap: "Minimum bp of a read's alignments within the region"
        spanning_only: "Keep only reads whose alignments cover the whole region"
        max_reads: "Maximum reads per haplotype (largest overlap first)"
        t2t_region: "T2T-CHM13 chrom:start-end to compare against instead of the lifted region (flank still added)"
        hg38_gzi: "Needed when hg38_fa is bgzipped"
        t2t_gzi: "Needed when t2t_fa is bgzipped"
        chain: "UCSC chain file hg38 -> T2T-CHM13 (e.g. grch38-chm13v2.chain or hg38ToHs1.over.chain.gz)"
        hg38_segdup_bed: "Optional SD BED for the hg38 plots"
        hg38_genes_gtf: "Optional gene GTF/GFF3 for the hg38 plots (protein-coding genes are drawn)"
        t2t_segdup_bed: "Optional SD BED for the T2T plots"
        t2t_genes_gtf: "Optional gene GTF/GFF3 for the T2T plots"
    }

    Array[Pair[String, Pair[File, File]]] haplotypes = [("hap1", (hap1_bam, hap1_bam_idx)), ("hap2", (hap2_bam, hap2_bam_idx))]

    call PrepareRegion {
        input:
            region = region,
            flank = flank,
            t2t_region = t2t_region,
            prefix = prefix,
            hg38_fa = hg38_fa,
            hg38_fai = hg38_fai,
            hg38_gzi = hg38_gzi,
            t2t_fa = t2t_fa,
            t2t_fai = t2t_fai,
            t2t_gzi = t2t_gzi,
            chain = chain,
            hg38_segdup_bed = hg38_segdup_bed,
            hg38_genes_gtf = hg38_genes_gtf,
            t2t_segdup_bed = t2t_segdup_bed,
            t2t_genes_gtf = t2t_genes_gtf,
            docker = recurrence_docker,
            runtime_attr_override = runtime_attr_prepare_region
    }

    scatter (hap in haplotypes) {
        call ListRegionReads {
            input:
                bam = hap.right.left,
                bam_idx = hap.right.right,
                region = region,
                min_overlap = min_overlap,
                spanning_only = spanning_only,
                max_reads = max_reads,
                prefix = hap.left,
                docker = recurrence_docker,
                runtime_attr_override = runtime_attr_list_region_reads
        }

        scatter (read_bam in ListRegionReads.read_bams) {
            call Recurrence as RecurrenceHg38 {
                input:
                    read_bam = read_bam,
                    region = region,
                    flank = flank,
                    window = window,
                    step = step,
                    max_ed_frac = max_ed_frac,
                    reference_name = "hg38",
                    ref_fa = PrepareRegion.hg38_segment_fa,
                    segdup_bed = PrepareRegion.hg38_segdup_region,
                    genes_gtf = PrepareRegion.hg38_genes_region,
                    threads = recurrence_threads,
                    prefix = hap.left,
                    docker = recurrence_docker,
                    runtime_attr_override = runtime_attr_recurrence
            }

            call Recurrence as RecurrenceT2T {
                input:
                    read_bam = read_bam,
                    region = region,
                    flank = flank,
                    window = window,
                    step = step,
                    max_ed_frac = max_ed_frac,
                    reference_name = "T2T",
                    ref_fa = PrepareRegion.t2t_segment_fa,
                    chain = chain,
                    t2t_segment = PrepareRegion.t2t_segment,
                    segdup_bed = PrepareRegion.t2t_segdup_region,
                    genes_gtf = PrepareRegion.t2t_genes_region,
                    threads = recurrence_threads,
                    prefix = hap.left,
                    docker = recurrence_docker,
                    runtime_attr_override = runtime_attr_recurrence
            }
        }
    }

    call BundleOutputs {
        input:
            pdfs = flatten([flatten(RecurrenceHg38.pdf), flatten(RecurrenceT2T.pdf)]),
            tables = flatten([flatten(flatten(RecurrenceHg38.tables)), flatten(flatten(RecurrenceT2T.tables))]),
            read_lists = ListRegionReads.reads_tsv,
            lift_tsv = PrepareRegion.lift_tsv,
            prefix = prefix,
            docker = recurrence_docker,
            runtime_attr_override = runtime_attr_bundle_outputs
    }

    output {
        File figures_tar = BundleOutputs.figures_tar
        File tables_tar = BundleOutputs.tables_tar
        File reads_tsv = BundleOutputs.reads_tsv
        String t2t_segment = PrepareRegion.t2t_segment
        File lift_tsv = PrepareRegion.lift_tsv
    }
}

task PrepareRegion {
    input {
        String region
        Int flank
        String? t2t_region
        String prefix
        File hg38_fa
        File hg38_fai
        File? hg38_gzi
        File t2t_fa
        File t2t_fai
        File? t2t_gzi
        File chain
        File? hg38_segdup_bed
        File? hg38_genes_gtf
        File? t2t_segdup_bed
        File? t2t_genes_gtf
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String hg38_ext = if sub(basename(hg38_fa), "\\.(gz|bgz)$", "") != basename(hg38_fa) then ".fa.gz" else ".fa"
    String t2t_ext = if sub(basename(t2t_fa), "\\.(gz|bgz)$", "") != basename(t2t_fa) then ".fa.gz" else ".fa"

    command <<<
        set -euo pipefail

        # FASTA, .fai and .gzi side by side under one name
        mkdir -p ref
        ln -s ~{hg38_fa} ref/hg38~{hg38_ext}
        ln -s ~{hg38_fai} ref/hg38~{hg38_ext}.fai
        ~{"ln -s " + hg38_gzi + " ref/hg38" + hg38_ext + ".gzi"}
        ln -s ~{t2t_fa} ref/t2t~{t2t_ext}
        ln -s ~{t2t_fai} ref/t2t~{t2t_ext}.fai
        ~{"ln -s " + t2t_gzi + " ref/t2t" + t2t_ext + ".gzi"}

        python3 /opt/read_recurrence/prepare_region.py \
            --region ~{region} \
            --flank ~{flank} \
            --hg38-ref ref/hg38~{hg38_ext} \
            --t2t-ref ref/t2t~{t2t_ext} \
            --chain ~{chain} \
            ~{"--t2t-region " + t2t_region} \
            --out ~{prefix}

        # SD and gene tracks cut to the two segments (empty files when a track is not given), so that the per-read jobs
        # get small, always-defined track files
        subset() {
            local in=$1 region=$2 kind=$3 out=$4
            : > "${out}"
            [[ -z "${in}" ]] && return 0
            local chrom=${region%%:*} rng=${region#*:}
            local s=${rng%-*} e=${rng#*-}
            if [[ "${in}" == *.gz ]]; then gzip -cd "${in}"; else cat "${in}"; fi \
                | awk -F'\t' -v c="${chrom}" -v s="${s}" -v e="${e}" -v k="${kind}" \
                    'k == "bed" && $1 == c && $3 >= s && $2 <= e; k == "gtf" && $1 == c && $3 == "gene" && $5 >= s && $4 <= e' \
                > "${out}"
            echo "${in}: $(wc -l < "${out}") features in ${region}" >&2
        }
        hg38_seg=$(head -n 1 ~{prefix}.hg38.fa | tr -d '>')
        t2t_seg=$(cat ~{prefix}.t2t_region.txt)
        subset "~{hg38_segdup_bed}" "${hg38_seg}" bed ~{prefix}.hg38.segdup.bed
        subset "~{hg38_genes_gtf}" "${hg38_seg}" gtf ~{prefix}.hg38.genes.gtf
        subset "~{t2t_segdup_bed}" "${t2t_seg}" bed ~{prefix}.t2t.segdup.bed
        subset "~{t2t_genes_gtf}" "${t2t_seg}" gtf ~{prefix}.t2t.genes.gtf
    >>>

    output {
        File hg38_segment_fa = "~{prefix}.hg38.fa"
        File t2t_segment_fa = "~{prefix}.t2t.fa"
        String t2t_segment = read_string("~{prefix}.t2t_region.txt")
        File lift_tsv = "~{prefix}.lift.tsv"
        File hg38_segdup_region = "~{prefix}.hg38.segdup.bed"
        File hg38_genes_region = "~{prefix}.hg38.genes.gtf"
        File t2t_segdup_region = "~{prefix}.t2t.segdup.bed"
        File t2t_genes_region = "~{prefix}.t2t.genes.gtf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size([hg38_fa, t2t_fa, chain], "GB") + size([hg38_segdup_bed, hg38_genes_gtf, t2t_segdup_bed, t2t_genes_gtf], "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task ListRegionReads {
    input {
        File bam
        File bam_idx
        String region
        Int min_overlap
        Boolean spanning_only
        Int max_reads
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p in reads
        ln -s ~{bam} in/input.bam
        ln -s ~{bam_idx} in/input.bam.bai

        python3 /opt/read_recurrence/list_region_reads.py \
            --bam in/input.bam \
            --region ~{region} \
            --min-overlap ~{min_overlap} \
            ~{if spanning_only then "--spanning-only" else ""} \
            --max-reads ~{max_reads} \
            --out ~{prefix} \
            --outdir reads
        rm -f reads/*.bai
    >>>

    output {
        File reads_tsv = "~{prefix}.reads.tsv"
        Array[File] read_bams = glob("reads/*.bam")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bam, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task Recurrence {
    input {
        File read_bam
        String region
        Int flank
        Int window
        Int step
        Float max_ed_frac
        String reference_name
        File ref_fa
        File? chain
        String? t2t_segment
        File segdup_bed
        File genes_gtf
        Int threads
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String file_id = sub(basename(read_bam, ".bam"), "^" + prefix + "\\.", "")
    String out = "~{prefix}__~{reference_name}__~{file_id}"

    command <<<
        set -euo pipefail

        mkdir -p in tables
        ln -s ~{read_bam} in/read.bam
        cp ~{ref_fa} in/ref.fa
        python3 -c "import pysam; pysam.index('in/read.bam'); pysam.faidx('in/ref.fa')"
        read_name=$(python3 -c "import pysam; print(next(pysam.AlignmentFile('in/read.bam').fetch(until_eof=True)).query_name)")

        if [[ "~{reference_name}" == "T2T" ]]; then
            python3 /opt/read_recurrence/read_t2t_recurrence.py \
                --bam in/read.bam \
                --read-name "$read_name" \
                --region ~{region} \
                --chain ~{chain} \
                --t2t-ref in/ref.fa \
                --t2t-region ~{t2t_segment} \
                --out tables/~{out} \
                --window ~{window} \
                --step ~{step} \
                --max-ed-frac ~{max_ed_frac} \
                --pad ~{flank} \
                --threads ~{threads} \
                --orient-to-ref
            label="T2T-CHM13"
        else
            python3 /opt/read_recurrence/read_ref_recurrence.py \
                --bam in/read.bam \
                --read-name "$read_name" \
                --region ~{region} \
                --ref in/ref.fa \
                --out tables/~{out} \
                --window ~{window} \
                --step ~{step} \
                --max-ed-frac ~{max_ed_frac} \
                --pad ~{flank} \
                --threads ~{threads} \
                --orient-to-ref
            label="hg38"
        fi
        gzip tables/~{out}.hits.tsv

        # tracks are empty files when not given
        tracks=()
        if [[ -s ~{segdup_bed} ]]; then tracks+=(--segdup ~{segdup_bed}); fi
        if [[ -s ~{genes_gtf} ]]; then tracks+=(--genes ~{genes_gtf}); fi
        Rscript /opt/read_recurrence/plot_read_ref_recurrence.R \
            --prefix tables/~{out} \
            --ref-label "$label" \
            --title "~{prefix} $read_name vs $label; ~{window} bp windows, step ~{step}" \
            ${tracks[@]+"${tracks[@]}"} \
            --out ~{out}.pdf
    >>>

    output {
        File pdf = "~{out}.pdf"
        Array[File] tables = glob("tables/~{out}.*")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: threads,
        mem_gb: 8,
        disk_gb: 2 * ceil(size([read_bam, ref_fa, segdup_bed, genes_gtf], "GB")) + ceil(size(chain, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task BundleOutputs {
    input {
        Array[File] pdfs
        Array[File] tables
        Array[File] read_lists
        File lift_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # files are named <haplotype>__<reference>__<read>.<ext>
        fig=~{prefix}_recurrence_figures
        tab=~{prefix}_recurrence_tables
        mkdir -p "$fig" "$tab"
        while read -r f; do
            IFS='/' read -r hap ref read <<< "$(basename "$f" | sed 's/__/\//; s/__/\//')"
            mkdir -p "$fig/$hap/$ref"
            cp "$f" "$fig/$hap/$ref/$read"
        done < ~{write_lines(pdfs)}
        while read -r f; do
            IFS='/' read -r hap ref read <<< "$(basename "$f" | sed 's/__/\//; s/__/\//')"
            mkdir -p "$tab/$hap/$ref"
            cp "$f" "$tab/$hap/$ref/$read"
        done < ~{write_lines(tables)}
        cp ~{lift_tsv} "$tab/"

        # read lists of both haplotypes in one table
        first=1
        while read -r f; do
            hap=$(basename "$f" .reads.tsv)
            if [[ $first == 1 ]]; then head -1 "$f" | sed 's/^/haplotype\t/'; first=0; fi
            tail -n +2 "$f" | sed "s/^/$hap\t/"
        done < ~{write_lines(read_lists)} > ~{prefix}.reads.tsv
        cp ~{prefix}.reads.tsv "$fig/"

        tar -czf ~{prefix}.recurrence_figures.tar.gz "$fig"
        tar -czf ~{prefix}.recurrence_tables.tar.gz "$tab"
    >>>

    output {
        File figures_tar = "~{prefix}.recurrence_figures.tar.gz"
        File tables_tar = "~{prefix}.recurrence_tables.tar.gz"
        File reads_tsv = "~{prefix}.reads.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 3 * ceil(size(pdfs, "GB") + size(tables, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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
