version 1.0

# Recurrence (dot) plots of the split-aligned reads or assembly contigs of two haplotype BAMs (aligned to hg38) in a
# region: reads with >= min_alignments primary/supplementary alignments (SA tag followed; secondary alignments ignored)
# on the region's chromosome, one of them overlapping the region, e.g. a contig aligned on both sides of a large
# deletion. Each read is plotted against all the reference segments it aligned to (each alignment's span +/- pad),
# side by side, so a deletion shows as one diagonal that jumps between panels.
# 1. ListSplitReads (per haplotype): the split reads, one bundle each (the read's alignments, its reference segments as a
#    small FASTA, SD and gene tracks cut to the segments) and a summary table. A haplotype without split reads gives no
#    bundle and no plot. With ref_segments, every read is compared against the same fixed reference segments (e.g. the
#    two SDs of a GD +/- 50 kb) instead of its own alignment spans; min_alignments = 1 plots every read in the region.
# 2. SplitRecurrence (one job per split read): read_split_recurrence.py (every window x bp of the read's aligned span
#    +/- pad matched against every segment, edit distance <= max_ed_frac x window, both strands; reads matching mostly on
#    the reverse strand reverse-complemented) and the plot: reference segments side by side on x (named by
#    segment_names), the read on y; under them the labelled regions (region_labels, e.g. the SDs), the reference block
#    layout (blocks_tsv, segment_structure.py sections), SDs and genes.
# 3. BundleOutputs: all figures in one tar.gz (<prefix>_split_recurrence_figures/<haplotype>/<read>.pdf, with the split
#    read summary), all tables in another.

import "Structs.wdl"

workflow SplitReadRecurrence {
    input {
        File hap1_bam
        File hap1_bam_idx
        File hap2_bam
        File hap2_bam_idx
        String region
        String prefix

        Int pad = 20000
        Int window = 100
        Int step = 100
        Float max_ed_frac = 0.1
        Int min_overlap = 1000
        Int min_alignments = 2
        Boolean all_chroms = false
        Int max_reads = 50
        String? ref_segments
        String? segment_names
        String? region_labels

        File hg38_fa
        File hg38_fai
        File? hg38_gzi
        File? segdup_bed
        File? genes_gtf
        File? blocks_tsv

        String recurrence_docker
        Int recurrence_threads = 4

        RuntimeAttr? runtime_attr_list_split_reads
        RuntimeAttr? runtime_attr_split_recurrence
        RuntimeAttr? runtime_attr_bundle_outputs
    }

    parameter_meta {
        hap1_bam: "Haplotype 1 reads or contigs aligned to hg38 (BAM)"
        hap2_bam: "Haplotype 2 reads or contigs aligned to hg38 (BAM)"
        region: "hg38 region, chrom:start-end; split reads with an alignment overlapping it are plotted"
        prefix: "Prefix for output files"
        pad: "bp added on each side of every alignment, on the reference segments and on the read"
        window: "Window size (bp)"
        step: "bp between window starts (= window: windows tile the read)"
        max_ed_frac: "Maximum edit distance of a match, as a fraction of the window"
        min_overlap: "Minimum bp of a read's alignments within the region"
        min_alignments: "Minimum primary/supplementary alignments for a read to count as split"
        all_chroms: "Count and plot alignments on any chromosome (default: only the region's chromosome)"
        max_reads: "Maximum split reads plotted per haplotype (most alignments first)"
        ref_segments: "Fixed reference segments for every read, comma-separated chrom:start-end (e.g. chr16:29274176-29688676,chr16:30138532-30393458 for the 16p11.2 BP4 and BP5 SDs +/- 50 kb); default: each read's alignment spans +/- pad"
        segment_names: "Names of the reference segments for the panel titles, comma-separated in the same order (e.g. BP4,BP5)"
        region_labels: "Named regions marked on the plots, NAME=chrom:start-end separated by ';' (e.g. BP4 SD=chr16:29324176-29638676;BP5 SD=chr16:30188532-30343458)"
        blocks_tsv: "Reference block layout drawn under the plots (segment_structure.py <prefix>.sections.tsv)"
        hg38_gzi: "Needed when hg38_fa is bgzipped"
        segdup_bed: "Optional SD BED for the SD track"
        genes_gtf: "Optional gene GTF/GFF3 for the gene track (protein-coding genes are drawn)"
    }

    Array[Pair[String, Pair[File, File]]] haplotypes = [("hap1", (hap1_bam, hap1_bam_idx)), ("hap2", (hap2_bam, hap2_bam_idx))]

    scatter (hap in haplotypes) {
        call ListSplitReads {
            input:
                bam = hap.right.left,
                bam_idx = hap.right.right,
                region = region,
                pad = pad,
                min_overlap = min_overlap,
                min_alignments = min_alignments,
                all_chroms = all_chroms,
                max_reads = max_reads,
                ref_segments = ref_segments,
                prefix = hap.left,
                hg38_fa = hg38_fa,
                hg38_fai = hg38_fai,
                hg38_gzi = hg38_gzi,
                segdup_bed = segdup_bed,
                genes_gtf = genes_gtf,
                docker = recurrence_docker,
                runtime_attr_override = runtime_attr_list_split_reads
        }

        Array[Pair[String, File]] hap_bundles = cross([hap.left], ListSplitReads.bundles)
    }

    # One flat list of (haplotype, read bundle) jobs, so that the per-read calls sit in a single (not nested) scatter
    Array[Pair[String, File]] read_jobs = flatten(hap_bundles)

    scatter (job in read_jobs) {
        call SplitRecurrence {
            input:
                bundle = job.right,
                region = region,
                pad = pad,
                window = window,
                step = step,
                max_ed_frac = max_ed_frac,
                all_chroms = all_chroms,
                segment_names = segment_names,
                region_labels = region_labels,
                blocks_tsv = blocks_tsv,
                threads = recurrence_threads,
                prefix = job.left,
                docker = recurrence_docker,
                runtime_attr_override = runtime_attr_split_recurrence
        }
    }

    call BundleOutputs {
        input:
            pdfs = SplitRecurrence.pdf,
            tables = flatten(SplitRecurrence.tables),
            summaries = ListSplitReads.summary,
            prefix = prefix,
            docker = recurrence_docker,
            runtime_attr_override = runtime_attr_bundle_outputs
    }

    output {
        File figures_tar = BundleOutputs.figures_tar
        File tables_tar = BundleOutputs.tables_tar
        File split_reads_tsv = BundleOutputs.split_reads_tsv
    }
}

task ListSplitReads {
    input {
        File bam
        File bam_idx
        String region
        Int pad
        Int min_overlap
        Int min_alignments
        Boolean all_chroms
        Int max_reads
        String? ref_segments
        String prefix
        File hg38_fa
        File hg38_fai
        File? hg38_gzi
        File? segdup_bed
        File? genes_gtf
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String ext = if sub(basename(hg38_fa), "\\.(gz|bgz)$", "") != basename(hg38_fa) then ".fa.gz" else ".fa"

    command <<<
        set -euo pipefail

        mkdir -p in ref bundles
        ln -s ~{bam} in/input.bam
        ln -s ~{bam_idx} in/input.bam.bai
        ln -s ~{hg38_fa} ref/hg38~{ext}
        ln -s ~{hg38_fai} ref/hg38~{ext}.fai
        ~{"ln -s " + hg38_gzi + " ref/hg38" + ext + ".gzi"}

        python3 /opt/read_recurrence/list_split_reads.py \
            --bam in/input.bam \
            --region ~{region} \
            --ref ref/hg38~{ext} \
            --out ~{prefix} \
            --outdir bundles \
            --pad ~{pad} \
            --min-overlap ~{min_overlap} \
            --min-alignments ~{min_alignments} \
            ~{if all_chroms then "--all-chroms" else ""} \
            --max-reads ~{max_reads} \
            ~{"--segments " + ref_segments} \
            ~{"--segdup " + segdup_bed} \
            ~{"--genes " + genes_gtf}
    >>>

    output {
        File summary = "~{prefix}.split_reads.tsv"
        Array[File] bundles = glob("bundles/*.tar.gz")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bam, "GB")) + ceil(size(hg38_fa, "GB") + size([segdup_bed, genes_gtf], "GB")) + 10,
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

task SplitRecurrence {
    input {
        File bundle
        String region
        Int pad
        Int window
        Int step
        Float max_ed_frac
        Boolean all_chroms
        String? segment_names
        String? region_labels
        File? blocks_tsv
        Int threads
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String bundle_id = basename(bundle, ".tar.gz")
    String file_id = sub(bundle_id, "^" + prefix + "\\.", "")
    String out = "~{prefix}__~{file_id}"

    command <<<
        set -euo pipefail

        tar -xzf ~{bundle}
        b=~{bundle_id}
        mkdir -p tables
        read_name=$(cat "${b}/read_name.txt")

        python3 /opt/read_recurrence/read_split_recurrence.py \
            --bam "${b}/read.bam" \
            --read-name "${read_name}" \
            --region ~{region} \
            --ref "${b}/segments.fa" \
            --segments "$(cat "${b}/segments.txt")" \
            ~{if all_chroms then "--all-chroms" else ""} \
            --out tables/~{out} \
            --window ~{window} \
            --step ~{step} \
            --max-ed-frac ~{max_ed_frac} \
            --pad ~{pad} \
            --threads ~{threads} \
            --orient-to-ref
        gzip tables/~{out}.hits.tsv

        # tracks are empty files when not given
        tracks=()
        if [[ -s "${b}/segdup.bed" ]]; then tracks+=(--segdup "${b}/segdup.bed"); fi
        if [[ -s "${b}/genes.gtf" ]]; then tracks+=(--genes "${b}/genes.gtf"); fi
        Rscript /opt/read_recurrence/plot_read_ref_recurrence.R \
            --prefix tables/~{out} \
            --ref-label hg38 \
            --title "~{prefix} ${read_name}; ~{window} bp windows, step ~{step}" \
            ~{"--segment-names " + segment_names} \
            ~{if defined(region_labels) then "--regions '" + select_first([region_labels]) + "'" else ""} \
            ~{"--sections " + blocks_tsv} \
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
        disk_gb: 4 * ceil(size(bundle, "GB")) + 10,
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
        Array[File] summaries
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # files are named <haplotype>__<read>.<ext>
        fig=~{prefix}_split_recurrence_figures
        tab=~{prefix}_split_recurrence_tables
        mkdir -p "${fig}" "${tab}"
        while read -r f; do
            [[ -z "${f}" ]] && continue
            hap=$(basename "${f}" | sed 's/__.*//')
            mkdir -p "${fig}/${hap}"
            cp "${f}" "${fig}/${hap}/$(basename "${f}" | sed 's/^[^_]*__//')"
        done < ~{write_lines(pdfs)}
        while read -r f; do
            [[ -z "${f}" ]] && continue
            hap=$(basename "${f}" | sed 's/__.*//')
            mkdir -p "${tab}/${hap}"
            cp "${f}" "${tab}/${hap}/$(basename "${f}" | sed 's/^[^_]*__//')"
        done < ~{write_lines(tables)}

        # split reads of both haplotypes in one table
        first=1
        while read -r f; do
            hap=$(basename "${f}" .split_reads.tsv)
            if [[ ${first} == 1 ]]; then head -n 1 "${f}" | sed 's/^/haplotype\t/'; first=0; fi
            tail -n +2 "${f}" | sed "s/^/${hap}\t/"
        done < ~{write_lines(summaries)} > ~{prefix}.split_reads.tsv
        cp ~{prefix}.split_reads.tsv "${fig}/"
        echo "split reads plotted: $(find "${fig}" -name '*.pdf' | wc -l)" >&2

        tar -czf ~{prefix}.split_recurrence_figures.tar.gz "${fig}"
        tar -czf ~{prefix}.split_recurrence_tables.tar.gz "${tab}"
    >>>

    output {
        File figures_tar = "~{prefix}.split_recurrence_figures.tar.gz"
        File tables_tar = "~{prefix}.split_recurrence_tables.tar.gz"
        File split_reads_tsv = "~{prefix}.split_reads.tsv"
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
