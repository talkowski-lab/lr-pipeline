version 1.0

# Breakpoint profiles of genomic-disorder (GD) loci: for each locus, the segmental duplication (SD) at its left and right
# breakpoints (+/- half the SD length) is profiled with read depth, BAF and LD, and plotted as one figure
# (left SD on the left, right SD on the right; DEL carriers red to orange, DUP carriers blue to green, arrows at the BAF
# breakpoint estimates of DEL carriers), under a full-width copy-number overview of the whole locus (left window start to
# right window end) with the SD tracks.
#
# loci_tsv: output of gd_sd_windows.py (header + one row per locus; windows 0-based half-open; default: the 39-locus
#   table built in the image from references/gd_calls.txt). Columns used:
#   1 locus, 2 chrom, 3 gd_start, 4 gd_end, 7 highlight_pattern, 8 left_sd_name, 12-13 left window,
#   14 right_sd_name, 18-19 right window.
# gvcfs: one VCF per contig (local or gs://, bgzipped + tabix index next to it); the VCF of a contig is the path
#   whose file name starts with "<contig>." or "<contig>_" (e.g. .../chr16.g.vcf.bgz). SelectLoci reads the loci table first,
#   lists the contigs that carry GD loci (gd_contigs) and matches each to its VCF; only those contigs' loci are run (loci
#   whose contig has no VCF are listed in skipped_loci), and only the windows are streamed. Run the workflow once with all
#   contig VCFs (e.g. on the contig set), not once per contig.
# Outputs: all profile PDFs in one tar.gz (<prefix>_GD_breakpoint_profiles/<chrom>/<locus>.breakpoint_profile.pdf) and
#   all tables in another (<prefix>_GD_breakpoint_tables/<chrom>/<locus>/...).
# rd_matrix: cohort read-depth matrix ('#Chr Start End <samples>', bgzipped + .tbi); only the windows are streamed.
# Samples in the output tables and figures are renamed with id_map (real ID -> new ID; unmapped samples are dropped).

import "Structs.wdl"

workflow GDBreakpointProfiles {
    input {
        File? loci_tsv
        Array[String] gvcfs
        String rd_matrix
        String prefix

        Int cn_bin = 1000
        Int page_size = 6
        Boolean order_by_edges = true
        File rd_medians
        File? samples
        File id_map
        String? requester_pays_project

        File? segdup_bed
        File? bp_sd_bed

        String gd_docker

        RuntimeAttr? runtime_attr_select_loci
        RuntimeAttr? runtime_attr_prepare_window
        RuntimeAttr? runtime_attr_prepare_locus_depth
        RuntimeAttr? runtime_attr_plot_locus
        RuntimeAttr? runtime_attr_bundle_outputs
    }

    parameter_meta {
        loci_tsv: "GD loci with left/right breakpoint-SD windows (gd_sd_windows.py output, with header; default: the 39-locus table in the image)"
        gvcfs: "One small-variant VCF per contig (gs:// allowed), each with a tabix index next to it"
        rd_matrix: "Cohort read-depth matrix ('#Chr Start End <samples>', bgzipped + tabix-indexed; gs:// allowed)"
        prefix: "Prefix for output files"
        cn_bin: "Copy-number bin size (bp) of all depth panels; 0 = automatic per window"
        page_size: "Carriers per page and type (N DEL + N DUP per page); 0 = all carriers on one page"
        order_by_edges: "Order carriers by breakpoint (DEL: BAF breakpoint estimates; DUP: copy-number edges) before paging"
        rd_medians: "2 rows: sample IDs / genome-wide median depth per sample (same scale as rd_matrix); real or new IDs"
        samples: "Optional real sample IDs (one per line) for BAF and genotypes; default: column 1 of id_map"
        id_map: "2 columns with a header line: real sample ID, new sample ID; output samples are renamed, unmapped ones dropped"
        requester_pays_project: "Google project billed for requester-pays buckets (gvcfs, rd_matrix)"
        segdup_bed: "Merged segmental-duplication BED for the SD track (default: hg38 copy in the image)"
        bp_sd_bed: "GD NAHR breakpoint-SD table for the SD breakpoint track (default: copy in the image)"
    }

    call SelectLoci {
        input:
            tsv = loci_tsv,
            gvcfs = gvcfs,
            prefix = prefix,
            docker = gd_docker,
            runtime_attr_override = runtime_attr_select_loci
    }

    scatter (row in read_tsv(SelectLoci.rows)) {

        call PrepareWindow as PrepareLeft {
            input:
                gvcf = row[19],
                rd_matrix = rd_matrix,
                chrom = row[1],
                win_start = row[11],
                win_end = row[12],
                prefix = "~{prefix}.~{row[0]}.left",
                samples = samples,
                id_map = id_map,
                requester_pays_project = requester_pays_project,
                docker = gd_docker,
                runtime_attr_override = runtime_attr_prepare_window
        }

        call PrepareWindow as PrepareRight {
            input:
                gvcf = row[19],
                rd_matrix = rd_matrix,
                chrom = row[1],
                win_start = row[17],
                win_end = row[18],
                prefix = "~{prefix}.~{row[0]}.right",
                samples = samples,
                id_map = id_map,
                requester_pays_project = requester_pays_project,
                docker = gd_docker,
                runtime_attr_override = runtime_attr_prepare_window
        }

        call PrepareLocusDepth {
            input:
                rd_matrix = rd_matrix,
                chrom = row[1],
                win_start = row[11],
                win_end = row[18],
                prefix = "~{prefix}.~{row[0]}.locus",
                rd_medians = rd_medians,
                id_map = id_map,
                requester_pays_project = requester_pays_project,
                docker = gd_docker,
                runtime_attr_override = runtime_attr_prepare_locus_depth
        }

        call PlotLocus {
            input:
                left_rd = PrepareLeft.rd,
                left_baf = PrepareLeft.baf,
                left_gt = PrepareLeft.gt,
                right_rd = PrepareRight.rd,
                right_baf = PrepareRight.baf,
                right_gt = PrepareRight.gt,
                locus_rd = PrepareLocusDepth.rd,
                highlight_pattern = row[6],
                gd_start = row[2],
                gd_end = row[3],
                left_sd_name = row[7],
                left_window = "~{row[11]},~{row[12]}",
                right_sd_name = row[13],
                right_window = "~{row[17]},~{row[18]}",
                prefix = "~{prefix}.~{row[0]}",
                cn_bin = cn_bin,
                page_size = page_size,
                order_by_edges = order_by_edges,
                rd_medians = PrepareLocusDepth.medians,
                segdup_bed = segdup_bed,
                bp_sd_bed = bp_sd_bed,
                docker = gd_docker,
                runtime_attr_override = runtime_attr_plot_locus
        }
    }

    call BundleOutputs {
        input:
            chroms = SelectLoci.chroms,
            loci = SelectLoci.loci,
            profile_pdfs = PlotLocus.profile_pdf,
            left_ranked = PlotLocus.left_ranked,
            right_ranked = PlotLocus.right_ranked,
            locus_tables = PlotLocus.tables_tar,
            left_window_rd = PrepareLeft.rd,
            left_window_baf = PrepareLeft.baf,
            left_window_gt = PrepareLeft.gt,
            right_window_rd = PrepareRight.rd,
            right_window_baf = PrepareRight.baf,
            right_window_gt = PrepareRight.gt,
            locus_rd = PrepareLocusDepth.rd,
            prefix = prefix,
            docker = gd_docker,
            runtime_attr_override = runtime_attr_bundle_outputs
    }

    output {
        File breakpoint_profiles_tar = BundleOutputs.profiles_tar
        File breakpoint_tables_tar = BundleOutputs.tables_tar
        File gd_contigs = SelectLoci.contigs
        File skipped_loci = SelectLoci.skipped
    }
}

task SelectLoci {
    input {
        File? tsv
        Array[String] gvcfs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Rows (header dropped) of the loci on contigs that have a VCF in gvcfs, with that VCF appended as column 20;
        # loci on other contigs are listed in <prefix>.loci.skipped.tsv
        loci=~{default="/opt/gd_scripts/references/gd_loci.windows.tsv" tsv}
        # The loci table must be gd_sd_windows.py output: header starting with "locus", 19 tab-separated columns
        if ! head -n 1 "${loci}" | grep -q "^locus"$'\t'"chrom"; then
            echo "loci_tsv must be the gd_sd_windows.py table (header 'locus chrom gd_start ...'); got: $(head -n 1 "${loci}")" >&2
            exit 1
        fi
        bad=$(awk -F'\t' 'NR > 1 && NF != 19' "${loci}" | wc -l)
        if [[ ${bad} -gt 0 ]]; then
            echo "loci_tsv: ${bad} rows do not have 19 tab-separated columns" >&2
            exit 1
        fi
        tail -n +2 "${loci}" | tr -d '\r' > loci.tsv

        # Contigs that carry GD loci, each matched to its VCF (file name starting with "<chrom>." or "<chrom>_")
        printf 'contig\tn_loci\tvcf\n' > ~{prefix}.gd_contigs.tsv
        : > ~{prefix}.loci.rows.tsv
        : > ~{prefix}.loci.skipped.tsv
        cut -f2 loci.tsv | awk '!seen[$0]++' > contigs.txt
        while read -r chrom; do
            n=$(awk -F'\t' -v c="${chrom}" '$2 == c' loci.tsv | wc -l | tr -d ' ')
            vcf=$(grep -E "(^|/)${chrom}[._][^/]*$" ~{write_lines(gvcfs)} | head -n 1 || true)
            if [[ -n "${vcf}" ]]; then
                printf '%s\t%s\t%s\n' "${chrom}" "${n}" "${vcf}" >> ~{prefix}.gd_contigs.tsv
                awk -F'\t' -v OFS='\t' -v c="${chrom}" -v v="${vcf}" '$2 == c {print $0, v}' loci.tsv >> ~{prefix}.loci.rows.tsv
            else
                printf '%s\t%s\tno VCF\n' "${chrom}" "${n}" >> ~{prefix}.gd_contigs.tsv
                awk -F'\t' -v OFS='\t' -v c="${chrom}" '$2 == c {print $1, $2, "no VCF for " c " in gvcfs"}' loci.tsv >> ~{prefix}.loci.skipped.tsv
            fi
        done < contigs.txt
        cut -f1 ~{prefix}.loci.rows.tsv > ~{prefix}.loci.names.txt
        cut -f2 ~{prefix}.loci.rows.tsv > ~{prefix}.loci.chroms.txt
        echo "contigs with GD loci: $(wc -l < contigs.txt) ($(paste -sd, contigs.txt)); loci kept: $(wc -l < ~{prefix}.loci.rows.tsv); skipped (no VCF for the contig): $(wc -l < ~{prefix}.loci.skipped.tsv)" >&2
    >>>

    output {
        File rows = "~{prefix}.loci.rows.tsv"
        File contigs = "~{prefix}.gd_contigs.tsv"
        File skipped = "~{prefix}.loci.skipped.tsv"
        Array[String] loci = read_lines("~{prefix}.loci.names.txt")
        Array[String] chroms = read_lines("~{prefix}.loci.chroms.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: 10,
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

task PrepareWindow {
    input {
        String gvcf
        String rd_matrix
        String chrom
        String win_start
        String win_end
        String prefix
        File? samples
        File id_map
        String? requester_pays_project
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # htslib (via pysam) reads gs:// with an OAuth token; the image has no gcloud, so use the VM metadata server
        export GCS_OAUTH_TOKEN=$(curl -sS -H "Metadata-Flavor: Google" \
            http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
            | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')

        bash /opt/gd_scripts/prepare_region_data.sh \
            --chr ~{chrom} \
            --start $(( ~{win_start} + 1 )) \
            --end ~{win_end} \
            --name ~{prefix} \
            --rd ~{rd_matrix} \
            --vcf ~{gvcf} \
            ~{"--samples " + samples} \
            --id-map ~{id_map} \
            ~{"--requester-pays-project " + requester_pays_project}
    >>>

    output {
        File rd = "~{prefix}.renamed.RD.txt.gz"
        File baf = "~{prefix}.renamed.baf.tsv.gz"
        File gt = "~{prefix}.renamed.GT.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 20,
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

task PrepareLocusDepth {
    input {
        String rd_matrix
        String chrom
        String win_start
        String win_end
        String prefix
        File rd_medians
        File id_map
        String? requester_pays_project
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(curl -sS -H "Metadata-Flavor: Google" \
            http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
            | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')

        bash /opt/gd_scripts/prepare_region_data.sh \
            --chr ~{chrom} \
            --start $(( ~{win_start} + 1 )) \
            --end ~{win_end} \
            --name ~{prefix} \
            --rd ~{rd_matrix} \
            --id-map ~{id_map} \
            --rd-only \
            ~{"--requester-pays-project " + requester_pays_project}

        # Per-sample medians under the new IDs (2 rows: IDs / medians); a sample with several new IDs gets each of them.
        # Medians already keyed by the new IDs are kept as they are.
        python3 - ~{rd_medians} ~{id_map} ~{prefix}.renamed.RD_median <<'PY'
import sys
med_path, map_path, out_path = sys.argv[1:4]
rows = [ln.rstrip("\n").split("\t") for ln in open(med_path) if ln.strip()]
med = dict(zip(rows[0], rows[1]))
pairs = [ln.split()[:2] for ln in open(map_path) if len(ln.split()) >= 2][1:]
new_ids = {n for _, n in pairs}
if sum(i in new_ids for i in med) >= sum(r in med for r, _ in pairs):
    out = med
else:
    out = {n: med[r] for r, n in pairs if r in med}
with open(out_path, "w") as fh:
    fh.write("\t".join(out) + "\n" + "\t".join(out.values()) + "\n")
print(f"medians: {len(out)} samples under the new IDs", file=sys.stderr)
PY
    >>>

    output {
        File rd = "~{prefix}.renamed.RD.txt.gz"
        File medians = "~{prefix}.renamed.RD_median"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 20,
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

task PlotLocus {
    input {
        File left_rd
        File left_baf
        File left_gt
        File right_rd
        File right_baf
        File right_gt
        File locus_rd
        String highlight_pattern
        String gd_start
        String gd_end
        String left_sd_name
        String left_window
        String right_sd_name
        String right_window
        String prefix
        Int cn_bin
        Int page_size
        Boolean order_by_edges
        File rd_medians
        File? segdup_bed
        File? bp_sd_bed
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # analyze_gd_locus.sh reads <prefix>.renamed.{RD.txt,baf.tsv,GT.tsv}.gz
        mkdir -p in out
        ln -s ~{left_rd} in/left.renamed.RD.txt.gz
        ln -s ~{left_baf} in/left.renamed.baf.tsv.gz
        ln -s ~{left_gt} in/left.renamed.GT.tsv.gz
        ln -s ~{right_rd} in/right.renamed.RD.txt.gz
        ln -s ~{right_baf} in/right.renamed.baf.tsv.gz
        ln -s ~{right_gt} in/right.renamed.GT.tsv.gz

        bash /opt/gd_scripts/analyze_gd_locus.sh \
            --name ~{prefix} \
            --highlight-pattern '~{highlight_pattern}' \
            --left-prefix in/left \
            --left-title '~{left_sd_name}' \
            --left-window ~{left_window} \
            --right-prefix in/right \
            --right-title '~{right_sd_name}' \
            --right-window ~{right_window} \
            --gd-start ~{gd_start} \
            --gd-end ~{gd_end} \
            --overview-rd ~{locus_rd} \
            --cn-bin ~{cn_bin} \
            --page-size ~{page_size} \
            ~{if order_by_edges then "--order-by-edges" else ""} \
            --medians ~{rd_medians} \
            --segdup ~{default="/opt/gd_scripts/references/hg38.SegDup.sorted.merged.bed.gz" segdup_bed} \
            --bp-sd ~{default="/opt/gd_scripts/references/GD_NAHR.breakpoint_SDs.bed.gz" bp_sd_bed} \
            --outdir out

        tar -czf ~{prefix}.tables.tar.gz -C out --exclude='*.pdf' .
        mv out/~{prefix}.breakpoint_profile.pdf .
        mv out/~{prefix}.left_loh.ranked.tsv out/~{prefix}.right_loh.ranked.tsv .
    >>>

    output {
        File profile_pdf = "~{prefix}.breakpoint_profile.pdf"
        File left_ranked = "~{prefix}.left_loh.ranked.tsv"
        File right_ranked = "~{prefix}.right_loh.ranked.tsv"
        File tables_tar = "~{prefix}.tables.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: 4 * ceil(size([left_rd, left_baf, left_gt, right_rd, right_baf, right_gt, locus_rd], "GB")) + 10,
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
        Array[String] chroms
        Array[String] loci
        Array[File] profile_pdfs
        Array[File] left_ranked
        Array[File] right_ranked
        Array[File] locus_tables
        Array[File] left_window_rd
        Array[File] left_window_baf
        Array[File] left_window_gt
        Array[File] right_window_rd
        Array[File] right_window_baf
        Array[File] right_window_gt
        Array[File] locus_rd
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # One line per locus (same order as the scatter): chrom, locus, then its files
        fig=~{prefix}_GD_breakpoint_profiles
        tab=~{prefix}_GD_breakpoint_tables
        mkdir -p "${fig}" "${tab}"
        paste ~{write_lines(chroms)} ~{write_lines(loci)} ~{write_lines(profile_pdfs)} > profiles.tsv
        while IFS=$'\t' read -r chrom locus pdf; do
            mkdir -p "${fig}/${chrom}"
            cp "${pdf}" "${fig}/${chrom}/"
        done < profiles.tsv
        paste ~{write_lines(chroms)} ~{write_lines(loci)} ~{write_lines(left_ranked)} ~{write_lines(right_ranked)} \
            ~{write_lines(locus_tables)} ~{write_lines(left_window_rd)} ~{write_lines(left_window_baf)} \
            ~{write_lines(left_window_gt)} ~{write_lines(right_window_rd)} ~{write_lines(right_window_baf)} \
            ~{write_lines(right_window_gt)} ~{write_lines(locus_rd)} > tables.tsv
        while IFS=$'\t' read -r chrom locus files; do
            mkdir -p "${tab}/${chrom}/${locus}"
            for f in ${files}; do
                cp "${f}" "${tab}/${chrom}/${locus}/"
            done
        done < <(awk -F'\t' -v OFS='\t' '{f = $3; for (i = 4; i <= NF; i++) f = f " " $i; print $1, $2, f}' tables.tsv)

        tar -czf ~{prefix}.GD_breakpoint_profiles.tar.gz "${fig}"
        tar -czf ~{prefix}.GD_breakpoint_tables.tar.gz "${tab}"
        echo "profiles: $(find "${fig}" -name '*.pdf' | wc -l) PDFs in $(ls "${fig}" | wc -l) chromosome folders" >&2
    >>>

    output {
        File profiles_tar = "~{prefix}.GD_breakpoint_profiles.tar.gz"
        File tables_tar = "~{prefix}.GD_breakpoint_tables.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 3 * ceil(size(profile_pdfs, "GB") + size(locus_tables, "GB") + size(left_window_baf, "GB") + size(right_window_baf, "GB") + size(left_window_rd, "GB") + size(right_window_rd, "GB") + size(left_window_gt, "GB") + size(right_window_gt, "GB") + size(locus_rd, "GB")) + 10,
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
