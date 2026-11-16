version 1.0

# Read-depth profiles of CNV sites, one PDF per input site BED, scattered so the groups run in parallel.
# Each site BED (plot_cnv_calls_depth.R format: #chrom start end SVID carriers SVLEN SVTYPE ...; carriers = comma-separated
# samples, all coloured) is plotted with plot_cnv_calls_depth.R on depth read on the VM: the windows (site +/- max(flank_mult x
# size, min_flank)) are sliced from the cohort RD matrix with pysam over gs:// (OAuth token from the VM metadata server), so the
# matrix is never localized. Normalised depth = RD / sample median coverage x 2 (median_coverage: 2-row sample / median file).

import "Structs.wdl"

workflow PlotCnvDepthProfiles {
    input {
        Array[File] site_beds
        String rd_matrix
        String prefix

        Float flank_mult = 1.0
        Int min_flank = 20000
        Int min_bin = 100
        String svtypes = "DEL,DUP"
        String extra_args = ""

        File median_coverage
        File segdup_bed
        File sr_bed
        File plot_script
        File font

        String docker

        RuntimeAttr? runtime_attr_plot_depth
    }

    parameter_meta {
        site_beds: "Site BEDs, one PDF each (#chrom start end SVID carriers SVLEN SVTYPE [sr_match_ids sr_match_detail ...])"
        rd_matrix: "Cohort RD matrix ('#Chr Start End <samples>', bgzipped + tabix-indexed), gs:// path, read in place"
        median_coverage: "Sample median coverage, 2 rows (sample IDs / medians)"
        sr_bed: "Short-read SV intervals (#chrom start end ID SVTYPE FILTER) for the srGS track"
        plot_script: "plot_cnv_calls_depth.R"
        font: "Arial TTF installed for cairo_pdf"
        extra_args: "Extra arguments appended to the plot_script call, e.g. '--mode fn'"
    }

    scatter (bed in site_beds) {
        call PlotDepth {
            input:
                sites = bed,
                rd_matrix = rd_matrix,
                prefix = prefix + "." + basename(bed, ".sites.bed.gz"),
                flank_mult = flank_mult,
                min_flank = min_flank,
                min_bin = min_bin,
                svtypes = svtypes,
                extra_args = extra_args,
                median_coverage = median_coverage,
                segdup_bed = segdup_bed,
                sr_bed = sr_bed,
                plot_script = plot_script,
                font = font,
                docker = docker,
                runtime_attr_override = runtime_attr_plot_depth
        }
    }

    output {
        Array[File] depth_profile_pdfs = PlotDepth.pdf
        Array[File] depth_slice_logs = PlotDepth.log
    }
}

task PlotDepth {
    input {
        File sites
        String rd_matrix
        String prefix
        Float flank_mult
        Int min_flank
        Int min_bin
        String svtypes
        String extra_args
        File median_coverage
        File segdup_bed
        File sr_bed
        File plot_script
        File font
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p ~/.fonts && cp ~{font} ~/.fonts/ && fc-cache -f > /dev/null

        export GCS_OAUTH_TOKEN=$(curl -sS -H "Metadata-Flavor: Google" \
            http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
            | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')

        # Plotting windows: site +/- max(flank_mult x size, min_flank), merged
        zcat ~{sites} | awk -F'\t' -v OFS='\t' -v m=~{flank_mult} -v mf=~{min_flank} \
            '!/^#/ {f = m * ($3 - $2); if (f < mf) f = mf; s = $2 - f; if (s < 0) s = 0; print $1, int(s), int($3 + f)}' \
            | sort -k1,1 -k2,2n \
            | awk -F'\t' -v OFS='\t' 'NR == 1 {c = $1; s = $2; e = $3; next}
                $1 == c && $2 <= e {if ($3 > e) e = $3; next}
                {print c, s, e; c = $1; s = $2; e = $3}
                END {if (NR) print c, s, e}' > windows.bed
        echo "sites: $(zcat ~{sites} | grep -vc '^#'); windows: $(wc -l < windows.bed)" | tee ~{prefix}.log

        # RD slice, one query per window through pysam (the image's tabix binary has no gs:// support), retried
        python3 - ~{rd_matrix} windows.bed rd.tsv <<'PY'
import sys
import time
import pysam
path, windows, out = sys.argv[1:4]
tbx = pysam.TabixFile(path)
rows = set()
with open(out, "w") as fo:
    fo.write("\n".join(tbx.header) + "\n")
    for line in open(windows):
        c, s, e = line.split()[:3]
        for attempt in range(5):
            try:
                got = list(tbx.fetch(c, int(s), int(e)))
                break
            except (OSError, ValueError):
                time.sleep(5 * (attempt + 1))
                tbx = pysam.TabixFile(path)
        else:
            sys.exit("failed to read %s:%s-%s" % (c, s, e))
        for r in got:
            f = r.split("\t", 2)
            k = (f[0], int(f[1]))
            if k not in rows:
                rows.add(k)
                fo.write(r + "\n")
print("RD bins:", len(rows))
PY
        (head -n 1 rd.tsv; tail -n +2 rd.tsv | sort -k1,1 -k2,2n) | bgzip -c > rd.txt.gz
        tabix -f -s1 -b2 -e3 rd.txt.gz
        echo "RD bins: $(($(wc -l < rd.tsv) - 1))" | tee -a ~{prefix}.log
        rm -f rd.tsv

        Rscript ~{plot_script} \
            --bed ~{sites} \
            --rd rd.txt.gz \
            --segdup ~{segdup_bed} \
            --medians ~{median_coverage} \
            --svtype ~{svtypes} \
            --min-size 0 \
            --flank-mult ~{flank_mult} \
            --min-flank ~{min_flank} \
            --bin-size 0 \
            --min-bin ~{min_bin} \
            --sr-bed ~{sr_bed} \
            ~{extra_args} \
            --out ~{prefix}.depth_profiles.pdf 2>&1 | tee -a ~{prefix}.log
    >>>

    output {
        File pdf = "~{prefix}.depth_profiles.pdf"
        File log = "~{prefix}.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: 50,
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
