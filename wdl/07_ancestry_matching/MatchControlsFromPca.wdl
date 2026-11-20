version 1.0

import "../utils/Structs.wdl"
import "MatchControlsByAncestry.wdl"

# Select ancestry-matched controls for each case from a finished PCA (e.g. the study_ancestry.tsv output of
# PcaAncestryTwoStudies: sample_id, cohort, assigned_ancestry, knn_label, knn_fraction, PC1..PC10).
# 1. SelectControls: for each sample whose cohort is case_name, the n_controls samples whose cohort is control_name that
#    are nearest in PC space (Euclidean distance over the first match_n_pcs PCs). With eigenvalues (the smartpca .eval
#    output of the same run), each PC is weighted by sqrt(eigenvalue), so distances follow the variance each PC
#    explains; without it, all PCs weigh the same. Matching is with replacement: a control can be matched to several
#    cases. Outputs: matched_controls.tsv (case sample, comma-separated matched controls, nearest first),
#    matched_pairs.tsv (one row per case-control pair: ancestries, rank, PC distance) and control_weights.tsv (times
#    each control was matched).
# 2. PlotMatches: PC1 vs PC2 and PC3 vs PC4 with each case as a filled dot, each matched control as an open circle and a
#    line from each case to each of its controls, coloured by assigned ancestry in the gnomAD palette (AFR, AMR, EAS,
#    EUR, SAS/CSA, ASJ, remaining/unassigned grey; MID and OCE, absent from that palette, get brown and black). PDF in
#    Arial when arial_ttf is given, otherwise the fontconfig Arial substitute.

workflow MatchControlsFromPca {
    input {
        File study_ancestry
        String prefix

        File? eigenvalues
        String case_name = "CDP"
        String control_name = "AoU"
        Int n_controls = 5
        Int match_n_pcs = 10
        File? arial_ttf

        String python_docker = "python:3.11-slim"
        String r_docker = "us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers/gd-breakpoint-profiles:2026-11-09.3"

        RuntimeAttr? runtime_attr_match_controls
        RuntimeAttr? runtime_attr_plot_matches
    }

    call MatchControlsByAncestry.MatchControls as SelectControls {
        input:
            study_ancestry = study_ancestry,
            eigenvalues = eigenvalues,
            case_name = case_name,
            control_name = control_name,
            n_controls = n_controls,
            match_n_pcs = match_n_pcs,
            prefix = prefix,
            docker = python_docker,
            runtime_attr_override = runtime_attr_match_controls
    }

    call PlotMatches {
        input:
            study_ancestry = study_ancestry,
            matched_pairs = SelectControls.matched_pairs,
            case_name = case_name,
            control_name = control_name,
            arial_ttf = arial_ttf,
            prefix = prefix,
            docker = r_docker,
            runtime_attr_override = runtime_attr_plot_matches
    }

    output {
        File matched_controls = SelectControls.matched_controls
        File matched_pairs = SelectControls.matched_pairs
        File control_weights = SelectControls.control_weights
        File matched_pcs_pdf = PlotMatches.matched_pcs_pdf
    }
}

task PlotMatches {
    input {
        File study_ancestry
        File matched_pairs
        String case_name
        String control_name
        File? arial_ttf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [[ -n "~{arial_ttf}" ]]; then
            mkdir -p ~/.fonts
            cp ~{arial_ttf} ~/.fonts/Arial.ttf
            fc-cache -f > /dev/null
        fi

        Rscript - ~{study_ancestry} ~{matched_pairs} "~{case_name}" "~{control_name}" ~{prefix} <<'RS'
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
a <- commandArgs(trailingOnly = TRUE)
anc <- fread(a[1], colClasses = c(sample_id = "character"))
pairs <- fread(a[2], colClasses = "character")
case_name <- a[3]; control_name <- a[4]; prefix <- a[5]
setnames(pairs, c(paste0(case_name, "_sample"), paste0(control_name, "_sample")), c("case", "control"))

fams <- unlist(strsplit(tryCatch(system("fc-list : family", intern = TRUE), error = function(e) "Arial"), ","))
font_family <- if ("Arial" %in% trimws(fams)) "Arial" else "sans"

# gnomAD genetic ancestry colours; CSA (HGDP central/south Asian) takes the SAS colour, unassigned the remaining grey;
# MID and OCE are not in that palette
pal <- c(AFR = "#941494", AMR = "#ED1E24", EAS = "#108C44", EUR = "#6AA5CD", SAS = "#FF9912", CSA = "#FF9912",
         ASJ = "#FF7F50", MID = "#A65628", OCE = "#000000", RMI = "#ABB9B9", unassigned = "#ABB9B9")
labs_anc <- unique(anc$assigned_ancestry)
extra <- setdiff(labs_anc, names(pal))
if (length(extra) > 0) pal <- c(pal, setNames(rep("#666666", length(extra)), extra))
anc[, ancestry := factor(assigned_ancestry, levels = intersect(names(pal), labs_anc))]

cases <- anc[sample_id %in% pairs$case]
controls <- anc[sample_id %in% pairs$control]
cases[, role := paste0(case_name, " case")]
controls[, role := paste0("matched ", control_name, " control")]
pts <- rbind(cases, controls)
seg <- merge(merge(pairs[, .(case, control)], anc[, .(case = sample_id, PC1, PC2, PC3, PC4, ancestry)], by = "case"),
             anc[, .(control = sample_id, cPC1 = PC1, cPC2 = PC2, cPC3 = PC3, cPC4 = PC4)], by = "control")

panel <- function(dt, x, y, label) data.table(dt, px = dt[[x]], py = dt[[y]], panel = label)
p12 <- "PC1 vs PC2"; p34 <- "PC3 vs PC4"
pts_l <- rbind(panel(pts, "PC1", "PC2", p12), panel(pts, "PC3", "PC4", p34))
seg_l <- rbind(data.table(seg[, .(ancestry)], x = seg$PC1, y = seg$PC2, xend = seg$cPC1, yend = seg$cPC2, panel = p12),
               data.table(seg[, .(ancestry)], x = seg$PC3, y = seg$PC4, xend = seg$cPC3, yend = seg$cPC4, panel = p34))
shapes <- setNames(c(16, 1), c(paste0(case_name, " case"), paste0("matched ", control_name, " control")))

p <- ggplot() +
  geom_segment(data = seg_l, aes(x = x, y = y, xend = xend, yend = yend, colour = ancestry), linewidth = 0.3, alpha = 0.35) +
  geom_point(data = pts_l[role == names(shapes)[2]], aes(px, py, colour = ancestry, shape = role), size = 4.5, stroke = 1.1) +
  geom_point(data = pts_l[role == names(shapes)[1]], aes(px, py, colour = ancestry, shape = role), size = 3.2, alpha = 0.9) +
  facet_wrap(~panel, scales = "free") +
  scale_colour_manual(values = pal, name = "Assigned ancestry", drop = TRUE) +
  scale_shape_manual(values = shapes, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 16, alpha = 1)),
         shape = guide_legend(override.aes = list(size = 5, colour = "grey30"))) +
  labs(x = "PC1 (left) / PC3 (right)", y = "PC2 (left) / PC4 (right)",
       title = sprintf("%s: %s cases and their matched %s controls", prefix, case_name, control_name),
       subtitle = sprintf("%s cases, %s distinct controls, %s case-control links",
                          format(nrow(cases), big.mark = ","), format(nrow(controls), big.mark = ","),
                          format(nrow(seg), big.mark = ","))) +
  theme_bw(base_size = 22, base_family = font_family) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"))
out <- paste0(prefix, ".matched_controls_pcs.pdf")
if (capabilities("cairo")) cairo_pdf(out, width = 26, height = 12, family = font_family) else pdf(out, width = 26, height = 12)
print(p)
invisible(dev.off())
cat(sprintf("%d cases, %d matched controls, %d links; ancestries: %s\n", nrow(cases), nrow(controls), nrow(seg),
            paste(levels(anc$ancestry), collapse = ", ")))
RS
    >>>

    output {
        File matched_pcs_pdf = "~{prefix}.matched_controls_pcs.pdf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
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
