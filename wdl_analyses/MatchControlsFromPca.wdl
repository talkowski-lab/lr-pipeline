version 1.0

import "Structs.wdl"
import "MatchControlsByAncestry.wdl"

# Select ancestry-matched controls for each case from a finished PCA (e.g. the study_ancestry.tsv output of
# PcaAncestryTwoStudies: sample_id, cohort, assigned_ancestry, knn_label, knn_fraction, PC1..PC10).
# For each sample whose cohort is case_name, the n_controls samples whose cohort is control_name that are nearest in PC
# space (Euclidean distance over the first match_n_pcs PCs). With eigenvalues (the smartpca .eval output of the same
# run), each PC is weighted by sqrt(eigenvalue), so distances follow the variance each PC explains; without it, all PCs
# weigh the same. Matching is with replacement: a control can be matched to several cases.
# Outputs: matched_controls.tsv (case sample, comma-separated matched controls, nearest first), matched_pairs.tsv (one
# row per case-control pair: ancestries, rank, PC distance) and control_weights.tsv (times each control was matched).

workflow MatchControlsFromPca {
    input {
        File study_ancestry
        String prefix

        File? eigenvalues
        String case_name = "CDP"
        String control_name = "AoU"
        Int n_controls = 5
        Int match_n_pcs = 10

        String python_docker = "python:3.11-slim"

        RuntimeAttr? runtime_attr_match_controls
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

    output {
        File matched_controls = SelectControls.matched_controls
        File matched_pairs = SelectControls.matched_pairs
        File control_weights = SelectControls.control_weights
    }
}
