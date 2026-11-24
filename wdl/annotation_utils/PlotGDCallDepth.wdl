version 1.0

import "../utils/Structs.wdl"
import "SummarizeGDCallDepth.wdl"

workflow PlotGDCallDepth {
    input {
        File depth_traces
        File call_windows
        File depth_bin_stats
        String prefix

        Boolean nahr_check = true
        Int outside_bins = 5
        Float max_dev = 0.5
        File? manual_labels
        File? keep_samples
        Int raster_px = 2000

        File? sd_pairs_bed
        File? segdup_bed
        File? genes_bed
        File? gene_info

        String plot_docker

        RuntimeAttr? runtime_attr_plot_region_depth
    }

    call SummarizeGDCallDepth.PlotRegionDepth {
        input:
            traces = depth_traces,
            windows = call_windows,
            bin_stats = depth_bin_stats,
            sd_pairs_bed = sd_pairs_bed,
            segdup_bed = segdup_bed,
            genes_bed = genes_bed,
            gene_info = gene_info,
            manual_labels = manual_labels,
            keep_samples = keep_samples,
            nahr_check = nahr_check,
            outside_bins = outside_bins,
            max_dev = max_dev,
            raster_px = raster_px,
            prefix = prefix,
            docker = plot_docker,
            runtime_attr_override = runtime_attr_plot_region_depth
    }

    output {
        File region_depth_pdf = PlotRegionDepth.region_depth_pdf
        File region_counts = PlotRegionDepth.region_counts
        File? carrier_flags = PlotRegionDepth.carrier_flags
    }
}
