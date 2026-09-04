version 1.0

## Scatter across chromosomes: for each chromosome, localize only the Hail
## Table partitions needed for that chromosome from a gnomAD-style coverage
## Hail Table (gs://.../*.ht -- a directory, not a single file), bin into
## fixed-size windows, and compute the mean of the table's "mean" column per
## bin. Combine into one genome-wide BED, and also report the per-chromosome
## BEDs individually.
##
## `binning_script` must point at a copy of extract_chrom_coverage_bins.py
## uploaded to a GCS path Terra can read (e.g. your workspace bucket). No
## custom Docker image is required: the script pip installs its one
## dependency (google-cloud-storage) at runtime on top of the official Hail
## Docker image.

workflow ExtractHailCoverageBins {
  input {
    String hail_table_gcs_path
    Array[String] chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    Int bin_size = 100
    File binning_script
    String output_basename = "coverage_bins"
    String docker = "hailgenetics/hail:0.2.135"
  }

  scatter (chrom in chroms) {
    call ExtractChromBins {
      input:
        hail_table_gcs_path = hail_table_gcs_path,
        chrom = chrom,
        bin_size = bin_size,
        binning_script = binning_script,
        docker = docker,
    }
  }

  call CombineBeds {
    input:
      beds = ExtractChromBins.out_bed,
      output_basename = output_basename,
      docker = docker,
  }

  output {
    File combined_bed = CombineBeds.combined_bed
    Array[File] per_chrom_beds = ExtractChromBins.out_bed
    Array[File] per_chrom_logs = ExtractChromBins.log
  }
}

task ExtractChromBins {
  input {
    String hail_table_gcs_path
    String chrom
    Int bin_size
    File binning_script
    String docker
  }

  command <<<
    set -euo pipefail
    pip install --quiet --no-cache-dir google-cloud-storage
    python3 ~{binning_script} \
      --hail-table-path ~{hail_table_gcs_path} \
      --chrom ~{chrom} \
      --bin-size ~{bin_size} \
      --out-bed ~{chrom}.~{bin_size}bp_bins.bed \
      2> ~{chrom}.log
  >>>

  output {
    File out_bed = "~{chrom}.~{bin_size}bp_bins.bed"
    File log = "~{chrom}.log"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "16 GB"
    disks: "local-disk 50 HDD"
    preemptible: 2
  }
}

task CombineBeds {
  input {
    Array[File] beds
    String output_basename
    String docker
  }

  command <<<
    set -euo pipefail
    cat ~{sep=" " beds} > ~{output_basename}.bed
    wc -l ~{output_basename}.bed
  >>>

  output {
    File combined_bed = "~{output_basename}.bed"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 20 HDD"
  }
}
