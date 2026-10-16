version 1.0

## Scatter across per-chromosome VCFs, extract singleton (AC=1) variants whose
## carrier is supported ONLY by assembly-based methods (dipcall, hapdiff -- no
## kanpig, no read-based caller), then combine into one genome-wide TSV.
##
## Output columns: ID, allele_type, allele_length, dbSNP_ID, gnomAD_V4_match_ID, FILTER
##
## `extraction_script` must point at a copy of extract_singleton_assembly_only.py
## uploaded to a GCS path Terra can read (e.g. your workspace bucket) -- pass
## its gs:// URI in the inputs JSON. No custom Docker image is required: the
## script is localized like any other File input and run with plain python3.

workflow SingletonAssemblyOnlyScan {
  input {
    Array[File] vcfs
    Array[File] vcf_indices
    File extraction_script
    String output_basename = "singleton_assembly_only"
    String docker = "python:3.11-slim"
  }

  scatter (pair in zip(vcfs, vcf_indices)) {
    call ExtractSingletonAssemblyOnly {
      input:
        vcf = pair.left,
        vcf_index = pair.right,
        extraction_script = extraction_script,
        docker = docker,
    }
  }

  call CombineTsvs {
    input:
      tsvs = ExtractSingletonAssemblyOnly.out_tsv,
      output_basename = output_basename,
      docker = docker,
  }

  output {
    File combined_tsv = CombineTsvs.combined_tsv
    Array[File] per_chrom_tsvs = ExtractSingletonAssemblyOnly.out_tsv
    Array[File] per_chrom_logs = ExtractSingletonAssemblyOnly.log
  }
}

task ExtractSingletonAssemblyOnly {
  input {
    File vcf
    File vcf_index
    File extraction_script
    String docker
  }

  String base = basename(vcf, ".vcf.gz")
  Int disk_gb = ceil(size(vcf, "GB") * 2) + 20

  command <<<
    set -euo pipefail
    pip install --quiet --no-cache-dir pysam
    python3 ~{extraction_script} ~{vcf} ~{base}.singleton_assembly_only.tsv 2> ~{base}.log
  >>>

  output {
    File out_tsv = "~{base}.singleton_assembly_only.tsv"
    File log = "~{base}.log"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task CombineTsvs {
  input {
    Array[File] tsvs
    String output_basename
    String docker
  }

  command <<<
    set -euo pipefail
    files=(~{sep=" " tsvs})
    head -n 1 "${files[0]}" > ~{output_basename}.tsv
    for f in "${files[@]}"; do
      tail -n +2 "$f" >> ~{output_basename}.tsv
    done
    wc -l ~{output_basename}.tsv
  >>>

  output {
    File combined_tsv = "~{output_basename}.tsv"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 20 HDD"
  }
}
