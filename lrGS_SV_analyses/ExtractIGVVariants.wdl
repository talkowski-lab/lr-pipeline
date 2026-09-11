version 1.0

## Given a list of variant IDs and one VCF per contig, scan each VCF for
## matching records and combine the hits into one genome-wide TSV in the
## same format as example.igv_variants.tsv.
##
## Output columns: chrom, start, end, ID, allele_type, samples
## `samples` lists every carrier (non-ref, non-missing GT) for that variant,
## comma-separated.

workflow ExtractIGVVariants {
  input {
    Array[String] variant_ids
    Array[File] vcfs
    String output_basename = "igv_variants"
    String docker = "python:3.11-slim"
  }

  File ids_file = write_lines(variant_ids)

  scatter (vcf in vcfs) {
    call ExtractFromVcf {
      input:
        vcf = vcf,
        ids_file = ids_file,
        docker = docker,
    }
  }

  call CombineAndSort {
    input:
      tsvs = ExtractFromVcf.out_tsv,
      ids_file = ids_file,
      output_basename = output_basename,
      docker = docker,
  }

  output {
    File igv_variants_tsv = CombineAndSort.combined_tsv
    File missing_ids = CombineAndSort.missing_ids
  }
}

task ExtractFromVcf {
  input {
    File vcf
    File ids_file
    String docker
  }

  String base = basename(vcf, ".vcf.gz")
  Int disk_gb = ceil(size(vcf, "GB") * 2) + 20

  command <<<
    set -euo pipefail

    python3 <<CODE
import gzip

with open("~{ids_file}") as f:
    wanted = {line.strip() for line in f if line.strip()}

def open_vcf(path):
    with open(path, "rb") as fh:
        is_gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gz else open(path, "rt")

samples = []
rows = []
with open_vcf("~{vcf}") as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.rstrip("\n").split("\t")[9:]
            continue
        fields = line.rstrip("\n").split("\t")
        vid = fields[2]
        if vid not in wanted:
            continue
        chrom, pos, _id, ref, alt = fields[0:5]
        info = {}
        for kv in fields[7].split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                info[k] = v
            else:
                info[kv] = "1"
        end = info.get("END", str(int(pos) + len(ref) - 1))
        allele_type = info.get("allele_type")
        if allele_type is None:
            alt0 = alt.split(",")[0]
            if len(ref) == len(alt0) == 1:
                allele_type = "snv"
            elif len(alt0) > len(ref):
                allele_type = "ins"
            elif len(alt0) < len(ref):
                allele_type = "del"
            else:
                allele_type = "."
        fmt = fields[8].split(":")
        gt_idx = fmt.index("GT")
        carriers = []
        for sample, sample_field in zip(samples, fields[9:]):
            gt = sample_field.split(":")[gt_idx]
            alleles = gt.replace("|", "/").split("/")
            if any(a not in ("0", ".", "") for a in alleles):
                carriers.append(sample)
        rows.append([chrom, pos, end, vid, allele_type, ",".join(carriers)])

with open("~{base}.igv_variants.tsv", "w") as out:
    for row in rows:
        out.write("\t".join(row) + "\n")
CODE
  >>>

  output {
    File out_tsv = "~{base}.igv_variants.tsv"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task CombineAndSort {
  input {
    Array[File] tsvs
    File ids_file
    String output_basename
    String docker
  }

  command <<<
    set -euo pipefail

    printf "#chrom\tstart\tend\tID\tallele_type\tsamples\n" > ~{output_basename}.tsv
    cat ~{sep=" " tsvs} | sort -k1,1 -k2,2n >> ~{output_basename}.tsv

    tail -n +2 ~{output_basename}.tsv | cut -f4 | sort -u > found_ids.txt
    sort -u ~{ids_file} > wanted_ids.txt
    comm -23 wanted_ids.txt found_ids.txt > ~{output_basename}.missing_ids.txt
  >>>

  output {
    File combined_tsv = "~{output_basename}.tsv"
    File missing_ids = "~{output_basename}.missing_ids.txt"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 20 HDD"
  }
}
