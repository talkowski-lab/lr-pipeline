#!/usr/bin/env bash
# Prepare depth (RD), BAF and genotype (GT) tables for one region, with samples subset and renamed, in one run:
#   1. RD:  tabix the region out of the cohort RD matrix (header kept)               -> <name>.RD.txt.gz
#   2. VCF: fetch the region from the per-chromosome small-variant VCF (local or gs://, requester pays via pysam)
#                                                                                     -> <name>.vcf.gz (+ .tbi)
#   3. BAF: per-sample BAF from FORMAT/AD (vcf_to_baf_table.py)                        -> <name>.baf.tsv.gz
#   4. GT:  per-sample genotypes (extract_gt.sh)                                       -> <name>.GT.tsv.gz
#   5. rename: keep and rename the samples of --id-map (subset_rename_samples.py)     -> <name>.renamed.{RD.txt,baf.tsv,GT.tsv}.gz
#
# Usage:
#   prepare_region_data.sh --chr chr16 --start 30088531 --end 30343458 --name 16p11.2_BP5 \
#     --rd aou_lrGS_GD_carriers.RD.txt.gz --vcf "gs://bucket/path/{chr}.g.vcf.bgz" \
#     --samples GD_carrier.samples --id-map AoU_srWGS_SV.v7OC.extractGDs.dropoverlapping.lrGS_GD_samples.tsv \
#     [--outdir .] [--scripts-dir <folder with the helper scripts; default: this script's folder>]
#     [--baf-script vcf_to_baf_table.py] [--requester-pays-project PROJECT] [--keep-intermediate] [--rd-only]
#   --rd-only  only the depth table (<name>.renamed.RD.txt.gz); --vcf and --samples are then not needed
#
#   --chr/--start/--end  region (1-based, inclusive); --name: output prefix (e.g. 16p11.2_BP5 or 15q11.2_BP1)
#   --rd       cohort RD matrix, bgzipped + tabix-indexed, header line '#Chr Start End <samples>'
#   --vcf      VCF path; '{chr}' is replaced by --chr (e.g. one VCF per chromosome). gs:// paths (VCF and --rd) are
#              read with pysam/htslib: the OAuth token comes from $GCS_OAUTH_TOKEN if set, else `gcloud auth
#              print-access-token`, else the VM metadata server; the requester-pays project from
#              --requester-pays-project, else $GCS_REQUESTER_PAYS_PROJECT, else $GOOGLE_PROJECT
#   --samples  optional: real sample IDs, one per line, for the BAF and GT tables; default: column 1 of --id-map without
#              its header line (--id-map-no-header if it has none); either way limited to samples present in the VCF
#   --id-map   2 columns: real sample ID, new sample ID (samples missing from a table are skipped)
#   --keep-intermediate  keep the non-renamed tables and the region VCF (removed by default)
# Needs: tabix, bgzip, bcftools, python3 with pysam and pandas; helper scripts vcf_to_baf_table.py (or --baf-script),
# extract_gt.sh, subset_rename_samples.py.
set -euo pipefail

usage() { sed -n '2,29p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }
chr=""; start=""; end=""; name=""; rd=""; vcf=""; samples=""; idmap=""; outdir="."; keep=0; rd_only=0; map_header=1
sdir="$(cd "$(dirname "$0")" && pwd)"; baf_script="vcf_to_baf_table.py"; rp_project=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --chr) chr=$2; shift 2 ;; --start) start=$2; shift 2 ;; --end) end=$2; shift 2 ;; --name) name=$2; shift 2 ;;
    --rd) rd=$2; shift 2 ;; --vcf) vcf=$2; shift 2 ;; --samples) samples=$2; shift 2 ;; --id-map) idmap=$2; shift 2 ;;
    --outdir) outdir=$2; shift 2 ;; --scripts-dir) sdir=$2; shift 2 ;; --baf-script) baf_script=$2; shift 2 ;;
    --requester-pays-project) rp_project=$2; shift 2 ;; --keep-intermediate) keep=1; shift ;;
    --rd-only) rd_only=1; shift ;; --id-map-no-header) map_header=0; shift ;; -h|--help) usage ;; *) echo "unknown option $1" >&2; usage ;;
  esac
done
required=(chr start end name rd idmap); [[ $rd_only -eq 0 ]] && required+=(vcf)
for v in "${required[@]}"; do [[ -z "${!v}" ]] && { echo "missing --${v/idmap/id-map}" >&2; usage; }; done
helpers=(subset_rename_samples.py); [[ $rd_only -eq 0 ]] && helpers+=("$baf_script" extract_gt.sh)
for s in "${helpers[@]}"; do
  [[ -f "$sdir/$s" ]] || { echo "helper script not found: $sdir/$s (set --scripts-dir)" >&2; exit 1; }
done
mkdir -p "$outdir"
region="${chr}:${start}-${end}"
p="$outdir/$name"
vcf=${vcf//\{chr\}/$chr}
echo "region $region -> $p.*" >&2

# Credentials for gs:// inputs (htslib reads GCS_OAUTH_TOKEN and GCS_REQUESTER_PAYS_PROJECT)
if [[ "$vcf" == gs://* || "$rd" == gs://* ]]; then
  if [[ -z "${GCS_OAUTH_TOKEN:-}" ]]; then
    if command -v gcloud >/dev/null; then
      GCS_OAUTH_TOKEN="$(gcloud auth print-access-token)"
    else
      GCS_OAUTH_TOKEN="$(curl -sS -H "Metadata-Flavor: Google" \
        http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
        | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')"
    fi
    export GCS_OAUTH_TOKEN
  fi
  rp_project="${rp_project:-${GCS_REQUESTER_PAYS_PROJECT:-${GOOGLE_PROJECT:-}}}"
  [[ -n "$rp_project" ]] && export GCS_REQUESTER_PAYS_PROJECT="$rp_project"
fi

# 1. RD matrix for the region, header kept (pysam: local or gs://, needs the .tbi next to the matrix)
python3 - "$rd" "$chr" "$start" "$end" <<'PY' | bgzip -c > "$p.RD.txt.gz"
import sys
import pysam
path, chrom, start, end = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
t = pysam.TabixFile(path)
for h in t.header:
    print(h)
n = 0
for line in t.fetch(chrom, start - 1, end):
    print(line)
    n += 1
print(f"RD: {n} bins", file=sys.stderr)
PY

outs=("$p.renamed.RD.txt.gz")
if [[ $rd_only -eq 0 ]]; then
  # 2. Region VCF (pysam handles gs:// with an OAuth token and requester pays)
  python3 - "$vcf" "$chr" "$start" "$end" "$p.vcf.gz" <<'PY'
import sys
import pysam
path, chrom, start, end, out_path = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
vcf = pysam.VariantFile(path)
if chrom not in vcf.header.contigs:
    sys.exit(f"{chrom} is not a contig of {path}")
out = pysam.VariantFile(out_path, "wz", header=vcf.header)
n = 0
for rec in vcf.fetch(chrom, start - 1, end):
    out.write(rec)
    n += 1
out.close()
print(f"VCF: {n} records, {len(vcf.header.samples)} samples", file=sys.stderr)
PY
  tabix -f -p vcf "$p.vcf.gz"

  # Samples for BAF/GT: --samples or column 1 of the ID map, limited to samples in the VCF
  src="$samples"
  if [[ -z "$src" ]]; then
    src="$p.idmap_samples.txt"
    tail -n +$((map_header + 1)) "$idmap" | awk 'NF{print $1}' > "$src"
  fi
  bcftools query -l "$p.vcf.gz" | grep -xFf - "$src" > "$p.samples.txt" || true
  echo "samples: $(wc -l < "$p.samples.txt" | tr -d ' ') of $(awk 'NF' "$src" | wc -l | tr -d ' ') requested are in the VCF" >&2
  [[ -s "$p.samples.txt" ]] || { echo "none of the requested samples is in the VCF" >&2; exit 1; }
  samples="$p.samples.txt"

  # 3-4. BAF and GT tables for the selected samples
  python3 "$sdir/$baf_script" "$p.vcf.gz" "$samples" "$p.baf.tsv.gz" --region "$region"
  bash "$sdir/extract_gt.sh" -v "$p.vcf.gz" -s "$samples" -o "$p.GT.tsv.gz" -r "$region"
  outs+=("$p.renamed.baf.tsv.gz" "$p.renamed.GT.tsv.gz")
fi

# 5. Subset + rename samples
python3 "$sdir/subset_rename_samples.py" --input "$p.RD.txt.gz" --id-map "$idmap" --out "$p.renamed.RD.txt.gz" \
  --fixed-cols 3 --allow-missing
if [[ $rd_only -eq 0 ]]; then
  python3 "$sdir/subset_rename_samples.py" --input "$p.baf.tsv.gz" --id-map "$idmap" --out "$p.renamed.baf.tsv.gz" \
    --fixed-cols 5 --allow-missing
  python3 "$sdir/subset_rename_samples.py" --input "$p.GT.tsv.gz" --id-map "$idmap" --out "$p.renamed.GT.tsv.gz" \
    --fixed-cols 5 --allow-missing
fi

# Check: every table is on the requested chromosome
for f in "${outs[@]}"; do
  bad=$(gzip -cd "$f" | awk -F'\t' -v c="$chr" 'NR > 1 && $1 != c' | wc -l | tr -d ' ')
  if [[ "$bad" -gt 0 ]]; then echo "WARNING: $f has $bad rows not on $chr" >&2; fi
done
if [[ $keep -eq 0 ]]; then
  rm -f "$p.RD.txt.gz" "$p.baf.tsv.gz" "$p.GT.tsv.gz" "$p.vcf.gz" "$p.vcf.gz.tbi" "$p.samples.txt" "$p.idmap_samples.txt"
fi
echo "done: ${outs[*]}" >&2
