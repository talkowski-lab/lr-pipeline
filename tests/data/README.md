# Test data


## hprc_hgsvc_v10.chr22.subset.vcf.gz
A 385-record subset of the open-access HPRC/HGSVC V10 cohort callset on chr22, restricted to five samples: HG00096, HG00097, HG00099, HG00126 and HG00128. Every record has `AC` above zero in those samples. It is the shared input of the workflow smoke tests.

Records per `INFO/allele_type`, `INFO/SOURCE` and `FILTER` stratum:

| Stratum | Records |
| --- | --- |
| snv, DeepVariant, PASS | 100 |
| trv, Vamos, PASS | 40 |
| trv, TRExplorer, PASS | 40 |
| ins, HPRC_SV_Integration, PASS | 40 |
| ins, DeepVariant, PASS | 40 |
| del, DeepVariant, PASS | 40 |
| del, HPRC_SV_Integration, PASS | 32 |
| del, DeepVariant, MONOALLELIC | 27 |
| dup, HPRC_SV_Integration, PASS | 19 |
| del, HPRC_SV_Integration, ASSEMBLY_ONLY_SINGLETON | 3 |
| ins, HPRC_SV_Integration, ASSEMBLY_ONLY_SINGLETON | 2 |
| ins, DeepVariant, MONOALLELIC | 2 |

Of these, 76 records are multiallelic.


### Regeneration
Set `SRC` to the full chr22 callset and `SAMPLES` to the five samples above, joined by commas. The steps below reproduce the committed records exactly.

```bash
# Subset the four windows to the five samples, keeping variable sites
bcftools view -s "$SAMPLES" --min-ac 1 \
    -r chr22:20000000-21000000,chr22:30000000-31000000,chr22:40000000-41000000,chr22:45000000-46000000 \
    -Oz -o windows.vcf.gz "$SRC"

# Keep every k-th record per stratum, with k = floor(count / target) and a target of 100 for snv and 40 otherwise
python sample_strata.py windows.vcf.gz windows.sampled.vcf.gz

# Select every chr22-wide record that is a dup or has an assembly-only or single-read filter
bcftools view -s "$SAMPLES" --min-ac 1 \
    -i 'INFO/allele_type="dup" || FILTER="ASSEMBLY_ONLY_SINGLETON" || FILTER="SINGLE_READ_SUPPORT"' \
    -Oz -o rare.vcf.gz "$SRC"

# Index both pieces, concatenate and sort them, then drop the VEP and histogram INFO fields
bcftools index -t windows.sampled.vcf.gz
bcftools index -t rare.vcf.gz
bcftools concat -a -D windows.sampled.vcf.gz rare.vcf.gz \
    | bcftools sort \
    | bcftools annotate \
        -x INFO/vep,INFO/gq_hist_all_bin_freq,INFO/gq_hist_alt_bin_freq,INFO/dp_hist_all_bin_freq,INFO/dp_hist_alt_bin_freq,INFO/sd_hist_all_bin_freq,INFO/sd_hist_alt_bin_freq,INFO/ab_hist_alt_bin_freq \
        -Oz -o hprc_hgsvc_v10.chr22.subset.vcf.gz
tabix -p vcf hprc_hgsvc_v10.chr22.subset.vcf.gz
```

`sample_strata.py` groups the window records on `allele_type|SOURCE|FILTER`:

```python
import collections
import sys

import pysam

src, dst = sys.argv[1:3]
groups = collections.defaultdict(list)
with pysam.VariantFile(src) as vcf:
    header = vcf.header
    for record in vcf:
        stratum = (record.info["allele_type"], record.info["SOURCE"], ",".join(record.filter.keys()))
        groups[stratum].append(record.copy())

selected = []
for stratum, records in groups.items():
    target = 100 if stratum[0] == "snv" else 40
    step = max(len(records) // target, 1)
    selected.extend(records[::step][:target])

with pysam.VariantFile(dst, "w", header=header) as out:
    for record in sorted(selected, key=lambda r: r.pos):
        out.write(record)
```


## GRCh38.chr22.fa.fai
A one-line FASTA index for chr22 (`chr22`, 50,818,468 bp). Workflows that read only the contig lengths take it in place of the reference.


## chr22.windows.fa
A 60 kb FASTA of 60 reference windows of 1,000 bp, joined into one contig named `chr22`, with its `.fai`. It lets `PreprocessVcfs` and `MergeVcfs` normalize and check REF alleles without a full reference. Window `i` holds the reference around the `i`-th record chosen by `window_records()` in `tests/conftest.py`, with the record's REF starting at the 401st base. The `window_vcf` fixture moves those records to position `i * 1000 + 401`, so their REF alleles match the FASTA.


### Regeneration
Run from `tests/`, with `SRC` set to an uppercase or soft-masked GRCh38 FASTA that has an index:

```bash
python - <<'PY'
import conftest as c
_, records = c.load_records(c.SUBSET_VCF)
c.write_window_fasta("$SRC", c.WINDOW_FASTA, records)
PY
```

Index the output with `samtools faidx chr22.windows.fa`.
