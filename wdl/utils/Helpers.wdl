version 1.0

import "Structs.wdl"

task AddFilter {
    input {
        File vcf
        File vcf_idx
        String filter_name
        String filter_description
        String filter_expression
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -h \
            ~{vcf} \
        | grep "^##" > header.txt

        echo '##FILTER=<ID=~{filter_name},Description="~{filter_description}">' >> header.txt

        bcftools view \
            -h \
            ~{vcf} \
        | grep "^#CHROM" >> header.txt

        bcftools reheader \
            -h header.txt \
            ~{vcf} \
        | bcftools filter \
            --mode + \
            -s ~{filter_name} \
            -e '~{filter_expression}' \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File flagged_vcf = "~{prefix}.vcf.gz"
        File flagged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task AddInfo {
    input {
        File vcf
        File vcf_idx
        String tag_id
        String tag_value
        String tag_description
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo '##INFO=<ID=~{tag_id},Number=1,Type=String,Description="~{tag_description}">' > header.lines

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t~{tag_value}\n' \
            ~{vcf} \
        | bgzip -c > annotations.txt.gz

        tabix -s1 -b2 -e2 annotations.txt.gz

        bcftools annotate -h header.lines -a annotations.txt.gz \
            -c CHROM,POS,REF,ALT,~ID,INFO/~{tag_id} \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task AddMissingInfoHeaderLines {
    input {
        File vcf
        File vcf_idx
        Array[String] info_fields
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        : > hdr_lines.txt
        while IFS= read -r field; do
            printf '##INFO=<ID=%s,Number=.,Type=String,Description="Auto-added missing header line">\n' "${field}" >> hdr_lines.txt
        done < ~{write_lines(info_fields)}

        bcftools annotate \
            -h hdr_lines.txt \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task AddTREndTag {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} \
            | awk 'BEGIN{FS=OFS="\t"; has_end_header=0}

                /^##INFO=<ID=END,/ {
                    has_end_header=1
                    print
                    next
                }
                /^#CHROM/ {
                    if (!has_end_header) {
                        print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"
                    }
                    print
                    next
                }
                /^#/ {
                    print
                    next
                }
                {
                    end_val=$2 + length($4) - 1

                    if ($8=="." || $8=="") {
                        $8="END=" end_val
                    } else if ($8 ~ /(^|;)END=/) {
                        n=split($8, info_parts, ";")
                        for (i=1; i<=n; i++) {
                            if (info_parts[i] ~ /^END=/) {
                                info_parts[i]="END=" end_val
                            }
                        }
                        $8=info_parts[1]
                        for (i=2; i<=n; i++) {
                            $8=$8 ";" info_parts[i]
                        }
                    } else {
                        $8=$8 ";END=" end_val
                    }

                    print
                }' \
            | bgzip -c > ~{prefix}.vcf.gz

        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File vcf_with_end = "~{prefix}.vcf.gz"
        File vcf_with_end_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size([vcf, vcf_idx], "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task AnnotateVariantAttributes {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        touch new_headers.txt
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=allele_length'; then
            echo '##INFO=<ID=allele_length,Number=1,Type=Integer,Description="Allele length">' >> new_headers.txt
        fi
        if ! bcftools view -h ~{vcf} | grep -q '##INFO=<ID=allele_type'; then
            echo '##INFO=<ID=allele_type,Number=1,Type=String,Description="Allele type">' >> new_headers.txt
        fi

        bcftools annotate \
            -h new_headers.txt \
            -Oz -o temp.vcf.gz \
            ~{vcf}

        tabix -p vcf temp.vcf.gz

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/allele_length\t%INFO/allele_type\n' \
            temp.vcf.gz \
        | awk -F'\t' '{
            ref_length = length($3)
            alt_len = length($4)
            calc_length = alt_len - ref_length
            calc_type = "snv"

            if (alt_len > ref_length) {
                calc_type = "ins"
            } else if (alt_len < ref_length) {
                calc_type = "del"
            }

            allele_length = ($6 == ".") ? calc_length : $6
            allele_type = ($7 == ".") ? calc_type : $7

            print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"allele_length"\t"allele_type
        }' \
            | bgzip -c > annot.txt.gz

        tabix -s1 -b2 -e2 annot.txt.gz

        bcftools annotate \
            -a annot.txt.gz \
            -c CHROM,POS,REF,ALT,~ID,INFO/allele_length,INFO/allele_type \
            -Oz -o ~{prefix}.vcf.gz \
            temp.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([vcf, vcf_idx], "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task AppendAnnotationsFromVcf {
    input {
        File annotation_tsv
        File truth_vcf
        File truth_vcf_idx
        Boolean is_sv_truth
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import subprocess
import re

annotation_tsv = "~{annotation_tsv}"
truth_vcf = "~{truth_vcf}"
is_sv_truth = ~{true="True" false="False" is_sv_truth}
prefix = "~{prefix}"

def get_ac_af_an_fields(vcf_path):
    cmd = f"bcftools view -h {vcf_path}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    fields = {'AC': {}, 'AF': {}, 'AN': {}}
    for line in proc.stdout:
        m = re.match(r'##INFO=<ID=([^,]+)', line)
        if m:
            fid = m.group(1)
            norm_id = fid.upper().replace('_REMAINING', '_RMI')
            for p in ['AC', 'AF', 'AN']:
                if norm_id == p or norm_id.startswith(p + '_'):
                    fields[p][norm_id] = fid
    proc.wait()
    return fields

vcf_fields = get_ac_af_an_fields(truth_vcf)
dyn_cols = sorted(vcf_fields['AC']) + sorted(vcf_fields['AF']) + sorted(vcf_fields['AN'])
norm_to_orig = {**vcf_fields['AC'], **vcf_fields['AF'], **vcf_fields['AN']}

if is_sv_truth:
    extra_fields = ['N_HOMREF', 'N_HET', 'N_HOMALT']
else:
    extra_fields = ['nhomalt']

query_field_pairs = [(c, norm_to_orig[c]) for c in dyn_cols] + [(f, f) for f in extra_fields]
fmt = '%ID\\t' + '\\t'.join(f'%INFO/{orig}' for _, orig in query_field_pairs) + '\\n'

# Collect the truth IDs the annotations actually reference so only those rows are retained below
wanted_ids = set()
with open(annotation_tsv) as fin:
    for line in fin:
        fields = line.rstrip('\n').split('\t')
        if len(fields) > 6:
            wanted_ids.add(fields[6])

cmd = f"bcftools query -f '{fmt}' {truth_vcf}"
proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
truth_info = {}
for line in proc.stdout:
    parts = line.rstrip('\n').split('\t')
    if len(parts) == len(query_field_pairs) + 1 and parts[0] in wanted_ids:
        truth_info[parts[0]] = {query_field_pairs[i][0]: parts[i + 1] for i in range(len(query_field_pairs))}
proc.wait()

def to_num(val):
    try:
        return float(val) if val and val != '.' else 0
    except Exception:
        return 0

def compute_genotype_counts(info):
    if is_sv_truth:
        return info.get('N_HOMREF', '.'), info.get('N_HET', '.'), info.get('N_HOMALT', '.')
    homalt = to_num(info.get('nhomalt', '.'))
    het = to_num(info.get('AC', '.')) - 2 * homalt
    homref = to_num(info.get('AN', '.')) / 2 - homalt - het
    return str(int(homref)), str(int(het)), info.get('nhomalt', '.')

extra_cols = ['match_type', 'truth_ID', 'source_tag', 'filter'] + dyn_cols + ['N_HOMREF', 'N_HET', 'N_HOMALT']
header_row = '\t'.join(['#CHROM', 'POS', 'REF', 'ALT', 'ID'] + extra_cols)

with open(annotation_tsv) as fin, open(f"{prefix}.tsv", 'w') as fout:
    fout.write(header_row + '\n')
    for line in fin:
        fields = line.rstrip('\n').split('\t')
        truth_id = fields[6]
        info = truth_info.get(truth_id, {})
        dyn_vals = [info.get(f, '.') for f in dyn_cols]
        homref, het, homalt = compute_genotype_counts(info)
        fout.write('\t'.join(fields + dyn_vals + [homref, het, homalt]) + '\n')

EOF
    >>>

    output {
        File annotated_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 6,
        disk_gb: 2 * ceil(size(annotation_tsv, "GB") + size(truth_vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task BedtoolsClosest {
    input {
        File bed_a
        File bed_b
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        paste <(head -1 ~{bed_a}) <(head -1 ~{bed_b}) \
            | sed -e "s/#//g" \
            > ~{prefix}.bed

        bedtools closest \
            -wo \
            -a <(sort -k1,1 -k2,2n ~{bed_a}) \
            -b <(sort -k1,1 -k2,2n ~{bed_b}) \
            | awk -F'\t' '$7 != "."' \
            >> ~{prefix}.bed
    >>>

    output {
        File output_bed = "~{prefix}.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bed_a, "GB") + size(bed_b, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task CheckSampleConsistency {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] sample_ids
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        printf '%s\n' ~{sep=' ' sample_ids} | sort > requested_samples.txt

        vcfs_array=(~{sep=' ' vcfs})

        for vcf in "${vcfs_array[@]}"; do
            bcftools query -l "$vcf" | sort > vcf_samples.txt

            comm -3 requested_samples.txt vcf_samples.txt > differences.txt

            if [ -s differences.txt ]; then
                echo "ERROR: Sample mismatch in $vcf"

                echo "--- Missing from VCF ---"
                comm -23 requested_samples.txt vcf_samples.txt

                echo "--- Extra in VCF (not requested) ---"
                comm -13 requested_samples.txt vcf_samples.txt

                exit 1
            fi
        done
    >>>

    output {
        String status = "success"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task CollapseMultiallelics {
    input {
        File tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import subprocess

groups = {}
with open('unsorted.tsv', 'w') as out:
    with open('~{tsv}', 'r') as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            vid = fields[4]
            if 'TRV' not in vid:
                out.write(line)
            else:
                if vid not in groups:
                    groups[vid] = fields[:]
                else:
                    groups[vid][3] += ',' + fields[3]
                    groups[vid][5] += ',' + fields[5]
    for row in groups.values():
        out.write('\t'.join(row) + '\n')

subprocess.run(['sort', '-k1,1', '-k2,2n', 'unsorted.tsv', '-o', '~{prefix}.tsv'], check=True)
CODE
    >>>

    output {
        File collapsed_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsv, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task ConcatAlignedTsvs {
    input {
        Array[File] tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import sys
import csv

input_files = "~{sep=',' tsvs}".split(',')
output_filename = "aligned_unsorted.tsv"
header_filename = "~{prefix}.header.txt"
fixed_cols = ["#CHROM", "POS", "REF", "ALT", "ID"]

all_keys = set()
for f in input_files:
    with open(f, 'r') as fh:
        line = fh.readline().strip()
        if not line: continue
        parts = line.split('\t')
        if len(parts) > 5:
            keys = parts[5:]
            all_keys.update(keys)

sorted_keys = sorted(list(all_keys))
master_header = fixed_cols + sorted_keys

with open(header_filename, 'w') as hout:
    for k in sorted_keys:
        hout.write(k + "\n")

with open(output_filename, 'w') as out:
    out.write("\t".join(master_header) + "\n")

    for f in input_files:
        with open(f, 'r') as fh:
            header_line = fh.readline().strip()
            if not header_line:
                continue

            file_cols = header_line.split('\t')
            col_map = {name: i for i, name in enumerate(file_cols)}

            for line in fh:
                parts = line.strip().split('\t')
                if not parts: continue

                out_row = []
                for target_col in master_header:
                    if target_col in col_map:
                        try:
                            val = parts[col_map[target_col]]
                            out_row.append(val)
                        except IndexError:
                            out_row.append(".")
                    else:
                        out_row.append(".")

                out.write("\t".join(out_row) + "\n")
CODE

        tail -n +2 aligned_unsorted.tsv | sort -k1,1 -k2,2n > ~{prefix}.tsv
    >>>

    output {
        File merged_tsv = "~{prefix}.tsv"
        File merged_header = "~{prefix}.header.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsvs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task ConcatTsvs {
    input {
        Array[File] tsvs
        Boolean sort_output
        Boolean preserve_header = false
        Boolean compressed_tsvs = false
        Boolean compressed_output = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Select the decompression command for the input TSVs
        if [ "~{compressed_tsvs}" == "true" ]; then
            CAT_CMD="gunzip -c"
        else
            CAT_CMD="cat"
        fi

        # Concatenate the TSVs, keeping a single header when requested
        if [ "~{preserve_header}" == "true" ]; then
            $CAT_CMD ~{tsvs[0]} | head -n 1 > combined_raw.tsv || true
            for file in ~{sep=' ' tsvs}; do
                $CAT_CMD "$file" | tail -n +2 >> combined_raw.tsv
            done
        else
            $CAT_CMD ~{sep=' ' tsvs} > combined_raw.tsv
        fi

        # Sort by coordinate when requested
        if [ "~{sort_output}" == "true" ]; then
            if [ "~{preserve_header}" == "true" ]; then
                head -n 1 combined_raw.tsv > ~{prefix}.tsv
                tail -n +2 combined_raw.tsv | sort -k1,1 -k2,2n >> ~{prefix}.tsv
            else
                sort -k1,1 -k2,2n combined_raw.tsv > ~{prefix}.tsv
            fi
        else
            mv combined_raw.tsv ~{prefix}.tsv
        fi

        # Compress the output when requested
        if [ "~{compressed_output}" == "true" ]; then
            gzip -1 "~{prefix}.tsv"
        fi
    >>>

    output {
        File concatenated_tsv = if compressed_output then "~{prefix}.tsv.gz" else "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsvs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task ConcatVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Boolean allow_overlaps
        Boolean naive
        Boolean sort_output = false
        Boolean no_version = false
        Boolean no_address = false
        Int? sort_mem_mb
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        VCFS_FILE="~{write_lines(vcfs)}"

        if [[ "~{sort_output}" == "true" ]]; then
            # Give bcftools sort its memory limit in megabytes because a bare --max-mem value is read as bytes
            bcftools concat \
                ~{if no_version then "--no-version" else ""} \
                ~{if allow_overlaps then "--allow-overlaps" else ""} \
                --file-list ${VCFS_FILE} \
                -Ou \
                | bcftools sort \
                    ~{if defined(sort_mem_mb) then "--max-mem " + select_first([sort_mem_mb]) + "M" else ""} \
                    -T ./bcftools-sort.XXXXXX \
                    -Oz -o "~{prefix}.vcf.gz"
        else
            bcftools concat \
                ~{if no_version then "--no-version" else ""} \
                ~{if allow_overlaps then "--allow-overlaps" else ""} \
                ~{if naive then "--naive" else ""} \
                --file-list ${VCFS_FILE} \
                -Oz -o "~{prefix}.vcf.gz"
        fi

        tabix -p vcf -f "~{prefix}.vcf.gz"
    >>>

    output {
        File concat_vcf = "~{prefix}.vcf.gz"
        File concat_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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
        noAddress: no_address
    }
}

task ConcatVcfsLR {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Boolean remove_dup = true
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools concat \
            --allow-overlaps \
            --file-list ~{write_lines(vcfs)} \
            -Oz -o merged.tmp.vcf.gz

        tabix -p vcf merged.tmp.vcf.gz

        if [[ ~{remove_dup} == "true" ]]; then
            bcftools norm \
                -d exact \
                -Oz -o ~{prefix}.vcf.gz \
                merged.tmp.vcf.gz
        else
            mv merged.tmp.vcf.gz ~{prefix}.vcf.gz
        fi

        tabix -p vcf ~{prefix}.vcf.gz

    >>>

    output {
        File concat_vcf = "~{prefix}.vcf.gz"
        File concat_vcf_idx =  "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 25,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task ConsolidateCollapsedSites {
    input {
        File vcf
        File vcf_idx
        Int breakpoint_window
        Float reciprocal_overlap
        Float sample_similarity
        Float sequence_similarity
        Float size_similarity
        Int size_min
        Int size_max
        String keep_strategy
        Boolean set_merge_annotations
        Boolean strip_format_to_gt
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Optionally strip FORMAT fields except GT
        if [[ "~{strip_format_to_gt}" == "true" ]]; then
            FMT_FIELDS=$(bcftools view -h ~{vcf} | grep '^##FORMAT' | grep -v 'ID=GT,' | \
                sed 's/.*ID=\([^,]*\).*/FORMAT\/\1/' | paste -sd',' - || true)
            if [[ -n "$FMT_FIELDS" ]]; then
                bcftools annotate \
                    -x "$FMT_FIELDS" \
                    -Oz -o truvari_input.vcf.gz \
                    ~{vcf}
            else
                cp ~{vcf} truvari_input.vcf.gz
            fi
        else
            cp ~{vcf} truvari_input.vcf.gz
        fi

        tabix -f -p vcf truvari_input.vcf.gz

        truvari collapse \
            -i truvari_input.vcf.gz \
            -o collapsed.vcf \
            -c removed.vcf \
            --keep ~{keep_strategy} \
            --pctovl ~{reciprocal_overlap} \
            --pctseq ~{sequence_similarity} \
            --pctsize ~{size_similarity} \
            --refdist ~{breakpoint_window} \
            --sizemin ~{size_min} \
            --sizemax ~{size_max}

        bgzip -f collapsed.vcf

        bgzip -f removed.vcf

        # Consolidate INFO and GT from the removed variants into retained records
        python3 <<CODE
import sys
import pysam
from collections import defaultdict

def get_nonref_samples(record):
    nonref = set()
    for sample in record.samples:
        gt = record.samples[sample]['GT']
        if gt is not None and any(a is not None and a != 0 for a in gt):
            nonref.add(sample)
    return nonref

def sample_overlap(set_a, set_b):
    if not set_a and not set_b:
        return 1.0
    if not set_a or not set_b:
        return 0.0
    return len(set_a & set_b) / len(set_a | set_b)

sample_sim_threshold = ~{sample_similarity}
set_merge_annot = ~{true="True" false="False" set_merge_annotations}

original_records = {}
orig_vcf = pysam.VariantFile("~{vcf}")
if set_merge_annot:
    if 'MERGE_COUNT' not in orig_vcf.header.info:
        orig_vcf.header.add_line(
            '##INFO=<ID=MERGE_COUNT,Number=1,Type=Integer,'
            'Description="Number of source records merged into this site">'
        )
    if 'MERGE_TYPE' not in orig_vcf.header.info:
        orig_vcf.header.add_line(
            '##INFO=<ID=MERGE_TYPE,Number=1,Type=String,'
            'Description="Merge strategy: EXACT, TRV_EXACT, TRUVARI, or UNIQUE">'
        )
orig_header = orig_vcf.header.copy()
for record in orig_vcf:
    original_records[record.id] = record
orig_vcf.close()

collapse_groups = defaultdict(list)
with pysam.VariantFile("removed.vcf.gz") as rm_vcf:
    for record in rm_vcf:
        match_id = record.info.get("MatchId", None)
        if match_id:
            key = match_id[0] if isinstance(match_id, (tuple, list)) else match_id
            collapse_groups[key].append(record.id)

kept_vcf = pysam.VariantFile("collapsed.vcf.gz")
out_vcf = pysam.VariantFile("~{prefix}.vcf.gz", 'w', header=orig_header)

for kept_record in kept_vcf:
    orig_kept = original_records.get(kept_record.id)
    if orig_kept is None:
        continue

    cluster_size = 1
    collapse_id = kept_record.info.get("CollapseId", None)

    if collapse_id and collapse_id in collapse_groups:
        kept_nonref = get_nonref_samples(orig_kept)
        for removed_id in collapse_groups[collapse_id]:
            orig_removed = original_records.get(removed_id)
            if orig_removed is None:
                continue
            removed_nonref = get_nonref_samples(orig_removed)

            if sample_sim_threshold > 0:
                overlap = sample_overlap(kept_nonref, removed_nonref)
                if overlap < sample_sim_threshold:
                    out_vcf.write(orig_removed)
                    continue

            # Pull INFO fields that exist on the removed record but not the kept record
            for info_key in orig_removed.info.keys():
                if info_key in orig_kept.info:
                    continue
                if info_key not in orig_kept.header.info:
                    continue
                try:
                    orig_kept.info[info_key] = orig_removed.info[info_key]
                except (TypeError, ValueError) as err:
                    hdr = orig_kept.header.info[info_key]
                    print(
                        f"[ConsolidateCollapsedSites] dropped INFO/{info_key} "
                        f"(header Number={hdr.number} Type={hdr.type}; "
                        f"kept {kept_record.id} nalts={len(orig_kept.alts or ())}, "
                        f"removed {removed_id} nalts={len(orig_removed.alts or ())}, "
                        f"value={orig_removed.info[info_key]!r}): {err}",
                        file=sys.stderr,
                    )
                    continue

            # Pull non-ref GTs from the removed record for samples missing on kept
            for sample in orig_kept.samples:
                kept_gt = orig_kept.samples[sample]['GT']
                rm_gt = orig_removed.samples[sample]['GT']
                if kept_gt is None or all(a is None or a == 0 for a in kept_gt):
                    if rm_gt is not None and any(a is not None and a != 0 for a in rm_gt):
                        orig_kept.samples[sample]['GT'] = rm_gt

            cluster_size += 1

    if set_merge_annot:
        orig_kept.info['MERGE_COUNT'] = cluster_size
        orig_kept.info['MERGE_TYPE'] = 'TRUVARI' if cluster_size > 1 else 'UNIQUE'

    out_vcf.write(orig_kept)

kept_vcf.close()
out_vcf.close()
CODE

        bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o sorted.vcf.gz \
            ~{prefix}.vcf.gz

        mv sorted.vcf.gz ~{prefix}.vcf.gz

        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File consolidated_vcf = "~{prefix}.vcf.gz"
        File consolidated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 25,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task ConvertPALMERToVcf {
    input {
        File palmer_calls
        File palmer_tsd_reads
        String mei_type
        String sample
        File ref_fa
        File ref_fai
        String haplotype
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python /opt/scripts/mei/PALMER_to_vcf.py \
            --palmer_calls ~{palmer_calls} \
            --palmer_tsd_reads ~{palmer_tsd_reads} \
            --mei_type ~{mei_type} \
            --sample ~{sample} \
            --ref_fa ~{ref_fa} \
            --ref_fai ~{ref_fai} \
            --haplotype "~{haplotype}" \
        | bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o ~{prefix}.palmer_calls.vcf.gz

        tabix -p vcf ~{prefix}.palmer_calls.vcf.gz
    >>>

    output {
        File vcf = "~{prefix}.palmer_calls.vcf.gz"
        File vcf_idx = "~{prefix}.palmer_calls.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(palmer_calls, "GB") + size(palmer_tsd_reads, "GB") + size(ref_fa, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

# Derived from broadinstitute/long-read-pipelines ConvertToHailMT.wdl.
task ConvertToHailMT {
    meta {
        description: "Convert a .vcf.bgz file to a Hail MatrixTable."
    }

    parameter_meta {
        gvcf: "VCF to convert to a MatrixTable."
        tbi: "Index for gvcf."
        reference: "Reference assembly label; only GRCh38 is supported."
        ref_fa: "Reference sequences FASTA file; downloaded from Hail if not provided."
        ref_fai: "Index for ref_fa; downloaded from Hail if not provided."
        prefix: "Prefix for the output MatrixTable."
        docker: "Docker image for Hail."
        runtime_attr_override: "Override runtime attributes for this task."
    }

    input {
        File gvcf
        File tbi
        String prefix = "out"

        String reference = "GRCh38"
        String? ref_fa
        String? ref_fai

        String docker

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 6 * ceil(size(gvcf, "GB"))

    command <<<
        set -euo pipefail

        python3 <<EOF

        import hail as hl
        hl.init(default_reference='GRCh38')

        if '~{defined(ref_fa)}' == 'true' and '~{defined(ref_fai)}' == 'true':
            ref = hl.ReferenceGenome.from_fasta_file('~{reference}', '~{ref_fa}', '~{ref_fai}')

        callset = hl.import_vcf(
            '~{gvcf}',
            array_elements_required=False,
            force_bgz=True,
            reference_genome='~{reference}'
        )

        callset.write('~{prefix}.mt')

        EOF

        echo "Created matrix table."
        echo "Tarring file now."
        tar -cf ~{prefix}.mt.tar ~{prefix}.mt

        touch completion_key_file
    >>>

    output {
        File mt_tar = "~{prefix}.mt.tar"
        File completion_file = "completion_key_file"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 64,
        disk_gb: disk_size,
        boot_disk_gb: 25,
        preemptible_tries: 0,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConvertToSymbolic {
    input {
        File vcf
        File vcf_idx
        Boolean move_dup_to_origin
        String type_field = "allele_type"
        String length_field = "allele_length"
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Collect the distinct allele types before converting records to symbolic representation
        bcftools query \
            -f '%INFO/~{type_field}\n' \
            ~{vcf} \
        | sort -u > raw_types.txt

        python3 <<CODE
import pysam
import re
import sys

move_dup = ~{true="True" false="False" move_dup_to_origin}

# Treat only an exact DUP as a duplication
def map_type(raw):
    t = raw.upper()
    if t == 'DUP':
        return 'DUP'
    elif 'DEL' in t:
        return 'DEL'
    elif 'INS' in t or 'DUP' in t or 'NUMT' in t:
        return 'INS'
    else:
        print(f"Error: unrecognized allele type '{raw}'", file=sys.stderr)
        sys.exit(1)

ORIGIN_RE = re.compile(r'(chr[^:]+):(\d+)-(\d+)')

def extract_origin_info(origin):
    if origin is None:
        return None, None, None
    origin_str = origin if isinstance(origin, str) else ",".join(origin)
    best, best_len = None, -1
    for val in origin_str.split(","):
        m = ORIGIN_RE.search(val.strip())
        if m:
            chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
            if abs(end - start) > best_len:
                best_len = abs(end - start)
                best = (chrom, start, end)
    return best if best is not None else (None, None, None)

# Collect the mapped allele_type values present in the VCF
with open("raw_types.txt") as f:
    present_types = {map_type(line.strip()) for line in f if line.strip()}

# Build updated header
vcf_in = pysam.VariantFile("~{vcf}")
header = vcf_in.header

# ORIGIN is only needed when DUPs are repositioned onto it, so a callset without the field still converts otherwise
if move_dup and 'ORIGIN' not in header.info:
    print("Error: move_dup_to_origin is true but the VCF header has no INFO/ORIGIN field", file=sys.stderr)
    sys.exit(1)

if 'END' not in header.info:
    header.add_line('##INFO=<ID=END,Number=.,Type=Integer,Description="End position of the variant">')
if 'SVTYPE' not in header.info:
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">')
if 'SVLEN' not in header.info:
    header.add_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Variant length">')
if move_dup and 'ORIGINAL_POS' not in header.info:
    header.add_line('##INFO=<ID=ORIGINAL_POS,Number=1,Type=Integer,Description="POS prior to DUP being repositioned to its source coordinate">')
if move_dup and 'ORIGINAL_CHROM' not in header.info:
    header.add_line('##INFO=<ID=ORIGINAL_CHROM,Number=1,Type=String,Description="CHROM prior to DUP being repositioned to its source coordinate">')
header.add_line('##ALT=<ID=N,Description="Baseline reference">')
for allele_type in present_types:
    if allele_type not in header.alts:
        header.add_line(f'##ALT=<ID={allele_type},Description="{allele_type} variant">')

# Convert each record to a symbolic-allele representation
vcf_out = pysam.VariantFile("unsorted.vcf.gz", 'w', header=header)
for record in vcf_in:
    raw = record.info["~{type_field}"]
    if isinstance(raw, (list, tuple)):
        raw = raw[0]
    allele_type = map_type(raw)

    record.ref = 'N'
    record.alts = (f'<{allele_type}>',)
    record.info['SVTYPE'] = allele_type

    allele_length = record.info["~{length_field}"]
    if isinstance(allele_length, (list, tuple)):
        allele_length = allele_length[0]
    svlen = abs(allele_length)
    record.info['SVLEN'] = svlen

    # Set END, repositioning a DUP onto its ORIGIN only when asked so a callset without ORIGIN still spans an interval
    if move_dup and allele_type == 'DUP':
        origin_chrom, origin_pos, origin_end = extract_origin_info(record.info.get('ORIGIN', None))
        if origin_chrom is None or origin_pos is None or origin_end is None:
            print(f"Error: cannot extract ORIGIN for DUP {record.id} at {record.chrom}:{record.pos} (ORIGIN={record.info.get('ORIGIN')})", file=sys.stderr)
            sys.exit(1)
        if origin_chrom != record.chrom:
            print(f"Error: ORIGIN contig {origin_chrom} differs from record contig for DUP {record.id} at {record.chrom}:{record.pos} (ORIGIN={record.info.get('ORIGIN')})", file=sys.stderr)
            sys.exit(1)
        record.info['ORIGINAL_POS'] = record.pos
        record.info['ORIGINAL_CHROM'] = record.chrom
        record.pos = origin_pos
        record.stop = origin_end
        record.info['SVLEN'] = origin_end - origin_pos
    elif allele_type == 'INS':
        record.stop = record.pos + 1
    else:
        record.stop = record.pos + svlen

    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
CODE

        bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o ~{prefix}.vcf.gz \
            unsorted.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File processed_vcf = "~{prefix}.vcf.gz"
        File processed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task ConvertTsvToParquet {
    input {
        File tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [[ -s "~{tsv}" ]]; then
            python3 <<CODE
import pandas as pd

df = pd.read_csv("~{tsv}", sep="\t", low_memory=False)
df.to_parquet("~{prefix}.parquet", index=False)
CODE
        fi
    >>>

    output {
        File? parquet = "~{prefix}.parquet"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4 * ceil(size(tsv, "GB")) + 8,
        disk_gb: 4 * ceil(size(tsv, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task CreateContigShards {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String contig
        Int shard_bin_size
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        vcfs_file="~{write_lines(vcfs)}"
        vcf_idxs_file="~{write_lines(vcf_idxs)}"
        paste "$vcfs_file" "$vcf_idxs_file" > vcf_pairs.tsv
        max_pos=0
        while IFS=$'\t' read -r vcf vcf_idx; do
            if [[ "$vcf_idx" != "$vcf.tbi" ]]; then
                ln -sf "$vcf_idx" "$vcf.tbi"
            fi
            pos=$(bcftools view -H -r ~{contig} "$vcf" | tail -n1 | cut -f2 || true)
            if [[ -n "$pos" ]] && (( pos > max_pos )); then
                max_pos=$pos
            fi
        done < vcf_pairs.tsv

        python3 - <<CODE
import math

contig = "~{contig}"
shard_bin_size = ~{shard_bin_size}
max_pos = int("$max_pos")

with open("~{prefix}.txt", "w") as out:
    if max_pos == 0:
        pass
    else:
        shard_count = int(math.ceil(max_pos / shard_bin_size))
        for shard_index in range(shard_count):
            start = shard_index * shard_bin_size + 1
            end = min((shard_index + 1) * shard_bin_size, max_pos)
            out.write(f"{contig}:{start}-{end}\n")
CODE
    >>>

    output {
        Array[String] shard_regions = read_lines("~{prefix}.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(2 * size(vcfs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task CreateGapAwareShards {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String contig
        Int shard_bin_size
        Int min_gap = 10000
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        vcfs_file="~{write_lines(vcfs)}"
        vcf_idxs_file="~{write_lines(vcf_idxs)}"
        paste "$vcfs_file" "$vcf_idxs_file" > vcf_pairs.tsv
        : > record_spans.txt
        while IFS=$'\t' read -r vcf vcf_idx; do
            if [[ "$vcf_idx" != "$vcf.tbi" ]]; then
                ln -sf "$vcf_idx" "$vcf.tbi"
            fi
            bcftools query -r ~{contig} -f '%POS\t%REF\n' "$vcf" \
                | awk 'BEGIN{OFS="\t"} {print $1, $1 + length($2) - 1}' \
                >> record_spans.txt
        done < vcf_pairs.tsv

        LC_ALL=C sort -k1,1n record_spans.txt > record_spans.sorted.txt

        python3 - <<CODE
contig = "~{contig}"
shard_bin_size = ~{shard_bin_size}
min_gap = ~{min_gap}

spans = []
with open("record_spans.sorted.txt") as f:
    for line in f:
        start, end = line.split()
        spans.append((int(start), int(end)))

# Record every gap wider than min_gap, the only positions a boundary may fall on because truvari chunks on such gaps
cut_points = []
max_end = spans[0][1] if spans else 0
for index in range(1, len(spans)):
    if spans[index][0] - max_end > min_gap:
        cut_points.append((max_end, spans[index][0]))
    max_end = max(max_end, spans[index][1])

# Grow each shard to at least shard_bin_size before snapping its end forward to the next safe gap
regions = []
start = 1
cut_index = 0
while spans and start <= max_end:
    target = start + shard_bin_size - 1
    if target >= max_end:
        regions.append((start, max_end))
        break
    while cut_index < len(cut_points) and cut_points[cut_index][0] < target:
        cut_index += 1
    if cut_index >= len(cut_points):
        regions.append((start, max_end))
        break
    regions.append((start, cut_points[cut_index][0]))
    start = cut_points[cut_index][1]
    cut_index += 1

with open("~{prefix}.txt", "w") as out:
    if not regions:
        out.write(f"{contig}\n")
    else:
        for region_start, region_end in regions:
            out.write(f"{contig}:{region_start}-{region_end}\n")
CODE
    >>>

    output {
        Array[String] shard_regions = read_lines("~{prefix}.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(2 * size(vcfs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task DropVcfFields {
    input {
        File vcf
        File vcf_idx
        String drop_fields
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            -x ~{drop_fields} \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File dropped_vcf = "~{prefix}.vcf.gz"
        File dropped_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task ExactMatch {
    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        String source_tag
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Intersect once without -n or -C, so a single pass writes the private set 0000 and the shared sets 0002/0003
        bcftools isec \
            -c none \
            -p isec \
            -Oz \
            ~{vcf} \
            ~{truth_snv_indel_vcf}

        mv isec/0003.vcf.gz ~{prefix}.matched_truth.vcf.gz
        tabix -p vcf -f ~{prefix}.matched_truth.vcf.gz

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            isec/0002.vcf.gz \
            > eval_matched.tsv

        bcftools query \
            -f '%ID\t%FILTER\n' \
            ~{prefix}.matched_truth.vcf.gz \
            | awk -F'\t' 'BEGIN{OFS="\t"} {
                n = split($2, parts, ";")
                out = ""
                for (i = 1; i <= n; i++) {
                    if (parts[i] != "." && parts[i] != "PASS") {
                        out = (out == "" ? parts[i] : out "," parts[i])
                    }
                }
                if (out == "") out = "."
                print $1, out
            }' > truth_matched.tsv

        # Fail loudly rather than shift every row, because the paste below pairs the two shared record sets positionally
        eval_count=$(wc -l < eval_matched.tsv)
        truth_count=$(wc -l < truth_matched.tsv)
        if [[ "${eval_count}" -ne "${truth_count}" ]]; then
            echo "Shared record counts differ: ${eval_count} callset vs ${truth_count} truth" >&2
            exit 1
        fi

        paste eval_matched.tsv truth_matched.tsv \
            | awk -v src="~{source_tag}" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,"EXACT",$6,src,$7}' \
            > ~{prefix}.tsv

        mv isec/0000.vcf.gz ~{prefix}.vcf.gz
        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File annotation_tsv = "~{prefix}.tsv"
        File matched_truth_vcf = "~{prefix}.matched_truth.vcf.gz"
        File matched_truth_vcf_idx = "~{prefix}.matched_truth.vcf.gz.tbi"
        File unmatched_vcf = "~{prefix}.vcf.gz"
        File unmatched_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf, "GB") + size(truth_snv_indel_vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task ExtractSample {
    input {
        File vcf
        File vcf_idx
        String sample
        String? extra_args
        Boolean normalize_output = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -s ~{sample} \
            --min-ac 1 \
            ~{if normalize_output then "--trim-alt-alleles" else ""} \
            ~{if defined(extra_args) then extra_args else ""} \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task ExtractVcfAnnotations {
    input {
        File vcf
        File vcf_idx
        File original_vcf
        File original_vcf_idx
        String prefix
        Boolean add_header_row = false
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam
import sys

vcf = pysam.VariantFile("~{vcf}")
orig = pysam.VariantFile("~{original_vcf}")

new_keys = sorted(list(set(vcf.header.info.keys()) - set(orig.header.info.keys())))

with open("~{prefix}.header.txt", "w") as out:
    for k in new_keys:
        out.write(k + "\n")

with open("~{prefix}.annotations.tsv", "w") as out:
    if "~{add_header_row}" == "true":
        fixed_cols = ["#CHROM", "POS", "REF", "ALT", "ID"]
        full_header = fixed_cols + new_keys
        out.write("\t".join(full_header) + "\n")

    for record in vcf:
        alts = ",".join(record.alts) if record.alts else "."
        rid = record.id if record.id else "."
        row = [record.chrom, str(record.pos), record.ref, alts, rid]
        for k in new_keys:
            if k in record.info:
                val = record.info[k]
                if isinstance(val, bool):
                    row.append("1" if val else "0")
                elif isinstance(val, (list, tuple)):
                    row.append(",".join(map(str, val)))
                else:
                    row.append(str(val))
            else:
                row.append(".")

        out.write("\t".join(row) + "\n")
CODE
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
        File annotations_header = "~{prefix}.header.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(original_vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task ExtractVcfCoords {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            ~{vcf} \
            > ~{prefix}.tsv
    >>>

    output {
        File coords_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task FilterDuplicateZeroDepthReferenceBlocks {
    meta {
        description: "Collapse duplicate zero-depth non-alt reference blocks while preserving gVCF coverage. Optionally process supplied ranges and create indexes."
    }

    input {
        File gvcf
        File gvcf_idx
        String prefix
        Array[String] ranges = []
        Boolean remove_duplicates = true
        Boolean create_indexes = false
        Int default_max_retries = 0
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2 * ceil(size(gvcf, "GB"))

    command <<<
        set -euo pipefail

        test -s ~{gvcf_idx}
        mkdir -p per_contig

        filter_records() {
            awk -F$'\t' '
                BEGIN { removed_count = 0 }
                function clear_group(    i) { for (i = 1; i <= row_count; i++) delete rows[i]; for (i in duplicate_count) delete duplicate_count[i]; for (i in removable) delete removable[i]; row_count = 0 }
                function record_end(line,    fields, info_fields, field_count, info_count, i) { field_count = split(line, fields, "\t"); info_count = split(fields[8], info_fields, ";"); for (i = 1; i <= info_count; i++) if (info_fields[i] ~ /^END=[0-9]+$/) return substr(info_fields[i], 5) + 0; return fields[2] + length(fields[4]) - 1 }
                function set_end(line, new_end,    fields, info_fields, field_count, info_count, i, output) { field_count = split(line, fields, "\t"); info_count = split(fields[8], info_fields, ";"); output = ""; for (i = 1; i <= info_count; i++) { if (info_fields[i] ~ /^END=[0-9]+$/) info_fields[i] = "END=" new_end; output = output (i == 1 ? "" : ";") info_fields[i] } fields[8] = output; output = fields[1]; for (i = 2; i <= field_count; i++) output = output "\t" fields[i]; return output }
                function update_coverage(chrom, end) { if (!(chrom in covered_until) || end > covered_until[chrom]) covered_until[chrom] = end }
                function emit_pending(next_chrom, next_pos,    line, end) { if (pending_line == "") return; line = pending_line; end = pending_end; if (pending_chrom == next_chrom && next_pos <= end) { end = next_pos - 1; line = set_end(line, end) } print line; update_coverage(pending_chrom, end); pending_line = "" }
                function set_pending(line, chrom, end) { pending_line = line; pending_chrom = chrom; pending_end = end; removed_count-- }
                function flush_group(    i, line, candidate_line, candidate_end, end, has_kept_row) {
                    candidate_line = ""; candidate_end = -1; has_kept_row = 0
                    for (i = 1; i <= row_count; i++) { line = rows[i]; if (removable[line] && duplicate_count[line] > 1) { if (!(line in seen_candidate)) { seen_candidate[line] = 1; end = record_end(line); if (end > candidate_end) { candidate_line = line; candidate_end = end } } removed_count++ } else has_kept_row = 1 }
                    if (has_kept_row) { emit_pending(current_chrom, current_pos); for (i = 1; i <= row_count; i++) { line = rows[i]; if (!(removable[line] && duplicate_count[line] > 1)) { print line; update_coverage(current_chrom, record_end(line)) } } }
                    else if (candidate_line != "") { if (pending_line != "") { if (pending_chrom != current_chrom || current_pos > pending_end) emit_pending("", 0); else if (candidate_end > pending_end) emit_pending(current_chrom, current_pos); else { for (i in seen_candidate) delete seen_candidate[i]; return } } if (!(current_chrom in covered_until) || current_pos > covered_until[current_chrom]) set_pending(candidate_line, current_chrom, candidate_end) }
                    for (i in seen_candidate) delete seen_candidate[i]
                }
                function is_removable_record(    format_fields, sample_fields, field_count, sample_count, i, gt_index, min_dp_index, gt) { field_count = split($9, format_fields, ":"); gt_index = 0; min_dp_index = 0; for (i = 1; i <= field_count; i++) { if (format_fields[i] == "GT") gt_index = i; if (format_fields[i] == "MIN_DP") min_dp_index = i } if (gt_index == 0 || min_dp_index == 0 || NF < 10) return 0; sample_count = split($10, sample_fields, ":"); if (sample_count < gt_index || sample_count < min_dp_index) return 0; gt = sample_fields[gt_index]; return sample_fields[min_dp_index] == "0" && gt ~ /^(0|\.)([\/|](0|\.))*$/ }
                /^#/ { print; next }
                { coordinate = $1 SUBSEP $2; if (row_count > 0 && coordinate != current_coordinate) { flush_group(); clear_group() } current_coordinate = coordinate; current_chrom = $1; current_pos = $2 + 0; rows[++row_count] = $0; duplicate_count[$0]++; if (is_removable_record()) removable[$0] = 1 }
                END { if (row_count > 0) flush_group(); emit_pending("", 0); print "Removed " removed_count " duplicate zero-depth non-alt gVCF records" > "/dev/stderr" }
            '
        }

        write_records() {
            local range="$1"
            local outfile="$2"
            if [[ "~{remove_duplicates}" == "true" ]]; then
                if [[ -n "$range" ]]; then bcftools view ~{gvcf} "$range"; else bcftools view ~{gvcf}; fi | filter_records | bgzip > "$outfile"
            else
                if [[ -n "$range" ]]; then bcftools view ~{gvcf} "$range"; else bcftools view ~{gvcf}; fi | bgzip > "$outfile"
            fi
            if [[ "~{create_indexes}" == "true" ]]; then tabix -p vcf "$outfile"; fi
        }

        if [[ ~{length(ranges)} -eq 0 ]]; then
            write_records "" "per_contig/~{prefix}.cleaned.g.vcf.gz"
        else
            index=0
            for range in ~{sep=' ' ranges}; do
                pindex=$(printf '%06d' "$index")
                frange=$(echo "$range" | sed 's/[:-]/___/g')
                write_records "$range" "per_contig/$pindex.~{basename(gvcf, ".g.vcf.gz")}.locus_$frange.g.vcf.gz"
                index=$((index + 1))
            done
        fi
    >>>

    output {
        Array[File] cleaned_gvcfs = glob("per_contig/*.g.vcf.gz")
        Array[File] cleaned_gvcf_idxs = glob("per_contig/*.g.vcf.gz.tbi")
    }

    RuntimeAttr default_attr = object { cpu_cores: 1, mem_gb: 1, disk_gb: disk_size, boot_disk_gb: 25, preemptible_tries: 1, max_retries: default_max_retries }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}

# Drop hand-listed (trid, motif) rows from an LPS table, e.g. one trgt-lps emitted from a spurious INFO/MOTIFS value
task FilterLpsTsvRows {
    input {
        File tsv
        Array[Array[String]] filter_trid_motif_pairs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        zcat -f ~{tsv} \
            | awk -F'\t' '
                # Every pair is read before the LPS table, so a pair matching no row can be reported at the end
                NR == FNR {
                    if (NF != 2) {
                        printf "Filter pair on line %d has %d field(s), expected a trid and a motif\n", FNR, NF > "/dev/stderr"
                        malformed = 1
                        exit 1
                    }
                    matched[$1 SUBSEP $2] = 0
                    next
                }

                FNR == 1 {
                    print
                    next
                }

                ($1 SUBSEP $2) in matched {
                    matched[$1 SUBSEP $2] += 1
                    next
                }

                {
                    print
                }

                END {
                    if (malformed) {
                        exit 1
                    }

                    for (key in matched) {
                        split(key, pair, SUBSEP)
                        printf "Filtered %d LPS row(s) for pair (%s, %s)\n", matched[key], pair[1], pair[2] > "/dev/stderr"
                        if (matched[key] == 0) {
                            unmatched = 1
                        }
                    }

                    if (unmatched) {
                        printf "Every filter pair must match at least one LPS row; the LPS table and the filter list disagree\n" > "/dev/stderr"
                        exit 1
                    }
                }
            ' ~{write_tsv(filter_trid_motif_pairs)} - \
            > ~{prefix}.tsv
    >>>

    output {
        File filtered_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10 * ceil(size(tsv, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task FilterTRGTVcf {
    input {
        File vcf
        File vcf_idx
        Int? min_repeat_unit
        Int? min_length_diff
        Int? max_catalog_length
        File? ref_fa
        Boolean normalize = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam


def motif_lengths(record):
    motifs = record.info.get("MOTIFS")
    if motifs is None:
        return []
    if not isinstance(motifs, (list, tuple)):
        motifs = [motifs]

    lengths = []
    for motif_group in motifs:
        for motif in str(motif_group).split(','):
            motif = motif.strip()
            if motif:
                lengths.append(len(motif))
    return lengths


def keep_record(record):
    min_repeat_unit = ~{if defined(min_repeat_unit) then min_repeat_unit else "None"}
    min_length_diff = ~{if defined(min_length_diff) then min_length_diff else "None"}
    max_catalog_length = ~{if defined(max_catalog_length) then max_catalog_length else "None"}

    if max_catalog_length is not None and len(record.ref) > max_catalog_length:
        return False

    if min_repeat_unit is not None:
        lengths = motif_lengths(record)
        if not lengths or min(lengths) < min_repeat_unit:
            return False

    if min_length_diff is not None:
        if not record.alts:
            return False
        length_diffs = [abs(len(record.ref) - len(alt)) for alt in record.alts if alt is not None]
        if not length_diffs or max(length_diffs) < min_length_diff:
            return False

    return True


vcf_in = pysam.VariantFile("~{vcf}")
vcf_out = pysam.VariantFile("preprocessed.vcf.gz", "wz", header=vcf_in.header)

for record in vcf_in:
    a_type = record.info.get("allele_type")
    if (a_type and a_type != "trv") or (keep_record(record)):
        vcf_out.write(record)

vcf_in.close()
vcf_out.close()
CODE

        if [[ "~{normalize}" == "true" ]]; then
            bcftools norm \
                -f ~{select_first([ref_fa, ""])} \
                preprocessed.vcf.gz \
            | awk -F'\t' 'BEGIN {OFS="\t"} /^#/ {print; next} { $4=toupper($4); $5=toupper($5); print }' \
            | bgzip -c > normalized.unsorted.vcf.gz

            bcftools sort \
                --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
                -T . \
                -Oz -o ~{prefix}.vcf.gz \
                normalized.unsorted.vcf.gz
        else
            mv preprocessed.vcf.gz ~{prefix}.vcf.gz
        fi

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File processed_vcf = "~{prefix}.vcf.gz"
        File processed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 4 * ceil(size(vcf, "GB") + size(select_first([ref_fa, vcf]), "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task FindTrios {
    input {
        File vcf
        File vcf_idx
        File? ped
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
import pysam

PED_PATH = "~{ped}"

trios = []
if PED_PATH:
    vcf = pysam.VariantFile("~{vcf}")
    vcf_samples = set(vcf.header.samples)
    vcf.close()

    with open(PED_PATH) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            sample, father, mother = fields[1], fields[2], fields[3]
            if father != "0" and mother != "0":
                if sample in vcf_samples and father in vcf_samples and mother in vcf_samples:
                    trios.append((sample, father, mother))

with open("~{prefix}.trio_definitions.tsv", "w") as out:
    for child, father, mother in trios:
        out.write(f"{child}\t{father}\t{mother}\n")

all_samples = set()
for child, father, mother in trios:
    all_samples.update([child, father, mother])

with open("~{prefix}.trio_sample_ids.txt", "w") as out:
    for sample in sorted(all_samples):
        out.write(sample + "\n")
PYCODE
    >>>

    output {
        File trio_definitions = "~{prefix}.trio_definitions.tsv"
        File trio_sample_ids_file = "~{prefix}.trio_sample_ids.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB") + size(ped, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task GetContigsFromTsv {
    input {
        File tsv
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cut -f1 ~{tsv} | sort -u > contigs.txt
    >>>

    output {
        Array[String] contigs = read_lines("contigs.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsv, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task GetHailMTSize {
    input {
        String mt_uri
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        tot_size=$(gsutil -m du -sh ~{mt_uri} | awk -F '    ' '{ print $1 }')

        python3 <<CODE > mt_size.txt
import sys

size = "$tot_size".split()[0]
unit = "$tot_size".split()[1]

def convert_to_gib(size, unit):
    size_dict = {"KiB": 2**10, "MiB": 2**20, "GB": 2**30, "TiB": 2**40}
    return float(size) * size_dict[unit] / size_dict["GB"]

size_in_gib = convert_to_gib(size, unit)
print(size_in_gib)
CODE
    >>>

    output {
        Float mt_size = read_lines('mt_size.txt')[0]
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 25,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task GetSamplesFromVcf {
    input {
        File vcf
        File vcf_idx
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf} > samples.txt
    >>>

    output {
        Array[String] samples = read_lines("samples.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task IndexVcf {
    input {
        File vcf
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String filename = basename(vcf)

    command <<<
        set -euo pipefail

        cp ~{vcf} ~{filename}

        tabix -p vcf ~{filename}
    >>>

    output {
        File vcf_idx = "~{filename}.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task MakeWindows {
    input {
        File ref_fai
        Array[String] contigs
        Int window_size
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        awk -v contigs="~{sep=',' contigs}" 'BEGIN {split(contigs, c, ","); for(i in c) req[c[i]]=1} {if(req[$1]) print $1 "\t" $2}' ~{ref_fai} > genome.txt

        bedtools makewindows \
            -g genome.txt \
            -w ~{window_size} \
            > windows.bed

        awk '{print $1":"$2"-"$3}' windows.bed > regions.txt
    >>>

    output {
        Array[String] regions = read_lines("regions.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(ref_fai, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task MergeAlignedTsvs {
    input {
        Array[File] tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
from collections import defaultdict

input_files = "~{sep=',' tsvs}".split(',')
header_filename = "~{prefix}.header.txt"
fixed_cols = {"#CHROM", "POS", "REF", "ALT", "ID"}

# Collect annotation keys (columns beyond the 5 fixed cols) in encounter order
all_keys = []
all_keys_seen = set()
file_col_maps = []

for f in input_files:
    with open(f, 'r') as fh:
        line = fh.readline().strip()
        parts = line.split('\t') if line else []
        col_map = {name: i for i, name in enumerate(parts)}
        file_col_maps.append(col_map)
        for k in parts[5:]:
            if k not in all_keys_seen:
                all_keys.append(k)
                all_keys_seen.add(k)

with open(header_filename, 'w') as hout:
    for k in all_keys:
        hout.write(k + "\n")

# Read all data indexed by 5-col variant key, joining across files
variant_order = []
variant_seen = set()
variant_data = defaultdict(dict)

for f, col_map in zip(input_files, file_col_maps):
    if not col_map:
        continue
    with open(f, 'r') as fh:
        # Skip the header row
        fh.readline()
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            key = tuple(parts[:5])
            if key not in variant_seen:
                variant_order.append(key)
                variant_seen.add(key)
            for col_name, col_idx in col_map.items():
                if col_name not in fixed_cols and col_idx < len(parts):
                    variant_data[key][col_name] = parts[col_idx]

with open("aligned_unsorted.tsv", 'w') as out:
    for key in variant_order:
        row = list(key) + [variant_data[key].get(col, '.') for col in all_keys]
        out.write('\t'.join(row) + '\n')
CODE

        sort -k1,1 -k2,2n aligned_unsorted.tsv > ~{prefix}.tsv
    >>>

    output {
        File merged_tsv = "~{prefix}.tsv"
        File merged_header = "~{prefix}.header.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsvs, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task MergeBams {
    input {
        Array[File] bams
        Array[File] bais
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        samtools merge \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            -f \
            -o ~{prefix}.bam \
            ~{sep=' ' bams}

        samtools index \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{prefix}.bam
    >>>

    output {
        File merged_bam = "~{prefix}.bam"
        File merged_bam_idx = "~{prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 4,
        disk_gb: ceil(2.5 * size(bams, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task MergeFastq {
    input {
        Array[File] fastq_gz_files
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{sep=" " fastq_gz_files} > ~{prefix}.merged.fastq.gz
    >>>

    output {
        File merged_fastq_gz = "~{prefix}.merged.fastq.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(fastq_gz_files, "GB")) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task MergeHeaderLines {
    input {
        Array[File] header_files
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{sep=' ' header_files} \
            | sort -u > ~{prefix}.txt
    >>>

    output {
        File merged_header = "~{prefix}.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(header_files, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task MergePALMEROutputs {
    input {
        Array[File] calls_shards
        Array[File] tsd_reads_shards
        String mei_type
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        head -n1 ~{calls_shards[0]} > ~{prefix}_~{mei_type}_calls.txt
        for f in ~{sep=' ' calls_shards}; do
            grep -v '^cluster_id' $f >> ~{prefix}_~{mei_type}_calls.txt || true
        done

        head -n1 ~{tsd_reads_shards[0]} > ~{prefix}_~{mei_type}_tsd_reads.txt
        for f in ~{sep=' ' tsd_reads_shards}; do
            grep -v '^cluster_id' $f >> ~{prefix}_~{mei_type}_tsd_reads.txt || true
        done
    >>>

    output {
        File calls = "~{prefix}_~{mei_type}_calls.txt"
        File tsd_reads = "~{prefix}_~{mei_type}_tsd_reads.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(calls_shards, "GB") + size(tsd_reads_shards, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task MergeVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String? contig
        String? extra_args
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools merge \
            -Oz -o ~{prefix}.vcf.gz \
            ~{if defined(contig) then "-r " + contig else ""} \
            ~{if defined(extra_args) then extra_args else ""} \
            -l ~{write_lines(vcfs)}

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task NormalizeTRGTHaploidGenotypes {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import pysam

with pysam.VariantFile("~{vcf}") as vcf_in, pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=vcf_in.header) as vcf_out:
    for record in vcf_in:
        for sample in record.samples.values():
            genotype = sample.get("GT")
            if genotype is None or len(genotype) != 2:
                continue

            first_allele, second_allele = genotype
            if first_allele is None:
                sample["GT"] = (second_allele,)
            elif second_allele is None:
                sample["GT"] = (first_allele,)

        vcf_out.write(record)
CODE

        tabix -f -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File normalized_vcf = "~{prefix}.vcf.gz"
        File normalized_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size([vcf, vcf_idx], "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task NormalizeVcf {
    input {
        File vcf
        File vcf_idx
        File ref_fa
        File ref_fai
        String? check_ref
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools norm \
            -m -any \
            -f ~{ref_fa} \
            ~{if defined(check_ref) then "-c " + check_ref else ""} \
            -Oz -o unsorted.vcf.gz \
            ~{vcf}

        bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o ~{prefix}.vcf.gz \
            unsorted.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File normalized_vcf = "~{prefix}.vcf.gz"
        File normalized_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10 * ceil(size(vcf, "GB") + size(ref_fa, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task RenameVariantIds {
    input {
        File vcf
        File vcf_idx
        String prefix
        String id_format = "%CHROM-%POS-%REF-%ALT"
        Boolean strip_chr = false
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            --set-id '~{id_format}' \
            -Oz -o temp_renamed.vcf.gz \
            ~{vcf}

        if [ "~{strip_chr}" == "true" ]; then
            bcftools view -h temp_renamed.vcf.gz > header.txt
            bcftools view -H temp_renamed.vcf.gz \
                | awk 'BEGIN{OFS="\t"} {gsub(/^chr/, "", $3); print}' \
                | cat header.txt - \
                | bgzip -c > ~{prefix}.vcf.gz
            rm temp_renamed.vcf.gz header.txt
        else
            mv temp_renamed.vcf.gz ~{prefix}.vcf.gz
        fi

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 8 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task ResetVcfFilters {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            -x FILTER \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File reset_vcf = "~{prefix}.vcf.gz"
        File reset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task RestoreOriginalAlleles {
    input {
        File tsv
        File coords_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
orig = {}
with open('~{coords_tsv}', 'r') as f:
    for line in f:
        fields = line.rstrip('\n').split('\t')
        vid = fields[4]
        if vid != '.':
            orig[vid] = fields

with open('~{tsv}', 'r') as f, open('~{prefix}.tsv', 'w') as out:
    for line in f:
        fields = line.rstrip('\n').split('\t')
        vid = fields[4]
        if vid in orig:
            fields[0] = orig[vid][0]
            fields[1] = orig[vid][1]
            fields[2] = orig[vid][2]
            fields[3] = orig[vid][3]
        out.write('\t'.join(fields) + '\n')
CODE
    >>>

    output {
        File restored_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsv, "GB") + size(coords_tsv, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task RevertSymbolicAlleles {
    input {
        File annotated_vcf
        File annotated_vcf_idx
        File original_vcf
        File original_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

# Index original VCF by variant ID
original_data = {}
with pysam.VariantFile("~{original_vcf}") as orig_vcf:
    has_svtype = 'SVTYPE' in orig_vcf.header.info
    has_svlen = 'SVLEN' in orig_vcf.header.info
    has_end = 'END' in orig_vcf.header.info
    for record in orig_vcf:
        original_data[record.id] = {
            'chrom': record.chrom,
            'pos': record.pos,
            'ref': record.ref,
            'alts': record.alts,
            'svtype': record.info.get('SVTYPE', None) if has_svtype else None,
            'svlen': record.info.get('SVLEN', None) if has_svlen else None,
            'end': record.info.get('END', None) if has_end else None,
        }

# Revert symbolic alleles in annotated VCF
vcf_in = pysam.VariantFile("~{annotated_vcf}")
vcf_out = pysam.VariantFile("unsorted.vcf.gz", 'w', header=vcf_in.header)
for record in vcf_in:
    if record.id in original_data:
        orig = original_data[record.id]
        record.chrom = orig['chrom']
        record.pos = orig['pos']
        record.ref = orig['ref']
        record.alts = orig['alts']

        if orig['svtype'] is not None:
            record.info['SVTYPE'] = orig['svtype']
        elif 'SVTYPE' in record.info:
            del record.info['SVTYPE']

        if orig['svlen'] is not None:
            record.info['SVLEN'] = orig['svlen']
        elif 'SVLEN' in record.info:
            del record.info['SVLEN']

        if orig['end'] is not None:
            record.info['END'] = orig['end']
        elif 'END' in record.info:
            del record.info['END']

    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
CODE

        bcftools annotate \
            -x INFO/ORIGINAL_POS,INFO/ORIGINAL_CHROM \
            -Ou unsorted.vcf.gz \
        | bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File reverted_vcf = "~{prefix}.vcf.gz"
        File reverted_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(annotated_vcf, "GB") + size(original_vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task RunPALMERShard {
    input {
        File bam
        File bai
        String mode
        String mei_type
        File ref_fa
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        dir=$(pwd)

        mv ~{bam} ./
        mv ~{bai} ./
        bam_base=$(basename ~{bam})
        chrom=$(echo $bam_base | sed 's/\.bam$//' | rev | cut -f1 -d '_' | rev)

        mkdir -p "${chrom}"
        /PALMER/PALMER \
            --input ${bam_base} \
            --ref_fa ~{ref_fa} \
            --ref_ver GRCh38 \
            --type ~{mei_type} \
            --mode ~{mode} \
            --output "~{prefix}" \
            --chr $chrom \
            --workdir "${dir}/${chrom}/"

        sed -i "s/$/\t~{mei_type}/" ${chrom}/~{prefix}_calls.txt
        sed -i "s/$/\t~{mei_type}/" ${chrom}/~{prefix}_TSD_reads.txt
        mv ${chrom}/~{prefix}_calls.txt ~{prefix}_calls_shard.txt
        mv ${chrom}/~{prefix}_TSD_reads.txt ~{prefix}_tsd_reads_shard.txt
    >>>

    output {
        File calls_shard = "~{prefix}_calls_shard.txt"
        File tsd_reads_shard = "~{prefix}_tsd_reads_shard.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task SetMissingFiltersToPass {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} \
            | awk 'BEGIN{OFS="\t"} /^#/ {print; next} $7=="." {$7="PASS"} {print}' \
            | bgzip -c > ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task ShardVcfByRecords {
    input {
        File vcf
        File vcf_idx
        Int records_per_shard
        Boolean use_ssd = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir scatter_output

        bcftools +scatter \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            -n ~{records_per_shard} \
            --prefix ~{prefix}. \
            -Oz -o scatter_output \
            ~{vcf}

        mkdir shards

        find scatter_output -maxdepth 1 -name "*.vcf.gz" | sort -k1,1V > vcfs.list

        i=0
        while read VCF; do
            if [[ -z "$VCF" ]]; then continue; fi
            shard_no=$(printf %06d $i)
            mv "$VCF" "shards/shard_${shard_no}.vcf.gz"
            tabix -p vcf "shards/shard_${shard_no}.vcf.gz"
            i=$((i+1))
        done < vcfs.list
    >>>

    output {
        Array[File] shards = glob("shards/shard_*.vcf.gz")
        Array[File] shard_idxs = glob("shards/shard_*.vcf.gz.tbi")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 25,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SortReadCounts {
    input {
        File read_counts
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        export LC_ALL=C

        zcat ~{read_counts} > counts.tsv
        grep "^@" counts.tsv > ~{prefix}.tsv
        grep -v -m 1 "^@" counts.tsv >> ~{prefix}.tsv
        grep "^@SQ" counts.tsv | sed -E 's/^@SQ\tSN:([^\t]+).*/\1/' > contig_order.txt
        grep -v "^@" counts.tsv | tail -n +2 \
            | awk -F'\t' 'NR==FNR{rank[$1]=NR; next} {print rank[$1]"\t"$0}' contig_order.txt - \
            | sort -s -k1,1n \
            | cut -f2- >> ~{prefix}.tsv

        bgzip ~{prefix}.tsv
    >>>

    output {
        File sorted_read_counts = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(read_counts, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task SplitBam {
    input {
        File bam
        File bai
        String prefix
        Array[String] contigs
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        for contig in ~{sep=" " contigs}
        do
            samtools view -b ~{bam} $contig > ~{prefix}_${contig}.bam
            samtools index ~{prefix}_${contig}.bam
        done
    >>>

    output {
        Array[File] bams = glob("~{prefix}_*bam")
        Array[File] bais = glob("~{prefix}_*bai")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(bam, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task SplitBamByContig {
    input {
        File bam
        File bai
        Array[String] contigs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        for contig in ~{sep=' ' contigs}
        do
            samtools view \
                -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
                -h \
                -b \
                -o ~{prefix}.$contig.bam \
                ~{bam} \
                $contig

            samtools index \
                -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
                ~{prefix}.$contig.bam
        done
    >>>

    output {
        Array[File] contig_bams = glob("*.bam")
        Array[File] contig_bais = glob("*.bam.bai")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 2,
        disk_gb: 2 * ceil(size(bam, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task StripGenotypes {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -G \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File stripped_vcf = "~{prefix}.vcf.gz"
        File stripped_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task StripInfoFields {
    input {
        File vcf
        File vcf_idx
        Array[String] info_fields
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            -x INFO/~{sep=",INFO/" info_fields} \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File stripped_vcf = "~{prefix}.vcf.gz"
        File stripped_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task SubsetBamToContig {
    input {
        File bam
        File bai
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: { localization_optional: true }
        bai: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            -h \
            -b \
            -o ~{prefix}.bam \
            ~{bam} \
            ~{contig}

        samtools index \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{prefix}.bam
    >>>

    output {
        File contig_bam = "~{prefix}.bam"
        File contig_bai = "~{prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 2,
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task SubsetBamToRegions {
    meta {
        description: "Subset a BAM to the regions listed in a file."
    }

    input {
        File bam
        File bai
        File region_file
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: { localization_optional: true }
        bai: { localization_optional: true }
    }

    Array[String] regions = read_lines(region_file)

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            -X \
            -M \
            -h \
            -b \
            -o ~{prefix}.bam \
            ~{bam} \
            ~{bai} \
            ~{sep=' ' regions}

        samtools index \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{prefix}.bam
    >>>

    output {
        File subset_bam = "~{prefix}.bam"
        File subset_bai = "~{prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 2,
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task SubsetBedToContig {
    input {
        File bed
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if gzip -t "~{bed}" 2>/dev/null; then
            gzip -cd "~{bed}" | awk -v contig="~{contig}" '$1 == contig' > "~{prefix}.bed"
        else
            awk -v contig="~{contig}" '$1 == contig' "~{bed}" > "~{prefix}.bed"
        fi
    >>>

    output {
        File subset_bed = "~{prefix}.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 2 * ceil(size(bed, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task SubsetTsvToContig {
    input {
        File tsv
        String contig
        Boolean sort_output = false
        Boolean compressed_tsv = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ "~{compressed_tsv}" == "true" ]; then
            gzip -dc ~{tsv} | awk -v contig="~{contig}" '$1 == contig' > ~{prefix}.tsv
        else
            awk -v contig="~{contig}" '$1 == contig' ~{tsv} > ~{prefix}.tsv
        fi

        if [ "~{sort_output}" == "true" ]; then
            sort -k1,1 -k2,2n ~{prefix}.tsv -o ~{prefix}.tsv
        fi
    >>>

    output {
        File subset_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsv, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task SubsetVcfByArgs {
    input {
        File vcf
        File vcf_idx
        String? include_args
        String? exclude_args
        String? extra_args
        Boolean use_ssd = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} \
            ~{if defined(include_args) then "-i '~{include_args}'" else ""} \
            ~{if defined(exclude_args) then "-e '~{exclude_args}'" else ""} \
            ~{if defined(extra_args) then extra_args else ""} \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf "~{prefix}.vcf.gz"
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size([vcf, vcf_idx], "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SubsetVcfByLength {
    input {
        File vcf
        File vcf_idx
        String length_field = "allele_length"
        Int? min_length
        Int? max_length
        String? extra_args
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String size_filter = if defined(min_length) && defined(max_length) then 'abs(INFO/~{length_field})>=~{min_length} && abs(INFO/~{length_field})<=~{max_length}' else if defined(min_length) then 'abs(INFO/~{length_field})>=~{min_length}' else if defined(max_length) then 'abs(INFO/~{length_field})<=~{max_length}' else '1==1'

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} \
            --include "~{size_filter}" \
            ~{if defined(extra_args) then extra_args else ""} \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf "~{prefix}.vcf.gz"
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size([vcf, vcf_idx], "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task SubsetVcfToContig {
    input {
        File vcf
        File? vcf_idx
        String contig
        String? extra_args
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if ~{!defined(vcf_idx)}; then
            tabix -p vcf ~{vcf}
        elif [[ "~{vcf_idx}" != "~{vcf}.tbi" ]]; then
            ln -sf "~{vcf_idx}" "~{vcf}.tbi"
        fi

        bcftools view \
            -r ~{contig} \
            ~{if defined(extra_args) then extra_args else ""} \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz
        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task SubsetVcfToRegion {
    input {
        File vcf
        File vcf_idx
        String region
        Boolean use_ssd = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [[ "~{vcf_idx}" != "~{vcf}.tbi" ]]; then
            ln -sf "~{vcf_idx}" "~{vcf}.tbi"
        fi

        bcftools view \
            -r ~{region} \
            -t ~{region} \
            --threads $(nproc) \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SubsetVcfToRegionStreaming {
    input {
        File vcf
        File vcf_idx
        String region
        String? include_args
        Boolean drop_genotypes = false
        Boolean use_ssd = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        vcf: { localization_optional: true }
        vcf_idx: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # Pair -r with -t so records are selected by POS alone, because -r alone also returns REF spans reaching in
        bcftools view \
            -r ~{region} \
            -t ~{region} \
            ~{if defined(include_args) then "-i '~{include_args}'" else ""} \
            ~{if drop_genotypes then "-G" else ""} \
            --threads $(nproc) \
            ~{vcf} \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SubsetVcfToSamples {
    input {
        File vcf
        File vcf_idx
        Array[String] samples
        Boolean keep_samples = true
        Boolean filter_to_sample = true
        String? extra_args
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat > samples.txt <<EOF
~{sep='\n' samples}
EOF

        bcftools annotate \
            -x INFO/AC,INFO/AN,INFO/AF \
            ~{vcf} \
        | bcftools view \
            --samples-file ~{if keep_samples then "samples.txt" else "^samples.txt"} \
            ~{if filter_to_sample then "--min-ac 1" else ""} \
            ~{if defined(extra_args) then extra_args else ""} \
            -Oz -o ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task SwapSampleIds {
    input {
        File vcf
        File vcf_idx
        File sample_swap_list
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf} > current_samples.txt

        awk 'FNR==NR {swap[$1]=$2; next} {if ($1 in swap) print swap[$1]; else print $1}' \
            ~{sample_swap_list} current_samples.txt > new_samples.txt

        bcftools reheader \
            --samples new_samples.txt \
            ~{vcf} \
        > ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File swapped_vcf = "~{prefix}.vcf.gz"
        File swapped_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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

task TransferAWSToGCS {
    input {
        String aws_path
        String output_gcs_path
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        filename=$(basename "~{aws_path}")

        aws s3 cp --no-sign-request "~{aws_path}" "$filename"

        gsutil cp "$filename" "~{output_gcs_path}"
    >>>

    output {
        String gcs_path = output_gcs_path
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 600,
        boot_disk_gb: 10,
        preemptible_tries: 1,
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
