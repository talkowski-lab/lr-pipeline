version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SplitVcfPerContig {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Boolean create_no_geno = false
        Boolean modify_snv_ids = false
        Boolean rename_dbsnp_contigs = false
        Boolean rename_dbvar_contigs = false
        Array[String]? missing_info_header_fields

        String utils_docker

        RuntimeAttr? runtime_attr_add_header_lines
        RuntimeAttr? runtime_attr_split_vcf
    }

    if (defined(missing_info_header_fields)) {
        call Helpers.AddMissingInfoHeaderLines {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                info_fields = select_first([missing_info_header_fields]),
                prefix = "~{prefix}.hdr_fixed",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_add_header_lines
        }
    }

    File base_vcf = select_first([AddMissingInfoHeaderLines.annotated_vcf, vcf])
    File base_vcf_idx = select_first([AddMissingInfoHeaderLines.annotated_vcf_idx, vcf_idx])

    call SplitByContig {
        input:
            vcf = base_vcf,
            vcf_idx = base_vcf_idx,
            contigs = contigs,
            prefix = prefix,
            create_no_geno = create_no_geno,
            modify_snv_ids = modify_snv_ids,
            rename_dbsnp_contigs = rename_dbsnp_contigs,
            rename_dbvar_contigs = rename_dbvar_contigs,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_split_vcf
    }

    output {
        Array[File] contig_vcfs = SplitByContig.contig_vcfs
        Array[File] contig_vcf_idxs = SplitByContig.contig_vcf_idxs
        Array[File] contig_no_geno_vcfs = SplitByContig.contig_no_geno_vcfs
        Array[File] contig_no_geno_vcf_idxs = SplitByContig.contig_no_geno_vcf_idxs
    }
}

task SplitByContig {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix
        Boolean create_no_geno
        Boolean modify_snv_ids
        Boolean rename_dbsnp_contigs
        Boolean rename_dbvar_contigs
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if ~{rename_dbsnp_contigs}; then
            cat > chr_name_mapping.txt <<EOF
NC_000001.11 chr1
NC_000002.12 chr2
NC_000003.12 chr3
NC_000004.12 chr4
NC_000005.10 chr5
NC_000006.12 chr6
NC_000007.14 chr7
NC_000008.11 chr8
NC_000009.12 chr9
NC_000010.11 chr10
NC_000011.10 chr11
NC_000012.12 chr12
NC_000013.11 chr13
NC_000014.9 chr14
NC_000015.10 chr15
NC_000016.10 chr16
NC_000017.11 chr17
NC_000018.10 chr18
NC_000019.10 chr19
NC_000020.11 chr20
NC_000021.9 chr21
NC_000022.11 chr22
NC_000023.11 chrX
NC_000024.10 chrY
EOF

            bcftools annotate \
                --rename-chrs chr_name_mapping.txt \
                -Oz -o ~{prefix}.renamed.vcf.gz \
                ~{vcf}

            tabix -p vcf ~{prefix}.renamed.vcf.gz

            input_vcf="~{prefix}.renamed.vcf.gz"
        elif ~{rename_dbvar_contigs}; then
            cat > chr_name_mapping.txt <<EOF
1 chr1
2 chr2
3 chr3
4 chr4
5 chr5
6 chr6
7 chr7
8 chr8
9 chr9
10 chr10
11 chr11
12 chr12
13 chr13
14 chr14
15 chr15
16 chr16
17 chr17
18 chr18
19 chr19
20 chr20
21 chr21
22 chr22
X chrX
Y chrY
MT chrM
EOF

            bcftools annotate \
                --rename-chrs chr_name_mapping.txt \
                -Oz -o ~{prefix}.renamed.vcf.gz \
                ~{vcf}

            tabix -p vcf ~{prefix}.renamed.vcf.gz

            input_vcf="~{prefix}.renamed.vcf.gz"
        else
            input_vcf="~{vcf}"
        fi

        while IFS= read -r contig; do
            if ~{modify_snv_ids}; then
                bcftools view \
                    -r "${contig}" \
                    "${input_vcf}" \
                | awk 'BEGIN{OFS="\t"} /^#/ {print; next} $8 ~ /(^|;)allele_type=snv(;|$)/ {$3=$1"-"$2"-"$4"-"$5} {print}' \
                | bgzip -c > "~{prefix}.${contig}.full.vcf.gz"
            else
                bcftools view \
                    -r "${contig}" \
                    "${input_vcf}" \
                    -Oz -o "~{prefix}.${contig}.full.vcf.gz"
            fi

            tabix -p vcf "~{prefix}.${contig}.full.vcf.gz"

            if ~{create_no_geno}; then
                bcftools view \
                    -G \
                    "~{prefix}.${contig}.full.vcf.gz" \
                    -Oz -o "~{prefix}.${contig}.no_geno.vcf.gz"
                
                tabix -p vcf "~{prefix}.${contig}.no_geno.vcf.gz"
            fi
        done < ~{write_lines(contigs)}
    >>>

    output {
        Array[File] contig_vcfs = glob("~{prefix}.*.full.vcf.gz")
        Array[File] contig_vcf_idxs = glob("~{prefix}.*.full.vcf.gz.tbi")
        Array[File] contig_no_geno_vcfs = glob("~{prefix}.*.no_geno.vcf.gz")
        Array[File] contig_no_geno_vcf_idxs = glob("~{prefix}.*.no_geno.vcf.gz.tbi")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 10,
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
