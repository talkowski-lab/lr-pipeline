version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateSmallInsertionOrigins {
    input {
        Array[String] contigs
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix

        Int min_length = 20
        Int max_length = 49

        File ref_fa
        File ref_fai
        File t2t_fa
        File tr_catalog_bed
        File t2t_genes

        String insertion_annotation_docker
        String blat_docker
        Int blat_threads = 4

        RuntimeAttr? runtime_attr_local
        RuntimeAttr? runtime_attr_blat
        RuntimeAttr? runtime_attr_final
        RuntimeAttr? runtime_attr_concat
    }

    scatter (i in range(length(contigs))) {
        call SmallInsertionsLocal {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                min_length = min_length,
                max_length = max_length,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                tr_catalog_bed = tr_catalog_bed,
                prefix = "~{prefix}.~{contigs[i]}",
                docker = insertion_annotation_docker,
                runtime_attr_override = runtime_attr_local
        }

        call BlatT2T {
            input:
                queries_fa = SmallInsertionsLocal.t2t_queries,
                t2t_fa = t2t_fa,
                threads = blat_threads,
                prefix = "~{prefix}.~{contigs[i]}",
                docker = blat_docker,
                runtime_attr_override = runtime_attr_blat
        }

        call SmallInsertionsFinal {
            input:
                remap_like_tsv = SmallInsertionsLocal.remap_like_tsv,
                tr_stage_tsv = SmallInsertionsLocal.tr_stage_tsv,
                local_dups_tsv = SmallInsertionsLocal.local_dups_tsv,
                no_local_copy_ids = SmallInsertionsLocal.no_local_copy_ids,
                t2t_psl = BlatT2T.psl,
                t2t_genes = t2t_genes,
                prefix = "~{prefix}.~{contigs[i]}",
                docker = insertion_annotation_docker,
                runtime_attr_override = runtime_attr_final
        }
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = SmallInsertionsFinal.final_annotation,
            sort_output = false,
            preserve_header = true,
            prefix = "~{prefix}.small_insertion_origins",
            docker = insertion_annotation_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File small_insertion_origins_tsv = ConcatTsvs.concatenated_tsv
    }
}

task SmallInsertionsLocal {
    input {
        File vcf
        File vcf_idx
        Int min_length
        Int max_length
        File ref_fa
        File ref_fai
        File tr_catalog_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        S=/opt/insertion_annotation
        bcftools query -i 'INFO/allele_type="ins" && INFO/allele_length>=~{min_length} && INFO/allele_length<=~{max_length}' \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ~{vcf} \
            | awk -F'\t' 'BEGIN{OFS="\t"} length($3)==1 && substr($4,1,1)==$3 {print $1,$2,$3,$4,$5,"not remapped",".",".","."}' \
            > ~{prefix}.remap_like.tsv
        (printf 'ID\tins_call\tremap_class\n'; awk -F'\t' 'BEGIN{OFS="\t"} {print $5,"","not remapped"}' ~{prefix}.remap_like.tsv) \
            | gzip > annotations.tsv.gz
        printf 'ID\tscore\n' | gzip > hit_repeats.empty.tsv.gz
        python3 $S/classify_tr_insertions.py --annotations annotations.tsv.gz --remap ~{prefix}.remap_like.tsv \
            --tr-catalog ~{tr_catalog_bed} --ref-fa ~{ref_fa} --hit-repeats hit_repeats.empty.tsv.gz --out ~{prefix}.tr_stage.tsv
        awk -F'\t' 'NR>1 && ($5=="not TR" || $5=="in TR, non-repeat sequence"){print $1}' ~{prefix}.tr_stage.tsv > nonTR.ids.txt
        python3 $S/small_insertion_local_dups.py --remap ~{prefix}.remap_like.tsv --ids nonTR.ids.txt --ref-fa ~{ref_fa} \
            --out ~{prefix}.local_dups.tsv
        awk -F'\t' 'NR>1 && $12=="no local copy"{print $1}' ~{prefix}.local_dups.tsv > ~{prefix}.no_local_copy.ids.txt

        python3 - ~{prefix}.remap_like.tsv ~{prefix}.no_local_copy.ids.txt ~{ref_fa} ~{prefix}.t2t_queries.fa <<'EOF'
import sys

import pysam

ids = {x.strip() for x in open(sys.argv[2]) if x.strip()}
fa = pysam.FastaFile(sys.argv[3])
with open(sys.argv[4], "w") as o:
    for line in open(sys.argv[1]):
        c = line.rstrip("\n").split("\t")
        if c[4] not in ids:
            continue
        pos = int(c[1])
        o.write(f">{c[4]}|ins\n{c[3][1:]}\n")
        o.write(f">{c[4]}|L\n{fa.fetch(c[0], max(0, pos - 200), pos).upper()}\n")
        o.write(f">{c[4]}|R\n{fa.fetch(c[0], pos, pos + 200).upper()}\n")
EOF
    >>>

    output {
        File remap_like_tsv = "~{prefix}.remap_like.tsv"
        File tr_stage_tsv = "~{prefix}.tr_stage.tsv"
        File local_dups_tsv = "~{prefix}.local_dups.tsv"
        File no_local_copy_ids = "~{prefix}.no_local_copy.ids.txt"
        File t2t_queries = "~{prefix}.t2t_queries.fa"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(ref_fa, "GB")) + 2 * ceil(size(vcf, "GB")) + 10,
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

task BlatT2T {
    input {
        File queries_fa
        File t2t_fa
        Int threads
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ ! -s ~{queries_fa} ]; then
            : > ~{prefix}.t2t.psl
            exit 0
        fi
        faToTwoBit ~{t2t_fa} t2t.2bit
        mkdir -p shards
        python3 - ~{queries_fa} ~{threads} <<'EOF'
import sys

n = int(sys.argv[2])
outs = [open(f"shards/q{i}.fa", "w") for i in range(n)]
k = -1
for line in open(sys.argv[1]):
    if line.startswith(">") and line.rstrip().endswith("|ins"):
        k += 1
    outs[k % n].write(line)
for o in outs:
    o.close()
EOF
        ls shards/q*.fa | xargs -P ~{threads} -I {} sh -c \
            'blat t2t.2bit {} -stepSize=5 -repMatch=2253 -minScore=15 -minIdentity=80 -noHead {}.psl'
        cat shards/q*.fa.psl > ~{prefix}.t2t.psl
    >>>

    output {
        File psl = "~{prefix}.t2t.psl"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: threads,
        mem_gb: 8 * threads,
        disk_gb: 3 * ceil(size(t2t_fa, "GB")) + 10,
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

task SmallInsertionsFinal {
    input {
        File remap_like_tsv
        File tr_stage_tsv
        File local_dups_tsv
        File no_local_copy_ids
        File t2t_psl
        File t2t_genes
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        S=/opt/insertion_annotation
        python3 $S/small_insertion_t2t_check.py --psl ~{t2t_psl} --ids ~{no_local_copy_ids} --out ~{prefix}.no_local_copy.t2t.tsv
        python3 $S/build_small_insertion_annotation.py --tr-stage ~{tr_stage_tsv} --remap ~{remap_like_tsv} \
            --local-dups ~{local_dups_tsv} --t2t ~{prefix}.no_local_copy.t2t.tsv --genes-table ~{t2t_genes} \
            --out ~{prefix}.small_insertion_origins.tsv
    >>>

    output {
        File final_annotation = "~{prefix}.small_insertion_origins.tsv"
        File t2t_check_tsv = "~{prefix}.no_local_copy.t2t.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([remap_like_tsv, tr_stage_tsv, local_dups_tsv, t2t_psl], "GB")) + 10,
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
