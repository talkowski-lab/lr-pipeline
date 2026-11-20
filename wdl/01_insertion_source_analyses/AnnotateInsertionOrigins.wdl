version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateInsertionOrigins {
    input {
        Array[String] contigs
        Array[File] remap_tsvs
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[File] svan_tsvs
        Array[File] svan_headers
        String prefix

        Int min_length = 50

        File ref_fa
        File ref_fai
        Array[File] ref_bwa_idx
        File ref_asm20_mmi
        File ref_sr_mmi
        File t2t_fa
        File repeats_bed
        File tr_catalog_bed
        File immune_loci_bed
        File known_cnv_loci_bed
        File cytoband
        File sd_gd_blocks_bed
        File sd_loci_bed
        File sd_copies_bed
        File genomic_superdups
        File t2t_genes

        String insertion_annotation_docker
        String repeatmasker_docker
        Int repeatmasker_threads = 8

        RuntimeAttr? runtime_attr_stage
        RuntimeAttr? runtime_attr_align_hg38
        RuntimeAttr? runtime_attr_align_t2t
        RuntimeAttr? runtime_attr_repeatmasker
        RuntimeAttr? runtime_attr_bwa
        RuntimeAttr? runtime_attr_iterative_hg38
        RuntimeAttr? runtime_attr_iterative_t2t
        RuntimeAttr? runtime_attr_concat
    }

    scatter (i in range(length(contigs))) {
        call RunStage as Prepare {
            input:
                stage = "prepare",
                contig = contigs[i],
                remap_tsv = remap_tsvs[i],
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                svan_tsv = svan_tsvs[i],
                svan_header = svan_headers[i],
                min_length = min_length,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                repeats_bed = repeats_bed,
                tr_catalog_bed = tr_catalog_bed,
                immune_loci_bed = immune_loci_bed,
                known_cnv_loci_bed = known_cnv_loci_bed,
                cytoband = cytoband,
                sd_gd_blocks_bed = sd_gd_blocks_bed,
                sd_loci_bed = sd_loci_bed,
                sd_copies_bed = sd_copies_bed,
                genomic_superdups = genomic_superdups,
                t2t_genes = t2t_genes,
                prefix = "~{prefix}.~{contigs[i]}.prepare",
                docker = insertion_annotation_docker,
                runtime_attr_override = runtime_attr_stage
        }

        call AlignMinimap2 as AlignHg38 {
            input:
                queries_fa = select_first([Prepare.next_queries]),
                ref_fa = ref_fa,
                prefix = "~{prefix}.~{contigs[i]}.hg38",
                docker = insertion_annotation_docker,
                runtime_attr_override = runtime_attr_align_hg38
        }

        call RepeatMasker {
            input:
                fa = select_first([Prepare.inserted_sequences]),
                threads = repeatmasker_threads,
                prefix = "~{prefix}.~{contigs[i]}",
                docker = repeatmasker_docker,
                runtime_attr_override = runtime_attr_repeatmasker
        }

        call RunStage as Stage1 {
            input:
                stage = "stage1",
                contig = contigs[i],
                workdir_tar = Prepare.workdir_tar,
                hg38_paf = AlignHg38.paf,
                rm_out = RepeatMasker.rm_out,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                repeats_bed = repeats_bed,
                tr_catalog_bed = tr_catalog_bed,
                immune_loci_bed = immune_loci_bed,
                known_cnv_loci_bed = known_cnv_loci_bed,
                cytoband = cytoband,
                sd_gd_blocks_bed = sd_gd_blocks_bed,
                sd_loci_bed = sd_loci_bed,
                sd_copies_bed = sd_copies_bed,
                genomic_superdups = genomic_superdups,
                t2t_genes = t2t_genes,
                prefix = "~{prefix}.~{contigs[i]}.stage1",
                docker = insertion_annotation_docker,
                runtime_attr_override = runtime_attr_stage
        }

        call AlignMinimap2 as AlignT2T {
            input:
                queries_fa = select_first([Stage1.next_queries]),
                ref_fa = t2t_fa,
                prefix = "~{prefix}.~{contigs[i]}.t2t",
                docker = insertion_annotation_docker,
                runtime_attr_override = runtime_attr_align_t2t
        }

        call RunStage as Stage2 {
            input:
                stage = "stage2",
                contig = contigs[i],
                workdir_tar = Stage1.workdir_tar,
                t2t_paf = AlignT2T.paf,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                repeats_bed = repeats_bed,
                tr_catalog_bed = tr_catalog_bed,
                immune_loci_bed = immune_loci_bed,
                known_cnv_loci_bed = known_cnv_loci_bed,
                cytoband = cytoband,
                sd_gd_blocks_bed = sd_gd_blocks_bed,
                sd_loci_bed = sd_loci_bed,
                sd_copies_bed = sd_copies_bed,
                genomic_superdups = genomic_superdups,
                t2t_genes = t2t_genes,
                prefix = "~{prefix}.~{contigs[i]}.stage2",
                docker = insertion_annotation_docker,
                runtime_attr_override = runtime_attr_stage
        }

        call BwaMem {
            input:
                queries_fa = select_first([Stage2.next_queries]),
                ref_fa = ref_fa,
                ref_bwa_idx = ref_bwa_idx,
                prefix = "~{prefix}.~{contigs[i]}",
                docker = insertion_annotation_docker,
                runtime_attr_override = runtime_attr_bwa
        }

        call RunStage as Stage3 {
            input:
                stage = "stage3",
                contig = contigs[i],
                workdir_tar = Stage2.workdir_tar,
                bwa_sam = BwaMem.sam,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                repeats_bed = repeats_bed,
                tr_catalog_bed = tr_catalog_bed,
                immune_loci_bed = immune_loci_bed,
                known_cnv_loci_bed = known_cnv_loci_bed,
                cytoband = cytoband,
                sd_gd_blocks_bed = sd_gd_blocks_bed,
                sd_loci_bed = sd_loci_bed,
                sd_copies_bed = sd_copies_bed,
                genomic_superdups = genomic_superdups,
                t2t_genes = t2t_genes,
                prefix = "~{prefix}.~{contigs[i]}.stage3",
                docker = insertion_annotation_docker,
                runtime_attr_override = runtime_attr_stage
        }

        call IterativeAlign as IterativeHg38 {
            input:
                queries_fa = select_first([Stage3.next_queries]),
                ref_fa = ref_fa,
                asm20_mmi = ref_asm20_mmi,
                sr_mmi = ref_sr_mmi,
                prefix = "~{prefix}.~{contigs[i]}.hg38",
                docker = insertion_annotation_docker,
                runtime_attr_override = runtime_attr_iterative_hg38
        }

        call IterativeAlign as IterativeT2T {
            input:
                queries_fa = select_first([Stage3.next_queries]),
                ref_fa = t2t_fa,
                prefix = "~{prefix}.~{contigs[i]}.t2t",
                docker = insertion_annotation_docker,
                runtime_attr_override = runtime_attr_iterative_t2t
        }

        call RunStage as Stage4 {
            input:
                stage = "stage4",
                contig = contigs[i],
                workdir_tar = Stage3.workdir_tar,
                hg38_iterative = IterativeHg38.segments_tsv,
                t2t_iterative = IterativeT2T.segments_tsv,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                repeats_bed = repeats_bed,
                tr_catalog_bed = tr_catalog_bed,
                immune_loci_bed = immune_loci_bed,
                known_cnv_loci_bed = known_cnv_loci_bed,
                cytoband = cytoband,
                sd_gd_blocks_bed = sd_gd_blocks_bed,
                sd_loci_bed = sd_loci_bed,
                sd_copies_bed = sd_copies_bed,
                genomic_superdups = genomic_superdups,
                t2t_genes = t2t_genes,
                prefix = "~{prefix}.~{contigs[i]}.stage4",
                docker = insertion_annotation_docker,
                runtime_attr_override = runtime_attr_stage
        }
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = select_all(Stage4.final_annotation),
            sort_output = false,
            preserve_header = true,
            prefix = "~{prefix}.insertion_origins",
            docker = insertion_annotation_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File insertion_origins_tsv = ConcatTsvs.concatenated_tsv
        Array[File] per_contig_workdirs = Stage4.workdir_tar
    }
}

task RunStage {
    input {
        String stage
        String contig
        File? workdir_tar
        File? remap_tsv
        File? vcf
        File? vcf_idx
        File? svan_tsv
        File? svan_header
        Int min_length = 50
        File? hg38_paf
        File? rm_out
        File? t2t_paf
        File? bwa_sam
        File? hg38_iterative
        File? t2t_iterative
        File ref_fa
        File ref_fai
        File repeats_bed
        File tr_catalog_bed
        File immune_loci_bed
        File known_cnv_loci_bed
        File cytoband
        File sd_gd_blocks_bed
        File sd_loci_bed
        File sd_copies_bed
        File genomic_superdups
        File t2t_genes
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p res work
        ln -s ~{repeats_bed} res/REPEATS_hg38.bed.gz
        ln -s ~{tr_catalog_bed} res/variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz
        ln -s ~{immune_loci_bed} res/immune_loci.hg38.bed
        ln -s ~{known_cnv_loci_bed} res/known_cnv_loci.hg38.bed
        ln -s ~{cytoband} res/hg38.cytoBand.txt
        ln -s ~{sd_gd_blocks_bed} res/hg38.sd_reference.gd_blocks.bed
        ln -s ~{sd_loci_bed} res/hg38.sd_reference.loci.bed
        ln -s ~{sd_copies_bed} res/hg38.sd_reference.copies.bed
        ln -s ~{genomic_superdups} res/genomicSuperDups.txt.gz
        ln -s ~{t2t_genes} res/hs1.ncbiRefSeqCurated.txt.gz
        ln -s ~{ref_fa} res/ref.fa
        ln -s ~{ref_fai} res/ref.fa.fai
        if [ -n "~{default="" workdir_tar}" ]; then
            tar -xzf ~{default="" workdir_tar} -C work
        fi

        S=/opt/insertion_annotation/annotate_insertions_stages.sh
        case ~{stage} in
            prepare)
                bash $S prepare -w work -c ~{contig} -r ~{default="" remap_tsv} -v ~{default="" vcf} -f res/ref.fa \
                    -s ~{default="" svan_tsv} -H ~{default="" svan_header} -m ~{min_length}
                cp work/output/~{contig}.queries.fa next_queries.fa ;;
            stage1)
                bash $S stage1 -w work -c ~{contig} -p ~{default="" hg38_paf} -o ~{default="" rm_out} -f res/ref.fa -R res
                cp work/output/~{contig}.t2t_queries.fa next_queries.fa ;;
            stage2)
                bash $S stage2 -w work -c ~{contig} -t ~{default="" t2t_paf} -f res/ref.fa -R res
                cp work/output/~{contig}.bwa_queries.fa next_queries.fa ;;
            stage3)
                bash $S stage3 -w work -c ~{contig} -b ~{default="" bwa_sam} -R res
                cp work/output/~{contig}.iterative_queries.fa next_queries.fa ;;
            stage4)
                bash $S stage4 -w work -c ~{contig} -i ~{default="" hg38_iterative} -j ~{default="" t2t_iterative} \
                    -f res/ref.fa -R res -t "$(nproc)"
                cp work/output/~{contig}.final_annotation.tsv final_annotation.tsv ;;
        esac
        if [ "~{stage}" == "prepare" ]; then
            cp work/output/~{contig}.inserted_sequences.fa inserted_sequences.fa
        fi
        rm -f work/scripts
        tar -czf ~{prefix}.workdir.tar.gz -C work .
    >>>

    output {
        File workdir_tar = "~{prefix}.workdir.tar.gz"
        File? next_queries = "next_queries.fa"
        File? inserted_sequences = "inserted_sequences.fa"
        File? final_annotation = "final_annotation.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: 3 * ceil(size(ref_fa, "GB")) + 10 * ceil(size(select_all([workdir_tar, remap_tsv, vcf, hg38_paf]), "GB")) + 30,
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

task AlignMinimap2 {
    input {
        File queries_fa
        File ref_fa
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ -s ~{queries_fa} ]; then
            minimap2 -x asm20 -c --max-chain-skip 50000 -N 500 -k 15 -w 5 -t "$(nproc)" ~{ref_fa} ~{queries_fa} \
                | gzip > ~{prefix}.alignments.paf.gz
        else
            : | gzip > ~{prefix}.alignments.paf.gz
        fi
    >>>

    output {
        File paf = "~{prefix}.alignments.paf.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 32,
        disk_gb: 2 * ceil(size(ref_fa, "GB")) + ceil(size(queries_fa, "GB")) + 20,
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

task RepeatMasker {
    input {
        File fa
        Int threads
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cp ~{fa} ~{prefix}.inserted_sequences.fa
        RepeatMasker -species human -pa ~{threads} -dir . ~{prefix}.inserted_sequences.fa
        if [ ! -s ~{prefix}.inserted_sequences.fa.out ]; then
            printf 'There were no repetitive sequences detected\n' > ~{prefix}.inserted_sequences.fa.out
        fi
    >>>

    output {
        File rm_out = "~{prefix}.inserted_sequences.fa.out"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: threads,
        mem_gb: 16,
        disk_gb: 10 * ceil(size(fa, "GB")) + 20,
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

task BwaMem {
    input {
        File queries_fa
        File ref_fa
        Array[File] ref_bwa_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir ref_files
        ln -s ~{ref_fa} ref_files/$(basename ~{ref_fa})
        for x in ~{sep=' ' ref_bwa_idx}; do
            ln -s "$x" ref_files/$(basename "$x")
        done
        if [ -s ~{queries_fa} ]; then
            bwa mem -t "$(nproc)" -a -Y ref_files/$(basename ~{ref_fa}) ~{queries_fa} | samtools view -h -F 4 - | gzip > ~{prefix}.bwa.sam.gz
        else
            : | gzip > ~{prefix}.bwa.sam.gz
        fi
    >>>

    output {
        File sam = "~{prefix}.bwa.sam.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: 2 * ceil(size(ref_fa, "GB") + size(ref_bwa_idx, "GB")) + 10,
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

task IterativeAlign {
    input {
        File queries_fa
        File ref_fa
        File? asm20_mmi
        File? sr_mmi
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ -s ~{queries_fa} ]; then
            A="~{default="" asm20_mmi}"; R="~{default="" sr_mmi}"
            if [ -z "$A" ]; then minimap2 -x asm20 -t "$(nproc)" -d asm20.mmi ~{ref_fa}; A=asm20.mmi; fi
            if [ -z "$R" ]; then minimap2 -x sr -t "$(nproc)" -d sr.mmi ~{ref_fa}; R=sr.mmi; fi
            python3 /opt/insertion_annotation/iterative_align.py --queries ~{queries_fa} --asm20-index "$A" --sr-index "$R" \
                --threads "$(nproc)" --out ~{prefix}.iterative_segments.tsv
        else
            printf 'ID\tq_start\tq_end\tref\tstrand\tidentity\tmapq\taligned_len\tpreset\tround\tstatus\n' > ~{prefix}.iterative_segments.tsv
        fi
    >>>

    output {
        File segments_tsv = "~{prefix}.iterative_segments.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 32,
        disk_gb: 6 * ceil(size(ref_fa, "GB")) + 20,
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
