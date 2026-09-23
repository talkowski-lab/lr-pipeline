version 1.0

import "../utils/Structs.wdl"

workflow AnnotateExternalVariants {
    meta {
        description: [
            "This utility matches an evaluation VCF against a truth VCF by structural variant type. Both callsets are converted to BED and split into deletions, duplications and insertions, compared with `bedtools` both within type and across the duplication and insertion types, and the per-type results are combined into one TSV of matched variants."
        ]
    }

    parameter_meta {
        vcf_eval: "VCF being evaluated."
        vcf_eval_idx: "Index for `vcf_eval`."
        vcf_truth: "Truth VCF to evaluate against."
        vcf_truth_idx: "Index for `vcf_truth`."
        population: "Population labels whose allele frequencies are carried across from the truth callset."
        matched_variants_tsv: "TSV pairing each evaluation variant with its matched truth variant."
    }


    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
        String prefix

        Array[String] population

        String sv_pipeline_docker
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_split_eval_vcf
        RuntimeAttr? runtime_attr_split_truth_vcf
        RuntimeAttr? runtime_attr_bedtools_closest
        RuntimeAttr? runtime_attr_select_matched_svs
        RuntimeAttr? runtime_attr_collect_matches
    }

    # Convert eval VCF to BED and split by SV type
    call SplitEvalVcf {
        input:
            vcf = vcf_eval,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_split_eval_vcf
    }

    # Convert truth VCF to BED and split by SV type
    call SplitTruthVcf {
        input:
            vcf = vcf_truth,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_split_truth_vcf
    }

    # DEL vs DEL
    call BedtoolsClosest as compare_del {
        input:
            bed_a = SplitEvalVcf.del,
            bed_b = SplitTruthVcf.del,
            svtype = "del",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }

    # DUP vs DUP
    call BedtoolsClosest as compare_dup {
        input:
            bed_a = SplitEvalVcf.dup,
            bed_b = SplitTruthVcf.dup,
            svtype = "dup",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }

    # INS vs INS
    call BedtoolsClosest as compare_ins {
        input:
            bed_a = SplitEvalVcf.ins,
            bed_b = SplitTruthVcf.ins,
            svtype = "ins",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }

    # DUP(eval) vs INS(truth)
    call BedtoolsClosest as compare_dup_as_ins {
        input:
            bed_a = SplitEvalVcf.dup,
            bed_b = SplitTruthVcf.ins,
            svtype = "dup_as_ins",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }

    # INS(eval) vs DUP(truth)
    call BedtoolsClosest as compare_ins_as_dup {
        input:
            bed_a = SplitEvalVcf.ins,
            bed_b = SplitTruthVcf.dup,
            svtype = "ins_as_dup",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }

    # Match DELs using reciprocal overlap (CNV-style)
    call SelectMatchedSVs as calcu_del {
        input:
            input_bed = compare_del.output_bed,
            svtype = "del",
            population = population,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    # Match DUPs using reciprocal overlap (CNV-style)
    call SelectMatchedSVs as calcu_dup {
        input:
            input_bed = compare_dup.output_bed,
            svtype = "dup",
            population = population,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    # Match INS vs INS using distance/ratio
    call SelectMatchedINSs as calcu_ins {
        input:
            input_bed = compare_ins.output_bed,
            svtype = "ins",
            population = population,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    # Match DUP(eval) vs INS(truth) using distance/ratio
    call SelectMatchedINSs as calcu_dup_as_ins {
        input:
            input_bed = compare_dup_as_ins.output_bed,
            svtype = "dup_as_ins",
            population = population,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    # Match INS(eval) vs DUP(truth) using distance/ratio
    call SelectMatchedINSs as calcu_ins_as_dup {
        input:
            input_bed = compare_ins_as_dup.output_bed,
            svtype = "ins_as_dup",
            population = population,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    # Collect all matches into a single TSV
    call CollectMatches {
        input:
            vcf_eval = vcf_eval,
            labeled_del = calcu_del.output_comp,
            labeled_dup = calcu_dup.output_comp,
            labeled_ins = calcu_ins.output_comp,
            labeled_dup_as_ins = calcu_dup_as_ins.output_comp,
            labeled_ins_as_dup = calcu_ins_as_dup.output_comp,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_collect_matches
    }

    output {
        File matched_variants_tsv = CollectMatches.matched_tsv
    }
}

task SplitEvalVcf {
    input {
        File vcf
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail

        python3 <<CODE
        import sys
        from pysam import VariantFile
        vcf_in = VariantFile("~{vcf}")
        beds = {}
        for svtype in ["DEL", "DUP", "INS"]:
            beds[svtype] = open(f"~{prefix}.{svtype}.bed", "w")
            beds[svtype].write("#chrom\tstart\tend\tname\tsvtype\tsvlen\n")
        for rec in vcf_in.fetch():
            allele_type = rec.info.get("allele_type", None)
            if allele_type is None:
                continue
            allele_type = allele_type.lower()
            allele_length = rec.info.get("allele_length", 0)
            if isinstance(allele_length, tuple):
                allele_length = allele_length[0]
            if "del" in allele_type:
                svtype = "DEL"
                svlen = abs(allele_length)
                start = rec.pos
                end = rec.pos + svlen
            elif "dup" in allele_type:
                svtype = "DUP"
                svlen = abs(allele_length)
                start = rec.pos
                end = rec.pos + svlen
            elif "ins" in allele_type:
                svtype = "INS"
                svlen = abs(allele_length)
                start = rec.pos
                end = rec.pos + 1
            else:
                svtype = "INS"
                svlen = abs(allele_length)
                start = rec.pos
                end = rec.pos + 1
            beds[svtype].write(f"{rec.chrom}\t{start}\t{end}\t{rec.id}\t{svtype}\t{svlen}\n")
        for f in beds.values():
            f.close()
        CODE
    >>>

    output {
        File del = "~{prefix}.DEL.bed"
        File dup = "~{prefix}.DUP.bed"
        File ins = "~{prefix}.INS.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SplitTruthVcf {
    input {
        File vcf
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail
        svtk vcf2bed -i SVTYPE -i SVLEN ~{vcf} tmp.bed
        cut -f1-4,7-8 tmp.bed > ~{prefix}.bed
        set +o pipefail
        head -1 ~{prefix}.bed > header
        set -o pipefail
        cat header <(awk '{if ($5=="DEL") print}' ~{prefix}.bed) > ~{prefix}.DEL.bed
        cat header <(awk '{if ($5=="DUP") print}' ~{prefix}.bed) > ~{prefix}.DUP.bed
        cat header <(awk '{if ($5=="INS" || $5=="INS:ME" || $5=="INS:ME:ALU" || $5=="INS:ME:LINE1" || $5=="INS:ME:SVA" || $5=="ALU" || $5=="LINE1" || $5=="SVA" || $5=="HERVK") print}' ~{prefix}.bed) > ~{prefix}.INS.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
        File del = "~{prefix}.DEL.bed"
        File dup = "~{prefix}.DUP.bed"
        File ins = "~{prefix}.INS.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3,
        disk_gb: 10,
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
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task BedtoolsClosest {
    input {
        File bed_a
        File bed_b
        String svtype
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eu
        paste <(head -1 ~{bed_a}) <(head -1 ~{bed_b}) | sed -e "s/#//g" > ~{svtype}.bed
        set -o pipefail
        bedtools closest -wo -a <(sort -k1,1 -k2,2n ~{bed_a}) -b <(sort -k1,1 -k2,2n ~{bed_b}) >> ~{svtype}.bed
    >>>

    output {
        File output_bed = "~{svtype}.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3,
        disk_gb: 5,
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
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SelectMatchedSVs {
    input {
        File input_bed
        String svtype
        Array[String] population
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    String matched_prefix = basename(input_bed, ".bed")
    File pop_list = write_lines(population)

    command <<<
        set -euo pipefail
        Rscript /opt/sv-pipeline/05_annotation/scripts/R1.bedtools_closest_CNV.R \
            -i ~{input_bed} \
            -o ~{matched_prefix}.comparison \
            -p ~{pop_list}
    >>>

    output {
        File output_comp = "~{matched_prefix}.comparison"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3,
        disk_gb: 5,
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
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SelectMatchedINSs {
    input {
        File input_bed
        String svtype
        Array[String] population
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    String matched_prefix = basename(input_bed, ".bed")
    File pop_list = write_lines(population)

    command <<<
        Rscript /opt/sv-pipeline/05_annotation/scripts/R2.bedtools_closest_INS.R \
            -i ~{input_bed} \
            -o ~{matched_prefix}.comparison \
            -p ~{pop_list}
    >>>

    output {
        File output_comp = "~{matched_prefix}.comparison"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3,
        disk_gb: 5,
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
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CollectMatches {
    input {
        File vcf_eval
        File labeled_del
        File labeled_dup
        File labeled_ins
        File labeled_dup_as_ins
        File labeled_ins_as_dup
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
        import csv
        from pysam import VariantFile
        matches = {}
        for comp_file in ["~{labeled_del}", "~{labeled_dup}", "~{labeled_ins}", "~{labeled_dup_as_ins}", "~{labeled_ins_as_dup}"]:
            with open(comp_file) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    qid = row["query_svid"]
                    rid = row["ref_svid"]
                    if qid not in matches:
                        matches[qid] = rid
        vcf_in = VariantFile("~{vcf_eval}")
        with open("~{prefix}.matched_variants.tsv", "w") as out:
            out.write("CHROM\tPOS\tREF\tALT\tID\tMATCHED_ID\n")
            for rec in vcf_in.fetch():
                if rec.id in matches:
                    alt = ",".join(str(a) for a in rec.alts) if rec.alts else "."
                    out.write(f"{rec.chrom}\t{rec.pos}\t{rec.ref}\t{alt}\t{rec.id}\t{matches[rec.id]}\n")
        CODE
    >>>

    output {
        File matched_tsv = "~{prefix}.matched_variants.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
