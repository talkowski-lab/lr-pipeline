version 1.0

import "../utils/Structs.wdl"

workflow AnnotateParaphase {
    input {
        File paraphase_json
        File? paraphase_vcfs
        Array[File] region_vcfs = []
        Array[String] region_names = []
        String sample_id
        String prefix

        String paraphase_version = "4.0.0"
        String genome_build = "GRCh38"
        Boolean gene1only
        Boolean targeted
        Int? x_copy_number
        String clinvar_release
        Int vep_cache_version
        String vep_executable = "vep"

        File rules
        File parser_script
        File ref_fa
        File ref_fai
        File ref_vep_cache
        File clinvar_vcf
        File clinvar_vcf_idx

        String python_docker
        String vep_bcftools_docker

        RuntimeAttr? runtime_attr_prepare
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_report
    }

    call PrepareParaphase {
        input:
            paraphase_json = paraphase_json,
            paraphase_vcfs = paraphase_vcfs,
            region_vcfs = region_vcfs,
            region_names = region_names,
            sample_id = sample_id,
            prefix = prefix,
            paraphase_version = paraphase_version,
            genome_build = genome_build,
            gene1only = gene1only,
            targeted = targeted,
            rules = rules,
            parser_script = parser_script,
            ref_fai = ref_fai,
            docker = python_docker,
            runtime_attr_override = runtime_attr_prepare
    }

    call AnnotateParaphaseSites {
        input:
            sites_vcf = PrepareParaphase.sites_vcf,
            prefix = prefix,
            clinvar_release = clinvar_release,
            vep_cache_version = vep_cache_version,
            vep_executable = vep_executable,
            parser_script = parser_script,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            ref_vep_cache = ref_vep_cache,
            clinvar_vcf = clinvar_vcf,
            clinvar_vcf_idx = clinvar_vcf_idx,
            docker = vep_bcftools_docker,
            runtime_attr_override = runtime_attr_annotate
    }

    call ReportParaphase {
        input:
            prepared = PrepareParaphase.prepared,
            annotated_vcf = AnnotateParaphaseSites.annotated_vcf,
            vep_json = AnnotateParaphaseSites.vep_json,
            resources = AnnotateParaphaseSites.resources,
            prefix = prefix,
            x_copy_number = x_copy_number,
            rules = rules,
            parser_script = parser_script,
            docker = python_docker,
            runtime_attr_override = runtime_attr_report
    }

    output {
        File findings_tsv = ReportParaphase.findings_tsv
        File gene_summary_tsv = ReportParaphase.gene_summary_tsv
        File variants_tsv = ReportParaphase.variants_tsv
        File haplotypes_tsv = ReportParaphase.haplotypes_tsv
        File qc_tsv = ReportParaphase.qc_tsv
        File report_json = ReportParaphase.report_json
        File provenance_json = ReportParaphase.provenance_json
        File prepared_json = PrepareParaphase.prepared
        File vep_json = AnnotateParaphaseSites.vep_json
        File annotated_vcf = AnnotateParaphaseSites.annotated_vcf
        File annotated_vcf_idx = AnnotateParaphaseSites.annotated_vcf_idx
    }
}

task PrepareParaphase {
    input {
        File paraphase_json
        File? paraphase_vcfs
        Array[File] region_vcfs
        Array[String] region_names
        String sample_id
        String prefix
        String paraphase_version
        String genome_build
        Boolean gene1only
        Boolean targeted
        File rules
        File parser_script
        File ref_fai
        String docker
        RuntimeAttr? runtime_attr_override
    }

    File request = write_json(object {
        paraphase_json: paraphase_json,
        paraphase_vcfs: paraphase_vcfs,
        region_vcfs: region_vcfs,
        region_names: region_names,
        sample_id: sample_id,
        prefix: prefix,
        paraphase_version: paraphase_version,
        genome_build: genome_build,
        gene1only: gene1only,
        targeted: targeted,
        rules: rules,
        container: docker,
        ref_fai: ref_fai
    })

    command <<<
        set -euo pipefail

        python3 '~{parser_script}' prepare --request '~{request}'
    >>>

    output {
        File prepared = prefix + ".prepared.json"
        File sites_vcf = prefix + ".sites.vcf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10 + ceil(3 * (size(paraphase_json, "GB") + size(region_vcfs, "GB")) + 10 * size(paraphase_vcfs, "GB")),
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

task AnnotateParaphaseSites {
    input {
        File sites_vcf
        String prefix
        String clinvar_release
        Int vep_cache_version
        String vep_executable
        File parser_script
        File ref_fa
        File ref_fai
        File ref_vep_cache
        File clinvar_vcf
        File clinvar_vcf_idx
        String docker
        RuntimeAttr? runtime_attr_override
    }

    File request = write_json(object {
        sites_vcf: sites_vcf,
        container: docker,
        prefix: prefix,
        clinvar_release: clinvar_release,
        vep_cache_version: vep_cache_version,
        vep_executable: vep_executable,
        ref_fa: ref_fa,
        ref_fai: ref_fai,
        ref_vep_cache: ref_vep_cache,
        clinvar_vcf: clinvar_vcf,
        clinvar_vcf_idx: clinvar_vcf_idx
    })

    command <<<
        set -euo pipefail

        python3 '~{parser_script}' annotate --request '~{request}'
    >>>

    output {
        File annotated_vcf = prefix + ".annotated.vcf.gz"
        File annotated_vcf_idx = prefix + ".annotated.vcf.gz.tbi"
        File vep_json = prefix + ".vep.jsonl"
        File resources = prefix + ".resources.json"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 20 + ceil(size(ref_fa, "GB") * 2 + size(ref_vep_cache, "GB") * 5 + size(clinvar_vcf, "GB") * 2 + size(sites_vcf, "GB") * 5),
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

task ReportParaphase {
    input {
        File prepared
        File annotated_vcf
        File vep_json
        File resources
        String prefix
        Int? x_copy_number
        File rules
        File parser_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    File request = write_json(object {
        prepared: prepared,
        annotated_vcf: annotated_vcf,
        vep_json: vep_json,
        resources: resources,
        prefix: prefix,
        x_copy_number: x_copy_number,
        rules: rules
    })

    command <<<
        set -euo pipefail

        python3 '~{parser_script}' report --request '~{request}'
    >>>

    output {
        File findings_tsv = prefix + ".findings.tsv"
        File gene_summary_tsv = prefix + ".gene_summary.tsv"
        File variants_tsv = prefix + ".variants.tsv"
        File haplotypes_tsv = prefix + ".haplotypes.tsv"
        File qc_tsv = prefix + ".qc.tsv"
        File report_json = prefix + ".report.json"
        File provenance_json = prefix + ".provenance.json"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10 + ceil(3 * (size(prepared, "GB") + size(annotated_vcf, "GB") + size(vep_json, "GB"))),
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
