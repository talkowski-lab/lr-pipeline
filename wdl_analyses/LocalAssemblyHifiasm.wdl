version 1.0

# Local hifiasm assembly of one sample's reads at one or more loci.
# Adapted from broadinstitute/long-read-pipelines (branch ym_testing) wdl/LocalAssemblyHifiasm.wdl, with
# the partially phased haplotype contigs (hap1/hap2) added to the outputs. hifiasm >= 0.15 is required for
# hap1/hap2 without parental data, so the default docker is lr-hifiasm:0.16.1 (was 0.13).
#
# Outputs:
#   primary_fa: primary contigs (hifiasm --primary; one pseudo-haplotype, not phased)
#   alt_fa:     alternate contigs (hifiasm --primary; the other haplotype where it diverges, fragmented)
#   hap1_fa:    haplotype 1 contigs (hifiasm default mode; phased within blocks only)
#   hap2_fa:    haplotype 2 contigs (hifiasm default mode; phased within blocks only)

import "Structs.wdl"

workflow LocalAssemblyHifiasm {
    input {
        Array[String]+ loci
        String aligned_bam
        File aligned_bai
        String prefix

        File? ref_fa
        File? ref_fai
        String? requester_pays_project
        Int hifiasm_cpus = 32

        String subset_docker = "quay.io/ymostovoy/lr-utils-basic:2.0"
        String utils_docker = "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
        String hifiasm_docker = "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.16.1"

        RuntimeAttr? runtime_attr_subset_bam
        RuntimeAttr? runtime_attr_merge_bams
        RuntimeAttr? runtime_attr_bam_to_fastq
        RuntimeAttr? runtime_attr_assemble
    }

    parameter_meta {
        loci: "Loci to assemble, e.g. [\"chr1:1000-2000\", \"chr1:5000-10000\"]; reads from all loci are merged for one assembly"
        aligned_bam: "aligned BAM/CRAM (gs:// path; streamed by samtools)"
        aligned_bai: "index of aligned_bam"
        prefix: "prefix for output files"
        ref_fa: "reference FASTA (required if aligned_bam is a CRAM)"
        requester_pays_project: "Google project to bill when aligned_bam is in a requester-pays bucket"
    }

    scatter (locus in loci) {
        call SubsetBam {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                locus = locus,
                prefix = prefix,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                requester_pays_project = requester_pays_project,
                docker = subset_docker,
                runtime_attr_override = runtime_attr_subset_bam
        }
    }

    if (length(loci) > 1) {
        call MergeBams {
            input:
                bams = SubsetBam.subset_bam,
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_bams
        }
    }

    File subset_bam = select_first([MergeBams.merged_bam, SubsetBam.subset_bam[0]])

    call BamToFastq {
        input:
            bam = subset_bam,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_bam_to_fastq
    }

    call Assemble {
        input:
            reads = BamToFastq.reads_fq,
            prefix = prefix,
            num_cpus = hifiasm_cpus,
            docker = hifiasm_docker,
            runtime_attr_override = runtime_attr_assemble
    }

    output {
        File primary_fa = Assemble.p_fa
        File alt_fa = Assemble.a_fa
        File hap1_fa = Assemble.h1_fa
        File hap2_fa = Assemble.h2_fa
    }
}

task SubsetBam {
    input {
        String bam
        File bai
        String locus
        String prefix
        File? ref_fa
        File? ref_fai
        String? requester_pays_project
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        # htslib >= 1.20 reads requester-pays buckets; the image has no gcloud, so get the token from the VM metadata server
        export GCS_OAUTH_TOKEN=$(curl -sS -H "Metadata-Flavor: Google" \
            http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
            | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')
        ~{if defined(requester_pays_project) then "export GCS_REQUESTER_PAYS_PROJECT=" + requester_pays_project else ""}

        samtools view -bhX ~{"-T " + ref_fa} ~{bam} ~{bai} ~{locus} > ~{prefix}.subset.bam
        samtools index ~{prefix}.subset.bam
    >>>

    output {
        File subset_bam = "~{prefix}.subset.bam"
        File subset_bai = "~{prefix}.subset.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10,
        disk_gb: 4 * ceil(size(bai, "GB")) + ceil(size(ref_fa, "GB")) + 20,
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

task MergeBams {
    input {
        Array[File] bams
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        samtools merge -p -c -@2 --no-PG ~{prefix}.merged.bam ~{sep=" " bams}
        samtools index ~{prefix}.merged.bam
    >>>

    output {
        File merged_bam = "~{prefix}.merged.bam"
        File merged_bai = "~{prefix}.merged.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 20,
        disk_gb: 4 * ceil(size(bams, "GB")) + 10,
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

task BamToFastq {
    input {
        File bam
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        samtools fastq ~{bam} | gzip > ~{prefix}.fq.gz
    >>>

    output {
        File reads_fq = "~{prefix}.fq.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(bam, "GB")) + 10,
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

task Assemble {
    input {
        File reads
        String prefix
        Int num_cpus = 32
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        # default mode (hifiasm >= 0.15): partially phased haplotype contigs, ~{prefix}.bp.hap{1,2}.p_ctg.gfa
        hifiasm -o ~{prefix} -t~{num_cpus} ~{reads}
        # --primary reuses the *.bin files of the first run and writes ~{prefix}.p_ctg.gfa and ~{prefix}.a_ctg.gfa
        hifiasm -o ~{prefix} -t~{num_cpus} --primary ~{reads}

        awk '/^S/{print ">"$2; print $3}' ~{prefix}.p_ctg.gfa > ~{prefix}.p_ctg.fa
        awk '/^S/{print ">"$2; print $3}' ~{prefix}.a_ctg.gfa > ~{prefix}.a_ctg.fa
        awk '/^S/{print ">"$2; print $3}' ~{prefix}.bp.hap1.p_ctg.gfa > ~{prefix}.bp.hap1.p_ctg.fa
        awk '/^S/{print ">"$2; print $3}' ~{prefix}.bp.hap2.p_ctg.gfa > ~{prefix}.bp.hap2.p_ctg.fa
    >>>

    output {
        File p_fa = "~{prefix}.p_ctg.fa"
        File a_fa = "~{prefix}.a_ctg.fa"
        File h1_fa = "~{prefix}.bp.hap1.p_ctg.fa"
        File h2_fa = "~{prefix}.bp.hap2.p_ctg.fa"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: num_cpus,
        mem_gb: 150,
        disk_gb: 10 * ceil(size(reads, "GB")) + 10,
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
