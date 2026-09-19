version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow RepeatMasker {
    input {
        File vcf
        File vcf_idx
        String prefix

        Int? min_length

        String repeatmasker_docker
        String utils_docker

        RuntimeAttr? runtime_attr_ins_to_fa
        RuntimeAttr? runtime_attr_repeat_masker
    }

    call INSToFa {
        input:
            vcf = vcf,
            min_length = min_length,
            prefix = "~{prefix}.ins",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_ins_to_fa
    }

    call RunRepeatMasker {
        input:
            fa = INSToFa.ins_fa,
            prefix = "~{prefix}.rm",
            docker = repeatmasker_docker,
            runtime_attr_override = runtime_attr_repeat_masker
    }

    output {
        File rm_out = RunRepeatMasker.rm_out
        File rm_fa = INSToFa.ins_fa
    }
}

task INSToFa {
    input {
        File vcf
        Int? min_length
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -i '~{if defined(min_length) then "abs(INFO/allele_length) >= " + min_length + " && " else ""}INFO/allele_type == "ins"' \
            -Oz -o ~{prefix}.subset.vcf.gz \
            ~{vcf}

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\n' \
            ~{prefix}.subset.vcf.gz \
        | awk 'length($3)==1 {print ">"$1":"$2";"$3"\n"$4}' > ~{prefix}.tmp.fa

        seqkit rename -N1 ~{prefix}.tmp.fa > ~{prefix}.fa
    >>>

    output {
        File ins_fa = '~{prefix}.fa'
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 20,
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

task RunRepeatMasker {
    input {
        File fa
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        RepeatMasker \
            -e rmblast \
            -species human \
            -pa 4 \
            -s \
            ~{fa}

        mv ~{fa}.out ./~{prefix}.out
    >>>

    output {
        File rm_out = '~{prefix}.out'
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(fa, "GB")) + 20,
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
