version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow BLAT {
    input {
        File queries_fa
        String prefix

        String blat_args = "-stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0"
        Int ooc_rep_match = 1024

        File ref_fa

        String blat_docker

        RuntimeAttr? runtime_attr_blat
    }

    call BlatQueries {
        input:
            queries_fa = queries_fa,
            blat_args = blat_args,
            ooc_rep_match = ooc_rep_match,
            ref_fa = ref_fa,
            prefix = prefix,
            docker = blat_docker,
            runtime_attr_override = runtime_attr_blat
    }

    output {
        File blat_psl = BlatQueries.blat_psl
    }
}

task BlatQueries {
    input {
        File queries_fa
        String blat_args
        Int ooc_rep_match
        File ref_fa
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        faToTwoBit ~{ref_fa} ref.2bit
        blat ref.2bit /dev/null /dev/null -tileSize=11 -makeOoc=11.ooc -repMatch=~{ooc_rep_match}

        n_shards=$(nproc)
        mkdir -p shards
        python3 - ~{queries_fa} "${n_shards}" <<'EOF'
import sys

n = int(sys.argv[2])
outs = [open(f"shards/q{i}.fa", "w") for i in range(n)]
k = -1
for line in open(sys.argv[1]):
    if line.startswith(">"):
        k += 1
    outs[k % n].write(line)
for o in outs:
    o.close()
EOF

        ls shards/q*.fa | xargs -P "${n_shards}" -I {} sh -c 'blat ref.2bit {} -ooc=11.ooc ~{blat_args} -noHead {}.psl'
        cat shards/q*.fa.psl | gzip > ~{prefix}.blat.psl.gz
    >>>

    output {
        File blat_psl = "~{prefix}.blat.psl.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 32,
        disk_gb: 3 * ceil(size(ref_fa, "GB")) + ceil(size(queries_fa, "GB")) + 10,
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
