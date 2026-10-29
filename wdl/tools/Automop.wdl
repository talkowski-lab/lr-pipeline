version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow Automop {
    meta {
        description: [
            "This tool runs `mop` (via FISS) to clean up unreferenced intermediate files in a Terra workspace, freeing storage. A dry-run mode reports what would be deleted without removing anything."
        ]
    }

    parameter_meta {
        workspace_namespace: "Terra workspace namespace to clean."
        workspace_name: "Terra workspace name to clean."
        dry_run: "Whether to report rather than perform deletions."
        fissfc_log: "Log of the cleanup run."
    }

    input {
        String workspace_namespace
        String workspace_name
        String prefix

        Boolean dry_run

        String automop_docker

        RuntimeAttr? runtime_attr_run_mop
    }

    call RunMop {
        input:
            workspace_namespace = workspace_namespace,
            workspace_name = workspace_name,
            dry_run = dry_run,
            prefix = prefix,
            docker = automop_docker,
            runtime_attr_override = runtime_attr_run_mop
    }

    output {
        File fissfc_log = RunMop.fissfc_log
    }
}

task RunMop {
    input {
        String workspace_namespace
        String workspace_name
        Boolean dry_run
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        fissfc --yes --verbose mop -w "~{workspace_name}" -p "~{workspace_namespace}"~{if dry_run then " --dry-run" else ""} > "~{prefix}.log"

        # fissfc exits 0 even when gsutil fails to delete files, so check its log for completion
        if grep -q '^CommandException' "~{prefix}.log" \
            || ! grep -q -e '^No files to mop in' -e '^~{if dry_run then "Total Size: " else "Operation completed over "}' "~{prefix}.log"; then
            echo "fissfc did not report a completed mop; see ~{prefix}.log" >&2
            exit 1
        fi
        grep -e '^No files to mop in' -e '^Total Size: ' -e '^Operation completed over ' "~{prefix}.log"
    >>>

    output {
        File fissfc_log = "~{prefix}.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 16,
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
