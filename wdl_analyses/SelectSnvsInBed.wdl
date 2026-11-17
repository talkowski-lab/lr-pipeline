version 1.0

import "Structs.wdl"

# Select SNVs from one or more VCF/BCF files that fall inside a BED of regions (e.g. high-confidence exome regions for
# PCA). Each input file (e.g. one per chromosome) is processed in parallel, then the per-file outputs are concatenated
# into one VCF. All input files must share the same samples in the same order. Records are kept when they overlap a BED
# interval and are SNVs. With biallelic_only, multiallelic sites are dropped (including those with spanning-deletion
# alleles); with pass_only, only FILTER=PASS records are kept (set false for pre-filtered files whose FILTER is ".").
# With min_af, only SNVs with AF > min_af are kept, where AF is recomputed from the genotypes in the file
# (bcftools +fill-tags, overwriting any INFO/AF). Genotypes are kept unchanged. Run the same workflow on case and
# control VCFs so both sides share the same site set. vcfs and vcf_idxs are gs:// paths read in place over the network
# and never localized: the BED is split into batches of regions_per_batch regions (a batch never spans contigs), each
# batch is fetched with a remote tabix query and trimmed (SNV, PASS, biallelic, AF) on the fly, n_parallel batches at a
# time, each with a fresh OAuth token from the VM metadata server, then the batches are concatenated in BED order. Set
# requester_pays_project for requester-pays buckets (e.g. All of Us). Set gvcf = true for gVCF input (*.g.vcf.bgz): every
# record carries a <NON_REF> (or <*>) allele, so alleles not seen in any genotype are trimmed first (bcftools view -a);
# reference blocks and RefCall sites then have no ALT and drop out, and A -> G,<NON_REF> becomes a biallelic A -> G SNV.
# With sites_tsv (tab-separated chr, pos, end, ref, alt, optionally gzipped, header skipped; pos is the 1-based VCF
# position, contig names as in the VCF), the output is further restricted to SNVs whose CHROM:POS:REF:ALT is in the
# list (e.g. SNVs common in other cohorts); kept records get ID = CHROM:POS:REF:ALT.

workflow SelectSnvsInBed {
    input {
        Array[String] vcfs
        Array[String] vcf_idxs
        File regions_bed
        String prefix

        Boolean pass_only = true
        Boolean biallelic_only = true
        Boolean gvcf = false
        Float? min_af
        File? sites_tsv
        String? requester_pays_project
        Int regions_per_batch = 500
        Int n_parallel = 16

        String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"

        RuntimeAttr? runtime_attr_select_snvs_in_bed
        RuntimeAttr? runtime_attr_concat_snvs
    }

    scatter (i in range(length(vcfs))) {
        call SelectSnvs {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                regions_bed = regions_bed,
                pass_only = pass_only,
                biallelic_only = biallelic_only,
                gvcf = gvcf,
                min_af = min_af,
                sites_tsv = sites_tsv,
                requester_pays_project = requester_pays_project,
                regions_per_batch = regions_per_batch,
                n_parallel = n_parallel,
                prefix = "~{prefix}.shard_~{i}",
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_select_snvs_in_bed
        }
    }

    call ConcatSnvs {
        input:
            vcfs = SelectSnvs.snv_vcf,
            vcf_idxs = SelectSnvs.snv_vcf_idx,
            prefix = prefix,
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat_snvs
    }

    output {
        File snv_vcf = ConcatSnvs.snv_vcf
        File snv_vcf_idx = ConcatSnvs.snv_vcf_idx
        File snv_count = ConcatSnvs.snv_count
    }
}

task SelectSnvs {
    input {
        String vcf
        String vcf_idx
        File regions_bed
        Boolean pass_only
        Boolean biallelic_only
        Boolean gvcf
        Float? min_af
        File? sites_tsv
        String? requester_pays_project
        Int regions_per_batch
        Int n_parallel
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String pass_arg = if pass_only then "-f PASS" else ""
    String biallelic_arg = if biallelic_only then "-m2 -M2" else ""
    String gvcf_cmd = if gvcf then "bcftools view -a -Ou | " else ""
    String sites_targets_arg = if defined(sites_tsv) then "-T sites.pos" else ""
    String sites_cmd = if defined(sites_tsv) then "bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Ou | bcftools view -i 'ID=@sites.ids' -Ou | " else ""
    String trim_cmd = if defined(min_af) then "bcftools +fill-tags - -Ou -- -t AF | bcftools view -i 'INFO/AF>" + select_first([min_af, 0.0]) + "' -Ob" else "bcftools view -Ob"

    command <<<
        set -euo pipefail

        cat > get_token.sh <<'EOF'
        wget -q -O - --header "Metadata-Flavor: Google" \
            "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token" \
            | sed -E 's/.*"access_token" *: *"([^"]*)".*/\1/'
        EOF
        if [ -n "~{requester_pays_project}" ]; then
            export GCS_REQUESTER_PAYS_PROJECT="~{requester_pays_project}"
        fi
        export INPUT="~{vcf}##idx##~{vcf_idx}"

        export GCS_OAUTH_TOKEN=$(bash get_token.sh)
        bcftools index -s "$INPUT" | cut -f1 > file_contigs.txt
        awk 'FNR == NR {keep[$1] = 1; next} ($1 in keep)' file_contigs.txt ~{regions_bed} > regions.bed
        if [ -n "~{sites_tsv}" ]; then
            gzip -cdf ~{sites_tsv} | awk -F'\t' '$2 ~ /^[0-9]+$/ {print $1 "\t" $2 > "sites.pos"; print $1 ":" $2 ":" $4 ":" $5 > "sites.ids"}'
            echo "$(wc -l < sites.ids) sites in sites_tsv"
        fi
        mkdir batches
        awk -v n=~{regions_per_batch} '
            $1 != prev || c == n {if (f) close(f); b++; c = 0; prev = $1; f = sprintf("batches/batch.%06d.bed", b)}
            {c++; print > f}' regions.bed
        echo "$(wc -l < regions.bed) regions on $(wc -l < file_contigs.txt) contigs in $(ls batches | wc -l) batches"

        cat > select_batch.sh <<'EOF'
        set -euo pipefail
        bed="$1"
        out="${bed%.bed}.bcf"
        for attempt in 1 2 3; do
            export GCS_OAUTH_TOKEN=$(bash get_token.sh)
            if bcftools view -R "$bed" ~{sites_targets_arg} ~{pass_arg} -Ou "$INPUT" | ~{gvcf_cmd}bcftools view -v snps ~{biallelic_arg} -Ou | ~{sites_cmd}~{trim_cmd} -o "$out"; then
                exit 0
            fi
            echo "retrying $bed (attempt $attempt failed)" >&2
            sleep $((attempt * 30))
        done
        exit 1
        EOF
        ls batches/*.bed | xargs -P ~{n_parallel} -n 1 bash select_batch.sh

        ls batches/*.bcf > parts.txt
        bcftools concat --file-list parts.txt --threads 2 -Oz -o ~{prefix}.snvs_in_bed.vcf.gz
        bcftools index -t ~{prefix}.snvs_in_bed.vcf.gz
        echo "$(bcftools index -n ~{prefix}.snvs_in_bed.vcf.gz) SNVs selected"
    >>>

    output {
        File snv_vcf = "~{prefix}.snvs_in_bed.vcf.gz"
        File snv_vcf_idx = "~{prefix}.snvs_in_bed.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 50,
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

task ConcatSnvs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools concat \
            -a \
            --file-list ~{write_lines(vcfs)} \
            --threads 2 \
            -Oz -o ~{prefix}.snvs_in_bed.vcf.gz
        bcftools index -t ~{prefix}.snvs_in_bed.vcf.gz
        bcftools index -n ~{prefix}.snvs_in_bed.vcf.gz > ~{prefix}.snvs_in_bed.count.txt
    >>>

    output {
        File snv_vcf = "~{prefix}.snvs_in_bed.vcf.gz"
        File snv_vcf_idx = "~{prefix}.snvs_in_bed.vcf.gz.tbi"
        File snv_count = "~{prefix}.snvs_in_bed.count.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 10,
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
