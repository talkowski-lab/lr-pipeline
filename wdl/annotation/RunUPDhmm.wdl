version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

# RunUPDhmm — detect uniparental disomy (UPD) across a cohort of families with the
# Bioconductor UPDhmm package.
#
# UPDhmm runs an HMM over the joint genotypes of a *complete trio* (proband + mother +
# father) and calls iso/heterodisomy segments. It cannot run on singletons or duos, so
# those are skipped rather than merged.
#
# Inputs are the per-sample phased SNV/indel VCFs (DeepVariant -> GLNexus -> split-by-
# sample -> HiPhase) referenced from a Terra data table, plus a 6-column PED describing
# family relationships. Sample VCF/index paths are passed as Array[String] (not File) so
# Cromwell localizes only the three VCFs needed by each scatter shard, not the whole cohort.
#
# For families with multiple children, one trio is emitted per child that has both parents
# genotyped. UPD is called on autosomes (chr1-22) only.
workflow RunUPDhmm {
    input {
        # 6-column PED: family_id, individual_id, father_id, mother_id, sex, phenotype.
        # father_id/mother_id of "0" denote a missing parent.
        File ped

        # Per-sample arrays, aligned by index, straight from the sample data table.
        Array[String] sample_ids
        Array[String] sample_vcfs
        Array[String] sample_vcf_indexes

        String prefix

        # Restrict trios to affected children (PED phenotype == 2).
        Boolean affected_only = false

        # Autosomes to call UPD on.
        Array[String] autosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                                   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                                   "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                                   "chr22"]

        # Optional per-genotype GQ floor; genotypes below it are set to missing before merge.
        Int? min_gq

        String bcftools_docker = "quay.io/ymostovoy/lr-utils-basic:latest"
        String updhmm_docker = "quay.io/ymostovoy/lr-updhmm:latest"

        RuntimeAttr? runtime_attr_prepare
        RuntimeAttr? runtime_attr_make_trio
        RuntimeAttr? runtime_attr_updhmm
        RuntimeAttr? runtime_attr_concat
    }

    call PrepareTrios {
        input:
            ped = ped,
            sample_ids = sample_ids,
            sample_vcfs = sample_vcfs,
            sample_vcf_indexes = sample_vcf_indexes,
            affected_only = affected_only,
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_prepare
    }

    scatter (trio in read_tsv(PrepareTrios.trios_data)) {
        call MakeTrioVcf {
            input:
                family_id = trio[0],
                proband_id = trio[1],
                father_id = trio[2],
                mother_id = trio[3],
                proband_vcf = trio[4],
                proband_vcf_idx = trio[5],
                father_vcf = trio[6],
                father_vcf_idx = trio[7],
                mother_vcf = trio[8],
                mother_vcf_idx = trio[9],
                autosomes = autosomes,
                min_gq = min_gq,
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_make_trio
        }

        call RunUPDhmmTask {
            input:
                trio_vcf = MakeTrioVcf.trio_vcf,
                trio_vcf_idx = MakeTrioVcf.trio_vcf_idx,
                family_id = trio[0],
                proband_id = trio[1],
                father_id = trio[2],
                mother_id = trio[3],
                docker = updhmm_docker,
                runtime_attr_override = runtime_attr_updhmm
        }
    }

    call Helpers.ConcatTsvs as ConcatEvents {
        input:
            tsvs = RunUPDhmmTask.events_tsv,
            prefix = "~{prefix}.upd_events",
            preserve_header = true,
            skip_sort = true,
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File cohort_upd_events = ConcatEvents.concatenated_tsv
        File trios_manifest = PrepareTrios.trios_manifest
        File skipped_manifest = PrepareTrios.skipped_manifest
    }
}

# Parse the PED + per-sample arrays into a per-trio manifest. Emits:
#  - trios.data.tsv    : headerless, one row per complete trio (consumed by read_tsv scatter)
#  - trios_manifest.tsv: same rows with a header (human-readable output)
#  - skipped_manifest.tsv: individuals that could not form a complete trio, with reason
# Fails if no complete trios are found. Takes only String arrays + the PED, so no cohort
# VCFs are localized here.
task PrepareTrios {
    input {
        File ped
        Array[String] sample_ids
        Array[String] sample_vcfs
        Array[String] sample_vcf_indexes
        Boolean affected_only
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail

        python3 <<'PYEOF'
import sys

ped_path      = "~{ped}"
affected_only = ~{if affected_only then "True" else "False"}

ids  = [l.strip() for l in open("~{write_lines(sample_ids)}")        if l.strip()]
vcfs = [l.strip() for l in open("~{write_lines(sample_vcfs)}")       if l.strip()]
idxs = [l.strip() for l in open("~{write_lines(sample_vcf_indexes)}") if l.strip()]

if not (len(ids) == len(vcfs) == len(idxs)):
    sys.exit("ERROR: sample_ids (%d), sample_vcfs (%d) and sample_vcf_indexes (%d) "
             "differ in length" % (len(ids), len(vcfs), len(idxs)))

smap = {s: (v, x) for s, v, x in zip(ids, vcfs, idxs)}

trios, skipped = [], []
with open(ped_path) as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 6:
            continue
        fid, iid, pat, mat, sex, pheno = f[0], f[1], f[2], f[3], f[4], f[5]

        if pat == "0" and mat == "0":
            skipped.append((fid, iid, "founder_or_singleton"))
            continue
        if pat == "0" or mat == "0":
            skipped.append((fid, iid, "duo_unsupported_by_updhmm"))
            continue
        if affected_only and pheno != "2":
            skipped.append((fid, iid, "child_not_affected"))
            continue

        missing = [m for m in (iid, pat, mat) if m not in smap]
        if missing:
            skipped.append((fid, iid, "missing_from_callset:" + ",".join(missing)))
            continue

        pv, px = smap[iid]
        fv, fx = smap[pat]
        mv, mx = smap[mat]
        trios.append((fid, iid, pat, mat, pv, px, fv, fx, mv, mx))

cols = ["family_id", "proband_id", "father_id", "mother_id",
        "proband_vcf", "proband_idx", "father_vcf", "father_idx",
        "mother_vcf", "mother_idx"]

with open("trios_manifest.tsv", "w") as out:
    out.write("\t".join(cols) + "\n")
    for t in trios:
        out.write("\t".join(t) + "\n")

with open("trios.data.tsv", "w") as out:   # headerless: consumed by read_tsv
    for t in trios:
        out.write("\t".join(t) + "\n")

with open("skipped_manifest.tsv", "w") as out:
    out.write("family_id\tindividual_id\treason\n")
    for s in skipped:
        out.write("\t".join(s) + "\n")

sys.stderr.write("Complete trios: %d\n" % len(trios))
sys.stderr.write("Skipped entries: %d\n" % len(skipped))
if not trios:
    sys.exit("ERROR: no complete trios in the PED are present in the callset.")
PYEOF
    >>>

    output {
        File trios_data = "trios.data.tsv"
        File trios_manifest = "trios_manifest.tsv"
        File skipped_manifest = "skipped_manifest.tsv"
    }

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

# Build the trio VCF UPDhmm needs: prefilter each single-sample VCF to autosomal biallelic
# PASS SNPs with GT only (optionally masking low-GQ genotypes), rename each sample column to
# its PED id, merge the three, and keep only biallelic SNP sites genotyped in all three.
task MakeTrioVcf {
    input {
        String family_id
        String proband_id
        String father_id
        String mother_id
        File proband_vcf
        File proband_vcf_idx
        File father_vcf
        File father_vcf_idx
        File mother_vcf
        File mother_vcf_idx
        Array[String] autosomes
        Int? min_gq
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(proband_vcf, "GB") + size(father_vcf, "GB") + size(mother_vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail

        REGIONS="~{sep=',' autosomes}"
        MINGQ=~{select_first([min_gq, 0])}

        # Filter one single-sample VCF to autosomal biallelic PASS SNPs, GT-only, and
        # force its sample column name to the given id. Uses -t (targets, streaming) so no
        # separate index co-location is required.
        filter_one () {
            local id="$1"; local invcf="$2"
            bcftools view -t "$REGIONS" -v snps -m2 -M2 -f PASS "$invcf" -Ou > "raw.$id.bcf"
            if [ "$MINGQ" -gt 0 ]; then
                bcftools +setGT "raw.$id.bcf" -Ou -- -t q -i "FMT/GQ<$MINGQ" -n . \
                    | bcftools annotate -x INFO -Ou \
                    | bcftools annotate -x '^FORMAT/GT' -Oz -o "filt.$id.vcf.gz"
            else
                bcftools annotate -x INFO "raw.$id.bcf" -Ou \
                    | bcftools annotate -x '^FORMAT/GT' -Oz -o "filt.$id.vcf.gz"
            fi
            printf '%s\n' "$id" > "sn.$id.txt"
            bcftools reheader -s "sn.$id.txt" "filt.$id.vcf.gz" -o "renamed.$id.vcf.gz"
            mv "renamed.$id.vcf.gz" "filt.$id.vcf.gz"
            tabix -p vcf "filt.$id.vcf.gz"
        }

        filter_one "~{proband_id}" "~{proband_vcf}"
        filter_one "~{father_id}"  "~{father_vcf}"
        filter_one "~{mother_id}"  "~{mother_vcf}"

        bcftools merge -m snps -Oz -o merged.vcf.gz \
            "filt.~{proband_id}.vcf.gz" "filt.~{father_id}.vcf.gz" "filt.~{mother_id}.vcf.gz"
        tabix -p vcf merged.vcf.gz

        # Keep only clean biallelic SNP sites genotyped in all three members.
        bcftools view -m2 -M2 -v snps -e 'GT[*]="mis"' merged.vcf.gz \
            -Oz -o "~{family_id}.~{proband_id}.trio.vcf.gz"
        tabix -p vcf "~{family_id}.~{proband_id}.trio.vcf.gz"
    >>>

    output {
        File trio_vcf = "~{family_id}.~{proband_id}.trio.vcf.gz"
        File trio_vcf_idx = "~{family_id}.~{proband_id}.trio.vcf.gz.tbi"
    }

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

# Run UPDhmm on one trio VCF and emit a UPD-events TSV (with a fixed header).
task RunUPDhmmTask {
    input {
        File trio_vcf
        File trio_vcf_idx
        String family_id
        String proband_id
        String father_id
        String mother_id
        String docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: 3 * ceil(size(trio_vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail

        Rscript /opt/gnomad-lr/scripts/updhmm/updhmm.R \
            "~{trio_vcf}" \
            "~{proband_id}" \
            "~{mother_id}" \
            "~{father_id}" \
            "~{family_id}" \
            "~{family_id}.~{proband_id}.upd_events.tsv" \
            ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])}
    >>>

    output {
        File events_tsv = "~{family_id}.~{proband_id}.upd_events.tsv"
    }

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
