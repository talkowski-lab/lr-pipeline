version 1.0

import "Structs.wdl"

task ConvertSymbolicAllelesToSequence {
    input {
        File vcf
        File vcf_idx
        Boolean drop_unsupported_symbolic_alleles
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools norm \
            -m -any \
            -Oz -o split.vcf.gz \
            ~{vcf}

        python3 <<'CODE'
import pysam


CONVERTIBLE_SYMBOLIC_ALTS = {"<DEL>", "<DUP>"}

drop_unsupported_symbolic_alleles = ~{true="True" false="False" drop_unsupported_symbolic_alleles}


def is_symbolic(alt):
    return alt.startswith("<") and alt.endswith(">")


def fail(record, message):
    raise ValueError(f"{record.chrom}:{record.pos}: {message}")


def get_end(record):
    end = record.stop
    if end > record.start + 1:
        return end

    svlen = record.info.get("SVLEN")
    if svlen is None:
        fail(record, "symbolic <DEL> requires END or SVLEN")
    if isinstance(svlen, (list, tuple)):
        if len(svlen) != 1:
            fail(record, "symbolic <DEL> requires one SVLEN value")
        svlen = svlen[0]
    try:
        svlen = int(svlen)
    except (TypeError, ValueError):
        fail(record, f"invalid SVLEN value {svlen!r}")
    if svlen == 0:
        fail(record, "symbolic <DEL> requires a nonzero END or SVLEN")
    return record.start + 1 + abs(svlen)


def get_symbolic_length(record):
    svlen = record.info.get("SVLEN")
    if svlen is not None:
        if isinstance(svlen, (list, tuple)):
            if len(svlen) != 1:
                fail(record, "symbolic <DUP> requires one SVLEN value")
            svlen = svlen[0]
        try:
            svlen = int(svlen)
        except (TypeError, ValueError):
            fail(record, f"invalid SVLEN value {svlen!r}")
        if svlen != 0:
            return abs(svlen)

    length = record.stop - record.pos
    if length > 0:
        return length
    fail(record, "symbolic <DUP> requires a nonzero SVLEN or END")


vcf_in = pysam.VariantFile("split.vcf.gz")
header = vcf_in.header.copy()
reference = pysam.FastaFile("~{ref_fa}", filepath_index="~{ref_fai}")
vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=header)

for record in vcf_in:
    record.translate(header)
    alt = record.alts[0]
    if not is_symbolic(alt):
        vcf_out.write(record)
        continue
    if alt not in CONVERTIBLE_SYMBOLIC_ALTS:
        if not drop_unsupported_symbolic_alleles:
            fail(record, f"unsupported symbolic ALT {alt}")
        continue

    symbolic_length = get_symbolic_length(record) if alt == "<DUP>" else None
    end = record.start + 1 + symbolic_length if alt == "<DUP>" else get_end(record)
    if end <= record.start + 1:
        fail(record, "symbolic <DEL> or <DUP> must span at least one non-anchor base")
    try:
        contig_length = reference.get_reference_length(record.chrom)
    except ValueError:
        fail(record, f"reference does not contain contig {record.chrom}")
    if end > contig_length:
        fail(record, f"END {end} exceeds reference contig length {contig_length}")

    sequence = reference.fetch(record.chrom, record.start, end).upper()
    if len(sequence) != end - record.start:
        fail(record, "reference sequence length does not match END")
    anchor = sequence[0]
    if alt == "<DEL>":
        record.alleles = (sequence, anchor)
    else:
        record.alleles = (anchor, anchor + sequence[1:])
    record.stop = record.start + len(record.ref)
    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
reference.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File converted_vcf = "~{prefix}.vcf.gz"
        File converted_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(ref_fa, "GB")) + 5,
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
