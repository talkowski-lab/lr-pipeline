version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow EvaluateOverlappingTRLoci {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_evaluate_overlapping_loci
        RuntimeAttr? runtime_attr_concat_tsvs
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call EvaluateOverlappingTRLociForContig {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_evaluate_overlapping_loci
        }
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = EvaluateOverlappingTRLociForContig.overlapping_loci_tsv,
            sort_output = false,
            preserve_header = true,
            prefix = "~{prefix}.overlapping_tr_loci",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_tsvs
    }

    output {
        File overlapping_loci_tsv = ConcatTsvs.concatenated_tsv
    }
}

task EvaluateOverlappingTRLociForContig {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import re

import edlib
import pysam

CIGAR_RE = re.compile(r'(\d+)([=XIDM])')
EMPTY_SEQ = '-'
ABSENT_HAP = '.'
HEADER = [
    'trid_a', 'trid_b', 'chrom', 'overlap_start', 'overlap_end', 'overlap_length',
    'hap1_seq_a', 'hap2_seq_a', 'hap1_seq_b', 'hap2_seq_b',
    'min_edit_distance', 'max_similarity',
]


def build_ref_cuts(ref, hap):
    """Map every reference base boundary of a TR record onto an offset in one of its haplotype sequences.

    Returns cuts[0..len(ref)], where hap[cuts[i]:cuts[j]] is the haplotype sequence spanning reference
    bases i through j - 1. Inserted bases are assigned to the reference base they follow, with a leading
    insertion assigned to the first reference base, so the per-base spans tile the haplotype exactly.
    """
    if hap == ref:
        return list(range(len(ref) + 1))

    cuts = [0] * (len(ref) + 1)
    ref_offset = 0
    hap_offset = 0
    for count, op in CIGAR_RE.findall(edlib.align(hap, ref, task='path')['cigar']):
        count = int(count)
        if op == 'I':
            hap_offset += count
            if ref_offset > 0:
                cuts[ref_offset] = hap_offset
        elif op == 'D':
            for _ in range(count):
                ref_offset += 1
                cuts[ref_offset] = hap_offset
        else:
            for _ in range(count):
                hap_offset += 1
                ref_offset += 1
                cuts[ref_offset] = hap_offset

    return cuts


def edit_distance(seq_a, seq_b):
    if seq_a == seq_b:
        return 0
    return edlib.align(seq_a, seq_b, task='distance')['editDistance']


def min_edit_distance(seqs_a, seqs_b):
    """Total edit distance summed over both haplotype pairs, minimized over the two unphased assignments.

    Also returns the number of haplotype pairs compared, which scales the similarity denominator.
    """
    if len(seqs_a) == 2 and len(seqs_b) == 2:
        as_called = edit_distance(seqs_a[0], seqs_b[0]) + edit_distance(seqs_a[1], seqs_b[1])
        swapped = edit_distance(seqs_a[0], seqs_b[1]) + edit_distance(seqs_a[1], seqs_b[0])
        return min(as_called, swapped), 2
    return min(edit_distance(seq_a, seq_b) for seq_a in seqs_a for seq_b in seqs_b), 1


def load_record(rec):
    """Collect the reference span, genotyped haplotype sequences and TRID of a TRGT record."""
    genotype = rec.samples[0]['GT']
    if not genotype or any(allele is None for allele in genotype):
        return None

    ref = rec.ref.upper()
    alleles = [ref] + [alt.upper() for alt in (rec.alts or [])]
    trid = rec.info['TRID']
    return {
        'chrom': rec.chrom,
        'start': rec.pos,
        'end': rec.pos + len(ref) - 1,
        'ref': ref,
        'haps': [alleles[allele] for allele in genotype],
        'nonref': any(allele > 0 for allele in genotype),
        'trid': trid if isinstance(trid, str) else ','.join(trid),
        'cuts': None,
    }


def hap_sequences(record, overlap_start, overlap_end):
    """Slice out the haplotype sequences of a record spanning an overlapping reference interval."""
    if record['cuts'] is None:
        record['cuts'] = [build_ref_cuts(record['ref'], hap) for hap in record['haps']]

    start_offset = overlap_start - record['start']
    end_offset = overlap_end - record['start'] + 1
    return [hap[cuts[start_offset]:cuts[end_offset]] for hap, cuts in zip(record['haps'], record['cuts'])]


def format_haps(seqs):
    formatted = [seq if seq else EMPTY_SEQ for seq in seqs]
    return formatted + [ABSENT_HAP] * (2 - len(formatted))


records = 0
pairs = 0
reference_pairs = 0
with pysam.VariantFile("~{vcf}") as vcf_in, open("~{prefix}.overlapping_tr_loci.tsv", 'w') as tsv_out:
    tsv_out.write('\t'.join(HEADER) + '\n')

    # Records arrive position-sorted, so every retained record still spanning the current position overlaps it.
    active = []
    for rec in vcf_in:
        record = load_record(rec)
        if record is None:
            continue
        records += 1

        active = [held for held in active if held['end'] >= record['start']]
        for held in active:
            if not (held['nonref'] or record['nonref']):
                reference_pairs += 1
                continue

            overlap_start = max(held['start'], record['start'])
            overlap_end = min(held['end'], record['end'])
            overlap_length = overlap_end - overlap_start + 1
            seqs_a = hap_sequences(held, overlap_start, overlap_end)
            seqs_b = hap_sequences(record, overlap_start, overlap_end)
            distance, ploidy = min_edit_distance(seqs_a, seqs_b)
            similarity = 1 - distance / (ploidy * overlap_length)
            tsv_out.write('\t'.join(str(field) for field in [
                held['trid'], record['trid'], record['chrom'], overlap_start, overlap_end, overlap_length,
                *format_haps(seqs_a), *format_haps(seqs_b), distance, f'{similarity:.4f}',
            ]) + '\n')
            pairs += 1

        active.append(record)

print(f'Evaluated {pairs} overlapping locus pairs with a non-reference call across {records} genotyped records, '
      f'skipping {reference_pairs} pairs where both loci were called reference')
CODE
    >>>

    output {
        File overlapping_loci_tsv = "~{prefix}.overlapping_tr_loci.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
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
