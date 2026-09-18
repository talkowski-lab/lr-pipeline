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
import edlib
import pysam

EMPTY_SEQ = '-'
ABSENT_HAP = '.'
HEADER = [
    'trid_a', 'trid_b', 'motifs_a', 'motifs_b', 'chrom',
    'locus_start_a', 'locus_end_a', 'locus_start_b', 'locus_end_b',
    'union_start', 'union_end', 'union_length', 'overlap_start', 'overlap_end', 'overlap_length',
    'hap1_seq_a', 'hap2_seq_a', 'hap1_seq_b', 'hap2_seq_b',
    'min_edit_distance', 'max_similarity',
]


def reference_slice(record, start, end):
    """Reference bases of an interval contained in a record's own span, read off its REF allele."""
    if start > end:
        return ''
    return record['ref'][start - record['start']:end - record['start'] + 1]


def union_sequences(record, other, union_start, union_end):
    """Haplotype sequences of a record extended to the union span of the two loci.

    Bases of the union outside the record's own span are filled with reference, taken from the other
    record's REF allele. Because the two spans overlap, the union is always covered by the two spans
    together, so every padded base is available without a reference FASTA. Both sequences of a pair then
    describe the same interval and can be compared directly, with no alignment of a haplotype back onto
    reference coordinates.
    """
    left = reference_slice(other, union_start, record['start'] - 1)
    right = reference_slice(other, record['end'] + 1, union_end)
    return [left + hap + right for hap in record['haps']]


def edit_distance(seq_a, seq_b):
    if seq_a == seq_b:
        return 0
    if not seq_a or not seq_b:
        return max(len(seq_a), len(seq_b))
    return edlib.align(seq_a, seq_b, task='distance')['editDistance']


def compare_loci(seqs_a, seqs_b):
    """Total edit distance over the best unphased haplotype pairing, and a length-normalized similarity.

    Normalizing by the compared sequence lengths rather than by the reference span keeps the similarity
    within [0, 1] however far a haplotype has expanded.
    """
    if len(seqs_a) == 2 and len(seqs_b) == 2:
        pairings = [((0, 0), (1, 1)), ((0, 1), (1, 0))]
    else:
        pairings = [((i, j),) for i in range(len(seqs_a)) for j in range(len(seqs_b))]

    distance = None
    denominator = 0
    for pairing in pairings:
        pairs = [(seqs_a[i], seqs_b[j]) for i, j in pairing]
        total = sum(edit_distance(seq_a, seq_b) for seq_a, seq_b in pairs)
        if distance is None or total < distance:
            distance = total
            denominator = sum(max(len(seq_a), len(seq_b)) for seq_a, seq_b in pairs)

    similarity = 1.0 if denominator == 0 else 1 - distance / denominator
    return distance, similarity


def info_text(rec, key):
    """Render a comma-separated INFO field exactly as written in the VCF, whichever way pysam splits it."""
    value = rec.info[key]
    return value if isinstance(value, str) else ','.join(value)


def load_record(rec):
    """Collect the reference span, genotyped haplotype sequences, TRID and motifs of a TRGT record."""
    genotype = rec.samples[0]['GT']
    if not genotype or any(allele is None for allele in genotype):
        return None

    ref = rec.ref.upper()
    alleles = [ref] + [alt.upper() for alt in (rec.alts or [])]
    return {
        'chrom': rec.chrom,
        'start': rec.pos,
        'end': rec.pos + len(ref) - 1,
        'ref': ref,
        'haps': [alleles[allele] for allele in genotype],
        'nonref': any(allele > 0 for allele in genotype),
        'trid': info_text(rec, 'TRID'),
        'motifs': info_text(rec, 'MOTIFS'),
    }


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
            union_start = min(held['start'], record['start'])
            union_end = max(held['end'], record['end'])
            seqs_a = union_sequences(held, record, union_start, union_end)
            seqs_b = union_sequences(record, held, union_start, union_end)
            distance, similarity = compare_loci(seqs_a, seqs_b)
            tsv_out.write('\t'.join(str(field) for field in [
                held['trid'], record['trid'], held['motifs'], record['motifs'], record['chrom'],
                held['start'], held['end'], record['start'], record['end'],
                union_start, union_end, union_end - union_start + 1,
                overlap_start, overlap_end, overlap_end - overlap_start + 1,
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
