#!/usr/bin/env python3

import argparse
import numpy
import pysam
import re


def parse_args():
    parser = argparse.ArgumentParser(description="Convert PALMER MEI calls to VCF format.")
    parser.add_argument("--palmer_calls", required=True, help="PALMER calls file")
    parser.add_argument("--palmer_tsd_reads", required=True, help="PALMER TSD reads file")
    parser.add_argument("--mei_type", required=True, help="Mobile element type (e.g., ALU, LINE1, SVA)")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--ref_fa", required=True, help="Reference genome FASTA file")
    parser.add_argument("--ref_fai", required=True, help="Reference genome FAI index file")
    parser.add_argument("--haplotype", required=True, help="Haplotype genotype (e.g., 1|0 or 0|1)")
    return parser.parse_args()


def write_vcf_header(ref_fai_path, sample):
    print("##fileformat=VCFv4.2")
    with open(ref_fai_path, "r") as ref_fai:
        for line in ref_fai:
            contig, length = line.strip().split("\t")[:2]
            print(f"##contig=<ID={contig},length={length}>")

    print('##ALT=<ID=INS,Description="Insertion of novel sequence relative to the reference">')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('##INFO=<ID=CONF_READS,Number=.,Type=Integer,Description="Number of confident supporting reads">')
    print('##INFO=<ID=ORI,Number=.,Type=String,Description="Orientation of insertion">')
    print('##INFO=<ID=ME_TYPE,Number=.,Type=String,Description="Type of mobile element">')
    print('##INFO=<ID=POLYA_LEN,Number=.,Type=Integer,Description="Length of polyA tail">')
    print('##INFO=<ID=TSD_5PRIME_LEN,Number=.,Type=Integer,Description="Length of 5prime TSD">')
    print('##INFO=<ID=TSD_3PRIME_LEN,Number=.,Type=Integer,Description="Length of 3prime TSD">')
    print('##INFO=<ID=TRANSD_LEN,Number=.,Type=Integer,Description="Length of transduction">')
    print('##INFO=<ID=INVERSION_5PRIME,Number=0,Type=Flag,Description="Whether the MEI has a 5 prime inversion">')
    print('##INFO=<ID=allele_type,Number=.,Type=String,Description="Type of allele">')
    print('##INFO=<ID=allele_length,Number=.,Type=Integer,Description="Length of allele">')
    print(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}")


def parse_tsd_reads(tsd_reads_path):
    insertion_seqs = {}
    with open(tsd_reads_path, "r") as tsd_file:
        for line in tsd_file:
            if line.startswith("cluster_id\t"):
                continue

            fields = re.split(r"\t+", line.strip())
            cluster_id = fields[0]
            if cluster_id not in insertion_seqs:
                whole_insertion_seq = fields[6]
                insertion_seqs[cluster_id] = whole_insertion_seq

    return insertion_seqs


def parse_palmer_calls(callfile_path, insertion_seqs, ref_fasta, min_conf=1):
    calls = {}
    seen_records = set()

    with open(callfile_path, "r") as callfile:
        for line in callfile:
            if line.startswith("cluster_id\t"):
                continue

            # Parse fields
            fields = re.split(r"\t+", line.strip())

            # Filter by confidence threshold
            conf = int(fields[10])
            if conf < min_conf:
                continue

            # Create a unique key from columns 2-10 and 14
            unique_key = tuple(fields[2:11] + [fields[14]])
            if unique_key in seen_records:
                continue
            seen_records.add(unique_key)

            # Parse call information
            cid = fields[0]
            chrom = fields[1]
            pos = round(
                numpy.median(
                    [int(fields[2]), int(fields[3]), int(fields[4]), int(fields[5])]
                )
            )

            # Derive REF and ALT
            insertion_seq = insertion_seqs.get(cid, "")
            ref = ref_fasta.fetch(chrom, pos - 1, pos).upper()
            alt = ref + insertion_seq.upper()

            # Set baseline VCF record fields
            call = {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "conf": fields[10],
                "ori": fields[13],
                "polyAlen": fields[14],
                "5TSDlen": fields[15],
                "3TSDlen": fields[16],
                "transDlen": fields[17],
                "5inv": int(fields[18]),
            }

            # Calculate MEI length
            tmp_start = round(numpy.median([int(fields[6]), int(fields[7])]))
            tmp_end = round(numpy.median([int(fields[8]), int(fields[9])]))
            if call["5inv"]:
                tmp_5invend = int(fields[19])
                tmp_5invstart = int(fields[20])
                call["length"] = (
                    tmp_end
                    - tmp_start
                    + int(call["transDlen"])
                    + tmp_5invend
                    - tmp_5invstart
                )
            else:
                call["length"] = tmp_end - tmp_start + int(call["transDlen"])

            calls[cid] = call

    return calls


def write_vcf_records(calls, mei_type, haplotype):
    for cid, call in calls.items():
        info_fields = (
            f"allele_length={call['length']};allele_type=ins;ME_TYPE={mei_type};"
            f"CONF_READS={call['conf']};ORI={call['ori']};POLYA_LEN={call['polyAlen']};"
            f"TSD_5PRIME_LEN={call['5TSDlen']};TSD_3PRIME_LEN={call['3TSDlen']};"
            f"TRANSD_LEN={call['transDlen']}"
        )
        if call["5inv"]:
            info_fields += ";INVERSION_5PRIME"

        print(
            f"{call['chrom']}\t{call['pos']}\t{cid}\t{call['ref']}\t{call['alt']}\t60\t.\t{info_fields}\tGT\t{haplotype}"
        )


def main():
    args = parse_args()

    # Load reference FASTA
    ref_fasta = pysam.FastaFile(args.ref_fa)

    # Parse TSD reads to get insertion sequences
    insertion_seqs = parse_tsd_reads(args.palmer_tsd_reads)

    # Parse calls to get VCF records
    calls = parse_palmer_calls(args.palmer_calls, insertion_seqs, ref_fasta)

    # Write VCF
    write_vcf_header(args.ref_fai, args.sample)
    write_vcf_records(calls, args.mei_type, args.haplotype)

    ref_fasta.close()


if __name__ == "__main__":
    main()
