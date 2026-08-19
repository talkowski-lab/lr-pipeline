#!/usr/bin/env python3


import hail as hl
import numpy as np
import argparse


parser = argparse.ArgumentParser(description="Parse arguments")
parser.add_argument("-i", dest="vcf_file", help="Input VCF file")
parser.add_argument("-o", dest="output_tsv", help="Output TSV file")
parser.add_argument("--cores", dest="cores", help="CPU cores")
parser.add_argument("--mem", dest="mem", help="Memory")
parser.add_argument("--build", dest="build", help="Genome build")

args = parser.parse_args()
vcf_file = args.vcf_file
output_tsv = args.output_tsv
cores = args.cores
mem = int(np.floor(float(args.mem)))
build = args.build

hl.init(
    min_block_size=128,
    spark_conf={
        "spark.executor.cores": cores,
        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
        "spark.driver.cores": cores,
        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g",
    },
    tmp_dir="tmp",
    local_tmpdir="tmp",
)

mt = hl.import_vcf(
    vcf_file,
    force_bgz=True,
    array_elements_required=False,
    call_fields=[],
    reference_genome=build,
)

mt = hl.vep(mt, config="vep_config.json", csq=True, tolerate_parse_error=True)

table = mt.rows()

table = table.select(
    CHROM=table.locus.contig,
    POS=table.locus.position,
    REF=table.alleles[0],
    ALT=hl.delimit(table.alleles[1:], ','),
    ID=hl.or_else(table.rsid, '.'),
    CSQ=hl.or_else(hl.delimit(table.vep, ','), '.')
)

table = table.key_by().select('CHROM', 'POS', 'REF', 'ALT', 'ID', 'CSQ')

table.export(output_tsv, header=False, delimiter='\t')
