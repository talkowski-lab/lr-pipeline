#!/usr/bin/env python3


import argparse
import os

import hail as hl
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Parse arguments")
    parser.add_argument("-i", dest="vcf_file", required=True, help="Input VCF file")
    parser.add_argument("-o", dest="output_tsv", required=True, help="Output TSV file")
    parser.add_argument("--cores", dest="cores", required=True, help="CPU cores")
    parser.add_argument("--mem", dest="mem", required=True, help="Memory")
    parser.add_argument("--build", dest="build", required=True, help="Genome build")

    args = parser.parse_args()
    vcf_file = args.vcf_file
    output_tsv = args.output_tsv
    cores = args.cores
    # Cap Spark memory at the VM's physical memory, since the task's runtime memory can be below --mem
    vm_mem_gb = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1024**3
    spark_mem_mb = int(np.floor(min(float(args.mem), vm_mem_gb) * 1024 * 0.4))
    build = args.build

    hl.init(
        min_block_size=128,
        spark_conf={
            "spark.executor.cores": cores,
            "spark.executor.memory": f"{spark_mem_mb}m",
            "spark.driver.cores": cores,
            "spark.driver.memory": f"{spark_mem_mb}m",
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
        ALT=hl.delimit(table.alleles[1:], ","),
        ID=hl.or_else(table.rsid, "."),
        CSQ=hl.or_else(hl.delimit(table.vep, ","), "."),
    )

    table = table.key_by().select("CHROM", "POS", "REF", "ALT", "ID", "CSQ")

    table.export(output_tsv, header=False, delimiter="\t")


if __name__ == "__main__":
    main()
