#!/usr/bin/env python3


import argparse

import hail as hl
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Shard a VCF or Hail MatrixTable into bgzipped VCF shards.")
    parser.add_argument("--vcf", required=True, help="Input VCF, or Hail MatrixTable ending in .mt")
    parser.add_argument("--n-shards", required=True, type=int, help="Number of shards, or 0 to keep Hail's partitioning")
    parser.add_argument("--records-per-shard", required=True, type=int, help="Target records per shard, or 0 to use --n-shards")
    parser.add_argument("--prefix", required=True, help="Output prefix; shards are written under <prefix>.vcf.bgz/")
    parser.add_argument("--cores", required=True, help="CPU cores")
    parser.add_argument("--mem", required=True, type=float, help="Memory in GB")
    parser.add_argument("--build", required=True, help="Reference genome build")
    args = parser.parse_args()

    n_shards = args.n_shards
    mem = int(np.floor(args.mem))

    hl.init(
        min_block_size=128,
        spark_conf={
            "spark.executor.cores": args.cores,
            "spark.executor.memory": f"{int(np.floor(mem * 0.4))}g",
            "spark.driver.cores": args.cores,
            "spark.driver.memory": f"{int(np.floor(mem * 0.4))}g",
        },
        tmp_dir="tmp",
        local_tmpdir="tmp",
    )

    is_vcf = args.vcf.split(".")[-1] != "mt"
    header = None
    if not is_vcf:
        mt = hl.read_matrix_table(args.vcf, n_partitions=n_shards if n_shards != 0 else None)
    else:
        mt = hl.import_vcf(
            args.vcf,
            force_bgz=True,
            array_elements_required=False,
            call_fields=[],
            reference_genome=args.build,
            find_replace=("nan", "."),
            n_partitions=n_shards if n_shards != 0 else None,
        )
        header = hl.get_vcf_metadata(args.vcf)

    # Sites-only inputs, such as AnnotateVEPHail's stripped VCF, have no GT
    if "GT" in mt.entry:
        # Normalize haploid calls (e.g. chrY) to diploid
        mt = mt.annotate_entries(GT=hl.if_else(mt.GT.ploidy == 1, hl.call(mt.GT[0], mt.GT[0]), mt.GT))

    if args.records_per_shard != 0:
        tot_num_records = mt.count_rows()
        n_shards = int(np.ceil(tot_num_records / args.records_per_shard))
        mt = mt.repartition(n_shards)

    tot_num_records = mt.count_rows()
    if tot_num_records > 0:
        hl.export_vcf(
            mt,
            output=args.prefix + ".vcf.bgz",
            parallel="header_per_shard",
            metadata=header,
        )


if __name__ == "__main__":
    main()
