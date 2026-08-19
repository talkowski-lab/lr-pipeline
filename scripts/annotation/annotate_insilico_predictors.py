import argparse
import numpy as np
import hail as hl

parser = argparse.ArgumentParser(description="Annotate VCF with GnomAD in silico predictors using Hail HTs.")

parser.add_argument('--build', required=True, help='Reference genome build (e.g., GRCh37 or GRCh38).')
parser.add_argument('--cores', required=True, help='CPU cores.')
parser.add_argument('--mem', required=True, help='Memory in GB.')

parser.add_argument('--cadd_ht', required=True, help='CADD Hail Table path.')
parser.add_argument('--pangolin_ht', required=True, help='Pangolin Hail Table path.')
parser.add_argument('--phylop_ht', required=True, help='PhyloP Hail Table path.')
parser.add_argument('--revel_ht', required=True, help='REVEL Hail Table path.')
parser.add_argument('--spliceai_ht', required=True, help='SpliceAI Hail Table path.')

parser.add_argument('--vcf', required=True, help='Input VCF path.')
parser.add_argument('--output_tsv', required=True, help='Output annotations TSV path.')

args = parser.parse_args()


# Helper to format output
def fmt(val):
    return hl.if_else(hl.is_defined(val), hl.format('%.6g', val), '.')


# Reference genome
build = args.build
cores = args.cores
mem = int(np.floor(float(args.mem)))

# Annotation HT paths
cadd_ht_uri = args.cadd_ht
pangolin_ht_uri = args.pangolin_ht
phylop_ht_uri = args.phylop_ht
revel_ht_uri = args.revel_ht
spliceai_ht_uri = args.spliceai_ht

# VCF paths
vcf_uri = args.vcf
output_tsv_uri = args.output_tsv

hl.init(
    default_reference=build,
    min_block_size=128,
    spark_conf={
        "spark.executor.cores": cores,
        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
        "spark.driver.cores": cores,
        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g",
    },
    tmp_dir="tmp", local_tmpdir="tmp"
)

# Load all predictor HTs
cadd_ht = hl.read_table(cadd_ht_uri)
pangolin_ht = hl.read_table(pangolin_ht_uri)
phylop_ht = hl.read_table(phylop_ht_uri)
revel_ht = hl.read_table(revel_ht_uri)
spliceai_ht = hl.read_table(spliceai_ht_uri)

# Import VCF
mt = hl.import_vcf(vcf_uri, force_bgz=True, array_elements_required=False, reference_genome=build)

# Annotate rows with predictor scores
mt = mt.annotate_rows(
    info=mt.info.annotate(
        CADD_raw_score=cadd_ht[mt.locus, mt.alleles].cadd.raw_score,
        CADD_PHRED_score=cadd_ht[mt.locus, mt.alleles].cadd.phred,
        pangolin_largest_ds=pangolin_ht[mt.locus, mt.alleles].pangolin_largest_ds,
        revel_max=revel_ht[mt.locus, mt.alleles].revel_max,
        phylop=phylop_ht[mt.locus].phylop,
        spliceai_ds_max=spliceai_ht[mt.locus, mt.alleles].spliceai_ds_max
    )
)

# Filter to variants with at least one predictor annotation
mt = mt.filter_rows(
    hl.is_defined(mt.info.CADD_raw_score) |
    hl.is_defined(mt.info.CADD_PHRED_score) |
    hl.is_defined(mt.info.pangolin_largest_ds) |
    hl.is_defined(mt.info.revel_max) |
    hl.is_defined(mt.info.phylop) |
    hl.is_defined(mt.info.spliceai_ds_max)
)


# Export annotations as TSV
ht = mt.rows()
ht = ht.select(
    CHROM=ht.locus.contig,
    POS=ht.locus.position,
    REF=ht.alleles[0],
    ALT=hl.delimit(ht.alleles[1:], ','),
    ID=hl.if_else(hl.is_defined(ht.rsid), ht.rsid, '.'),
    CADD_raw_score=fmt(ht.info.CADD_raw_score),
    CADD_PHRED_score=fmt(ht.info.CADD_PHRED_score),
    pangolin_largest_ds=fmt(ht.info.pangolin_largest_ds),
    revel_max=fmt(ht.info.revel_max),
    phylop=fmt(ht.info.phylop),
    spliceai_ds_max=fmt(ht.info.spliceai_ds_max)
)
ht = ht.key_by()
ht = ht.drop('locus', 'alleles')
ht.export(output_tsv_uri, delimiter='\t', header=False)
