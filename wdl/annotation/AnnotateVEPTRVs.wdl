version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "AnnotateVEPHail.wdl"

workflow AnnotateVEPTRVs {
    meta {
        description: [
            "This workflow annotates the alleles of tandem repeat variants (INFO/allele_type='trv') in a single-contig VCF with the Ensembl Variant Effect Predictor (VEP) (https://useast.ensembl.org/info/docs/tools/vep/index.html), producing one row per ALT allele that links its original and normalized representation to its allele-specific AC, AF and VEP annotations.",
            "Genotypes are stripped and the TRV sites are sharded by record before any splitting, so each shard is split into biallelic records, left-aligned and trimmed with `bcftools norm`, uppercased (the soft-masked reference otherwise leaves lowercase bases that VEP cannot parse) and run through VEP independently. Each split record carries the original ID suffixed with its allele index, so annotations are joined back by allele rather than by normalized coordinates, which can collide across records.",
            "AC and AF are read per allele from the Number=A INFO fields, so they must be consistent with the genotypes. Shards are concatenated in input order, so the output stays sorted by the original CHROM and POS, then by allele index."
        ]
    }

    parameter_meta {
        vcf: "Single-contig VCF to annotate."
        vcf_idx: "Index for VCF to annotate."
        records_per_shard: "Number of TRV records per shard, counted before splitting into biallelic records."
        vep_annotate_hail_python_script: "Path to the Hail script used to run VEP (defaults to this repository's copy on `main`)."
        genome_build: "Genome build to annotate against."
        vep_json_schema: "Hail type schema describing the structure of VEP's JSON output."
        normalize_check_ref: "`bcftools norm` `--check-ref` mode used when normalizing."
        ref_fa: "From references."
        ref_fai: "From references."
        ref_fa_gz: "bgzipped `ref_fa`, from references."
        ref_fai_gz: "Index for `ref_fa_gz`, from references."
        ref_vep_cache: "From references."
        annotations_tsv_vep_trv: "TSV with a header and one row per TRV ALT allele, with columns CHROM, POS, REF, ALT, ID, NORM_POS, NORM_REF, NORM_ALT, AC, AF and VEP, where VEP holds the comma-separated VEP annotations for that allele."
    }

    input {
        File vcf
        File vcf_idx
        String prefix

        Int records_per_shard = 2000
        String vep_annotate_hail_python_script = "https://raw.githubusercontent.com/talkowski-lab/lr-pipeline/main/scripts/helper/vep_annotate_hail.py"
        String genome_build = "GRCh38"
        String vep_json_schema = "Struct{allele_string:String,colocated_variants:Array[Struct{allele_string:String,clin_sig:Array[String],clin_sig_allele:String,end:Int32,id:String,phenotype_or_disease:Int32,pubmed:Array[Int32],somatic:Int32,start:Int32,strand:Int32}],context:String,end:Int32,id:String,input:String,intergenic_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],impact:String,minimised:Int32,variant_allele:String}],most_severe_consequence:String,motif_feature_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],high_inf_pos:String,impact:String,minimised:Int32,motif_feature_id:String,motif_name:String,motif_pos:Int32,motif_score_change:Float64,transcription_factors:Array[String],strand:Int32,variant_allele:String}],regulatory_feature_consequences:Array[Struct{allele_num:Int32,biotype:String,consequence_terms:Array[String],impact:String,minimised:Int32,regulatory_feature_id:String,variant_allele:String}],seq_region_name:String,start:Int32,strand:Int32,transcript_consequences:Array[Struct{allele_num:Int32,amino_acids:String,appris:String,biotype:String,canonical:Int32,ccds:String,cdna_start:Int32,cdna_end:Int32,cds_end:Int32,cds_start:Int32,codons:String,consequence_terms:Array[String],distance:Int32,domains:Array[Struct{db:String,name:String}],exon:String,flags:String,gene_id:String,gene_pheno:Int32,gene_symbol:String,gene_symbol_source:String,hgnc_id:String,hgvsc:String,hgvsp:String,hgvs_offset:Int32,impact:String,intron:String,lof:String,lof_flags:String,lof_filter:String,lof_info:String,mane_select:String,mane_plus_clinical:String,minimised:Int32,pick:Int32,mirna:Array[String],polyphen_prediction:String,polyphen_score:Float64,protein_end:Int32,protein_start:Int32,protein_id:String,sift_prediction:String,sift_score:Float64,source:String,strand:Int32,swissprot:String,transcript_id:String,trembl:String,tsl:Int32,uniparc:String,uniprot_isoform:Array[String],variant_allele:String}],variant_class:String}"
        String normalize_check_ref = "w"

        File ref_fa
        File ref_fai
        File ref_fa_gz
        File ref_fai_gz
        File ref_vep_cache

        String utils_docker
        String vep_hail_docker

        RuntimeAttr? runtime_attr_subset_trv
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_normalize_trv_alleles
        RuntimeAttr? runtime_attr_vep_annotate
        RuntimeAttr? runtime_attr_create_trv_annotation_table
        RuntimeAttr? runtime_attr_concat_shards
    }

    call Helpers.SubsetVcfByArgs {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            include_args = "INFO/allele_type=\"trv\"",
            extra_args = "-G",
            prefix = "~{prefix}.trv_sites",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_trv
    }

    call Helpers.ShardVcfByRecords {
        input:
            vcf = SubsetVcfByArgs.subset_vcf,
            vcf_idx = SubsetVcfByArgs.subset_vcf_idx,
            records_per_shard = records_per_shard,
            prefix = "~{prefix}.trv_sites",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_shard
    }

    scatter (i in range(length(ShardVcfByRecords.shards))) {
        call NormalizeTRVAlleles {
            input:
                vcf = ShardVcfByRecords.shards[i],
                vcf_idx = ShardVcfByRecords.shard_idxs[i],
                check_ref = normalize_check_ref,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.trv.shard_~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_normalize_trv_alleles
        }

        call AnnotateVEPHail.VepAnnotate {
            input:
                vcf = NormalizeTRVAlleles.normalized_vcf,
                ref_vep_cache = ref_vep_cache,
                ref_fa_gz = ref_fa_gz,
                ref_fai_gz = ref_fai_gz,
                vep_annotate_hail_python_script = vep_annotate_hail_python_script,
                genome_build = genome_build,
                vep_json_schema = vep_json_schema,
                prefix = "~{prefix}.trv.shard_~{i}",
                docker = vep_hail_docker,
                runtime_attr_override = runtime_attr_vep_annotate
        }

        call CreateTRVAnnotationTable {
            input:
                alleles_tsv = NormalizeTRVAlleles.alleles_tsv,
                vep_tsv = VepAnnotate.vep_tsv_file,
                prefix = "~{prefix}.trv.shard_~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_trv_annotation_table
        }
    }

    call Helpers.ConcatTsvs as ConcatShards {
        input:
            tsvs = CreateTRVAnnotationTable.annotations_tsv,
            sort_output = false,
            preserve_header = true,
            prefix = "~{prefix}.trv.vep_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_shards
    }

    output {
        File annotations_tsv_vep_trv = ConcatShards.concatenated_tsv
    }
}

task NormalizeTRVAlleles {
    input {
        File vcf
        File vcf_idx
        String check_ref
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Split each record into biallelic records tagged with their allele index, recording per-allele AC and AF
        python3 <<CODE
import pysam

with pysam.VariantFile("~{vcf}") as vcf_in:
    header = pysam.VariantHeader()
    for contig in vcf_in.header.contigs.values():
        header.contigs.add(contig.name, length=contig.length)

    with pysam.VariantFile("split.vcf", "w", header=header) as vcf_out, open("~{prefix}.alleles.tsv", "w") as tsv_out:
        for record in vcf_in:
            ac = record.info["AC"]
            af = record.info["AF"]
            for i, alt in enumerate(record.alts):
                allele_id = f"{record.id}_{i + 1}"
                tsv_out.write(f"{record.chrom}\t{record.pos}\t{record.ref}\t{alt}\t{record.id}\t{allele_id}\t{ac[i]}\t{af[i]:.6g}\n")
                vcf_out.write(vcf_out.new_record(
                    contig=record.chrom,
                    start=record.start,
                    stop=record.start + len(record.ref),
                    alleles=(record.ref, alt),
                    id=allele_id,
                ))
CODE

        # Left-align and trim the biallelic records, uppercasing alleles since the soft-masked reference leaks lowercase bases
        bcftools norm \
            -f ~{ref_fa} \
            -c ~{check_ref} \
            -Ou \
            split.vcf \
        | bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Ov \
        | awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {$4 = toupper($4); $5 = toupper($5); print}' \
        | bgzip -c > ~{prefix}.normalized.vcf.gz

        tabix -p vcf ~{prefix}.normalized.vcf.gz
    >>>

    output {
        File alleles_tsv = "~{prefix}.alleles.tsv"
        File normalized_vcf = "~{prefix}.normalized.vcf.gz"
        File normalized_vcf_idx = "~{prefix}.normalized.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 20 * ceil(size(vcf, "GB")) + ceil(size(ref_fa, "GB")) + 10,
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

task CreateTRVAnnotationTable {
    input {
        File alleles_tsv
        File vep_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Join normalized alleles and VEP annotations to the original alleles by allele ID, keeping input order
        python3 <<CODE
vep = {}
with open("~{vep_tsv}") as f:
    for line in f:
        _, norm_pos, norm_ref, norm_alt, allele_id, csq = line.rstrip("\n").split("\t")
        if allele_id in vep:
            raise ValueError(f"Duplicate allele ID in VEP output: {allele_id}")
        vep[allele_id] = (norm_pos, norm_ref, norm_alt, csq)

with open("~{alleles_tsv}") as f, open("~{prefix}.annotations.tsv", "w") as out:
    out.write("CHROM\tPOS\tREF\tALT\tID\tNORM_POS\tNORM_REF\tNORM_ALT\tAC\tAF\tVEP\n")
    n_rows = 0
    for line in f:
        chrom, pos, ref, alt, variant_id, allele_id, ac, af = line.rstrip("\n").split("\t")
        norm_pos, norm_ref, norm_alt, csq = vep.pop(allele_id)
        out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{variant_id}\t{norm_pos}\t{norm_ref}\t{norm_alt}\t{ac}\t{af}\t{csq}\n")
        n_rows += 1

if vep:
    raise ValueError(f"{len(vep)} VEP rows did not match an input allele")
print(f"Wrote {n_rows} allele rows")
CODE
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size([alleles_tsv, vep_tsv], "GB")) + 10,
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
