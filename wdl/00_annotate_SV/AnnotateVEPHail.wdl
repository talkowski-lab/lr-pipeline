version 1.0

import "../utils/Helpers.wdl"
import "../utils/ScatterVcf.wdl"
import "../utils/Structs.wdl"

workflow AnnotateVEPHail {
    input {
        File vcf
        File vcf_idx
        File ref_fa
        File ref_fai
        File ref_fa_gz
        File ref_fai_gz
        File ref_vep_cache
        String prefix

        String? subset_vcf_string
        String split_vcf_hail_script = "https://raw.githubusercontent.com/talkowski-lab/lr-annotation/main/scripts/vep/split_vcf_hail.py"
        String vep_annotate_hail_python_script = "https://raw.githubusercontent.com/talkowski-lab/lr-annotation/main/scripts/vep/vep_annotate_hail.py"
        String genome_build = "GRCh38"
        String vep_json_schema = "Struct{allele_string:String,colocated_variants:Array[Struct{allele_string:String,clin_sig:Array[String],clin_sig_allele:String,end:Int32,id:String,phenotype_or_disease:Int32,pubmed:Array[Int32],somatic:Int32,start:Int32,strand:Int32}],context:String,end:Int32,id:String,input:String,intergenic_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],impact:String,minimised:Int32,variant_allele:String}],most_severe_consequence:String,motif_feature_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],high_inf_pos:String,impact:String,minimised:Int32,motif_feature_id:String,motif_name:String,motif_pos:Int32,motif_score_change:Float64,transcription_factors:Array[String],strand:Int32,variant_allele:String}],regulatory_feature_consequences:Array[Struct{allele_num:Int32,biotype:String,consequence_terms:Array[String],impact:String,minimised:Int32,regulatory_feature_id:String,variant_allele:String}],seq_region_name:String,start:Int32,strand:Int32,transcript_consequences:Array[Struct{allele_num:Int32,amino_acids:String,appris:String,biotype:String,canonical:Int32,ccds:String,cdna_start:Int32,cdna_end:Int32,cds_end:Int32,cds_start:Int32,codons:String,consequence_terms:Array[String],distance:Int32,domains:Array[Struct{db:String,name:String}],exon:String,flags:String,gene_id:String,gene_pheno:Int32,gene_symbol:String,gene_symbol_source:String,hgnc_id:String,hgvsc:String,hgvsp:String,hgvs_offset:Int32,impact:String,intron:String,lof:String,lof_flags:String,lof_filter:String,lof_info:String,mane_select:String,mane_plus_clinical:String,minimised:Int32,pick:Int32,mirna:Array[String],polyphen_prediction:String,polyphen_score:Float64,protein_end:Int32,protein_start:Int32,protein_id:String,sift_prediction:String,sift_score:Float64,source:String,strand:Int32,swissprot:String,transcript_id:String,trembl:String,tsl:Int32,uniparc:String,uniprot_isoform:Array[String],variant_allele:String}],variant_class:String}"
        String normalize_check_ref = "w"
        Boolean normalize_vcf = false
        Boolean localize_vcf = true
        Boolean has_index = true
        Boolean get_chromosome_sizes = false
        Boolean split_by_chromosome = false
        Boolean split_into_shards = false

        String vep_hail_docker
        String hail_docker
        String sv_base_mini_docker
        String utils_docker

        RuntimeAttr? runtime_attr_strip_genotypes
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_index_shard
        RuntimeAttr? runtime_attr_extract_coords
        RuntimeAttr? runtime_attr_normalize
        RuntimeAttr? runtime_attr_vep_annotate
        RuntimeAttr? runtime_attr_collapse_multiallelics
        RuntimeAttr? runtime_attr_restore_alleles
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_split_by_chr
        RuntimeAttr? runtime_attr_split_into_shards
    }

    call Helpers.StripGenotypes {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = "~{prefix}.stripped",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_strip_genotypes
    }

    if (defined(subset_vcf_string)) {
        call Helpers.SubsetVcfByArgs {
            input:
                vcf = StripGenotypes.stripped_vcf,
                vcf_idx = StripGenotypes.stripped_vcf_idx,
                extra_args = subset_vcf_string,
                prefix = "~{prefix}.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }
    }

    call ScatterVcf.ScatterVcf {
        input:
            file = select_first([SubsetVcfByArgs.subset_vcf, StripGenotypes.stripped_vcf]),
            split_vcf_hail_script = split_vcf_hail_script,
            prefix = "~{prefix}.scattered",
            genome_build = genome_build,
            localize_vcf = localize_vcf,
            get_chromosome_sizes = get_chromosome_sizes,
            split_by_chromosome = split_by_chromosome,
            split_into_shards = split_into_shards,
            has_index = has_index,
            hail_docker = hail_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_split_by_chr = runtime_attr_split_by_chr,
            runtime_attr_split_into_shards = runtime_attr_split_into_shards
    }

    scatter (vcf_shard in ScatterVcf.vcf_shards) {
        String shard_prefix = sub(sub(basename(vcf_shard), "\\.vcf\\.bgz$", ""), "\\.vcf\\.gz$", "")

        if (normalize_vcf) {
            call Helpers.IndexVcf as IndexShard {
                input:
                    vcf = vcf_shard,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_index_shard
            }

            call Helpers.ExtractVcfCoords as ExtractShardCoords {
                input:
                    vcf = vcf_shard,
                    vcf_idx = select_first([IndexShard.vcf_idx]),
                    prefix = "~{shard_prefix}.coords",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_coords
            }

            call Helpers.NormalizeVcf as NormalizeVcfShard {
                input:
                    vcf = vcf_shard,
                    vcf_idx = select_first([IndexShard.vcf_idx]),
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    check_ref = normalize_check_ref,
                    prefix = "~{shard_prefix}.normalized",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_normalize
            }
        }

        call VepAnnotate {
            input:
                vcf = select_first([NormalizeVcfShard.normalized_vcf, vcf_shard]),
                ref_vep_cache = ref_vep_cache,
                ref_fa_gz = ref_fa_gz,
                ref_fai_gz = ref_fai_gz,
                vep_annotate_hail_python_script = vep_annotate_hail_python_script,
                genome_build = genome_build,
                vep_json_schema = vep_json_schema,
                docker = vep_hail_docker,
                runtime_attr_override = runtime_attr_vep_annotate
        }

        if (normalize_vcf) {
            call Helpers.CollapseMultiallelics as CollapseShardMultiallelics {
                input:
                    tsv = VepAnnotate.vep_tsv_file,
                    prefix = "~{shard_prefix}.collapsed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_collapse_multiallelics
            }

            call Helpers.RestoreOriginalAlleles as RestoreShardAlleles {
                input:
                    tsv = CollapseShardMultiallelics.collapsed_tsv,
                    coords_tsv = select_first([ExtractShardCoords.coords_tsv]),
                    prefix = "~{shard_prefix}.restored",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_restore_alleles
            }
        }
    }

    call Helpers.ConcatTsvs as ConcatShards {
        input:
            tsvs = if normalize_vcf then select_all(RestoreShardAlleles.restored_tsv) else VepAnnotate.vep_tsv_file,
            sort_output = false,
            prefix = "~{prefix}.vep_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_shards
    }

    output {
        File annotations_tsv_vep = ConcatShards.concatenated_tsv
    }
}

task VepAnnotate {
    input {
        File vcf
        File ref_vep_cache
        File ref_fa_gz
        File ref_fai_gz
        String vep_annotate_hail_python_script
        String genome_build
        String vep_json_schema
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String filename = basename(vcf)
    String prefix = if (sub(filename, "\\.gz", "")!=filename) then basename(vcf, ".vcf.gz") else basename(vcf, ".vcf.bgz")

    command <<<
        set -euo pipefail

        mkdir -p cache
        tar xzf ~{ref_vep_cache} -C ./cache

        dir_cache_path=$(pwd)/cache

        echo '{"command": [
        "/opt/vep/src/ensembl-vep/vep",
            "--format", "vcf",
            "__OUTPUT_FORMAT_FLAG__",
            "--everything",
            "--allele_number",
            "--no_stats",
            "--cache",
            "--offline",
            "--minimal",
            "--assembly", "~{genome_build}",
            "--merged",
            "--fasta", "~{ref_fa_gz}",
            "--dir_cache", "'$dir_cache_path'",
            "--flag_pick",
            "-o", "STDOUT"
        ],
        "env": {},
        "vep_json_schema": "~{vep_json_schema}"
        }' > vep_config.json

        curl ~{vep_annotate_hail_python_script} > vep_annotate.py

        python3 vep_annotate.py \
            -i ~{vcf} \
            -o ~{prefix}.vep.tsv \
            --cores ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} \
            --build ~{genome_build}

        cp $(ls . | grep hail*.log) hail_log.txt
    >>>

    output {
        File vep_tsv_file = "~{prefix}.vep.tsv"
        File hail_log = "hail_log.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(ref_vep_cache, "GB")) + 25,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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
