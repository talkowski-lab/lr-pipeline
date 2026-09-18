version 1.0

import "Helpers.wdl"
import "Structs.wdl"

workflow ScatterVcf {
    input {
        File file
        String prefix

        Int n_shards = 0
        Int records_per_shard = 0

        String split_vcf_hail_script = "https://raw.githubusercontent.com/talkowski-lab/annotations/refs/heads/main/scripts/split_vcf_hail.py"
        String genome_build = 'GRCh38'
        Boolean localize_vcf
        Boolean get_chromosome_sizes
        Boolean split_by_chromosome
        Boolean split_into_shards
        Boolean has_index

        String hail_docker
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_split_by_chr
        RuntimeAttr? runtime_attr_split_into_shards
    }
    
    if (split_by_chromosome) {
        if (!localize_vcf) {
            String vcf_uri = file

            if (get_chromosome_sizes) {
                call GetChromosomeSizes {
                    input:
                        vcf_file = vcf_uri,
                        has_index = select_first([has_index]),
                        docker = sv_base_mini_docker
                }
            }
        }

        Map[String, Array[String]] chromosomes_dict = {
            'GRCh38': ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"],
            'GRCh37': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y']
        }
        Array[String] chromosomes = chromosomes_dict[genome_build]

        scatter (chromosome in chromosomes) {
            if (localize_vcf) {
                call SplitByChromosome {
                    input:
                        vcf_file = file,
                        chromosome = chromosome,
                        prefix = prefix,
                        docker = sv_base_mini_docker,
                        runtime_attr_override = runtime_attr_split_by_chr
                }
            }

            if (!localize_vcf) {
                String vcf_uri = file
                # Estimate remote input size from contig length scaled by sample count
                Float input_size_ = if (get_chromosome_sizes) then select_first([GetChromosomeSizes.contig_lengths])[chromosome] * ceil(select_first([GetChromosomeSizes.n_samples])*0.001) / 1000000 else size(vcf_uri, 'GB')
                call SplitByChromosomeRemote {
                    input:
                        vcf_file = vcf_uri,
                        chromosome = chromosome,
                        prefix = prefix,
                        input_size = input_size_,
                        has_index = select_first([has_index]),
                        docker = sv_base_mini_docker,
                        runtime_attr_override = runtime_attr_split_by_chr
                }
            }

            File splitChromosomeShards = select_first([SplitByChromosome.shards, SplitByChromosomeRemote.shards])
            Float splitChromosomeContigLengths = select_first([SplitByChromosome.contig_lengths, SplitByChromosomeRemote.contig_lengths])
            Pair[File, Float] split_chromosomes = (splitChromosomeShards, splitChromosomeContigLengths)
        }
    }

    if (split_into_shards) {
        if (defined(split_chromosomes)) {
            scatter (chrom_pair in select_first([split_chromosomes])) {
                File chrom_shard = select_first([chrom_pair.left])
                Float chrom_n_records = select_first([chrom_pair.right])
                Int chrom_n_shards = ceil(chrom_n_records / select_first([records_per_shard, 0]))
                String chrom_shard_prefix = basename(chrom_shard, ".vcf.gz")

                call ExecuteScattering as scatterChromosomes {
                    input:
                        vcf_file = chrom_shard,
                        split_vcf_hail_script = split_vcf_hail_script,
                        n_shards = chrom_n_shards,
                        records_per_shard = 0,
                        prefix = chrom_shard_prefix,
                        genome_build = genome_build,
                        docker = hail_docker,
                        runtime_attr_override = runtime_attr_split_into_shards
                }
            }
            Array[File] chromosome_shards = flatten(scatterChromosomes.shards)
        }
        
        if (!defined(split_chromosomes)) {
            if (localize_vcf) {
                call ExecuteScattering {
                    input:
                        vcf_file = file,
                        split_vcf_hail_script = split_vcf_hail_script,
                        n_shards = select_first([n_shards]),
                        records_per_shard = select_first([records_per_shard, 0]),
                        prefix = prefix,
                        genome_build = genome_build,
                        docker = hail_docker,
                        runtime_attr_override = runtime_attr_split_into_shards
                    }
            }

            if (!localize_vcf) {
                String mt_uri = file

                call Helpers.GetHailMTSize as getHailMTSize {
                    input:
                        mt_uri = mt_uri,
                        docker = hail_docker
                }
                
                call ScatterVcfRemote {
                    input:
                        vcf_file = mt_uri,
                        input_size = getHailMTSize.mt_size,
                        split_vcf_hail_script = split_vcf_hail_script,
                        n_shards = select_first([n_shards]),
                        records_per_shard = select_first([records_per_shard, 0]),
                        prefix = prefix,
                        genome_build = genome_build,
                        docker = hail_docker,
                        runtime_attr_override = runtime_attr_split_into_shards
                }
            }
        }
    }    

    output {
        Array[File] vcf_shards = select_first([ExecuteScattering.shards, ScatterVcfRemote.shards, chromosome_shards, splitChromosomeShards, [file]])
    }
}   

task GetChromosomeSizes {
    input {
        String vcf_file
        Boolean has_index
        String docker
        RuntimeAttr? runtime_attr_override
    }
    
    Float base_disk_gb = 10.0

    command <<<
        set -euo pipefail

        if [[ "~{has_index}" == "false" ]]; then
            mkfifo /tmp/token_fifo
            ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
            HTS_AUTH_LOCATION=/tmp/token_fifo tabix --verbosity 3 ~{vcf_file}
        fi;
        
        export GCS_OAUTH_TOKEN=`/google-cloud-sdk/bin/gcloud auth application-default print-access-token`
        
        bcftools index -s ~{vcf_file} | cut -f1,3 > contig_lengths.txt
        bcftools query -l ~{vcf_file} | wc -l > n_samples.txt
    >>>

    output {
        Float n_samples = read_lines('n_samples.txt')[0]
        Map[String, Float] contig_lengths = read_map('contig_lengths.txt')
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb),
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

task SplitByChromosomeRemote {
    input {
        String vcf_file
        String chromosome
        String prefix
        Float input_size
        Boolean has_index
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    command <<<
        set -euo pipefail
        
        mkfifo /tmp/token_fifo
        ( while true ; do curl -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/token_fifo ; done ) &
        
        HTS_AUTH_LOCATION=/tmp/token_fifo tabix --verbosity 3 -h ~{vcf_file} ~{chromosome} | bgzip -c > ~{prefix}."~{chromosome}".vcf.gz
        
        tabix -p vcf ~{prefix}."~{chromosome}".vcf.gz
        # Count records in the contig
        HTS_AUTH_LOCATION=/tmp/token_fifo bcftools index -n ~{prefix}."~{chromosome}".vcf.gz > contig_length.txt
    >>>

    output {
        File shards = "~{prefix}.~{chromosome}.vcf.gz"
        File shards_idx = "~{prefix}.~{chromosome}.vcf.gz.tbi"
        Float contig_lengths = read_lines('contig_length.txt')[0]
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
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

task SplitByChromosome {
    input {
        File vcf_file
        String chromosome
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 2.0

    command <<<
        set -euo pipefail

        tabix --verbosity 3 ~{vcf_file}
        
        tabix --verbosity 3 -h ~{vcf_file} ~{chromosome} | bgzip -c > ~{prefix}."~{chromosome}".vcf.gz
        
        tabix -p vcf ~{prefix}."~{chromosome}".vcf.gz
        
        # Count records in the contig
        HTS_AUTH_LOCATION=/tmp/token_fifo bcftools index -n ~{prefix}."~{chromosome}".vcf.gz > contig_length.txt
    >>>

    output {
        File shards = "~{prefix}.~{chromosome}.vcf.gz"
        File shards_idx = "~{prefix}.~{chromosome}.vcf.gz.tbi"
        Float contig_lengths = read_lines('contig_length.txt')[0]
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
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

task ExecuteScattering {
    input {
        File vcf_file
        Int n_shards
        Int records_per_shard
        String prefix
        String split_vcf_hail_script
        String genome_build
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    command <<<
        set -euo pipefail
        
        curl  ~{split_vcf_hail_script} > split_vcf.py
        
        python3 split_vcf.py ~{vcf_file} ~{n_shards} ~{records_per_shard} ~{prefix} ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} ~{genome_build}
        
        for file in $(ls ~{prefix}.vcf.bgz | grep '.bgz'); do
            shard_num=$(echo $file | cut -d '-' -f2);
            mv ~{prefix}.vcf.bgz/$file ~{prefix}.shard_"$shard_num".vcf.bgz
        done
    >>>

    output {
        Array[File] shards = glob("~{prefix}.shard_*.vcf.bgz")
        Array[String] shards_string = glob("~{prefix}.shard_*.vcf.bgz")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
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

task ScatterVcfRemote {
    input {
        String vcf_file
        Float input_size
        Int n_shards
        Int records_per_shard
        String prefix
        String split_vcf_hail_script
        String genome_build
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    command <<<
        set -euo pipefail
        
        curl  ~{split_vcf_hail_script} > split_vcf.py
        
        python3 split_vcf.py ~{vcf_file} ~{n_shards} ~{records_per_shard} ~{prefix} ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} ~{genome_build}
        
        for file in $(ls ~{prefix}.vcf.bgz | grep '.bgz'); do
            shard_num=$(echo $file | cut -d '-' -f2);
            mv ~{prefix}.vcf.bgz/$file ~{prefix}.shard_"$shard_num".vcf.bgz
        done
    >>>

    output {
        Array[File] shards = glob("~{prefix}.shard_*.vcf.bgz")
        Array[String] shards_string = glob("~{prefix}.shard_*.vcf.bgz")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
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
