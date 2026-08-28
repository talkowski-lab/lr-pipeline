version 1.0

import "../utils/Structs.wdl"

workflow AnnotateExternalAFs {
    input {
        File vcf
        File vcf_index
        String prefix

        Array[File] ref_beds
        Array[String] ref_prefixes

        Array[String] contigs

        String pipeline_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_annotate_external
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_preprocess
        RuntimeAttr? runtime_attr_postprocess
    }

    scatter (contig in contigs) {
        call SubsetVcf {
            input:
                vcf = vcf,
                vcf_index = vcf_index,
                prefix = "~{prefix}.~{contig}",
                locus = contig,
                runtime_attr_override = runtime_attr_subset_vcf
          }

        call PreprocessMergedVcf {
            input:
                vcf = SubsetVcf.subset_vcf,
                vcf_index = SubsetVcf.subset_tbi,
                prefix = "~{prefix}.~{contig}.preprocessed",
                runtime_attr_override = runtime_attr_preprocess
        }

        call AnnotateExternalAFs {
            input:
                vcf = PreprocessMergedVcf.processed_vcf,
                vcf_index = PreprocessMergedVcf.processed_tbi,
                ref_beds = ref_beds,
                ref_prefixes = ref_prefixes,
                contig = contig,
                prefix = "~{prefix}.~{contig}.anno_func_AFs",
                pipeline_docker = pipeline_docker,
                runtime_attr_override = runtime_attr_annotate_external,
        }
    }

    call ConcatVcfs {
        input:
            vcfs=AnnotateExternalAFs.annotated_vcf,
            vcfs_idx=AnnotateExternalAFs.annotated_tbi,
            allow_overlaps=true,
            prefix="~{prefix}.concat",
            pipeline_docker=pipeline_docker,
            runtime_attr_override=runtime_attr_concat
    }

    call PostprocessVcf {
        input:
            vcf = ConcatVcfs.concat_vcf,
            vcf_index = ConcatVcfs.concat_vcf_index,
            prefix = "~{prefix}.annotated",
            runtime_attr_override = runtime_attr_postprocess
    }

    output {
        File external_af_annotated_vcf = PostprocessVcf.final_vcf
        File external_af_annotated_vcf_index = PostprocessVcf.final_tbi
    }
}

task SubsetVcf {
    input {
        File vcf
        File vcf_index
        String locus
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view ~{vcf} --regions ~{locus} | bgzip > ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_tbi = "~{prefix}.vcf.gz.tbi"
    }

    Int disk_size = 2*ceil(size([vcf, vcf_index], "GB")) + 1
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: disk_size,
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
        docker: "quay.io/ymostovoy/lr-utils-basic:latest"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task PreprocessMergedVcf {
    input {
        File vcf
        File vcf_index
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # DUP TO INS, retaining original info in tag
        python /dup_to_ins_LR.py ~{vcf} | \
            bcftools view -Oz > ~{prefix}.vcf.gz

        tabix ~{prefix}.vcf.gz
    >>>

    output {
        File processed_vcf = "~{prefix}.vcf.gz"
        File processed_tbi = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(10 + size(vcf, "GB")),
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "quay.io/ymostovoy/lr-process-mendelian:latest"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task AnnotateFunctionalConsequences {
    input {
        File vcf
        File vcf_index
        File noncoding_bed
        File coding_gtf
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        gatk --java-options "-Xmx~{java_mem_mb}m" SVAnnotate -V ~{vcf} --non-coding-bed ~{noncoding_bed} --protein-coding-gtf ~{coding_gtf} -O ~{prefix}.vcf

        bcftools view -Oz ~{prefix}.vcf > ~{prefix}.vcf.gz

        tabix ~{prefix}.vcf.gz
    >>>

    output {
        File anno_vcf = "~{prefix}.vcf.gz"
        File anno_tbi = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(10 + size(vcf, "GB") * 5),
        cpu_cores: 2,
        preemptible_tries: 1,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int java_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 1000 * 0.7)
    runtime {
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "quay.io/ymostovoy/lr-svannotate:latest"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CollapseDoubledDups {
    input {
        File vcf
        File vcf_index
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python /ins_to_dup_collapse.py \
            ~{vcf} \
        | bcftools view -Oz > ~{prefix}.vcf.gz

        tabix ~{prefix}.vcf.gz

    >>>

    output {
        File collapsed_vcf = "~{prefix}.vcf.gz"
        File collapsed_tbi = "~{prefix}.vcf.gz"
    }

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(10 + size(vcf, "GB") * 2),
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "quay.io/ymostovoy/lr-process-mendelian:latest"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task AnnotateExternalAFs {
    input {
        File vcf
        File vcf_index
        Array[File] ref_beds
        Array[String] ref_prefixes
        String prefix
        String contig
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eu
        
        #### split VCF by SV type into beds (SplitQueryVcf)
        svtk vcf2bed -i SVTYPE -i SVLEN ~{vcf} tmp.bed
        cut -f1-4,7-8 tmp.bed > invcf.bed
        head -1 invcf.bed > header
        cat header <(awk '{if ($5=="DEL") print}' invcf.bed )> invcf.DEL.bed
        cat header <(awk '{if ($5=="DUP") print}' invcf.bed )> invcf.DUP.bed
        cat header <(awk '{if ($5=="INS" || $5=="INS:ME" || $5=="INS:ME:ALU" || $5=="INS:ME:LINE1" || $5=="INS:ME:SVA" || $5=="ALU" || $5=="LINE1" || $5=="SVA" || $5=="HERVK" ) print}' invcf.bed )> invcf.INS.bed
        cat header <(awk '{if ($5=="INV" || $5=="CPX") print}' invcf.bed )> invcf.INV.bed
        cat header <(awk '{if ($5=="BND" || $5=="CTX") print}' invcf.bed )> invcf.BND.bed
        rm header

        #### parallelize by ref bed + ref prefix
        for x in ~{sep=" " ref_beds}; do echo $x >> tmp; done
        for x in ~{sep=" " ref_prefixes}; do echo $x >> tmp2; done
        paste tmp tmp2 > beds_prefixes

        # start with the original VCF
        cur_vcf=~{vcf}

        while read bed bed_prefix; do
            #### split ref_bed by SV type and filter by contig (SplitRefBed)
            zcat ${bed} | head -1 > header
            cat header <(zcat ${bed} | awk '{if ($1=="~{contig}" && $6=="DEL") print}') > ${bed_prefix}.~{contig}.DEL.bed
            cat header <(zcat ${bed} | awk '{if ($1=="~{contig}" && $6=="DUP") print}') > ${bed_prefix}.~{contig}.DUP.bed
            cat header <(zcat ${bed} | awk '{if ($1=="~{contig}" && $6=="INS" || $6=="INS:ME" || $6=="INS:ME:ALU" || $6=="INS:ME:LINE1" || $6=="INS:ME:SVA" || $6=="ALU" || $6=="LINE1" || $6=="SVA" || $6=="HERVK" ) print}') > ${bed_prefix}.~{contig}.INS.bed
            cat header <(zcat ${bed} | awk '{if ($1=="~{contig}" && $6=="INV" || $6=="CPX") print}' ) > ${bed_prefix}.~{contig}.INV.bed
            cat header <(zcat ${bed} | awk '{if ($1=="~{contig}" && $6=="BND" || $6=="CTX") print}' ) > ${bed_prefix}.~{contig}.BND.bed

            #### define annotations based on BED header
            echo "ALL" > pop
            sed 's/\t/\n/g' header | grep '_AF' | sed 's/_AF//g' >> pop

            #### process by SVTYPE
            for svtype in DEL DUP INS INV BND; do
                # BedtoolsClosest
                paste <(head -1 invcf.${svtype}.bed) <(head -1 ${bed_prefix}.~{contig}.${svtype}.bed) | sed -e "s/#//g" > ${svtype}.bed
                bedtools closest -wo -a <(sort -k1,1 -k2,2n invcf.${svtype}.bed) -b <(sort -k1,1 -k2,2n ${bed_prefix}.~{contig}.${svtype}.bed) >> ${svtype}.bed

                # SelectMatchedSVs
                if [ $svtype = "DEL" ] || [ $svtype = "DUP" ] || [ $svtype = "INV" ]; then
                    Rscript /opt/sv-pipeline/05_annotation/scripts/R1.bedtools_closest_CNV.R \
                        -i ${svtype}.bed \
                        -o ${bed_prefix}.~{contig}.${svtype}.comparison \
                        -p pop
                elif [ $svtype = "INS" ] || [ $svtype = "BND" ]; then
                    Rscript /opt/sv-pipeline/05_annotation/scripts/R2.bedtools_closest_INS.R \
                        -i ${svtype}.bed \
                        -o ${bed_prefix}.~{contig}.${svtype}.comparison \
                        -p pop
                fi
                cat ${bed_prefix}.~{contig}.${svtype}.comparison >> ${bed_prefix}.labeled.bed
            done

            #### ModifyVcf
        python <<CODE
        import os
        fin=os.popen(r'''zcat %s'''%("$cur_vcf"))
        header = []
        body = {}
        SVID_key = []
        for line in fin:
            pin=line.strip().split()
            if pin[0][:2]=='##':
                header.append(pin)
            else:
                body[pin[2]]=pin
                SVID_key.append(pin[2])
        header.append(['##INFO=<ID='+"${bed_prefix}"+'_SVID'+',Number=1,Type=String,Description="SVID of an overlapping event used for external allele frequency annotation.">'])

        fin.close()
        fin=open("${bed_prefix}.labeled.bed")
        colname = fin.readline().strip().split()

        for j in range(len(colname)-1):
            if j>1:
                header.append(['##INFO=<ID='+"${bed_prefix}"+'_'+colname[j]+',Number=1,Type=Float,Description="Allele frequency (for biallelic sites) or copy-state frequency (for multiallelic sites) of an overlapping event in the reference callset.">'])

        for line in fin:
            pin=line.strip().split()
            if pin[0]=='query_svid': continue
            if not pin[0] in body.keys(): continue
            info_add = ["${bed_prefix}"+'_SVID'+'='+pin[1]]
            for j in range(len(colname)-1):
                if j>1:
                    info_add.append("${bed_prefix}"+'_'+colname[j]+'='+pin[j])
            body[pin[0]][7]+=';'+';'.join(info_add)
        fin.close()

        fo=open('~{prefix}.annotated.tmp.vcf','w')
        for i in header:
            print(' '.join(i), file=fo)
        for i in SVID_key:
            print('\t'.join(body[i]), file=fo)
        fo.close()
        CODE

            bcftools view -Oz ~{prefix}.annotated.tmp.vcf > ~{prefix}.annotated.tmp.vcf.gz
            cur_vcf=~{prefix}.annotated.tmp.vcf.gz
        done < beds_prefixes
    
        bcftools sort -O z ~{prefix}.annotated.tmp.vcf.gz > ~{prefix}.annotated.vcf.gz
        tabix ~{prefix}.annotated.vcf.gz
    >>>
 
    output {
        File annotated_vcf = "~{prefix}.annotated.vcf.gz"
        File annotated_tbi = "~{prefix}.annotated.vcf.gz.tbi"
    }        

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(10 + 5*(size(vcf, "GB") + size(ref_beds, "GB"))),
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatVcfs {
    input {
        Array[File] vcfs
        Array[File]? vcfs_idx
        Boolean allow_overlaps = false
        Boolean naive = false
        Boolean sites_only = false
        Boolean sort_vcf_list = false
        String? prefix
        String pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
    # be held in memory or disk while working, potentially in a form that takes up more space)
    String outfile_name = prefix + ".vcf.gz"
    String allow_overlaps_flag = if allow_overlaps then "--allow-overlaps" else ""
    String naive_flag = if naive then "--naive" else ""
    String sites_only_command = if sites_only then "| bcftools view --no-version -G -Oz" else ""

    command <<<
        set -euo pipefail
        if ~{sort_vcf_list}; then
            VCFS=vcfs.list
            awk -F '/' '{print $NF"\t"$0}' ~{write_lines(vcfs)} | sort -k1,1V | awk '{print $2}' > $VCFS
        else
            VCFS="~{write_lines(vcfs)}"
        fi
        bcftools concat --no-version ~{allow_overlaps_flag} ~{naive_flag} -Oz --file-list $VCFS \
            ~{sites_only_command} > ~{outfile_name}
        tabix ~{outfile_name}
    >>>

    output {
      File concat_vcf = outfile_name
      File concat_vcf_index = outfile_name + ".tbi"
    }

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(10 + size(vcfs, "GB") * 2),
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task PostprocessVcf {
    input {
        File vcf
        File vcf_index
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # set variants originally called as DUP back to DUP
        python /orig_dup_to_dup.py ~{vcf} | bcftools view -Oz > ~{prefix}.vcf.gz

        tabix ~{prefix}.vcf.gz
    >>>

    output {
        File final_vcf = "~{prefix}.vcf.gz"
        File final_tbi = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        mem_gb: 4,
        disk_gb: ceil(10 + size(vcf, "GB") * 3),
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 0,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "quay.io/ymostovoy/lr-process-mendelian:latest"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
