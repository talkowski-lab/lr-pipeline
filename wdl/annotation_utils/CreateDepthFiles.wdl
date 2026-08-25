version 1.0

import "../utils/Structs.wdl"

workflow CreateDepthFiles {
    input {
        Array[String] sample_ids
        Array[File] mosdepth_bed_files
        String prefix

        Int bin_size = 100
        Int ploidy_bin_size = 1000000
        Int random_seed = 42

        String utils_docker

        RuntimeAttr? runtime_attr_bin_mosdepth
        RuntimeAttr? runtime_attr_set_bins
        RuntimeAttr? runtime_attr_make_bincov_columns
        RuntimeAttr? runtime_attr_zpaste
        RuntimeAttr? runtime_attr_build_ploidy_matrix
        RuntimeAttr? runtime_attr_ploidy_score
        RuntimeAttr? runtime_attr_median_cov
    }

    scatter (i in range(length(sample_ids))) {
        call BinMosdepth {
            input:
                mosdepth_bed_file = mosdepth_bed_files[i],
                bin_size = bin_size,
                prefix = sample_ids[i],
                docker = utils_docker,
                runtime_attr_override = runtime_attr_bin_mosdepth
        }
    }

    call SetBins {
        input:
            count_file = BinMosdepth.binned_bed[0],
            binsize = bin_size,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_set_bins
    }

    scatter (i in range(length(sample_ids))) {
        call MakeBincovMatrixColumns {
            input:
                count_file = BinMosdepth.binned_bed[i],
                sample = sample_ids[i],
                binsize = bin_size,
                bin_locs = SetBins.bin_locs,
                prefix = sample_ids[i],
                docker = utils_docker,
                runtime_attr_override = runtime_attr_make_bincov_columns
        }
    }

    call ZPaste {
        input:
            column_files = flatten([[SetBins.bin_locs], MakeBincovMatrixColumns.bincov_bed]),
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_zpaste
    }

    call BuildPloidyMatrix {
        input:
            bincov_matrix = ZPaste.bincov_matrix,
            bin_size = ploidy_bin_size,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_build_ploidy_matrix
    }

    call PloidyScore {
        input:
            ploidy_matrix = BuildPloidyMatrix.ploidy_matrix,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_ploidy_score
    }

    call MedianCov {
        input:
            bincov_matrix = ZPaste.bincov_matrix,
            prefix = prefix,
            random_seed = random_seed,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_median_cov
    }

    output {
        File binned_coverage = ZPaste.bincov_matrix
        File binned_coverage_idx = ZPaste.bincov_matrix_idx
        File median_coverage = MedianCov.median_cov
        File binned_estimated_ecn = PloidyScore.binwise_estimated_copy_numbers
        File estimated_cn = PloidyScore.estimated_copy_numbers
    }
}

task BinMosdepth {
    input {
        File mosdepth_bed_file
        Int bin_size
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import gzip
from statistics import median

bin_size = ~{bin_size}

def flush(contig, bins, out):
    for bin_start in sorted(bins.keys()):
        if len(bins[bin_start]) == bin_size:
            out.write(f"{contig}\t{bin_start}\t{bin_start + bin_size}\t{int(median(bins[bin_start]))}\n")

current_contig = None
bins = {}
with gzip.open("~{mosdepth_bed_file}", 'rt') as f, open("~{prefix}.bed", 'w') as out:
    for line in f:
        chrom, start, end, coverage = line.strip().split('\t')
        start = int(start)
        end = int(end)
        coverage = int(float(coverage))

        if chrom != current_contig:
            if current_contig is not None:
                flush(current_contig, bins, out)
            current_contig = chrom
            bins = {}

        start_bin = (start // bin_size) * bin_size
        end_bin = ((end - 1) // bin_size) * bin_size
        for bs in range(start_bin, end_bin + 1, bin_size):
            be = bs + bin_size
            overlap = min(end, be) - max(start, bs)
            bins.setdefault(bs, []).extend([coverage] * overlap)

    flush(current_contig, bins, out)
CODE

        bgzip ~{prefix}.bed
    >>>

    output {
        File binned_bed = "~{prefix}.bed.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size(mosdepth_bed_file, "GB")) + 10,
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

task SetBins {
    input {
        File count_file
        Int binsize
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        zcat ~{count_file} \
            | awk -v FS="\t" -v OFS="\t" -v b=~{binsize} 'BEGIN{ print "#Chr\tStart\tEnd" } { if ($3-$2==b) print $1,$2,$3 }' \
            | bgzip -c \
            > ~{prefix}.locs.bed.gz
    >>>

    output {
        File bin_locs = "~{prefix}.locs.bed.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10 * ceil(size(count_file, "GB")) + 10,
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

task MakeBincovMatrixColumns {
    input {
        File count_file
        String sample
        Int binsize
        File bin_locs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        TMP_BED="~{prefix}.tmp.bed"
        printf "#Chr\tStart\tEnd\t%s\n" "~{sample}" > "$TMP_BED"
        zcat ~{count_file} \
            | awk -v FS="\t" -v OFS="\t" -v b=~{binsize} '{ if ($3-$2==b) print $0 }' \
            >> "$TMP_BED"

        if ! cut -f1-3 "$TMP_BED" | cmp <(bgzip -cd ~{bin_locs}); then
            echo "~{count_file} has different intervals than ~{bin_locs}"
            exit 1
        fi
        cut -f4- "$TMP_BED" | bgzip -c > ~{prefix}.RD.txt.gz
    >>>

    output {
        File bincov_bed = "~{prefix}.RD.txt.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10 * ceil(size(count_file, "GB") + size(bin_locs, "GB")) + 10,
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

task ZPaste {
    input {
        Array[File]+ column_files
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p column_file_fifos
        FILE_NUM=0
        while read -r COLUMN_FILE; do
            FIFO=$(printf "column_file_fifos/%08d" $FILE_NUM)
            mkfifo "$FIFO"
            bgzip -@$(nproc) -cd "$COLUMN_FILE" > "$FIFO" &
            ((++FILE_NUM))
        done < ~{write_lines(column_files)}

        paste column_file_fifos/* | bgzip -@$(nproc) -c > ~{prefix}.RD.txt.gz
        tabix -p bed ~{prefix}.RD.txt.gz
    >>>

    output {
        File bincov_matrix = "~{prefix}.RD.txt.gz"
        File bincov_matrix_idx = "~{prefix}.RD.txt.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 1 + 0.003 * length(column_files),
        disk_gb: 3 * ceil(size(column_files, "GB")) + 10,
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

task BuildPloidyMatrix {
    input {
        File bincov_matrix
        Int bin_size
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        zcat ~{bincov_matrix} \
            | awk ' \
                function printRow() \
                    {printf "%s\t%d\t%d",chr,start,stop; \
                     for(i=4;i<=nf;++i) {printf "\t%d",vals[i]; vals[i]=0}; \
                     print ""} \
                BEGIN {binSize=~{bin_size}} \
                NR==1 {print substr($0,2)} \
                NR==2 {chr=$1; start=$2; stop=start+binSize; nf=NF; for(i=4;i<=nf;++i) {vals[i]=$i}} \
                NR>2  {if($1!=chr){printRow(); chr=$1; start=$2; stop=start+binSize} \
                       else if($2>=stop) {printRow(); while($2>=stop) {start=stop; stop=start+binSize}} \
                       for(i=4;i<=nf;++i) {vals[i]+=$i}} \
                END   {if(nf!=0)printRow()}' \
            | bgzip > ~{prefix}.ploidy_matrix.bed.gz
    >>>

    output {
        File ploidy_matrix = "~{prefix}.ploidy_matrix.bed.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(bincov_matrix, "GB")) + 10,
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

task PloidyScore {
    input {
        File ploidy_matrix
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir ploidy_est
        Rscript /opt/scripts/helper/estimatePloidy.R -z -O ./ploidy_est ~{ploidy_matrix}

        mv ploidy_est/estimated_copy_numbers.txt.gz ~{prefix}.estimated_copy_numbers.txt.gz
        mv ploidy_est/binwise_estimated_copy_numbers.bed.gz ~{prefix}.binwise_estimated_copy_numbers.bed.gz
    >>>

    output {
        File estimated_copy_numbers = "~{prefix}.estimated_copy_numbers.txt.gz"
        File binwise_estimated_copy_numbers = "~{prefix}.binwise_estimated_copy_numbers.bed.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(ploidy_matrix, "GB")) + 10,
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

task MedianCov {
    input {
        File bincov_matrix
        String prefix
        Int max_bins = 1000000
        Int random_seed = 42
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Reservoir-sample bins down to max_bins before loading into R so peak
        # memory stays bounded regardless of matrix size. Per-sample median
        # coverage is unchanged because medianCoverage.R already downsamples to
        # this many bins internally after reading the full matrix.
        zcat ~{bincov_matrix} \
            | awk -v N=~{max_bins} -v SEED=~{random_seed} '
                BEGIN { srand(SEED) }
                NR==1 { print; next }
                { n++
                  if (n <= N) { res[n] = $0 }
                  else { r = int(rand() * n) + 1; if (r <= N) res[r] = $0 } }
                END { m = (n < N ? n : N); for (i = 1; i <= m; i++) print res[i] }' \
            > ~{prefix}_fixed.bed

        Rscript /opt/scripts/helper/medianCoverage.R ~{prefix}_fixed.bed -H ~{prefix}_medianCov.bed
        
        Rscript -e "x <- read.table(\"~{prefix}_medianCov.bed\",check.names=FALSE); xtransposed <- t(x[,c(1,2)]); write.table(xtransposed,file=\"~{prefix}_medianCov.transposed.bed\",sep=\"\\t\",row.names=F,col.names=F,quote=F)"
    >>>

    output {
        File median_cov = "~{prefix}_medianCov.transposed.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 16,
        disk_gb: 3 * ceil(size(bincov_matrix, "GB")) + 10,
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
