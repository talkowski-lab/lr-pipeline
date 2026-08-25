version 1.0

workflow MakeDepthMetrics {
  input {
    Array[String] sample_ids
    Array[File] binned_read_counts
    Array[File] mosdepth_per_base

    File duckdb
    String output_prefix
    String unzip_docker
    String sv_base_mini_docker
  }

  call MergeBinnedReadCounts {
    input:
      sample_ids = sample_ids,
      binned_read_counts = binned_read_counts,
      output_prefix = output_prefix,
      sv_base_mini_docker = sv_base_mini_docker
  }

  scatter (i in range(length(sample_ids))) {
    call MedianCov {
      input:
        sample_id = sample_ids[i],
        mosdepth_per_base = mosdepth_per_base[i],
        duckdb = duckdb,
        unzip_docker = unzip_docker
    }
  }

  call MergeMedianCov {
    input:
      sample_ids = sample_ids,
      median_covs = MedianCov.median_cov,
      output_prefix = output_prefix,
      unzip_docker = unzip_docker
  }

  output {
    File merged_bincov = MergeBinnedReadCounts.merged_bincov
    File merged_bincov_index = MergeBinnedReadCounts.merged_bincov_index
    File median_cov = MergeMedianCov.merged_median_cov
  }
}

task MergeBinnedReadCounts {
  input {
    Array[String] sample_ids
    Array[File] binned_read_counts
    String output_prefix

    String sv_base_mini_docker
    Float? mem_gib
    Int? disk_gb
    Int? cpu
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
  }

  Int default_disk_gb = ceil(size(binned_read_counts, "GB") * 3) + 32

  runtime {
    cpu: select_first([cpu, 1])
    memory: select_first([mem_gib, 8]) + " GiB"
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    docker: sv_base_mini_docker
    preemptible: select_first([preemptible_tries, 3])
    maxRetries: select_first([max_retries, 1])
    noAddress: true
  }

  command <<<
    set -euo pipefail

    paste '~{write_lines(sample_ids)}' '~{write_lines(binned_read_counts)}' > manifest.tsv 
    awk '
      function do_read(cmd,    rc) {
        rc = cmd | getline Line
        while (Line ~ /^(@|(CONTIG))/ && rc == 1) {
          rc = cmd | getline Line 
        }

        return rc
      }
      BEGIN { FS = "\t"; OFS = "\t"; Stderr = "/dev/stderr"; printf "#Chr\tStart\tEnd" }
      { samples[NR] = $1; gsub(/\047/, "\047\\\047", $2); depth_paths[NR] = $2; printf "\t%s", $1 }
      END {
        printf "\n"
        while (1) {
          read_count = 0
          for (i = 1; i <= NR; ++i) {
            cmd = sprintf("bgzip -cd \047%s\047", depth_paths[i])
            rc = do_read(cmd)
            if (rc == -1) {
              printf "error reading file %s\n", depth_paths[i] > Stderr
              exit 1
            } else if (rc == 0) {
              continue
            } else {
              nf = split(Line, fields, /\t/)
              contig = fields[1]
              start = fields[2] - 1
              end = fields[3]
              depth = fields[4]
              ++read_count
            }

            if (i == 1) {
              contig0 = contig
              start0 = start
              end0 = end
              printf "%s\t%d\t%d", contig0, start0, end0
            }

            if (contig == contig0 && start == start0 && end == end0) {
              printf "\t%d", depth
            } else {
              printf "mismatched bin coordinates for %s\n", depth_paths[i] > Stderr
              exit 1
            }
          }
          printf "\n"

          if (read_count == 0) {
            break
          } else if (read_count != NR) {
            printf "depth files have differing number of bins\n" > Stderr
            exit 1
          }
        }
      }
    ' manifest.tsv \
      | bgzip -c > '~{output_prefix}.RD.txt.gz'

    tabix --zero-based --begin 2 --end 3 --sequence 1 --comment '#' \
      '~{output_prefix}.RD.txt.gz'
  >>>

  output {
    File merged_bincov = '~{output_prefix}.RD.txt.gz'
    File merged_bincov_index ='~{output_prefix}.RD.txt.gz.tbi'
  }
}

task MedianCov {
  input {
    String sample_id
    File mosdepth_per_base
    File duckdb

    String unzip_docker
    Float? mem_gib
    Int? disk_gb
    Int? cpu
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
  }

  Int default_disk_gb = ceil(size(mosdepth_per_base, "GB") * 3) + 32

  runtime {
    cpu: select_first([cpu, 1])
    memory: select_first([mem_gib, 2]) + " GiB"
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    docker: unzip_docker
    preemptible: select_first([preemptible_tries, 3])
    maxRetries: select_first([max_retries, 1])
    noAddress: true
  }

  command <<<
    set -euo pipefail

    unzip '~{duckdb}'
    gzip -cd '~{mosdepth_per_base}' \
      | awk -F'\t' '{d=$3-$2; for(i=0; i<d; ++i){print $4}}' \
      | gzip -c > temp.tsv.gz

    cat > commands.sql <<EOF
      SET threads = '~{select_first([cpu, 1])}';
      SET memory_limit = '~{ceil(select_first([mem_gib, 8]) / 2)}GiB';
      SET preserve_insertion_order = false;
      CREATE TABLE depths AS SELECT d FROM read_csv('temp.tsv.gz', header = false, columns = {'d': 'INTEGER'});
      COPY (SELECT CAST(floor(approx_quantile(d, 0.5)) AS INTEGER) FROM depths) TO 'median_cov.txt' (FORMAT csv, HEADER false);
EOF

    ./duckdb -bail -f commands.sql temp.duckdb
  >>>

  output {
    Int median_cov = read_int("median_cov.txt")
  }
}

task MergeMedianCov {
  input {
    Array[String] sample_ids
    Array[Int] median_covs
    String output_prefix

    String unzip_docker
    Float? mem_gib
    Int? disk_gb
    Int? cpu
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
  }

  runtime {
    cpu: select_first([cpu, 1])
    memory: select_first([mem_gib, 2]) + " GiB"
    disks: "local-disk " + select_first([disk_gb, 32]) + " HDD"
    bootDiskSizeGb: select_first([boot_disk_gb, 8])
    docker: unzip_docker
    preemptible: select_first([preemptible_tries, 3])
    maxRetries: select_first([max_retries, 1])
    noAddress: true
  }

  command <<<
    set -euo pipefail

    paste -s '~{write_lines(sample_ids)}' '~{write_lines(median_covs)}' > '~{output_prefix}_medianCov.transposed.bed'
  >>>

  output {
    File merged_median_cov = "~{output_prefix}_medianCov.transposed.bed"
  }
}
