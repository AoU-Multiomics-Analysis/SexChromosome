version 1.0

struct XYalignSample {
  String sample_id
  File bam
  File bam_index
  String? y_status
}

workflow XYalign {
  input {
    Array[XYalignSample]+ samples
    File reference_fasta
    File? reference_fasta_index
    File? reference_mask
    Array[String]+ chromosomes = ["chr19", "chrX", "chrY"]
    String x_chromosome = "chrX"
    String y_chromosome = "chrY"
    Boolean characterize = false
    Boolean call_variants = false
    Boolean run_remapping = false
    Boolean single_end = false
    Int window_size = 5000
    Int num_permutations = 10000
    Int num_bootstraps = 10000
    Int analysis_cpu = 2
    Int analysis_memory_gb = 16
    Int analysis_disk_gb = 100
    Int reference_cpu = 4
    Int reference_memory_gb = 16
    Int reference_disk_gb = 200
    Int remap_cpu = 8
    Int remap_memory_gb = 32
    Int remap_disk_gb = 300
    String java_memory = "8g"
    Int preemptible_attempts = 0
    String docker_image = "ghcr.io/aou-multiomics-analysis/sexchromosome:main"
  }

  if (run_remapping) {
    call PrepareXYReferences {
      input:
        reference = reference_fasta,
        reference_index = reference_fasta_index,
        reference_mask = reference_mask,
        y_chromosome = y_chromosome,
        cpus = reference_cpu,
        memory_gb = reference_memory_gb,
        disk_gb = reference_disk_gb,
        preemptible_attempts = preemptible_attempts,
        docker_image = docker_image
    }
  }

  scatter (sample in samples) {
    call AnalyzeXYChromosomes {
      input:
        sample_id = sample.sample_id,
        bam = sample.bam,
        bam_index = sample.bam_index,
        reference = reference_fasta,
        reference_index = reference_fasta_index,
        chromosomes = chromosomes,
        x_chromosome = x_chromosome,
        y_chromosome = y_chromosome,
        y_status = select_first([sample.y_status, "unknown"]),
        remapping_requested = run_remapping,
        characterize = characterize,
        call_variants = call_variants,
        window_size = window_size,
        num_permutations = num_permutations,
        num_bootstraps = num_bootstraps,
        cpus = analysis_cpu,
        memory_gb = analysis_memory_gb,
        disk_gb = analysis_disk_gb,
        preemptible_attempts = preemptible_attempts,
        docker_image = docker_image
    }
    if (run_remapping) {
      call RemapXYChromosomes {
        input:
          sample_id = sample.sample_id,
          bam = sample.bam,
          bam_index = sample.bam_index,
          reference = reference_fasta,
          reference_index = reference_fasta_index,
          prepared_references = select_first([PrepareXYReferences.references]),
          analysis_log = AnalyzeXYChromosomes.log,
          x_chromosome = x_chromosome,
          y_chromosome = y_chromosome,
          y_status = select_first([sample.y_status, "unknown"]),
          single_end = single_end,
          java_memory = java_memory,
          cpus = remap_cpu,
          memory_gb = remap_memory_gb,
          disk_gb = remap_disk_gb,
          preemptible_attempts = preemptible_attempts,
          docker_image = docker_image
      }
    }
  }

  output {
    Array[String] sample_ids = AnalyzeXYChromosomes.sample_id_out
    Array[File] chromosome_depth = AnalyzeXYChromosomes.depth
    Array[File] chromosome_mapq = AnalyzeXYChromosomes.mapq
    Array[File] analysis_reports = AnalyzeXYChromosomes.reports
    Array[File] analysis_logs = AnalyzeXYChromosomes.log
    File? prepared_references = PrepareXYReferences.references
    File? reference_log = PrepareXYReferences.log
    Array[File?] remapped_bams = RemapXYChromosomes.bam_out
    Array[File?] remapped_bam_indexes = RemapXYChromosomes.bam_index_out
    Array[File?] remapping_reports = RemapXYChromosomes.reports
    Array[File?] remapping_logs = RemapXYChromosomes.log
  }

  meta {
    description: "XYalign chromosome analysis and optional read remapping for Terra."
  }
}

task AnalyzeXYChromosomes {
  input {
    String sample_id
    File bam
    File bam_index
    File reference
    File? reference_index
    Array[String]+ chromosomes
    String x_chromosome
    String y_chromosome
    String y_status
    Boolean remapping_requested
    Boolean characterize
    Boolean call_variants
    Int window_size
    Int num_permutations
    Int num_bootstraps
    Int cpus
    Int memory_gb
    Int disk_gb
    Int preemptible_attempts
    String docker_image
  }

  File chromosome_list = write_lines(chromosomes)

  command <<<
    set -euo pipefail
    echo "[XYalign] Validate localized inputs and start analyze"
    optional_args=()
    reference_index='~{sub(select_first([reference_index, ""]), "'", "'\"'\"'")}'
    if [[ -n "$reference_index" ]]; then optional_args+=(--reference-index "$reference_index"); fi
    python /opt/sexchromosome/scripts/xyalign_runner.py analyze \
    --bam '~{sub(bam, "'", "'\"'\"'")}' \
    --bam-index '~{sub(bam_index, "'", "'\"'\"'")}' \
    --reference '~{sub(reference, "'", "'\"'\"'")}' \
    --sample-id '~{sub(sample_id, "'", "'\"'\"'")}' \
    --chromosomes-file '~{sub(chromosome_list, "'", "'\"'\"'")}' \
    --x-chromosome '~{sub(x_chromosome, "'", "'\"'\"'")}' \
    --y-chromosome '~{sub(y_chromosome, "'", "'\"'\"'")}' \
    --y-status '~{sub(y_status, "'", "'\"'\"'")}' \
    ~{if characterize then "--characterize" else ""} \
    ~{if call_variants then "--call-variants" else ""} \
    ~{if remapping_requested then "--remapping-requested" else ""} \
    --window-size ~{window_size} --num-permutations ~{num_permutations} --num-bootstraps ~{num_bootstraps} \
    --cpus ~{cpus} --output-dir result "${optional_args[@]}"
    echo "[XYalign] Completed analyze"
  >>>

  output {
    String sample_id_out = sample_id
    File depth = "result/chromosome_depth.tsv"
    File mapq = "result/chromosome_mapq.tsv"
    File reports = "result/reports.tar.gz"
    File log = "result/runner.log"
  }

  runtime {
    docker: docker_image
    cpu: cpus
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_attempts
    maxRetries: 1
  }
}

task PrepareXYReferences {
  input {
    File reference
    File? reference_index
    File? reference_mask
    String y_chromosome
    Int cpus
    Int memory_gb
    Int disk_gb
    Int preemptible_attempts
    String docker_image
  }

  command <<<
    set -euo pipefail
    echo "[XYalign] Validate localized inputs and start prepare-reference"
    optional_args=()
    reference_index='~{sub(select_first([reference_index, ""]), "'", "'\"'\"'")}'
    if [[ -n "$reference_index" ]]; then optional_args+=(--reference-index "$reference_index"); fi
    reference_mask='~{sub(select_first([reference_mask, ""]), "'", "'\"'\"'")}'
    if [[ -n "$reference_mask" ]]; then optional_args+=(--reference-mask "$reference_mask"); fi
    python /opt/sexchromosome/scripts/xyalign_runner.py prepare-reference \
    --reference '~{sub(reference, "'", "'\"'\"'")}' \
    --y-chromosome '~{sub(y_chromosome, "'", "'\"'\"'")}' \
    --cpus ~{cpus} --output-dir result "${optional_args[@]}"
    echo "[XYalign] Completed prepare-reference"
  >>>

  output {
    File references = "result/references.tar"
    File log = "result/runner.log"
  }

  runtime {
    docker: docker_image
    cpu: cpus
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_attempts
    maxRetries: 1
  }
}

task RemapXYChromosomes {
  input {
    String sample_id
    File bam
    File bam_index
    File reference
    File? reference_index
    File prepared_references
    File analysis_log
    String x_chromosome
    String y_chromosome
    String y_status
    Boolean single_end
    String java_memory
    Int cpus
    Int memory_gb
    Int disk_gb
    Int preemptible_attempts
    String docker_image
  }

  command <<<
    set -euo pipefail
    echo "[XYalign] Validate localized inputs and start remap"
    optional_args=()
    reference_index='~{sub(select_first([reference_index, ""]), "'", "'\"'\"'")}'
    if [[ -n "$reference_index" ]]; then optional_args+=(--reference-index "$reference_index"); fi
    test -r '~{sub(analysis_log, "'", "'\"'\"'")}'
    python /opt/sexchromosome/scripts/xyalign_runner.py remap \
    --bam '~{sub(bam, "'", "'\"'\"'")}' \
    --bam-index '~{sub(bam_index, "'", "'\"'\"'")}' \
    --reference '~{sub(reference, "'", "'\"'\"'")}' \
    --prepared-references '~{sub(prepared_references, "'", "'\"'\"'")}' \
    --sample-id '~{sub(sample_id, "'", "'\"'\"'")}' \
    --x-chromosome '~{sub(x_chromosome, "'", "'\"'\"'")}' \
    --y-chromosome '~{sub(y_chromosome, "'", "'\"'\"'")}' \
    --y-status '~{sub(y_status, "'", "'\"'\"'")}' \
    --java-memory '~{sub(java_memory, "'", "'\"'\"'")}' \
    ~{if single_end then "--single-end" else ""} \
    --cpus ~{cpus} --output-dir result "${optional_args[@]}"
    echo "[XYalign] Completed remap"
  >>>

  output {
    File bam_out = "result/remapped.bam"
    File bam_index_out = "result/remapped.bam.bai"
    File reports = "result/reports.tar.gz"
    File log = "result/runner.log"
  }

  runtime {
    docker: docker_image
    cpu: cpus
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_attempts
    maxRetries: 1
  }
}
