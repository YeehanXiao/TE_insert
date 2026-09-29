version 1.0

workflow TEiShortRead {
  input {
    File host_bam
    File te_reference_fasta
    String output_prefix = "tei"
    Int threads = 4
    Int min_mapq = 20
    Int min_clip_length = 30
    Int tsd_length = 10
    Float assignment_ratio = 0.1
  }

  call ExtractSoftClips {
    input:
      host_bam = host_bam,
      output_prefix = output_prefix,
      min_mapq = min_mapq,
      min_clip_length = min_clip_length,
      tsd_length = tsd_length
  }

  call AlignToTeReference {
    input:
      soft_clips_fastq = ExtractSoftClips.soft_clips_fastq,
      te_reference_fasta = te_reference_fasta,
      output_prefix = output_prefix,
      threads = threads
  }

  call CallInsertions {
    input:
      te_alignment_bam = AlignToTeReference.te_alignment_bam,
      output_prefix = output_prefix,
      assignment_ratio = assignment_ratio
  }

  output {
    File soft_clips_fastq = ExtractSoftClips.soft_clips_fastq
    File te_alignment_bam = AlignToTeReference.te_alignment_bam
    File insertion_breakpoints = CallInsertions.insertion_breakpoints
  }
}

task ExtractSoftClips {
  input {
    File host_bam
    String output_prefix
    Int min_mapq
    Int min_clip_length
    Int tsd_length
  }

  command <<<
    Rscript -e 'args <- commandArgs(trailingOnly = TRUE); TEi::extractSoftClip(args[[1]], args[[2]], mapq = as.integer(args[[3]]), length = as.integer(args[[4]]), tsd = as.integer(args[[5]]), overwrite = TRUE)' "~{host_bam}" "~{output_prefix}.softclips.fastq" "~{min_mapq}" "~{min_clip_length}" "~{tsd_length}"
  >>>

  output {
    File soft_clips_fastq = "~{output_prefix}.softclips.fastq"
  }
}

task AlignToTeReference {
  input {
    File soft_clips_fastq
    File te_reference_fasta
    String output_prefix
    Int threads
  }

  command <<<
    Rscript -e 'args <- commandArgs(trailingOnly = TRUE); TEi::buildIndex(args[[1]], "tei_reference", overwrite = TRUE); TEi::alignment("tei_reference", args[[2]], args[[3]], threads = as.integer(args[[4]]), overwrite = TRUE)' "~{te_reference_fasta}" "~{soft_clips_fastq}" "~{output_prefix}.te.bam" "~{threads}"
  >>>

  output {
    File te_alignment_bam = "~{output_prefix}.te.sort.bam"
  }
}

task CallInsertions {
  input {
    File te_alignment_bam
    String output_prefix
    Float assignment_ratio
  }

  command <<<
    Rscript -e 'args <- commandArgs(trailingOnly = TRUE); TEi::insertLocation(args[[1]], args[[2]], ratio = as.numeric(args[[3]]), overwrite = TRUE)' "~{te_alignment_bam}" "~{output_prefix}.insertions.bed" "~{assignment_ratio}"
  >>>

  output {
    File insertion_breakpoints = "~{output_prefix}.insertions.bed"
  }
}
