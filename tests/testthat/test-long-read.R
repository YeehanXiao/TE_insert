make_long_read_bam <- function() {
  sam <- tempfile(fileext = ".sam")
  prefix <- tempfile()
  record <- function(qname, flag, rname, pos, mapq, cigar, sequence, tag = NULL,
                     quality = NULL) {
    if (is.null(quality)) {
      quality <- if (sequence == "*") "*" else paste(rep("I", nchar(sequence)), collapse = "")
    }
    paste(c(qname, flag, rname, pos, mapq, cigar, "*", 0, 0, sequence, quality, tag),
          collapse = "\t")
  }
  insertion <- paste0(strrep("A", 100), strrep("C", 300), strrep("G", 100))
  split <- paste0(strrep("A", 100), strrep("C", 300), strrep("G", 100))
  lines <- c(
    "@HD\tVN:1.6\tSO:unsorted",
    "@SQ\tSN:chr1\tLN:5000",
    record("r_ins", 0, "chr1", 101, 60, "100M300I100M", insertion),
    record("r_left", 0, "chr1", 501, 60, "150S100M",
           paste0(strrep("T", 150), strrep("A", 100))),
    record("r_right", 0, "chr1", 701, 60, "100M150S",
           paste0(strrep("A", 100), strrep("G", 150))),
    record("r_split", 0, "chr1", 1, 60, "100M400S", split,
           "SA:Z:chr1,101,+,400S100M,60,0;"),
    record("r_reverse", 16, "chr1", 901, 60, "50M120I50M",
           paste0(strrep("A", 50), strrep("T", 120), strrep("G", 50))),
    record("r_noqual", 0, "chr1", 1101, 60, "50M110I50M", strrep("A", 210),
           quality = "*"),
    record("r_low", 0, "chr1", 1201, 5, "100M200I100M", strrep("A", 400)),
    record("r_secondary", 256, "chr1", 1401, 60, "100M200I100M", strrep("A", 400)),
    record("r_duplicate", 1024, "chr1", 1601, 60, "100M200I100M", strrep("A", 400)),
    record("r_qcfail", 512, "chr1", 1801, 60, "100M200I100M", strrep("A", 400)),
    record("r_supplementary", 2048, "chr1", 2001, 60, "100M200I100M", strrep("A", 400)),
    record("r_unmapped", 4, "*", 0, 0, "*", "*"),
    record("r_hard", 0, "chr1", 2201, 60, "50H100M", strrep("A", 100))
  )
  writeLines(lines, sam, useBytes = TRUE)
  bam <- Rsamtools::asBam(sam, destination = prefix, overwrite = TRUE)
  unlink(sam)
  bam
}

test_that("long-read BAM extraction streams core evidence types", {
  bam <- make_long_read_bam()
  on.exit(unlink(bam), add = TRUE)
  expect_warning(
    evidence <- extractLongReadEvidence(bam, minLength = 100L, chunkSize = 2L),
    "hard clipping"
  )

  expect_equal(nrow(evidence), 7L)
  expect_setequal(
    unique(evidence$read_id),
    c("r_ins", "r_left", "r_right", "r_split", "r_reverse", "r_noqual")
  )
  expect_true(all(grepl("^TEI_E[0-9]{9}$", evidence$evidence_id)))
  expect_false(any(grepl("r_", evidence$evidence_id)))

  insertion <- evidence[evidence$read_id == "r_ins", ]
  expect_equal(insertion$evidence_type, "cigar_insertion")
  expect_equal(c(insertion$left_breakpoint_0, insertion$right_breakpoint_0), c(200, 200))
  expect_equal(c(insertion$query_start_0, insertion$query_end_0), c(100, 400))
  expect_equal(insertion$sequence, strrep("C", 300))
  expect_equal(insertion$mean_baseq, 40)

  reverse <- evidence[evidence$read_id == "r_reverse", ]
  expect_equal(reverse$host_strand, "-")
  expect_equal(reverse$sequence, strrep("A", 120))

  noqual <- evidence[evidence$read_id == "r_noqual", ]
  expect_true(is.na(noqual$mean_baseq))

  split <- evidence[evidence$read_id == "r_split", ]
  expect_setequal(split$evidence_type, c("split_gap", "right_soft_clip"))
  gap <- split[split$evidence_type == "split_gap", ]
  expect_equal(c(gap$left_breakpoint_0, gap$right_breakpoint_0), c(100, 100))
  expect_equal(c(gap$query_start_0, gap$query_end_0), c(100, 400))
  expect_equal(gap$sequence, strrep("C", 300))
  expect_equal(attr(evidence, "skipped_hard_clipped"), 1L)
})

test_that("soft clips can be disabled without disabling insertion or SA evidence", {
  bam <- make_long_read_bam()
  on.exit(unlink(bam), add = TRUE)
  expect_warning(
    evidence <- extractLongReadEvidence(
      bam, minLength = 100L, includeSoftClips = FALSE, chunkSize = 3L
    ),
    "hard clipping"
  )
  expect_setequal(evidence$evidence_type, c("cigar_insertion", "split_gap"))
  expect_equal(nrow(evidence), 4L)
})

test_that("reverse-strand SA gaps use original read and genomic coordinates", {
  original <- paste0(strrep("A", 100), strrep("C", 300), strrep("G", 100))
  stored <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(original)))
  record <- .tei_lr_record_events(
    "reverse_split", "chr1", 401, 16L, 60L, "400S100M", stored, "*",
    "chr1,301,-,100M400S,60,0;", 20L, 100L, 50L, FALSE, TRUE
  )

  expect_length(record$events, 1L)
  event <- record$events[[1L]]
  expect_equal(event$evidence_type, "split_gap")
  expect_equal(c(event$left_breakpoint_0, event$right_breakpoint_0), c(400, 400))
  expect_equal(c(event$query_start_0, event$query_end_0), c(100, 400))
  expect_equal(event$sequence, strrep("C", 300))
})

test_that("PAF parsing preserves alignment metrics", {
  hits <- .tei_lr_parse_paf_lines(c(
    "TEI_E000000001\t300\t0\t290\t+\tAluY\t300\t0\t290\t280\t290\t60\ttp:A:P\tAS:i:560\tcg:Z:290M",
    "TEI_E000000001\t300\t2\t288\t-\tAluS\t300\t1\t287\t270\t286\t20\ttp:A:S\tAS:i:540"
  ))
  expect_equal(nrow(hits), 2L)
  expect_equal(hits$te_aligned_bp, c(290, 286))
  expect_equal(hits$te_query_fraction, c(290 / 300, 286 / 300))
  expect_equal(hits$te_identity, c(280 / 290, 270 / 286))
  expect_equal(hits$alignment_score, c(560, 540))
  expect_equal(hits$alignment_type, c("P", "S"))
})

test_that("classification reports near-tied TE families as ambiguous", {
  evidence <- data.frame(
    evidence_id = c("E1", "E2", "E3"), sequence = rep(strrep("A", 300), 3),
    stringsAsFactors = FALSE
  )
  hits <- .tei_lr_parse_paf_lines(c(
    "E1\t300\t0\t300\t+\tA1\t300\t0\t300\t290\t300\t60\ttp:A:P\tAS:i:100",
    "E1\t300\t0\t300\t+\tB1\t300\t0\t300\t288\t300\t50\ttp:A:S\tAS:i:98",
    "E2\t300\t0\t300\t-\tA1\t300\t0\t300\t290\t300\t60\ttp:A:P\tAS:i:100"
  ))
  metadata <- data.frame(
    te_id = c("A1", "B1"), te_family = c("Alu", "LINE1"),
    te_subfamily = c("AluY", "L1HS"), stringsAsFactors = FALSE
  )
  classification <- .tei_lr_select_classifications(
    evidence, hits, metadata, minAlignedBp = 80, minQueryFraction = 0.5,
    minIdentity = 0.7, minTeMapQ = 0, ambiguityFraction = 0.95
  )
  expect_equal(classification$classification_status, c("ambiguous", "unique", "unclassified"))
  expect_match(classification$alternatives[1], "Alu/AluY")
  expect_match(classification$alternatives[1], "LINE1/L1HS")
  expect_equal(classification$te_strand[2], "-")
})

test_that("classification executes an external PAF-producing command", {
  executable <- tempfile()
  reference <- tempfile(fileext = ".fa")
  on.exit(unlink(c(executable, reference)), add = TRUE)
  writeLines(c(
    "#!/bin/sh",
    "printf 'E1\\t300\\t0\\t300\\t+\\tTE1\\t300\\t0\\t300\\t295\\t300\\t60\\ttp:A:P\\tAS:i:590\\n'"
  ), executable)
  Sys.chmod(executable, mode = "0755")
  writeLines(c(">TE1", strrep("A", 300)), reference)
  evidence <- data.frame(evidence_id = "E1", sequence = strrep("A", 300),
                         stringsAsFactors = FALSE)
  result <- classifyLongReadEvidence(evidence, reference, minimap2 = executable)
  expect_equal(result$classifications$classification_status, "unique")
  expect_equal(result$classifications$te_id, "TE1")
  expect_equal(result$classifications$te_identity, 295 / 300)
})

test_that("summarization counts unique reads and retains low-support calls", {
  evidence <- data.frame(
    evidence_id = c("E1", "E2", "E3"), read_id = c("r1", "r2", "r3"),
    evidence_type = c("cigar_insertion", "split_gap", "left_soft_clip"),
    seqname = c("chr1", "chr1", "chr2"),
    host_strand = rep("+", 3),
    left_breakpoint_0 = c(100, 105, 500), right_breakpoint_0 = c(100, 105, 500),
    evidence_length = c(300, 310, 250), stringsAsFactors = FALSE
  )
  classification <- data.frame(
    evidence_id = c("E1", "E2", "E3"),
    classification_status = rep("unique", 3), te_id = rep("AluY", 3),
    te_family = rep("Alu", 3), te_subfamily = rep("AluY", 3),
    te_strand = rep("+", 3), te_identity = c(0.98, 0.97, 0.95),
    stringsAsFactors = FALSE
  )
  calls <- summarizeLongReadInsertions(evidence, classification,
                                       clusterWindow = 10L, minSupport = 2L)
  expect_equal(nrow(calls), 2L)
  expect_equal(calls$supporting_reads, c(2L, 1L))
  expect_equal(calls$filter, c("PASS", "LOW_SUPPORT"))
  expect_equal(calls$cigar_support, c(1L, 0L))
  expect_equal(calls$split_support, c(1L, 0L))
  expect_equal(calls$clip_support, c(0L, 1L))
})

test_that("opposite host strands support one genomic TE orientation", {
  evidence <- data.frame(
    evidence_id = c("E1", "E2"), read_id = c("forward", "reverse"),
    evidence_type = rep("cigar_insertion", 2), seqname = rep("chr1", 2),
    host_strand = c("+", "-"), left_breakpoint_0 = c(100, 104),
    right_breakpoint_0 = c(100, 104), evidence_length = c(300, 305),
    stringsAsFactors = FALSE
  )
  classification <- data.frame(
    evidence_id = c("E1", "E2"), classification_status = rep("unique", 2),
    te_id = rep("L1HS", 2), te_family = rep("LINE1", 2),
    te_subfamily = rep("L1HS", 2), te_strand = c("+", "-"),
    te_identity = c(0.98, 0.97), stringsAsFactors = FALSE
  )

  calls <- summarizeLongReadInsertions(evidence, classification, minSupport = 2L)

  expect_equal(nrow(calls), 1L)
  expect_equal(calls$supporting_reads, 2L)
  expect_equal(calls$te_strand, "+")
  expect_equal(calls$filter, "PASS")
})

test_that("minimap2 classifies an exact long candidate when available", {
  skip_if(Sys.which("minimap2") == "", "minimap2 is not installed")
  set.seed(11)
  target <- paste(sample(c("A", "C", "G", "T"), 2000, replace = TRUE), collapse = "")
  reference <- tempfile(fileext = ".fa")
  on.exit(unlink(reference), add = TRUE)
  writeLines(c(">synthetic_TE", target), reference)
  evidence <- data.frame(
    evidence_id = "TEI_E000000001", sequence = substr(target, 501, 1500),
    stringsAsFactors = FALSE
  )
  result <- classifyLongReadEvidence(
    evidence, reference, platform = "ont", minAlignedBp = 500L,
    minQueryFraction = 0.8, minIdentity = 0.8
  )
  expect_equal(result$classifications$classification_status, "unique")
  expect_equal(result$classifications$te_id, "synthetic_TE")
  expect_gte(result$classifications$te_identity, 0.99)
})
