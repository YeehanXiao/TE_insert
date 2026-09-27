test_that("terminal soft clips are extracted with stable metadata", {
  sam <- write_test_sam(c(
    paste(
      "left", 0, "chr1", 101, 60, "5S10M", "*", 0, 0,
      "AAAAACCCCCGGGGG", "*", sep = "\t"
    ),
    paste(
      "right", 0, "chr1", 201, 60, "10M5S", "*", 0, 0,
      "CCCCCGGGGGTTTTT", "IIIIIIIIIIIIIII", sep = "\t"
    ),
    paste(
      "duplicate", 1024, "chr1", 301, 60, "5S10M", "*", 0, 0,
      "AAAAACCCCCGGGGG", "IIIIIIIIIIIIIII", sep = "\t"
    )
  ))
  fastq <- tempfile(fileext = ".fastq")

  expect_invisible(extractSoftClip(sam, fastq, length = 5L, tsd = 3L))
  lines <- readLines(fastq)

  expect_length(lines, 8L)
  expect_match(lines[1], "^@chr1:100:CCC:1:0$")
  expect_identical(lines[2], "AAAAA")
  expect_identical(lines[4], "!!!!!")
  expect_match(lines[5], "^@chr1:210:GGG:2:1$")
  expect_identical(lines[6], "TTTTT")
})

test_that("multimapping assignments use integer target counts", {
  sam <- write_test_sam(
    c(
      "chr1:100:AC:1:0\t0\tL1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
      "chr1:100:AC:1:0\t256\tL1\t2\t40\t10M\t*\t0\t0\t*\t*",
      "chr1:100:AC:1:0\t256\tAlu\t1\t30\t10M\t*\t0\t0\t*\t*",
      "chr2:200:TT:2:1\t0\tL1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
      "chr2:200:TT:2:1\t256\tAlu\t1\t60\t10M\t*\t0\t0\t*\t*"
    ),
    references = c(L1 = 6000L, Alu = 500L)
  )
  bed <- tempfile(fileext = ".bed")

  expect_invisible(insertLocation(sam, bed, ratio = 0.5))
  result <- utils::read.delim(bed, header = FALSE, stringsAsFactors = FALSE)

  expect_equal(nrow(result), 1L)
  expect_identical(result[1, 1:6], data.frame(
    V1 = "chr1", V2 = 100L, V3 = 101L, V4 = "L1", V5 = "AC", V6 = 1L
  ))
  expect_equal(result$V7, 1)
  expect_equal(result$V8, 0.67)
})

test_that("nearby breakpoint evidence is summarized", {
  bed <- tempfile(fileext = ".bed")
  writeLines(c(
    "chr1\t100\t101\tL1\tACGT\t1\t1\t0.90",
    "chr1\t99\t100\tL1\tTTACGT\t2\t1\t0.80"
  ), bed)

  result <- processInsertion(bed, max.gapwidth = 2L)

  expect_s4_class(result, "GRanges")
  expect_length(result, 1L)
  expect_identical(as.character(result$name), "L1")
  expect_identical(as.character(result$te_name), "L1")
  expect_equal(result$left, 1)
  expect_equal(result$right, 1)
  expect_equal(result$support, 2)
  expect_identical(as.character(result$tsd), "ACGT")
})

test_that("an empty breakpoint file returns an empty result", {
  bed <- tempfile(fileext = ".bed")
  file.create(bed)

  expect_s4_class(result <- processInsertion(bed), "GRanges")
  expect_length(result, 0L)
})

test_that("breakpoint clustering preserves the max gap boundary", {
  bed <- tempfile(fileext = ".bed")
  writeLines(c(
    "chr1\t100\t101\tL1\tAAAA\t1\t1\t1",
    "chr1\t110\t111\tL1\tAAAA\t1\t1\t1",
    "chr1\t121\t122\tL1\tAAAA\t1\t1\t1"
  ), bed)

  result <- processInsertion(bed, max.gapwidth = 10L)

  expect_length(result, 2L)
  expect_equal(result$support, c(2, 1))
})
