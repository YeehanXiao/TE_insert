write_test_sam <- function(records, references = c(chr1 = 100000L)) {
  path <- tempfile(fileext = ".sam")
  header <- c(
    "@HD\tVN:1.6\tSO:queryname",
    sprintf("@SQ\tSN:%s\tLN:%d", names(references), references)
  )
  writeLines(c(header, records), path)
  path
}
