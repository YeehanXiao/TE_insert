args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(
    "Usage: validate_long_read.R host.bam te_reference.fa PLATFORM OUTPUT_DIR te_metadata.tsv|- THREADS",
    call. = FALSE
  )
}

runValidation <- function() {
  bam <- normalizePath(args[1L], mustWork = TRUE)
  teReference <- normalizePath(args[2L], mustWork = TRUE)
  platform <- match.arg(args[3L], c("ont", "hifi", "clr"))
  outputDir <- normalizePath(args[4L], mustWork = TRUE)
  teMetadata <- if (identical(args[5L], "-")) NULL else normalizePath(args[5L], mustWork = TRUE)
  threads <- suppressWarnings(as.integer(args[6L]))
  if (is.na(threads) || threads < 1L || as.character(threads) != args[6L]) {
    stop("THREADS must be a positive integer.", call. = FALSE)
  }

  message("Extracting long-read insertion evidence")
  evidence <- TEi::extractLongReadEvidence(bam)

  message("Classifying candidate sequences against the TE reference")
  classification <- TEi::classifyLongReadEvidence(
    evidence,
    teReference,
    teMetadata = teMetadata,
    platform = platform,
    threads = threads
  )

  message("Summarizing insertion calls")
  calls <- TEi::summarizeLongReadInsertions(evidence, classification)
  parameters <- list(
    platform = platform,
    threads = threads,
    minMapQ = 20L,
    minLength = 100L,
    maxReferenceGap = 50L,
    includeSoftClips = TRUE,
    includeSplitGaps = TRUE,
    chunkSize = 10000L,
    minimap2 = "minimap2",
    minAlignedBp = 80L,
    minQueryFraction = 0.5,
    minIdentity = 0.7,
    minTeMapQ = 0L,
    ambiguityFraction = 0.95,
    clusterWindow = 50L,
    minSupport = 2L,
    includeAmbiguous = FALSE
  )
  result <- structure(
    list(
      calls = calls,
      evidence = evidence,
      classifications = classification$classifications,
      hits = classification$hits,
      parameters = parameters,
      limitations = attr(evidence, "experimental_limitations")
    ),
    class = "TEiLongReadResult"
  )

  writeTsv <- function(x, name) {
    utils::write.table(
      x,
      file.path(outputDir, name),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "NA"
    )
  }
  countRows <- function(prefix, x) {
    counts <- table(x, useNA = "ifany")
    labels <- names(counts)
    labels[is.na(labels) | !nzchar(labels)] <- "NA"
    data.frame(
      metric = paste(prefix, labels, sep = "."),
      value = as.character(as.integer(counts)),
      stringsAsFactors = FALSE
    )
  }

  writeTsv(evidence, "evidence.tsv")
  writeTsv(classification$classifications, "classifications.tsv")
  writeTsv(classification$hits, "te_hits.tsv")
  writeTsv(calls, "insertion_calls.tsv")
  writeTsv(
    rbind(
      data.frame(
        metric = c("TEi.version", "platform", "evidence.total", "classifications.total", "te_hits.total", "calls.total", "calls.pass"),
        value = c(
          as.character(utils::packageVersion("TEi")),
          platform,
          as.character(c(
            nrow(evidence),
            nrow(classification$classifications),
            nrow(classification$hits),
            nrow(calls),
            sum(calls$filter == "PASS")
          ))
        ),
        stringsAsFactors = FALSE
      ),
      countRows("evidence", evidence$evidence_type),
      countRows("classification", classification$classifications$classification_status),
      countRows("call_filter", calls$filter)
    ),
    "validation_summary.tsv"
  )
  saveRDS(result, file.path(outputDir, "long_read_result.rds"), compress = "xz")
  writeTsv(
    data.frame(
      parameter = names(parameters),
      value = vapply(parameters, as.character, character(1L)),
      stringsAsFactors = FALSE
    ),
    "parameters.tsv"
  )
  writeLines(result$limitations, file.path(outputDir, "limitations.txt"))
  writeLines(capture.output(sessionInfo()), file.path(outputDir, "sessionInfo.txt"))

  message(
    "Completed with ", nrow(evidence), " evidence records, ",
    nrow(classification$classifications), " classifications, and ",
    nrow(calls), " insertion calls."
  )
}

tryCatch(
  runValidation(),
  error = function(error) {
    message("ERROR: ", conditionMessage(error))
    quit(save = "no", status = 1L)
  }
)
