.tei_lr_limitations <- c(
  "Experimental support is limited to genomic long reads aligned in BAM format.",
  "Evidence is limited to CIGAR insertions, terminal soft clips, and colinear same-strand SA gaps.",
  "Consensus assembly, genotyping, mosaic calibration, methylation, phasing, TSD reconstruction, and transduction annotation are not implemented."
)

.tei_lr_check_scalar <- function(x, name, lower, integer = FALSE) {
  if (length(x) != 1L || is.na(x) || !is.numeric(x) || x < lower ||
      (integer && x != as.integer(x))) {
    stop(name, " must be a single ", if (integer) "integer " else "numeric ",
         "value >= ", lower, ".", call. = FALSE)
  }
  if (integer) as.integer(x) else as.numeric(x)
}

.tei_lr_empty_evidence <- function() {
  data.frame(
    evidence_id = character(), read_id = character(), evidence_type = character(),
    seqname = character(), left_breakpoint_0 = numeric(),
    right_breakpoint_0 = numeric(), reference_gap = numeric(),
    host_strand = character(), host_mapq = integer(), query_start_0 = numeric(),
    query_end_0 = numeric(), evidence_length = numeric(), mean_baseq = numeric(),
    sequence = character(), source_cigar = character(), sa_entry = character(),
    stringsAsFactors = FALSE
  )
}

.tei_lr_empty_hits <- function() {
  data.frame(
    evidence_id = character(), query_length = numeric(), query_start_0 = numeric(),
    query_end_0 = numeric(), te_strand = character(), te_id = character(),
    target_length = numeric(), target_start_0 = numeric(), target_end_0 = numeric(),
    matches = numeric(), block_length = numeric(), te_mapq = integer(),
    te_aligned_bp = numeric(), te_query_fraction = numeric(), te_identity = numeric(),
    alignment_score = numeric(), alignment_type = character(),
    stringsAsFactors = FALSE
  )
}

.tei_lr_empty_classifications <- function() {
  data.frame(
    evidence_id = character(), classification_status = character(), te_id = character(),
    te_family = character(), te_subfamily = character(), te_strand = character(),
    te_mapq = integer(), te_aligned_bp = numeric(), te_query_fraction = numeric(),
    te_identity = numeric(), alignment_score = numeric(), alternatives = character(),
    stringsAsFactors = FALSE
  )
}

.tei_lr_parse_cigar <- function(cigar) {
  if (length(cigar) != 1L || is.na(cigar) || cigar == "*") {
    return(data.frame(length = numeric(), op = character(), stringsAsFactors = FALSE))
  }
  match <- gregexpr("[0-9]+[MIDNSHP=X]", cigar, perl = TRUE)[[1L]]
  if (match[1L] == -1L) stop("Invalid CIGAR: ", cigar, call. = FALSE)
  token <- regmatches(cigar, list(match))[[1L]]
  if (paste0(token, collapse = "") != cigar) stop("Invalid CIGAR: ", cigar, call. = FALSE)
  data.frame(
    length = as.numeric(sub("[MIDNSHP=X]$", "", token)),
    op = sub("^[0-9]+", "", token),
    stringsAsFactors = FALSE
  )
}

.tei_lr_consumes_query <- function(op) op %in% c("M", "I", "S", "=", "X")
.tei_lr_consumes_reference <- function(op) op %in% c("M", "D", "N", "=", "X")

.tei_lr_mean_quality <- function(quality, start, end) {
  if (is.na(quality) || quality == "*" || end <= start || nchar(quality) < end) return(NA_real_)
  value <- utf8ToInt(substr(quality, start + 1, end))
  if (any(value < 33L | value > 126L)) return(NA_real_)
  mean(value - 33L)
}

.tei_lr_event <- function(readId, type, seqname, left, right, referenceGap,
                          strand, mapq, queryStart, queryEnd, sequence, quality,
                          cigar, saEntry = NA_character_) {
  list(
    evidence_id = NA_character_, read_id = readId, evidence_type = type,
    seqname = seqname, left_breakpoint_0 = as.numeric(left),
    right_breakpoint_0 = as.numeric(right), reference_gap = as.numeric(referenceGap),
    host_strand = strand, host_mapq = as.integer(mapq),
    query_start_0 = as.numeric(queryStart), query_end_0 = as.numeric(queryEnd),
    evidence_length = as.numeric(queryEnd - queryStart),
    mean_baseq = .tei_lr_mean_quality(quality, queryStart, queryEnd),
    sequence = substr(sequence, queryStart + 1, queryEnd),
    source_cigar = cigar, sa_entry = saEntry
  )
}

.tei_lr_segment <- function(seqname, pos, strand, cigar, mapq, readLength,
                            saEntry = NA_character_) {
  parsed <- .tei_lr_parse_cigar(cigar)
  if (!nrow(parsed)) return(NULL)
  fullLength <- sum(parsed$length[parsed$op %in% c("M", "I", "S", "H", "=", "X")])
  if (fullLength != readLength) return(NULL)
  left <- 1L
  while (left <= nrow(parsed) && parsed$op[left] %in% c("H", "S")) left <- left + 1L
  right <- nrow(parsed)
  while (right >= 1L && parsed$op[right] %in% c("H", "S")) right <- right - 1L
  if (left > right) return(NULL)
  leading <- if (left > 1L) sum(parsed$length[seq_len(left - 1L)]) else 0
  trailing <- if (right < nrow(parsed)) sum(parsed$length[seq.int(right + 1L, nrow(parsed))]) else 0
  queryStart <- if (strand == "+") leading else trailing
  queryEnd <- if (strand == "+") fullLength - trailing else fullLength - leading
  referenceStart <- as.numeric(pos) - 1
  referenceEnd <- referenceStart + sum(parsed$length[.tei_lr_consumes_reference(parsed$op)])
  list(
    seqname = seqname, strand = strand, cigar = cigar, mapq = as.integer(mapq),
    query_start = queryStart, query_end = queryEnd,
    reference_start = referenceStart, reference_end = referenceEnd,
    sa_entry = saEntry
  )
}

.tei_lr_parse_sa <- function(sa) {
  if (length(sa) != 1L || is.na(sa) || !nzchar(sa)) return(list())
  entries <- strsplit(sa, ";", fixed = TRUE)[[1L]]
  entries <- entries[nzchar(entries)]
  lapply(entries, function(entry) {
    field <- strsplit(entry, ",", fixed = TRUE)[[1L]]
    if (length(field) != 6L || !field[3L] %in% c("+", "-")) return(NULL)
    pos <- suppressWarnings(as.numeric(field[2L]))
    mapq <- suppressWarnings(as.integer(field[5L]))
    if (is.na(pos) || is.na(mapq)) return(NULL)
    list(seqname = field[1L], pos = pos, strand = field[3L], cigar = field[4L],
         mapq = mapq, entry = entry)
  })
}

.tei_lr_record_events <- function(readId, seqname, pos, flag, mapq, cigar,
                                  sequence, quality, sa, minMapQ, minLength,
                                  maxReferenceGap, includeSoftClips,
                                  includeSplitGaps) {
  parsed <- .tei_lr_parse_cigar(cigar)
  if (!nrow(parsed)) return(list(events = list(), hardClipped = FALSE))
  if (any(parsed$op == "H")) return(list(events = list(), hardClipped = TRUE))
  queryLength <- sum(parsed$length[.tei_lr_consumes_query(parsed$op)])
  if (queryLength != nchar(sequence)) {
    stop("CIGAR and SEQ lengths disagree for read ", readId, ".", call. = FALSE)
  }
  strand <- if (bitwAnd(as.integer(flag), 16L)) "-" else "+"
  originalSequence <- sequence
  originalQuality <- quality
  if (strand == "-") {
    originalSequence <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAString(sequence))
    )
    if (!is.na(quality) && quality != "*") {
      originalQuality <- paste(rev(strsplit(quality, "", fixed = TRUE)[[1L]]),
                               collapse = "")
    }
  }
  originalInterval <- function(start, end) {
    if (strand == "+") c(start, end) else c(queryLength - end, queryLength - start)
  }
  query <- 0
  reference <- as.numeric(pos) - 1
  events <- list()
  add <- function(event) events[[length(events) + 1L]] <<- event
  nonHard <- which(parsed$op != "H")
  first <- nonHard[1L]
  last <- nonHard[length(nonHard)]
  for (i in seq_len(nrow(parsed))) {
    length <- parsed$length[i]
    op <- parsed$op[i]
    if (op == "I" && length >= minLength) {
      interval <- originalInterval(query, query + length)
      add(.tei_lr_event(readId, "cigar_insertion", seqname, reference, reference, 0,
                        strand, mapq, interval[1L], interval[2L], originalSequence,
                        originalQuality, cigar))
    }
    if (includeSoftClips && op == "S" && length >= minLength && i %in% c(first, last)) {
      type <- if (i == first) "left_soft_clip" else "right_soft_clip"
      interval <- originalInterval(query, query + length)
      add(.tei_lr_event(readId, type, seqname, reference, reference, 0,
                        strand, mapq, interval[1L], interval[2L], originalSequence,
                        originalQuality, cigar))
    }
    if (.tei_lr_consumes_query(op)) query <- query + length
    if (.tei_lr_consumes_reference(op)) reference <- reference + length
  }
  if (!includeSplitGaps || is.na(sa) || !nzchar(sa)) {
    return(list(events = events, hardClipped = FALSE))
  }
  segments <- list(.tei_lr_segment(seqname, pos, strand, cigar, mapq, queryLength))
  for (entry in .tei_lr_parse_sa(sa)) {
    if (is.null(entry) || entry$mapq < minMapQ || entry$seqname != seqname ||
        entry$strand != strand) next
    segment <- .tei_lr_segment(entry$seqname, entry$pos, entry$strand, entry$cigar,
                               entry$mapq, queryLength, entry$entry)
    if (!is.null(segment)) segments[[length(segments) + 1L]] <- segment
  }
  segments <- Filter(Negate(is.null), segments)
  if (length(segments) < 2L) return(list(events = events, hardClipped = FALSE))
  key <- vapply(segments, function(x) paste(x$seqname, x$reference_start,
                                             x$query_start, x$cigar, x$strand, sep = "\r"), "")
  segments <- segments[!duplicated(key)]
  segments <- segments[order(vapply(segments, `[[`, numeric(1), "query_start"),
                             vapply(segments, `[[`, numeric(1), "query_end"))]
  if (length(segments) < 2L) return(list(events = events, hardClipped = FALSE))
  for (i in seq_len(length(segments) - 1L)) {
    leftSegment <- segments[[i]]
    rightSegment <- segments[[i + 1L]]
    queryGap <- rightSegment$query_start - leftSegment$query_end
    referenceGap <- if (strand == "+") {
      rightSegment$reference_start - leftSegment$reference_end
    } else {
      leftSegment$reference_start - rightSegment$reference_end
    }
    if (queryGap < minLength || abs(referenceGap) > maxReferenceGap) next
    breakpointPair <- if (strand == "+") {
      c(leftSegment$reference_end, rightSegment$reference_start)
    } else {
      c(rightSegment$reference_end, leftSegment$reference_start)
    }
    leftBreakpoint <- min(breakpointPair)
    rightBreakpoint <- max(breakpointPair)
    add(.tei_lr_event(
      readId, "split_gap", seqname, leftBreakpoint, rightBreakpoint, referenceGap,
      strand, min(leftSegment$mapq, rightSegment$mapq), leftSegment$query_end,
      rightSegment$query_start, originalSequence, originalQuality,
      paste(leftSegment$cigar, rightSegment$cigar, sep = "|"),
      paste(stats::na.omit(c(leftSegment$sa_entry, rightSegment$sa_entry)), collapse = ";")
    ))
  }
  list(events = events, hardClipped = FALSE)
}

#' Extract candidate TE insertion sequences from long-read alignments
#'
#' Streams primary BAM alignments and extracts large CIGAR insertions, terminal
#' soft clips, and reliable colinear same-strand gaps described by SA tags.
#'
#' @param alignment Path to a BAM file.
#' @param minMapQ Minimum host mapping quality.
#' @param minLength Minimum candidate sequence length.
#' @param maxReferenceGap Maximum absolute host gap or overlap for SA evidence.
#' @param includeSoftClips Include terminal soft-clip evidence.
#' @param includeSplitGaps Include colinear same-strand SA-gap evidence.
#' @param chunkSize Number of BAM records read per chunk.
#' @return A data frame with one row per candidate sequence.
#' @export
extractLongReadEvidence <- function(alignment, minMapQ = 20L, minLength = 100L,
                                    maxReferenceGap = 50L, includeSoftClips = TRUE,
                                    includeSplitGaps = TRUE, chunkSize = 10000L) {
  if (length(alignment) != 1L || !file.exists(alignment)) {
    stop("alignment must be an existing BAM file.", call. = FALSE)
  }
  if (!grepl("\\.bam$", alignment, ignore.case = TRUE)) {
    stop("Experimental long-read extraction currently supports BAM input only.", call. = FALSE)
  }
  minMapQ <- .tei_lr_check_scalar(minMapQ, "minMapQ", 0, TRUE)
  minLength <- .tei_lr_check_scalar(minLength, "minLength", 1, TRUE)
  maxReferenceGap <- .tei_lr_check_scalar(maxReferenceGap, "maxReferenceGap", 0, TRUE)
  chunkSize <- .tei_lr_check_scalar(chunkSize, "chunkSize", 1, TRUE)
  if (length(includeSoftClips) != 1L || is.na(includeSoftClips) ||
      length(includeSplitGaps) != 1L || is.na(includeSplitGaps)) {
    stop("includeSoftClips and includeSplitGaps must be TRUE or FALSE.", call. = FALSE)
  }
  bam <- Rsamtools::BamFile(alignment, yieldSize = chunkSize)
  Rsamtools::open.BamFile(bam)
  on.exit(Rsamtools::close.BamFile(bam), add = TRUE)
  param <- Rsamtools::ScanBamParam(
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar", "seq", "qual"),
    tag = "SA"
  )
  eventChunks <- list()
  hardClipped <- 0L
  repeat {
    chunk <- Rsamtools::scanBam(bam, param = param)[[1L]]
    n <- length(chunk$qname)
    if (!n) break
    chunkEvents <- vector("list", n)
    sa <- chunk$tag[["SA"]]
    if (is.null(sa)) sa <- rep(NA_character_, n)
    sa <- as.character(sa)
    if (length(sa) != n) sa <- rep(NA_character_, n)
    sequence <- as.character(chunk$seq)
    quality <- as.character(chunk$qual)
    for (i in seq_len(n)) {
      flag <- as.integer(chunk$flag[i])
      if (is.na(flag) || bitwAnd(flag, 3844L) != 0L || is.na(chunk$mapq[i]) ||
          chunk$mapq[i] < minMapQ) next
      record <- .tei_lr_record_events(
        as.character(chunk$qname[i]), as.character(chunk$rname[i]), chunk$pos[i], flag,
        chunk$mapq[i], as.character(chunk$cigar[i]), sequence[i], quality[i], sa[i],
        minMapQ, minLength, maxReferenceGap, isTRUE(includeSoftClips),
        isTRUE(includeSplitGaps)
      )
      hardClipped <- hardClipped + as.integer(record$hardClipped)
      if (length(record$events)) chunkEvents[[i]] <- record$events
    }
    chunkEvents <- unlist(chunkEvents, recursive = FALSE, use.names = FALSE)
    if (length(chunkEvents)) eventChunks[[length(eventChunks) + 1L]] <- chunkEvents
  }
  if (hardClipped) {
    warning("Skipped ", hardClipped,
            " primary alignment(s) with hard clipping because the full read sequence was unavailable.",
            call. = FALSE)
  }
  events <- unlist(eventChunks, recursive = FALSE, use.names = FALSE)
  if (!length(events)) {
    output <- .tei_lr_empty_evidence()
  } else {
    output <- do.call(rbind.data.frame, c(events, list(stringsAsFactors = FALSE)))
    rownames(output) <- NULL
    priority <- match(output$evidence_type,
                      c("cigar_insertion", "split_gap", "left_soft_clip", "right_soft_clip"))
    output <- output[order(priority), , drop = FALSE]
    key <- paste(output$read_id, output$seqname, output$left_breakpoint_0,
                 output$right_breakpoint_0, output$query_start_0,
                 output$query_end_0, sep = "\r")
    output <- output[!duplicated(key), , drop = FALSE]
    output <- output[order(output$seqname, output$left_breakpoint_0, output$read_id,
                           output$query_start_0), , drop = FALSE]
    output$evidence_id <- sprintf("TEI_E%09d", seq_len(nrow(output)))
    rownames(output) <- NULL
  }
  attr(output, "experimental_limitations") <- .tei_lr_limitations
  attr(output, "skipped_hard_clipped") <- hardClipped
  output
}

.tei_lr_parse_paf_lines <- function(lines) {
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(.tei_lr_empty_hits())
  parsed <- lapply(lines, function(line) {
    field <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(field) < 12L) stop("Malformed minimap2 PAF output.", call. = FALSE)
    tag <- field[-seq_len(12L)]
    value <- function(prefix) {
      hit <- tag[startsWith(tag, prefix)]
      if (length(hit)) sub(prefix, "", hit[1L], fixed = TRUE) else NA_character_
    }
    queryLength <- as.numeric(field[2L])
    queryStart <- as.numeric(field[3L])
    queryEnd <- as.numeric(field[4L])
    matches <- as.numeric(field[10L])
    blockLength <- as.numeric(field[11L])
    data.frame(
      evidence_id = field[1L], query_length = queryLength,
      query_start_0 = queryStart, query_end_0 = queryEnd, te_strand = field[5L],
      te_id = field[6L], target_length = as.numeric(field[7L]),
      target_start_0 = as.numeric(field[8L]), target_end_0 = as.numeric(field[9L]),
      matches = matches, block_length = blockLength, te_mapq = as.integer(field[12L]),
      te_aligned_bp = queryEnd - queryStart,
      te_query_fraction = if (queryLength > 0) (queryEnd - queryStart) / queryLength else NA_real_,
      te_identity = if (blockLength > 0) matches / blockLength else NA_real_,
      alignment_score = suppressWarnings(as.numeric(value("AS:i:"))),
      alignment_type = value("tp:A:"), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, parsed)
}

.tei_lr_read_metadata <- function(teMetadata) {
  if (is.null(teMetadata)) return(NULL)
  if (is.character(teMetadata) && length(teMetadata) == 1L) {
    if (!file.exists(teMetadata)) stop("teMetadata file does not exist.", call. = FALSE)
    teMetadata <- utils::read.delim(teMetadata, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!is.data.frame(teMetadata) || !"te_id" %in% names(teMetadata)) {
    stop("teMetadata must be a data frame or TSV with a te_id column.", call. = FALSE)
  }
  if (anyDuplicated(teMetadata$te_id)) stop("teMetadata te_id values must be unique.", call. = FALSE)
  if (!"te_family" %in% names(teMetadata)) teMetadata$te_family <- NA_character_
  if (!"te_subfamily" %in% names(teMetadata)) teMetadata$te_subfamily <- NA_character_
  teMetadata[, c("te_id", "te_family", "te_subfamily"), drop = FALSE]
}

.tei_lr_select_classifications <- function(evidence, hits, metadata, minAlignedBp,
                                            minQueryFraction, minIdentity, minTeMapQ,
                                            ambiguityFraction) {
  if (!nrow(evidence)) return(.tei_lr_empty_classifications())
  if (nrow(hits)) {
    hits$te_family <- NA_character_
    hits$te_subfamily <- NA_character_
    if (!is.null(metadata)) {
      index <- match(hits$te_id, metadata$te_id)
      hits$te_family <- metadata$te_family[index]
      hits$te_subfamily <- metadata$te_subfamily[index]
    }
  }
  hitIndex <- if (nrow(hits)) {
    split(seq_len(nrow(hits)), hits$evidence_id)
  } else {
    list()
  }
  rows <- lapply(evidence$evidence_id, function(id) {
    index <- hitIndex[[id]]
    candidate <- if (is.null(index)) hits[0L, , drop = FALSE] else hits[index, , drop = FALSE]
    passing <- candidate[
      candidate$te_aligned_bp >= minAlignedBp &
        candidate$te_query_fraction >= minQueryFraction &
        candidate$te_identity >= minIdentity & candidate$te_mapq >= minTeMapQ,
      , drop = FALSE
    ]
    pool <- if (nrow(passing)) passing else candidate
    if (!nrow(pool)) {
      return(data.frame(
        evidence_id = id, classification_status = "unclassified", te_id = NA_character_,
        te_family = NA_character_, te_subfamily = NA_character_, te_strand = NA_character_,
        te_mapq = NA_integer_, te_aligned_bp = NA_real_, te_query_fraction = NA_real_,
        te_identity = NA_real_, alignment_score = NA_real_, alternatives = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    score <- pool$alignment_score
    score[is.na(score)] <- pool$matches[is.na(score)]
    orderIndex <- order(-score, -pool$te_aligned_bp, -pool$te_identity, -pool$te_mapq,
                        pool$te_id, na.last = TRUE)
    pool <- pool[orderIndex, , drop = FALSE]
    score <- score[orderIndex]
    best <- pool[1L, , drop = FALSE]
    if (!nrow(passing)) {
      status <- "unclassified"
      alternatives <- NA_character_
    } else {
      close <- if (is.na(score[1L])) rep(FALSE, length(score)) else if (score[1L] > 0) {
        !is.na(score) & score >= score[1L] * ambiguityFraction
      } else {
        !is.na(score) & score == score[1L]
      }
      label <- ifelse(!is.na(pool$te_subfamily) & nzchar(pool$te_subfamily),
                      paste(pool$te_family, pool$te_subfamily, sep = "/"),
                      ifelse(!is.na(pool$te_family) & nzchar(pool$te_family),
                             pool$te_family, pool$te_id))
      alternatives <- paste(sort(unique(label[close])), collapse = ";")
      status <- if (length(unique(label[close])) > 1L) "ambiguous" else "unique"
    }
    data.frame(
      evidence_id = id, classification_status = status, te_id = best$te_id,
      te_family = best$te_family, te_subfamily = best$te_subfamily,
      te_strand = best$te_strand, te_mapq = best$te_mapq,
      te_aligned_bp = best$te_aligned_bp,
      te_query_fraction = best$te_query_fraction, te_identity = best$te_identity,
      alignment_score = best$alignment_score, alternatives = alternatives,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Classify long-read insertion evidence against a TE reference
#'
#' Runs minimap2 against a TE FASTA and retains ambiguity rather than forcing
#' near-tied hits into one TE family.
#'
#' @param evidence Output from [extractLongReadEvidence()].
#' @param teReference TE reference FASTA.
#' @param teMetadata Optional data frame or TSV with te_id, te_family, and te_subfamily.
#' @param platform Long-read platform: ont, hifi, or clr.
#' @param minimap2 Path or command name for minimap2.
#' @param threads Number of minimap2 threads.
#' @param minAlignedBp Minimum aligned query bases.
#' @param minQueryFraction Minimum fraction of the candidate aligned to a TE.
#' @param minIdentity Minimum PAF match identity.
#' @param minTeMapQ Minimum TE alignment mapping quality.
#' @param ambiguityFraction Hits within this fraction of the best score are compared.
#' @return A TEiLongReadClassification list containing classifications and all
#'   PAF hits. Classification `te_strand` values are relative to the extracted
#'   read sequence.
#' @export
classifyLongReadEvidence <- function(evidence, teReference, teMetadata = NULL,
                                     platform = c("ont", "hifi", "clr"),
                                     minimap2 = "minimap2", threads = 1L,
                                     minAlignedBp = 80L, minQueryFraction = 0.5,
                                     minIdentity = 0.7, minTeMapQ = 0L,
                                     ambiguityFraction = 0.95) {
  required <- c("evidence_id", "sequence")
  if (!is.data.frame(evidence) || !all(required %in% names(evidence))) {
    stop("evidence must contain evidence_id and sequence columns.", call. = FALSE)
  }
  if (anyDuplicated(evidence$evidence_id) || any(!nzchar(evidence$evidence_id))) {
    stop("evidence_id values must be non-empty and unique.", call. = FALSE)
  }
  if (length(teReference) != 1L || !file.exists(teReference)) {
    stop("teReference must be an existing FASTA file.", call. = FALSE)
  }
  platform <- match.arg(platform)
  threads <- .tei_lr_check_scalar(threads, "threads", 1, TRUE)
  minAlignedBp <- .tei_lr_check_scalar(minAlignedBp, "minAlignedBp", 1, TRUE)
  minQueryFraction <- .tei_lr_check_scalar(minQueryFraction, "minQueryFraction", 0)
  minIdentity <- .tei_lr_check_scalar(minIdentity, "minIdentity", 0)
  minTeMapQ <- .tei_lr_check_scalar(minTeMapQ, "minTeMapQ", 0, TRUE)
  ambiguityFraction <- .tei_lr_check_scalar(ambiguityFraction, "ambiguityFraction", 0)
  if (minQueryFraction > 1 || minIdentity > 1 || ambiguityFraction > 1) {
    stop("Fraction thresholds must not exceed 1.", call. = FALSE)
  }
  metadata <- .tei_lr_read_metadata(teMetadata)
  if (!nrow(evidence)) {
    output <- list(classifications = .tei_lr_empty_classifications(),
                   hits = .tei_lr_empty_hits(), parameters = list(platform = platform))
    class(output) <- "TEiLongReadClassification"
    return(output)
  }
  executable <- if (grepl("[/\\\\]", minimap2)) minimap2 else Sys.which(minimap2)
  if (!nzchar(executable) || !file.exists(executable)) {
    stop("minimap2 was not found. Install it or provide its executable path.", call. = FALSE)
  }
  queryFile <- tempfile(fileext = ".fa")
  pafFile <- tempfile(fileext = ".paf")
  stderrFile <- tempfile(fileext = ".log")
  on.exit(unlink(c(queryFile, pafFile, stderrFile)), add = TRUE)
  fasta <- as.vector(rbind(paste0(">", evidence$evidence_id), toupper(evidence$sequence)))
  writeLines(fasta, queryFile, useBytes = TRUE)
  preset <- c(ont = "map-ont", hifi = "map-hifi", clr = "map-pb")[[platform]]
  arguments <- c("-c", "-x", preset, "--secondary=yes", "-N", "20", "-t",
                 as.character(threads), shQuote(normalizePath(teReference)), shQuote(queryFile))
  status <- system2(executable, arguments, stdout = pafFile, stderr = stderrFile)
  if (!identical(status, 0L)) {
    message <- paste(readLines(stderrFile, warn = FALSE), collapse = "\n")
    stop("minimap2 failed", if (nzchar(message)) paste0(": ", message) else ".", call. = FALSE)
  }
  hits <- .tei_lr_parse_paf_lines(readLines(pafFile, warn = FALSE))
  classifications <- .tei_lr_select_classifications(
    evidence, hits, metadata, minAlignedBp, minQueryFraction, minIdentity,
    minTeMapQ, ambiguityFraction
  )
  if (nrow(hits) && !is.null(metadata)) {
    index <- match(hits$te_id, metadata$te_id)
    hits$te_family <- metadata$te_family[index]
    hits$te_subfamily <- metadata$te_subfamily[index]
  }
  output <- list(
    classifications = classifications, hits = hits,
    parameters = list(platform = platform, minAlignedBp = minAlignedBp,
                      minQueryFraction = minQueryFraction, minIdentity = minIdentity,
                      minTeMapQ = minTeMapQ, ambiguityFraction = ambiguityFraction)
  )
  class(output) <- "TEiLongReadClassification"
  output
}

.tei_lr_empty_calls <- function() {
  data.frame(
    call_id = character(), seqname = character(), left_breakpoint_0 = numeric(),
    right_breakpoint_0 = numeric(), te_id = character(), te_family = character(),
    te_subfamily = character(), te_strand = character(), supporting_reads = integer(),
    cigar_support = integer(), split_support = integer(), clip_support = integer(),
    median_inserted_length = numeric(), median_te_identity = numeric(),
    read_ids = character(), filter = character(), stringsAsFactors = FALSE
  )
}

#' Summarize classified long-read evidence into insertion calls
#'
#' @param evidence Output from [extractLongReadEvidence()].
#' @param classification Output from [classifyLongReadEvidence()] or its classifications table.
#' @param clusterWindow Maximum breakpoint distance within a call cluster.
#' @param minSupport Unique supporting reads required for PASS.
#' @param includeAmbiguous Include evidence with ambiguous TE classifications.
#' @return A data frame with one row per insertion call. Call `te_strand`
#'   values are oriented relative to the host reference genome.
#' @export
summarizeLongReadInsertions <- function(evidence, classification, clusterWindow = 50L,
                                        minSupport = 2L, includeAmbiguous = FALSE) {
  clusterWindow <- .tei_lr_check_scalar(clusterWindow, "clusterWindow", 0, TRUE)
  minSupport <- .tei_lr_check_scalar(minSupport, "minSupport", 1, TRUE)
  if (inherits(classification, "TEiLongReadClassification") ||
      (is.list(classification) && "classifications" %in% names(classification))) {
    classification <- classification$classifications
  }
  if (!is.data.frame(evidence) || !is.data.frame(classification) ||
      !all(c("evidence_id", "read_id", "evidence_type", "seqname", "host_strand",
             "left_breakpoint_0", "right_breakpoint_0", "evidence_length") %in% names(evidence)) ||
      !all(c("evidence_id", "classification_status", "te_id", "te_family",
             "te_subfamily", "te_strand", "te_identity") %in% names(classification))) {
    stop("evidence or classification has an invalid schema.", call. = FALSE)
  }
  if (!nrow(evidence) || !nrow(classification)) return(.tei_lr_empty_calls())
  index <- match(evidence$evidence_id, classification$evidence_id)
  data <- cbind(evidence, classification[index, setdiff(names(classification), "evidence_id"),
                                         drop = FALSE])
  accepted <- data$classification_status == "unique"
  if (isTRUE(includeAmbiguous)) accepted <- accepted | data$classification_status == "ambiguous"
  data <- data[!is.na(accepted) & accepted & !is.na(data$te_id), , drop = FALSE]
  if (!nrow(data)) return(.tei_lr_empty_calls())
  if (any(!data$host_strand %in% c("+", "-")) || any(!data$te_strand %in% c("+", "-"))) {
    stop("host_strand and te_strand must contain only + or -.", call. = FALSE)
  }
  data$te_strand <- ifelse(
    data$host_strand == "-",
    ifelse(data$te_strand == "+", "-", "+"),
    data$te_strand
  )
  family <- ifelse(!is.na(data$te_family) & nzchar(data$te_family), data$te_family, data$te_id)
  midpoint <- (data$left_breakpoint_0 + data$right_breakpoint_0) / 2
  group <- paste(data$seqname, family, data$te_strand, sep = "\r")
  rows <- list()
  for (groupName in unique(group)) {
    selected <- which(group == groupName)
    selected <- selected[order(midpoint[selected], data$read_id[selected])]
    cluster <- cumsum(c(TRUE, diff(midpoint[selected]) > clusterWindow))
    for (clusterId in unique(cluster)) {
      x <- data[selected[cluster == clusterId], , drop = FALSE]
      uniqueValue <- function(value) {
        value <- unique(value[!is.na(value) & nzchar(value)])
        if (length(value) == 1L) value else NA_character_
      }
      readIds <- sort(unique(x$read_id))
      support <- length(readIds)
      rows[[length(rows) + 1L]] <- data.frame(
        call_id = NA_character_, seqname = x$seqname[1L],
        left_breakpoint_0 = round(stats::median(x$left_breakpoint_0)),
        right_breakpoint_0 = round(stats::median(x$right_breakpoint_0)),
        te_id = uniqueValue(x$te_id), te_family = uniqueValue(x$te_family),
        te_subfamily = uniqueValue(x$te_subfamily), te_strand = uniqueValue(x$te_strand),
        supporting_reads = as.integer(support),
        cigar_support = as.integer(length(unique(x$read_id[x$evidence_type == "cigar_insertion"]))),
        split_support = as.integer(length(unique(x$read_id[x$evidence_type == "split_gap"]))),
        clip_support = as.integer(length(unique(x$read_id[x$evidence_type %in%
                                                             c("left_soft_clip", "right_soft_clip")]))),
        median_inserted_length = stats::median(x$evidence_length),
        median_te_identity = stats::median(x$te_identity, na.rm = TRUE),
        read_ids = paste(readIds, collapse = ","),
        filter = if (support >= minSupport) "PASS" else "LOW_SUPPORT",
        stringsAsFactors = FALSE
      )
    }
  }
  output <- do.call(rbind, rows)
  output <- output[order(output$seqname, output$left_breakpoint_0, output$te_family,
                         output$te_id), , drop = FALSE]
  output$call_id <- sprintf("TEI_CALL_%06d", seq_len(nrow(output)))
  rownames(output) <- NULL
  output
}

#' Detect candidate TE insertions from long-read BAM alignments
#'
#' @inheritParams extractLongReadEvidence
#' @inheritParams classifyLongReadEvidence
#' @inheritParams summarizeLongReadInsertions
#' @return A TEiLongReadResult list containing calls, evidence, TE hits, and limitations.
#' @export
detectLongReadInsertions <- function(alignment, teReference, teMetadata = NULL,
                                     platform = c("ont", "hifi", "clr"),
                                     minimap2 = "minimap2", threads = 1L,
                                     minMapQ = 20L, minLength = 100L,
                                     maxReferenceGap = 50L, includeSoftClips = TRUE,
                                     includeSplitGaps = TRUE, chunkSize = 10000L,
                                     minAlignedBp = 80L, minQueryFraction = 0.5,
                                     minIdentity = 0.7, minTeMapQ = 0L,
                                     ambiguityFraction = 0.95, clusterWindow = 50L,
                                     minSupport = 2L, includeAmbiguous = FALSE) {
  platform <- match.arg(platform)
  evidence <- extractLongReadEvidence(
    alignment, minMapQ, minLength, maxReferenceGap, includeSoftClips,
    includeSplitGaps, chunkSize
  )
  classification <- classifyLongReadEvidence(
    evidence, teReference, teMetadata, platform, minimap2, threads, minAlignedBp,
    minQueryFraction, minIdentity, minTeMapQ, ambiguityFraction
  )
  calls <- summarizeLongReadInsertions(
    evidence, classification, clusterWindow, minSupport, includeAmbiguous
  )
  output <- list(
    calls = calls, evidence = evidence,
    classifications = classification$classifications, hits = classification$hits,
    parameters = c(classification$parameters,
                   list(minMapQ = minMapQ, minLength = minLength,
                        maxReferenceGap = maxReferenceGap, clusterWindow = clusterWindow,
                        minSupport = minSupport)),
    limitations = .tei_lr_limitations
  )
  class(output) <- "TEiLongReadResult"
  output
}
