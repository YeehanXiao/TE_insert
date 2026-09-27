# TEi 1.1.0.9000

- Restored the complete upstream history and joint authorship of Tao Chen and
  Yihan Xiao.
- Reworked native BAM/SAM handling for portable builds and deterministic error
  reporting.
- Fixed transposable element multimap counting, final query-group handling,
  terminal soft-clip coordinates, missing qualities, and native resource
  cleanup.
- Added unique clipped-segment identifiers while retaining compatibility with
  the historical record-name format.
- Added deterministic short-read regression tests and macOS/Linux package
  checks.
- Added an experimental PacBio/ONT workflow for CIGAR insertion, terminal
  soft-clip, and same-strand split-alignment evidence.
- Added minimap2-based long-read TE classification with explicit ambiguity and
  low-support reporting.
- Normalized long-read TE orientation to the host reference strand and removed
  quadratic event-collection and classification lookups for large datasets.
- Preserved the legacy short-read `name` metadata field alongside the clearer
  `te_name` alias and corrected breakpoint-clustering boundary handling.
- Removed compiled objects, personal paths, stale dependencies, and operating
  system metadata from the maintained source tree.
