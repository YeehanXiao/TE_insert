#include "teseq.h"

#include <Rcpp.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct HtsFileCloser {
  void operator()(htsFile* file) const {
    if (file != NULL) hts_close(file);
  }
};

struct BamHeaderCloser {
  void operator()(bam_hdr_t* header) const {
    if (header != NULL) bam_hdr_destroy(header);
  }
};

struct BamRecordCloser {
  void operator()(bam1_t* record) const {
    if (record != NULL) bam_destroy1(record);
  }
};

typedef std::unique_ptr<htsFile, HtsFileCloser> HtsFilePtr;
typedef std::unique_ptr<bam_hdr_t, BamHeaderCloser> BamHeaderPtr;
typedef std::unique_ptr<bam1_t, BamRecordCloser> BamRecordPtr;

const char kBaseLookup[] = "NACNGNNNTNNNNNNN";

bool isAlignedBaseOperation(const int op) {
  return op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF;
}

bool validRecord(const bam1_t* record, const bam_hdr_t* header) {
  const bam1_core_t& core = record->core;
  return !(core.flag & BAM_FUNMAP) && core.tid >= 0 &&
         core.tid < header->n_targets && core.pos >= 0;
}

bool validCigar(const bam1_t* record) {
  const bam1_core_t& core = record->core;
  const uint32_t* cigar = bam_get_cigar(record);
  int64_t query_length = 0;
  for (uint32_t i = 0; i < core.n_cigar; ++i) {
    const int op = bam_cigar_op(cigar[i]);
    if (op < BAM_CMATCH || op > BAM_CBACK) return false;
    if (bam_cigar_type(op) & 1) {
      query_length += static_cast<int64_t>(bam_cigar_oplen(cigar[i]));
      if (query_length > core.l_qseq) return false;
    }
  }
  return query_length == core.l_qseq;
}

char baseAt(const bam1_t* record, const int64_t index) {
  return kBaseLookup[bam_seqi(bam_get_seq(record), index)];
}

char qualityAt(const bam1_t* record, const int64_t index) {
  const uint8_t quality = bam_get_qual(record)[index];
  return quality == 0xff ? '!' : static_cast<char>(std::min<int>(quality, 93) + 33);
}

void writeFastqRecord(std::ostream& output, const bam1_t* record,
                      const bam_hdr_t* header, const int64_t query_offset,
                      const uint32_t clip_length, const int64_t position,
                      const std::string& tsd, const int type,
                      const uint64_t segment_id) {
  const char* reference = header->target_name[record->core.tid];
  if (reference == NULL || reference[0] == '\0') return;

  output << '@' << reference << ':' << position << ':' << tsd << ':' << type
         << ':' << segment_id << '\n';
  for (uint32_t i = 0; i < clip_length; ++i) {
    output.put(baseAt(record, query_offset + i));
  }
  output << "\n+\n";
  for (uint32_t i = 0; i < clip_length; ++i) {
    output.put(qualityAt(record, query_offset + i));
  }
  output.put('\n');
}

void extractSoftClips(const bam1_t* record, const bam_hdr_t* header,
                      std::ostream& output, const int minimum_length,
                      const int maximum_tsd, uint64_t& segment_id) {
  const bam1_core_t& core = record->core;
  const uint32_t* cigar = bam_get_cigar(record);
  int64_t query_offset = 0;
  int64_t reference_position = core.pos;

  for (uint32_t i = 0; i < core.n_cigar; ++i) {
    const int op = bam_cigar_op(cigar[i]);
    const uint32_t op_length = bam_cigar_oplen(cigar[i]);

    if (op == BAM_CSOFT_CLIP && op_length >= static_cast<uint32_t>(minimum_length)) {
      int type = 0;
      uint32_t adjacent_length = 0;
      if (i > 0 && isAlignedBaseOperation(bam_cigar_op(cigar[i - 1]))) {
        type = 2;
        adjacent_length = bam_cigar_oplen(cigar[i - 1]);
      } else if (i + 1 < core.n_cigar &&
                 isAlignedBaseOperation(bam_cigar_op(cigar[i + 1]))) {
        type = 1;
        adjacent_length = bam_cigar_oplen(cigar[i + 1]);
      }

      if (type != 0 && query_offset + op_length <= core.l_qseq) {
        const uint32_t tsd_length = std::min<uint32_t>(maximum_tsd, adjacent_length);
        std::string tsd;
        tsd.reserve(tsd_length);
        if (type == 1 && query_offset + op_length + tsd_length <= core.l_qseq) {
          for (uint32_t j = 0; j < tsd_length; ++j) {
            tsd.push_back(baseAt(record, query_offset + op_length + j));
          }
          writeFastqRecord(output, record, header, query_offset, op_length,
                           reference_position, tsd, type, segment_id++);
        } else if (type == 2 && query_offset >= tsd_length) {
          for (uint32_t j = tsd_length; j > 0; --j) {
            tsd.push_back(baseAt(record, query_offset - j));
          }
          writeFastqRecord(output, record, header, query_offset, op_length,
                           reference_position, tsd, type, segment_id++);
        }
      }
    }

    if (bam_cigar_type(op) & 1) query_offset += op_length;
    if (bam_cigar_type(op) & 2) reference_position += op_length;
  }
}

struct AlignmentGroup {
  std::string query_name;
  std::unordered_map<int32_t, std::size_t> target_counts;
  std::size_t total;

  AlignmentGroup() : total(0) {}

  void clear(const std::string& name) {
    query_name = name;
    target_counts.clear();
    total = 0;
  }
};

std::vector<std::string> splitQueryName(const std::string& query_name) {
  std::vector<std::string> fields;
  std::size_t start = 0;
  for (;;) {
    const std::size_t separator = query_name.find(':', start);
    fields.push_back(query_name.substr(start, separator - start));
    if (separator == std::string::npos) break;
    start = separator + 1;
  }
  return fields;
}

bool parseInteger(const std::string& value, int64_t& output) {
  if (value.empty()) return false;
  char* end = NULL;
  errno = 0;
  const long long parsed = std::strtoll(value.c_str(), &end, 10);
  if (errno == ERANGE || end == value.c_str() || *end != '\0') return false;
  output = static_cast<int64_t>(parsed);
  return true;
}

bool safeField(const std::string& value) {
  return value.find_first_of("\t\r\n") == std::string::npos;
}

void writeAlignmentGroup(const AlignmentGroup& group, const bam_hdr_t* header,
                         std::ostream& output, const double threshold) {
  if (group.total == 0) return;

  int32_t selected_target = -1;
  double selected_ratio = 0.0;
  std::size_t targets_above_threshold = 0;
  for (std::unordered_map<int32_t, std::size_t>::const_iterator it =
           group.target_counts.begin();
       it != group.target_counts.end(); ++it) {
    const double fraction = static_cast<double>(it->second) / group.total;
    if (fraction > threshold) {
      selected_target = it->first;
      selected_ratio = fraction;
      ++targets_above_threshold;
    }
  }
  if (targets_above_threshold != 1 || selected_target < 0 ||
      selected_target >= header->n_targets) {
    return;
  }

  const std::vector<std::string> fields = splitQueryName(group.query_name);
  if ((fields.size() != 4 && fields.size() != 5) || fields[0].empty() ||
      !safeField(fields[0]) ||
      !safeField(fields[2])) {
    return;
  }

  int64_t position = 0;
  int64_t type = 0;
  if (!parseInteger(fields[1], position) || !parseInteger(fields[3], type) ||
      position < 0 || (type != 1 && type != 2)) {
    return;
  }

  int64_t start = position;
  int64_t end = position;
  if (type == 1) {
    if (position == std::numeric_limits<int64_t>::max()) return;
    end = position + 1;
  } else {
    if (position == 0) return;
    start = position - 1;
  }

  const char* target = header->target_name[selected_target];
  if (target == NULL || target[0] == '\0' ||
      !safeField(std::string(target))) {
    return;
  }

  output << fields[0] << '\t' << start << '\t' << end << '\t' << target << '\t'
         << fields[2] << '\t' << type << "\t1\t" << std::fixed
         << std::setprecision(2) << selected_ratio << '\n';
}

}  // namespace

// [[Rcpp::export]]
void processbam(std::string bamfile, std::string outfq, double quantile,
                int length, int tsd) {
  if (bamfile.empty() || outfq.empty()) Rcpp::stop("Input and output paths must not be empty.");
  if (!std::isfinite(quantile) || quantile < 0.0 || quantile > 255.0) {
    Rcpp::stop("quantile must be between 0 and 255.");
  }
  if (length < 1) Rcpp::stop("length must be at least 1.");
  if (tsd < 0) Rcpp::stop("tsd must be non-negative.");

  HtsFilePtr input(sam_open(bamfile.c_str(), "r"));
  if (!input) Rcpp::stop("Cannot open alignment file: %s", bamfile.c_str());
  BamHeaderPtr header(sam_hdr_read(input.get()));
  if (!header) Rcpp::stop("Cannot read alignment header: %s", bamfile.c_str());
  if (header->n_targets <= 0 || header->target_name == NULL) {
    Rcpp::stop("Alignment header contains no reference sequences: %s", bamfile.c_str());
  }

  std::ofstream output(outfq.c_str(), std::ios::out | std::ios::binary);
  if (!output) Rcpp::stop("Cannot open FASTQ output: %s", outfq.c_str());
  BamRecordPtr record(bam_init1());
  if (!record) Rcpp::stop("Cannot allocate an alignment record.");

  uint64_t segment_id = 0;
  int status = 0;
  while ((status = sam_read1(input.get(), header.get(), record.get())) >= 0) {
    if (!validRecord(record.get(), header.get()) ||
        record->core.n_cigar == 0 || record->core.l_qseq <= 0 ||
        !validCigar(record.get()) ||
        record->core.flag &
            (BAM_FDUP | BAM_FQCFAIL | BAM_FSECONDARY | BAM_FSUPPLEMENTARY) ||
        static_cast<double>(record->core.qual) < quantile) {
      continue;
    }
    extractSoftClips(record.get(), header.get(), output, length, tsd,
                     segment_id);
  }
  if (status < -1) Rcpp::stop("Failed while reading alignment file: %s", bamfile.c_str());
  output.close();
  if (!output) Rcpp::stop("Failed while writing FASTQ output: %s", outfq.c_str());
}

// [[Rcpp::export]]
void processSam2bed(std::string alignmentfile, std::string outbedfile,
                    float ratio) {
  if (alignmentfile.empty() || outbedfile.empty()) {
    Rcpp::stop("Input and output paths must not be empty.");
  }
  if (!std::isfinite(ratio) || ratio < 0.0f || ratio > 1.0f) {
    Rcpp::stop("ratio must be between 0 and 1.");
  }

  HtsFilePtr input(sam_open(alignmentfile.c_str(), "r"));
  if (!input) Rcpp::stop("Cannot open alignment file: %s", alignmentfile.c_str());
  BamHeaderPtr header(sam_hdr_read(input.get()));
  if (!header) Rcpp::stop("Cannot read alignment header: %s", alignmentfile.c_str());
  if (header->n_targets <= 0 || header->target_name == NULL) {
    Rcpp::stop("Alignment header contains no reference sequences: %s",
               alignmentfile.c_str());
  }

  std::ofstream output(outbedfile.c_str(), std::ios::out | std::ios::binary);
  if (!output) Rcpp::stop("Cannot open BED output: %s", outbedfile.c_str());
  BamRecordPtr record(bam_init1());
  if (!record) Rcpp::stop("Cannot allocate an alignment record.");

  AlignmentGroup group;
  int status = 0;
  while ((status = sam_read1(input.get(), header.get(), record.get())) >= 0) {
    if (!validRecord(record.get(), header.get())) continue;
    const char* query_name = bam_get_qname(record.get());
    if (query_name == NULL || query_name[0] == '\0') continue;

    if (group.total == 0) {
      group.clear(query_name);
    } else if (group.query_name != query_name) {
      writeAlignmentGroup(group, header.get(), output, ratio);
      group.clear(query_name);
    }
    ++group.target_counts[record->core.tid];
    ++group.total;
  }
  if (status < -1) {
    Rcpp::stop("Failed while reading alignment file: %s", alignmentfile.c_str());
  }
  writeAlignmentGroup(group, header.get(), output, ratio);
  output.close();
  if (!output) Rcpp::stop("Failed while writing BED output: %s", outbedfile.c_str());
}
