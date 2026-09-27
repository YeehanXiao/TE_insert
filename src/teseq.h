#ifndef TEI_TESEQ_H
#define TEI_TESEQ_H

#include <string>

void processbam(std::string bamfile, std::string outfq, double quantile,
                int length, int tsd);
void processSam2bed(std::string alignmentfile, std::string outbedfile,
                    float ratio);

#endif
