#ifndef KmerCounter_HPP
#define KmerCounter_HPP
#include <stdint.h>
#include <string>
#include <gqf.h>
#include "../HashUtils/hashutil.h"


void loadIntoMQF(std::string sequenceFilename, int k,int noThreads,Hasher *hasher,QF * memoryMQF);

void dumpMQF(QF * memoryMQF,int ksize,std::string outputFilename);

void estimateMemRequirement(std::string ntcardFilename,
   uint64_t numHashBits,uint64_t tagSize,
   uint64_t *res_noSlots,uint64_t *res_fixedSizeCounter
   , uint64_t *res_memory);
#endif
