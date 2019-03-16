#ifndef KmerCounter_HPP
#define KmerCounter_HPP
#include <stdint.h>
#include <string>
#include <gqf.hpp>
#include "HashUtils/hashutil.h"
#include <math.h>


void loadIntoMQF(std::string sequenceFilename, int k,int noThreads,Hasher *hasher,QF * memoryMQF,QF * diskMQF=NULL);

void dumpMQF(QF * memoryMQF,int ksize,std::string outputFilename);

void estimateMemRequirement(std::string ntcardFilename,
   uint64_t numHashBits,uint64_t tagSize,
   uint64_t *res_noSlots,uint64_t *res_fixedSizeCounter
   , uint64_t *res_memory);

inline uint64_t estimateMemory(uint64_t nslots,uint64_t slotSize, uint64_t fcounter, uint64_t tagSize)
   {
     uint64_t SLOTS_PER_BLOCK=64;
     uint64_t xnslots = nslots + 10*sqrt((double)nslots);
   	uint64_t nblocks = (xnslots + SLOTS_PER_BLOCK - 1) / SLOTS_PER_BLOCK;
     uint64_t blocksize=17;

     return ((nblocks)*(blocksize+8*(slotSize+fcounter+tagSize)))/1024;

   }
#endif
