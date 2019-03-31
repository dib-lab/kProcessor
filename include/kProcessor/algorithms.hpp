/**
 * @file algorithms.hpp
 *
 * Algorithms to help the calculations on the kDataframes using the Functional programming paradigm..
 *
 */

#ifndef KmerCounter_HPP
#define KmerCounter_HPP
#include <stdint.h>
#include <string>
#include <gqf.hpp>
#include "HashUtils/hashutil.h"
#include "kDataFrame.hpp"
#include <math.h>
#include <vector>



namespace kProcessor{
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

/// Load the kmers in the input file into the output kDataframe. Input File can be of formats: fastq,fasta, sam, and bam.
void parseSequences(string seqFileName,int nThreads,kDataFrame* output);

/// Applies a function on all the kmers in the input kDataframe. The output is another kDataframe with the transformed kmers.
kDataFrame* transform(kDataFrame* input,kmerRow (*fn)(kmerRow i));

/*! Merge the a list of kDataframes into a one.
 *
 * The Input is a list of kDataframes and a function to decide how to merge the kmers. The merged kmers will be inserted into the result kDataframe.
 * The function to decide the merging behavior takes a list of kmerRows,all have the same kmer hashedKmer, as an input and output one kmer row to be inserted in the result kDataFrame. If the returned kmer row has count equal zero nothing to be inserted.
 * The function will be called repeatedly through merging. This function is used to implement the set functions.
 */
void merge(const vector<kDataFrame*>& input,kDataFrame* result,kmerRow (*fn)(vector<kmerRow>& i));

/// Calculate the union of the kDataFrames. The result kDataframe will have all the kmers in the input list of kDataframes. The count of the kmers equals to the sum of the kmer count in the input list.
kDataFrame* kFrameUnion(const vector<kDataFrame*>& input);

/// Calculate the intersect of the kDataFrames. The result kDataframe will have only kmers that exists in all the kDataframes. The count of the kmers equals to the minimum of the kmer count in the input list.
kDataFrame* kFrameIntersect(const vector<kDataFrame*>& input);

/// Calculate the difference of the kDataframes. The result kDataframe will have only kmers that exists in the first kDataframe and not in any of the rest input kDataframes. The count of the kmers equals to the count in the first kDataframe.
kDataFrame* kFrameDiff(const vector<kDataFrame*>& input);


}
#endif
