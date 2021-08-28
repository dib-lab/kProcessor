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
#include <gqf.h>
#include "kDataFrame.hpp"
#include <math.h>
#include <vector>
#include <map>
#include <any>
#include <functional>


namespace kProcessor {

    // TO BE REMOVED TODO V2
    // void loadIntoMQF(std::string sequenceFilename, int k, int noThreads, Hasher *hasher, QF *memoryMQF, QF *diskMQF = NULL);

    void dumpMQF(QF *memoryMQF, int ksize, std::string outputFilename);

    void estimateMemRequirement(std::string ntcardFilename,
                                uint64_t numHashBits, uint64_t tagSize,
                                uint64_t *res_noSlots, uint64_t *res_fixedSizeCounter, uint64_t *res_memory);

    inline uint64_t estimateMemory(uint64_t nslots, uint64_t slotSize, uint64_t fcounter, uint64_t tagSize) {
        uint64_t SLOTS_PER_BLOCK_t = 64;
        uint64_t xnslots = nslots + 10 * sqrt((double) nslots);
        uint64_t nblocks = (xnslots + SLOTS_PER_BLOCK_t - 1) / SLOTS_PER_BLOCK_t;
        uint64_t blocksize = 17;

        return ((nblocks) * (blocksize + 8 * (slotSize + fcounter + tagSize))) / 1024;

    }

    vector<uint64_t> estimateKmersHistogram(string fileName, int kSize, int threads);

// Load the kmers in the input file into the output kDataframe. Input File can be of formats: fastq,fasta, sam, and bam.
// void parseSequences(string seqFileName, int nThreads, kDataFrame *output);

/// Load the kmers in the input file into the output kDataframe. Input File can be of formats: fastq,fasta.
    void parseSequences(kmerDecoder *KD, kDataFrame *output);

/// Load the kmers in the input file into the output kDataframe. Input File can be of formats: fastq,fasta.
    void countKmersFromFile(kDataFrame *kframe, std::map<std::string, int> parse_params, string filename,
                            int chunk_size = 1000);

/// Load the kmers from KMC DB.
    void loadFromKMC(kDataFrame *kframe, std::string KMC_DB_filename);

/// Load the kmers in the input string into the output kDataframe.
    void countKmersFromString(kmerDecoder *KD, string sequence, kDataFrame *output);

/// Load the kmers in the input string into the output kDataframe.
    void countKmersFromString(kDataFrame *frame, std::map<std::string, int> parse_params, string sequence);

/// Applies a function on all the kmers in the input kDataframe. The output is another kDataframe with the transformed kmers.
    kDataFrame *transform(kDataFrame *input,  function<kmerRow (kmerRow i)>);

    void transformInPlace(kDataFrame *input,  function<void (kDataFrameIterator& i)>);

/// filter the kmers in the kdataframe. The output is another kDataframe with the filtered kmers.
    kDataFrame *filter(kDataFrame *input,  function<bool (kmerRow& i)>);
    kDataFrame *filter(kDataFrame *input,  function<bool (kDataFrameIterator& i)>);

/// aggregate all the kmers in the kdataframe into a single value. The output is one value.
    any aggregate(kDataFrame *input, any initial,  function<any (kmerRow i, any v)>);
    any aggregate(kDataFrame *input, any initial,  function<any (kDataFrameIterator& i, any v)>);

/*! Merge the a list of kDataframes into a one.
 *
 * The Input is a list of kDataframes and a function to decide how to merge the kmers. The merged kmers will be inserted into the result kDataframe.
 * The function to decide the merging behavior takes a list of kmerRows,all have the same kmer hashedKmer, as an input and output one kmer row to be inserted in the result kDataFrame. If the returned kmer row has count equal zero nothing to be inserted.
 * The function will be called repeatedly through merging. This function is used to implement the set functions.
 */
    void merge(const vector<kDataFrame *> &input, kDataFrame *result, function<kmerRow (vector<kDataFrameIterator*> &i)>);

/// Calculate the union of the kDataFrames. The result kDataframe will have all the kmers in the input list of kDataframes. The count of the kmers equals to the sum of the kmer count in the input list.
    kDataFrame *kFrameUnion(const vector<kDataFrame *> &input);

/// Calculate the intersect of the kDataFrames. The result kDataframe will have only kmers that exists in all the kDataframes. The count of the kmers equals to the minimum of the kmer count in the input list.
    kDataFrame *kFrameIntersect(const vector<kDataFrame *> &input);

/// Calculate the difference of the kDataframes. The result kDataframe will have only kmers that exists in the first kDataframe and not in any of the rest input kDataframes. The count of the kmers equals to the count in the first kDataframe.
    kDataFrame *kFrameDiff(const vector<kDataFrame *> &input);

/// Initialize kmerDecoder to decode kmers from a FASTA/Q file with predefined mode.
    kmerDecoder *initialize_kmerDecoder(std::string filename, int chunkSize, std::string mode,
                                        std::map<std::string, int> parse_params);

/// Initialize kmerDecoder to decode kmers from a sequence string with predefined mode.
    kmerDecoder *initialize_kmerDecoder(std::string mode, std::map<std::string, int> parse_params);

// Initialize kmerDecoder to hash and Ihash single kmer.
    kmerDecoder *initialize_kmerDecoder(int kSize, hashingModes HM);

/* set hashing mode for kmerDecoder object
 *
 * Mode 0: Murmar Hashing | Irreversible
 * Mode 1: Integer Hashing | Reversible | Full Hashing
 * Mode 2: TwoBitsHashing | Not considered hashing, just store the two bits representation
*/
    void kmerDecoder_setHashing(kmerDecoder *KD, hashingModes hash_mode);

    void kmerDecoder_setHashing(kDataFrame *KF, hashingModes hash_mode);

/// Perform indexing to a sequences file with predefined kmers decoding mode, returns a colored kDataframe.
    void index(kmerDecoder *KD, string names_fileName, kDataFrame *frame);
    void indexMega(string fastaFileName, string tmpFolder, kDataFrame *frame);

/// Index function without needing the kmerDecoder
    void
    index(kDataFrame *frame, std::map<std::string, int> parse_params, string filename, int chunks,
          string names_fileName);

/*
 * Create an index for the kDataframes in the vector<input>. the index is in the form of kDataframe and it has default column of type mix Vector
 * mix vector split the colors into num_vectors vectors where each has maximum size vector_size.
 */
    void indexPriorityQueue(vector<kDataFrame *> &input, string tmpFolder, kDataFrame *output,uint32_t num_vectors=20,uint32_t vector_size=1000000);

    void mergeIndexes(vector<kDataFrame *> &input, string tmpFolder, kDataFrame *output);

    kDataFrame* innerJoin(vector<kDataFrame *> input, vector<uint32_t> kmersToKeep);
    //colored_kDataFrame * indexPriorityQueue(kmerDecoder *KD, string names_fileName, kDataFrame *frame);
    //colored_kDataFrame * indexPriorityQueue2(kmerDecoder *KD, string names_fileName, kDataFrame *frame);

    
}
#endif
