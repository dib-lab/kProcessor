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

/// Load the kmers in the input file into the output kDataframe. Input File can be of formats: fastq,fasta.
    void countKmersFromFile(kDataFrame *kframe, std::map<std::string, int> parse_params, string filename,
                            int chunk_size = 1000);

/// Load the kmers from KMC DB.
    void loadFromKMC(kDataFrame *kframe, std::string KMC_DB_filename);

/// Load the kmers in the input string into the output kDataframe.
    void countKmersFromString(kmerDecoder *KD, string sequence, kDataFrame *output);

/// Load the kmers in the input string into the output kDataframe.
    void countKmersFromString(kDataFrame *frame, std::map<std::string, int> parse_params, string sequence);

/// Applies a function on each row in the input kDataframe. The output is another kDataframe with the same kmers and new columns defined based on the user defined function.
    kDataFrame *transform(kDataFrame *input,  function<void (kDataFrameIterator& i)>);

/// Applies a function on each row in the input kDataframe. the function updates one or more columns in the input kdataframe based on the user defined function.
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

    /* create a prefix forest from an kdataframe index. kdataframe should have columns of type deduplicatedColumn<prefixTrie,phmap::flat_hash_map<uint32_t,uint32_t>.
     * the name of the columns are in the form of color<id>. I am using the id to sort the prefix tries in the forest.
     * the index is updated by the new forest and the prefix tries are deleted.
    */
    void createPrefixForest(kDataFrame* index, string tmpFolder,uint32_t num_vectors,uint32_t vector_size);

    /* Joins the input kdataframes into one kDataframe. All the columns in the input kdataframes will be copied to the output kdataframe and they will have new name in the format of "<inputColumnName><index in the input>".
     * kmersToKeep range is from 0-input.size() and it should be unique. kmers from kDataframes whose index exists in kmersToKeep will be inserted in the output dataframe.
     *  please note that the kdataframe are expected to be ordered. If you want to join unordered kdataframe then use parallelJoin
     *  example:
     *  kdataframe 0:                           kdataframe 1:               kdataframe 2:
     *  +-------+-------+-------+------------+  +-------+-------+-------+  +-------+----------+----------+
     *  |       | count | color | foldChange |  |       | count | color |  |       | isCoding | isRepeat |
     *  +-------+-------+-------+------------+  +-------+-------+-------+  +-------+----------+----------+
     *  | ACGAT | 10    | [1]   | 0.2        |  | ACGAT | 20    | [1]   |  | ACGAT | 0        | 1        |
     *  +-------+-------+-------+------------+  +-------+-------+-------+  +-------+----------+----------+
     *  | AGCAC | 15    | [1,3] | 0.1        |  | ATTTC | 117   | [2]   |  | ACCCC | 0        | 0        |
     *  +-------+-------+-------+------------+  +-------+-------+-------+  +-------+----------+----------+
     *  | ATCAC | 14    | [2]   | 0.5        |                             | AGCAC | 1        | 1        |
     *  +-------+-------+-------+------------+                             +-------+----------+----------+
     *                                                                     | ATTTC | 0        | 0        |
     *                                                                     +-------+----------+----------+
     *
     *  kmersToKeep = [0,1]
     *
     *  Output KDataframe
     * +-------+--------+--------+-------------+--------+--------+-----------+----------+
     * |       | count0 | color0 | foldChange0 | count1 | color1 | isCoding2 | isRepeat |
     * +-------+--------+--------+-------------+--------+--------+-----------+----------+
     * | ACGAT | 10     | [1]    | 0.2         | 20     | [1]    | 0         | 1        |
     * +-------+--------+--------+-------------+--------+--------+-----------+----------+
     * | AGCAC | 15     | [1,3]  | 0.1         | 0      | []     | 1         | 1        |
     * +-------+--------+--------+-------------+--------+--------+-----------+----------+
     * | ATCAC | 14     | [2]    | 0.5         | 0      | []     | 0         | 0        |
     * +-------+--------+--------+-------------+--------+--------+-----------+----------+
     * | ATTTC | 0      | []     | 0.0         | 117    | [2]    | 0         | 0        |
     * +-------+--------+--------+-------------+--------+--------+-----------+----------+
    */
    kDataFrame* innerJoin(vector<kDataFrame *> input, vector<uint32_t> kmersToKeep);
    //colored_kDataFrame * indexPriorityQueue(kmerDecoder *KD, string names_fileName, kDataFrame *frame);
    //colored_kDataFrame * indexPriorityQueue2(kmerDecoder *KD, string names_fileName, kDataFrame *frame);

    /*
     * Exact same behavior as innerJoin with two modifications: input kDataframes are not required to be sorted, and multithreaded implementation is provided.
     */
    kDataFrame* parallelJoin(vector<string>& kdataframeFileNames, vector<uint32_t> kmersToKeep,uint64_t numThreads=1);
    
}
#endif
