/**
 * @file algorithms.hpp
 *
 * Algorithms to help the calculations on the kDataframes using the Functional programming paradigm..
 *
 */

namespace kProcessor {
    void dumpMQF(QF* memoryMQF, int ksize, std::string outputFilename);

    void estimateMemRequirement(std::string ntcardFilename,
        uint64_t numHashBits, uint64_t tagSize,
        uint64_t* res_noSlots, uint64_t* res_fixedSizeCounter, uint64_t* res_memory);

    inline uint64_t estimateMemory(uint64_t nslots, uint64_t slotSize, uint64_t fcounter, uint64_t tagSize) {
        uint64_t SLOTS_PER_BLOCK_t = 64;
        uint64_t xnslots = nslots + 10 * sqrt((double)nslots);
        uint64_t nblocks = (xnslots + SLOTS_PER_BLOCK_t - 1) / SLOTS_PER_BLOCK_t;
        uint64_t blocksize = 17;

        return ((nblocks) * (blocksize + 8 * (slotSize + fcounter + tagSize))) / 1024;

    }

    vector<uint64_t> estimateKmersHistogram(string fileName, int kSize, int threads);

    // Load the kmers in the input file into the output kDataframe. Input File can be of formats: fastq,fasta, sam, and bam.
    // void parseSequences(string seqFileName, int nThreads, kDataFrame *output);

    /// Load the kmers in the input file into the output kDataframe. Input File can be of formats: fastq,fasta.
    void parseSequences(kmerDecoder* KD, kDataFrame* output);

    /// Load the kmers in the input file into the output kDataframe. Input File can be of formats: fastq,fasta.
    void countKmersFromFile(kDataFrame* kframe, std::map<std::string, int> parse_params, string filename,
        int chunk_size = 1000);

    /// Load the kmers from KMC DB.
    void loadFromKMC(kDataFrame* kframe, std::string KMC_DB_filename);

    /// Load the kmers in the input string into the output kDataframe.
    void countKmersFromString(kmerDecoder* KD, string sequence, kDataFrame* output);

    /// Load the kmers in the input string into the output kDataframe.
    void countKmersFromString(kDataFrame* frame, std::map<std::string, int> parse_params, string sequence);


    /// Calculate the union of the kDataFrames. The result kDataframe will have all the kmers in the input list of kDataframes. The count of the kmers equals to the sum of the kmer count in the input list.
    kDataFrame* kFrameUnion(const vector<kDataFrame*>& input);

    /// Calculate the intersect of the kDataFrames. The result kDataframe will have only kmers that exists in all the kDataframes. The count of the kmers equals to the minimum of the kmer count in the input list.
    kDataFrame* kFrameIntersect(const vector<kDataFrame*>& input);

    /// Calculate the difference of the kDataframes. The result kDataframe will have only kmers that exists in the first kDataframe and not in any of the rest input kDataframes. The count of the kmers equals to the count in the first kDataframe.
    kDataFrame* kFrameDiff(const vector<kDataFrame*>& input);

    /// Initialize kmerDecoder to decode kmers from a FASTA/Q file with predefined mode.
    kmerDecoder* initialize_kmerDecoder(std::string filename, int chunkSize, std::string mode,
        std::map<std::string, int> parse_params);

    /// Initialize kmerDecoder to decode kmers from a sequence string with predefined mode.
    kmerDecoder* initialize_kmerDecoder(std::string mode, std::map<std::string, int> parse_params);

    // Initialize kmerDecoder to hash and Ihash single kmer.
    kmerDecoder* initialize_kmerDecoder(int kSize, hashingModes HM);

    /* set hashing mode for kmerDecoder object
     *
     * Mode 0: Murmar Hashing | Irreversible
     * Mode 1: Integer Hashing | Reversible | Full Hashing
     * Mode 2: TwoBitsHashing | Not considered hashing, just store the two bits representation
    */
    void kmerDecoder_setHashing(kmerDecoder* KD, hashingModes hash_mode);

    void kmerDecoder_setHashing(kDataFrame* KF, hashingModes hash_mode);

    /// Perform indexing to a sequences file with predefined kmers decoding mode, returns a colored kDataframe.
    void index(kmerDecoder* KD, string names_fileName, kDataFrame* frame);
    void indexMega(string fastaFileName, string tmpFolder, kDataFrame* frame);

    /// Index function without needing the kmerDecoder
    void
        index(kDataFrame* frame, std::map<std::string, int> parse_params, string filename, int chunks,
            string names_fileName);

    /*
     * Create an index for the kDataframes in the vector<input>. the index is in the form of kDataframe and it has default column of type mix Vector
     * mix vector split the colors into num_vectors vectors where each has maximum size vector_size.
     */
    void indexPriorityQueue(vector<kDataFrame*>& input, string tmpFolder, kDataFrame* output, uint32_t num_vectors = 20, uint32_t vector_size = 1000000);

    void mergeIndexes(vector<kDataFrame*>& input, string tmpFolder, kDataFrame* output);

    kDataFrame* innerJoin(vector<kDataFrame*> input, vector<uint32_t> kmersToKeep);


    // Extended functions
    uint64_t aggregate_count(kDataFrame* kf, const string& countColName);
    void transform_normalize(kDataFrame* kf, const string& countColName, uint64_t totalCount);
    kDataFrame* filter_zeroCounts(kDataFrame* res, uint32_t allDatasets);
    void transform_foldchange(kDataFrame* res, uint32_t nSamples, uint32_t nControl, uint32_t allDatasets, const string& foldChangeColName);
    void aggregate_foldChangeByGene(kDataFrame* res, unordered_map<string, vector<double>>* foldChangeByGene, const string& foldChangeColName, string& colorColumnName);
    

}
