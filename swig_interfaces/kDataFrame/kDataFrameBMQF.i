class kDataFrameBMQF : public kDataFrame {

public:
        kDataFrameBMQF();
        kDataFrameBMQF(std::uint64_t kSize, string path);
        kDataFrameBMQF(std::uint64_t kSize, uint64_t nKmers, string path);
        kDataFrameBMQF(std::uint64_t ksize, uint8_t q, uint8_t fixedCounterSize, uint8_t tagSize, double falsePositiveRate, string path);
        kDataFrameBMQF(bufferedMQF* bufferedmqf, std::uint64_t ksize, double falsePositiveRate);
        //count histogram is array where count of kmers repeated n times is found at index n. index 0 holds number of distinct kmers.
        kDataFrameBMQF(std::uint64_t ksize, vector<std::uint64_t> countHistogram, uint8_t tagSize
                , double falsePositiveRate);

        void _reserve(std::uint64_t n);
        void _reserve(vector<std::uint64_t> countHistogram);

        kDataFrame* getTwin();

        static std::uint64_t estimateMemory(std::uint64_t nslots, std::uint64_t slotSize,
                std::uint64_t fcounter, std::uint64_t tagSize);


        static void estimateParameters(vector<std::uint64_t> countHistogram,
                std::uint64_t numHashBits, std::uint64_t tagSize,
                std::uint64_t* res_noSlots, std::uint64_t* res_fixedSizeCounter, std::uint64_t* res_memory);



        bool insert(const string& kmer);
        bool insert(std::uint64_t kmer);
        bool setOrder(string& kmer, std::uint64_t count);
        bool setOrder(std::uint64_t kmer, std::uint64_t count);


        bool setCount(const string& kmer, std::uint64_t N) override;
        bool setCount(std::uint64_t kmer, std::uint64_t N) override;




        std::uint64_t getkmerOrder(const string& kmer);
        std::uint64_t getkmerOrder(std::uint64_t kmer);


        bool erase(const string& kmer);
        bool erase(std::uint64_t kmer);

        std::uint64_t size();
        /// max_size function returns the estimated maximum number of kmers that the kDataframeBMQF can hold.
        /*! The number of kmers is estimated as if all the kmers repeated 2^(fixed counter size)-1 times.*/
        std::uint64_t max_size();
        float load_factor();
        float max_load_factor();


        bufferedMQF* getBMQF() {
                return bufferedmqf;
        }

        void serialize(string filePath);
        static kDataFrame* load(string filePath);

        kDataFrameIterator begin();
        // kDataFrameIterator end();
        kDataFrameIterator find(const string& kmer);
        kDataFrameIterator find(uint64_t kmer);
        string getFilename();

        void deleteMemoryBuffer();
        bool kmerExist(string kmer);
        bool kmerExist(uint64_t kmer);
};