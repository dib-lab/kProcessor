class kDataFrameMAP : public kDataFrame{

    public:
        kDataFrameMAP();
        kDataFrameMAP(std::uint64_t ksize);
//        kDataFrameMAP(std::uint64_t kSize, vector<std::uint64_t> kmersHistogram);
        kDataFrame* getTwin();
        void reserve (std::uint64_t n);
//        void reserve (vector<std::uint64_t> countHistogram);

        bool kmerExist(string kmer);

        bool setCount(string kmer, std::uint64_t count);
        bool setCount(std::uint64_t kmer, std::uint64_t count);
        bool insert(string kmer);
        bool insert(string kmer, std::uint64_t count);
        bool insert(std::uint64_t kmer, std::uint64_t count);
        bool insert(std::uint64_t kmer);
        std::uint64_t getCount(string kmer);
        std::uint64_t getCount(std::uint64_t kmerS);
        bool erase(string kmer);
        bool erase(std::uint64_t kmer);

        std::uint64_t size();
        std::uint64_t max_size();
        float load_factor();
        float max_load_factor();
        kDataFrameIterator begin();
         kDataFrameIterator end();
        kDataFrameIterator find(string kmer);
        kDataFrameIterator find(uint64_t kmer);
        std::uint64_t bucket(string kmer);
        void serialize(string filePath);
        static kDataFrame *load(string filePath);

        ~kDataFrameMAP() {
            this->MAP.clear();
        }
};