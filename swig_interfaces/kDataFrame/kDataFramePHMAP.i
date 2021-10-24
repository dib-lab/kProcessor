class kDataFramePHMAP : public kDataFrame {
        public:
        kDataFramePHMAP(uint64_t ksize);
        kDataFramePHMAP(readingModes RM, hashingModes hash_mode, map<string, int> params);
        kDataFrame* getTwin();
        void reserve(uint64_t n);

        bool setCount(string kmer, uint64_t count);
        bool insert(string kmer);
        bool insert(string kmer, uint64_t count);
        uint64_t getCount(string kmer);
        bool erase(string kmer);

        uint64_t size();
        uint64_t max_size();
        float load_factor();
        float max_load_factor();
        kDataFrameIterator begin();
        kDataFrameIterator end();

        uint64_t bucket(string kmer);
        void save(string filePath);
        static kDataFrame* load(string filePath);

        ~kDataFramePHMAP() { this->MAP.clear(); }
};