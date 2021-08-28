
class kDataFramePHMAP : public kDataFrame {

public:
    kDataFramePHMAP();

    explicit kDataFramePHMAP(uint64_t ksize);
    kDataFramePHMAP(std::uint64_t kSize, uint64_t nKmers);
    kDataFramePHMAP(readingModes RM, hashingModes hash_mode, map<string, int> params);
    kDataFramePHMAP(uint64_t ksize, hashingModes hash_mode);
    kDataFramePHMAP(uint64_t kSize, vector<uint64_t> kmersHistogram);


    kDataFrame* getTwin();

    void _reserve(std::uint64_t n);
    void _reserve(vector<std::uint64_t> countHistogram);

    bool kmerExist(string kmer);
    bool kmerExist(uint64_t kmer);

    bool setOrder(string& kmer, std::uint64_t count);
    bool setOrder(std::uint64_t kmer, std::uint64_t count);

    bool insert(const string& kmer);
    bool insert(std::uint64_t kmer);

    std::uint64_t getkmerOrder(const string& kmer);
    std::uint64_t getkmerOrder(std::uint64_t kmer);

    bool erase(const string& kmer);
    bool erase(std::uint64_t kmer);

    std::uint64_t size();

    std::uint64_t max_size();

    float load_factor();

    float max_load_factor();

    kDataFrameIterator begin();

    kDataFrameIterator end();
    kDataFrameIterator find(const string& kmer);
    kDataFrameIterator find(uint64_t kmer);

    std::uint64_t bucket(string kmer);

    void serialize(string filePath);

    static kDataFrame* load(string filePath);

};