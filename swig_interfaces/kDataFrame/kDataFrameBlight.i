class kDataFrameBlight : public kDataFrame {
public:
    kDataFrameBlight();

    kDataFrameBlight(uint64_t ksize)
    {
        blight_index = new kmer_Set_Light(ksize);
        kSize = ksize;
    }

    kDataFrameBlight(std::uint64_t ksize, string input_fasta_file);
    kDataFrame* getTwin();

    void _reserve(std::uint64_t n);
    void _reserve(vector<std::uint64_t> countHistogram);

    bool kmerExist(string kmer);
    bool kmerExist(uint64_t kmer);

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

    // kDataFrameIterator end();
    kDataFrameIterator find(const string& kmer);
    kDataFrameIterator find(uint64_t kmer);


    void serialize(string filePath);

    static kDataFrame* load(string filePath);

    template<typename T, typename Container>
    T getKmerColumnValue(const string& columnName, string kmer);

    template<typename T, typename Container>
    void setKmerColumnValue(const string& columnName, string kmer, T value);


};