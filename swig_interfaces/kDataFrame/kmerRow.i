class kmerRow {
public:
    string kmer;
    uint64_t hashedKmer;
    uint64_t count;
    uint64_t order;
    kDataFrame* origin;
    kmerRow() {
        kmer = "";
        hashedKmer = 0;
        count = 0;
        order = 0;
        origin = nullptr;
    }
    kmerRow(string kmer, std::uint64_t hashedKmer, std::uint64_t count, std::uint64_t order, kDataFrame* o)
    {
        this->kmer = kmer;
        this->hashedKmer = hashedKmer;
        this->count = count;
        this->order = order;
        this->origin = o;
    }
    kmerRow(const kmerRow& other)
    {
        kmer = other.kmer;
        hashedKmer = other.hashedKmer;
        count = other.count;
        origin = other.origin;
        order = 0;
    }

    static kmerRow copy(const kmerRow& other)
    {
        return *(new kmerRow(other));
    }

    template<typename T, typename Container>
    void getColumnValue(const string& colName, T& res);

    template<typename T, typename Container>
    void setColumnValue(const string& colName, T value);

};