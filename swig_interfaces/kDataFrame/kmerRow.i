namespace kProcessor {
    class kmerRow {
    public:
        uint64_t hashedKmer;
        uint64_t count;
        uint64_t order;
        kDataFrame* origin;
        kmerRow() {
            hashedKmer = 0;
            count = 0;
            order = 0;
            origin = nullptr;
        }
        kmerRow(std::uint64_t hashedKmer, std::uint64_t count, std::uint64_t order, kDataFrame* o)
        {
            this->hashedKmer = hashedKmer;
            this->count = count;
            this->order = order;
            this->origin = o;
        }
        kmerRow(const kmerRow& other)
        {
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

        bool operator==(kmerRow& other) const
        {
            return hashedKmer == other.hashedKmer;
        }

        bool operator < (kmerRow& other) const
        {
            return hashedKmer < other.hashedKmer;
        }
        bool operator > (kmerRow& other) const
        {
            return hashedKmer > other.hashedKmer;
        }
        bool operator < (const kmerRow& other) const
        {
            return hashedKmer < other.hashedKmer;
        }
        bool operator > (const kmerRow& other) const
        {
            return hashedKmer > other.hashedKmer;
        }

    };
}