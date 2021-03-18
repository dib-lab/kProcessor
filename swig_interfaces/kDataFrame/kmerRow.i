class kmerRow{
    public:
        string kmer;
        uint64_t hashedKmer;
        uint64_t count;
        kDataFrame * origin;
        kmerRow(){
            kmer="";
            hashedKmer=0;
            count=0;
        }
        kmerRow(string kmer,std::uint64_t hashedKmer,std::uint64_t count,kDataFrame* o)
        {
            this->kmer=kmer;
            this->hashedKmer=hashedKmer;
            this->count=count;
            this->origin = o;
        }
        kmerRow(const kmerRow& other)
        {
            kmer=other.kmer;
            hashedKmer=other.hashedKmer;
            count=other.count;
            origin = other.origin ;
        }

        kmerRow copy(const kmerRow& other)
        {
            return * (new kmerRow(other));
        }

        template<typename T,typename Container>
        void getColumnValue(string colName, T& res);

        template<typename T,typename Container>
        void setColumnValue(string colName, T value);

        bool operator==(kmerRow &other)
        {
            return hashedKmer == other.hashedKmer;
        }

        bool operator < (kmerRow& other)
        {
            return hashedKmer<other.hashedKmer;
        }
        bool operator > ( kmerRow& other)
        {
            return hashedKmer>other.hashedKmer;
        }
        bool operator < (const kmerRow& other)
        {
            return hashedKmer<other.hashedKmer;
        }
        bool operator > (const kmerRow& other)
        {
            return hashedKmer>other.hashedKmer;
        }

};
