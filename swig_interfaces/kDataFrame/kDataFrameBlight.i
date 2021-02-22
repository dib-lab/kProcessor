class kDataFrameBlight : public kDataFrame {
    public:
        kDataFrameBlight();

        kDataFrameBlight(uint64_t ksize)
        {
            blight_index=new kmer_Set_Light(ksize);
            kSize=ksize;
        }

        kDataFrameBlight(std::uint64_t ksize,string input_fasta_file);
        kDataFrame *getTwin();

        void reserve(std::uint64_t n);
        void reserve (vector<std::uint64_t> countHistogram);

        bool kmerExist(string kmer);

        bool setCount(string kmer, std::uint64_t count);
        bool setCount(std::uint64_t kmer, std::uint64_t count);

        bool insert(string kmer);

        bool insert(string kmer, std::uint64_t count);

        bool insert(std::uint64_t kmer, std::uint64_t count);

        bool insert(std::uint64_t kmer);

        std::uint64_t getCount(string kmer);
        std::uint64_t getCount(std::uint64_t kmer);

        bool erase(string kmer);
        bool erase(std::uint64_t kmer);

        std::uint64_t size();

        std::uint64_t max_size();

        float load_factor();

        float max_load_factor();

        kDataFrameIterator begin();

        // kDataFrameIterator end();
        kDataFrameIterator find(string kmer);
        kDataFrameIterator find(uint64_t kmer);


        void serialize(string filePath);

        static kDataFrame *load(string filePath);

        ~kDataFrameBlight() {

        }


        template<typename T,typename Container>
        T getKmerColumnValue(string columnName,string kmer);

        template<typename T,typename Container>
        void setKmerColumnValue(string columnName,string kmer, T value);



//        void changeDefaultColumnType(Column*);

        Column* getDefaultColumn(){
            return defaultColumn;
        }

        template<typename T,typename Container>
        T getKmerDefaultColumnValue(string kmer);

        template<typename T,typename Container>
        T getKmerDefaultColumnValue(std::uint64_t kmer);

        template<typename T,typename Container>
        void setKmerDefaultColumnValue(string kmer, T value);

        template<typename T,typename Container>
        void setKmerDefaultColumnValue(std::uint64_t kmer, T value);


};