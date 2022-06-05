class kDataFrame {
public:
    uint64_t lastKmerOrder;
    unordered_map<string, Column*> columns;
    typedef kDataFrameIterator iterator;
    kmerDecoder* KD;
    virtual string get_class_name() { return class_name; }  // Temporary until resolving #17
    kDataFrame();
    explicit kDataFrame(uint8_t kSize);


    virtual kDataFrame* clone() = 0;
    virtual ~kDataFrame();
    /// creates a new kDataframe using the same parameters as the current kDataFrame.
    /*! It is like clone but without copying the data */
    virtual kDataFrame* getTwin() = 0;

    /// request a capacity change so that the kDataFrame can approximately hold at least n kmers
    void reserve(std::uint64_t n);
    /// request a capacity change so that the kDataFrame can approximately hold at least n kmers
    virtual void _reserve(std::uint64_t n) = 0;
    /// request a capacity change so that the kDataFrame can approximately hold kmers with countHistogram distribution
    virtual void _reserve(vector<std::uint64_t> countHistogram) = 0;
    /// insert the kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
    /*! Returns bool value indicating whether the kmer is inserted or not*/
    virtual uint32_t insert(const string& kmer) = 0;
    /// insert the hashed kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
    /*! Returns bool value indicating whether the kmer is inserted or not*/
    virtual uint32_t insert(std::uint64_t kmer) = 0;

    /// set the kmer's count to N time in the kDataFrame
   /*! Returns bool value indicating whether the kmer is inserted or not.
   The difference between setCount and insert is that setCount set the count to N no matter the previous kmer count was*/
    void addCountColumn();
    virtual bool setCount(const string& kmer, std::uint64_t N);
    virtual bool setCount(std::uint64_t kmer, std::uint64_t N);
    std::uint64_t getCount(const string& kmer);
    std::uint64_t getCount(std::uint64_t kmer);
    void incrementCount(std::uint64_t kmer);
    void incrementCount(const string kmer);
    /// returns the count of the kmer in the kDataFrame, i.e. the number of times the kmer is inserted in the kdataFrame.
    virtual std::uint64_t getkmerOrder(const string& kmer) = 0;
    virtual std::uint64_t getkmerOrder(std::uint64_t kmer) = 0;
    virtual bool setOrder(const string& kmer, std::uint64_t N) = 0;
    virtual  bool setOrder(std::uint64_t kmer, std::uint64_t N) = 0;
    // Removes  a kmer from the kDataFrame
    /*! Returns bool value indicating whether the kmer is erased or not*/
    virtual bool erase(const string& kmer) = 0;
    virtual bool erase(std::uint64_t kmer) = 0;

    /// Returns the number of kmers in the kDataframe.
    virtual std::uint64_t size() = 0;
    /// Returns the maximum number of kmers that the kDataframe can hold.
    virtual std::uint64_t max_size() = 0;
    /// Test whether the kDataFrame is empty.
    /*! Returns a bool value indicating whether the kDataFrame is empty, i.e. whether its size is 0.*/
    bool empty();
    /// Returns the current load factor in the kDataFrame.
    virtual float load_factor() = 0;
    /// Returns the current maximum load factor for the kDataFrame.
    virtual float max_load_factor() = 0;

    ///Returns an iterator at the beginning of the kDataFrame.
    virtual kDataFrameIterator begin() = 0;
    ///Returns an iterator at the end of the kDataFrame.
    virtual kDataFrameIterator end();
    ///Returns an iterator at the specific kmer.
    virtual kDataFrameIterator find(const string& kmer) = 0;
    virtual kDataFrameIterator find(uint64_t kmer) = 0;


    dbgIterator getDBGIterator(const string& kmer);
    virtual void serialize(string filePath) = 0;
    /// Returns the kmerDecoder used by kDataframe.
    kmerDecoder* getkmerDecoder() const {
        return KD;
    };



    static kDataFrame* load(string filePath);
    void save(string filePath);

    vector<string> getColumnNames();
    std::uint64_t getkSize() const { return kSize; }

    // duplicate for easier name in python, getkSize won't be wrapped
    std::uint64_t ksize() const { return kSize; }

    void setkSize(std::uint64_t k) {
        kSize = k;
        if (KD != nullptr)
        {
            kmerDecoder* newKD = new Kmers(kSize, KD->hash_mode);
            delete KD;
            KD = newKD;
        }        
else {
            KD = new Kmers(kSize);
        }
    }


    void addColumn(string columnName, Column*);
    void removeColumn(string columnName);
    template<typename Container>
    typename Container::dataType getKmerColumnValue(const string& columnName, string kmer);

    template<typename Container>
    void setKmerColumnValue(const string& columnName, string kmer, typename Container::dataType value);


    template<typename Container>
    typename Container::dataType getKmerColumnValue(const string& columnName, uint64_t kmer);

    template<typename Container>
    void setKmerColumnValue(const string& columnName, uint64_t kmer, typename Container::dataType value);


    template<typename Container>
    typename Container::dataType getKmerColumnValueByOrder(const string& columnName, uint64_t kmerOrder);

    template<typename Container>
    void setKmerColumnValueByOrder(const string& columnName, uint64_t kmerOrder, typename Container::dataType value);




    virtual bool kmerExist(string kmer) = 0;
    virtual bool kmerExist(uint64_t kmer) = 0;

    void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName, std::uint64_t kmer);
    void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName, string kmer);

};