class kDataFrame{
public:
//    unordered_map<string, Column*> columns;
//    typedef kDataFrameIterator iterator;
//    typedef kmerRow value_type;
//    kmerDecoder * KD;
    virtual string get_class_name(){ return class_name;}  // Temporary until resolving #17
    kDataFrame();
    kDataFrame(uint8_t kSize);


    virtual ~kDataFrame();
/// creates a new kDataframe using the same parameters as the current kDataFrame.
/*! It is like clone but without copying the data */
    virtual kDataFrame* getTwin()=0;
/// request a capacity change so that the kDataFrame can approximately hold at least n kmers
    virtual void reserve (std::uint64_t n )=0;
/// request a capacity change so that the kDataFrame can approximately hold kmers with countHistogram distribution
//    virtual void reserve (vector<std::uint64_t> countHistogram)=0;
/// insert the kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
    virtual bool insert(string kmer)=0;
/// insert the kmer N time in the kDataFrame, or increment the kmer count with N if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
    virtual bool insert(string kmer,std::uint64_t N)=0;
/// insert the hashed kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
    virtual bool insert(std::uint64_t kmer)=0;
/// insert the hashed kmer N time in the kDataFrame, or increment the kmer count with N if it is already exists.
/*! Returns bool value indicating whether the hashed kmer is inserted or not*/
    virtual bool insert(std::uint64_t kmer,std::uint64_t N)=0;
    /// insert the kmer in the kmer row time in the kDataFrame, or increment the kmer count with the count in the row if it is already exists.
    /*! Returns bool value indicating whether the kmer is inserted or not*/
    bool insert(kmerRow k);
    kDataFrame::iterator insert(kDataFrame::iterator& it,kmerRow k);
/// set the kmer's count to N time in the kDataFrame
/*! Returns bool value indicating whether the kmer is inserted or not.
The difference between setCount and insert is that setCount set the count to N no matter the previous kmer count was*/
    virtual bool setCount(string kmer,std::uint64_t N)=0;
    virtual bool setCount(std::uint64_t kmer,std::uint64_t N)=0;
/// returns the count of the kmer in the kDataFrame, i.e. the number of times the kmer is inserted in the kdataFrame.
    virtual std::uint64_t getCount(string kmer)=0;
    virtual std::uint64_t getCount(std::uint64_t kmer)=0;
// Removes  a kmer from the kDataFrame
/*! Returns bool value indicating whether the kmer is erased or not*/
    virtual bool erase(string kmer)=0;
    virtual bool erase(std::uint64_t kmer)=0;

/// Returns the number of kmers in the kDataframe.
    virtual std::uint64_t size()=0;
/// Returns the maximum number of kmers that the kDataframe can hold.
    virtual std::uint64_t max_size()=0;
/// Test whether the kDataFrame is empty.
/*! Returns a bool value indicating whether the kDataFrame is empty, i.e. whether its size is 0.*/
    bool empty();
/// Returns the current load factor in the kDataFrame.
    virtual float load_factor()=0;
/// Returns the current maximum load factor for the kDataFrame.
    virtual float max_load_factor()=0;

///Returns an iterator at the beginning of the kDataFrame.
    virtual kDataFrameIterator begin()=0;
///Returns an iterator at the end of the kDataFrame.
    virtual kDataFrameIterator end();
///Returns an iterator at the specific kmer.
    virtual kDataFrameIterator find(string kmer)=0;
    virtual kDataFrameIterator find(uint64_t kmer)=0;


    dbgIterator getDBGIterator(string kmer);
    virtual void serialize(string filePath)=0;
/// Returns the kmerDecoder used by kDataframe.
    kmerDecoder* getkmerDecoder(){
        return KD;
    };



    static kDataFrame* load(string filePath);
    void save(string filePath);


    std::uint64_t getkSize(){return kSize;}

    // duplicate for easier name in python, getkSize won't be wrapped
    std::uint64_t ksize(){return kSize;}

    void setkSize(std::uint64_t k){

        kSize=k;
        if(KD!= nullptr)
            delete KD;
        KD= new Kmers(kSize);
    }


    void addColumn(string columnName, Column*);

    template<typename T,typename Container>
    T getKmerColumnValue(string columnName,string kmer);

    template<typename T,typename Container>
    void setKmerColumnValue(string columnName,string kmer, T value);


    template<typename T,typename Container>
    T getKmerColumnValue(string columnName,uint64_t kmer);

    template<typename T,typename Container>
    void setKmerColumnValue(string columnName,uint64_t kmer, T value);



    void changeDefaultColumnType(Column*);
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

    virtual bool kmerExist(string kmer)=0;

    void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName, std::uint64_t kmer);
    void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName, string kmer);

};

// TODO
//%extend kDataFrame {
//    %template(addColumn_int) addColumn<int>;
//    %template(addColumn_bool) addColumn<bool>;
//    %template(addColumn_double) addColumn<double>;
//
//    %template(getKmerColumnValue_int) getKmerColumnValue<int>;
//    %template(getKmerColumnValue_bool) getKmerColumnValue<bool>;
//    %template(getKmerColumnValue_double) getKmerColumnValue<double>;
//
//    %template(setKmerColumnValue_int) setKmerColumnValue<int>;
//    %template(setKmerColumnValue_bool) setKmerColumnValue<bool>;
//    %template(setKmerColumnValue_double) setKmerColumnValue<double>;
//};