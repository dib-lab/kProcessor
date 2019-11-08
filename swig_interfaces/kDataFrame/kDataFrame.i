class kDataFrame{


public:
  kDataFrame();
  kDataFrame(uint8_t kSize);
  virtual ~kDataFrame(){
  }

/// creates a new kDataframe using the same parameters as the current kDataFrame.
/*! It is like clone but without copying the data */
  virtual kDataFrame* getTwin()=0;

/// request a capacity change so that the kDataFrame can approximately at least n kmers
  virtual void reserve (uint64_t n )=0;
  
/// insert the kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/  
  virtual bool insert(string kmer)=0;
  
/// insert the kmer N time in the kDataFrame, or increment the kmer count with N if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
  virtual bool insert(string kmer,uint64_t N)=0;  
  
  /// insert the kmer in the kmer row time in the kDataFrame, or increment the kmer count with the count in the row if it is already exists.
  /*! Returns bool value indicating whether the kmer is inserted or not*/
  bool insert(kmerRow k);
  
/// set the kmer's count to N time in the kDataFrame
/*! Returns bool value indicating whether the kmer is inserted or not.
The difference between setCount and insert is that setCount set the count to N no matter the previous kmer count was*/
  virtual bool setCount(string kmer,uint64_t N)=0;
  
/// returns the count of the kmer in the kDataFrame, i.e. the number of times the kmer is inserted in the kdataFrame.
  virtual uint64_t getCount(string kmer)=0;
  
// Removes  a kmer from the kDataFrame
/*! Returns bool value indicating whether the kmer is erased or not*/
  virtual bool erase(string kmer)=0;
  
/// Returns the number of kmers in the kDataframe.
  virtual uint64_t size()=0;
  
/// Returns the maximum number of kmers that the kDataframe can hold.
  virtual uint64_t max_size()=0;
  
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
  virtual kDataFrameIterator end()=0;
  
/// Export the kDataFrame in a file.
  virtual void save(string filePath)=0;
  
/// Returns the  hash function used by kDataframe.
    kmerDecoder* getkmerDecoder(){
        return KD;
    };
  
/// Load the kDataFrame from a file.
  static kDataFrame* load(string filePath);
  
/// Get kmer size used in the kDataFrame.  
    uint64_t ksize(){return kSize;}
  
/// Set the kmer size that will be used in the kDataFrame.  
  void setkSize(uint64_t k){kSize=k;}


  template<typename T>
  void addColumn(string columnName);

  template<typename T>
  T getKmerColumnValue(string columnName,string kmer);

  template<typename T>
  void setKmerColumnValue(string columnName,string kmer, T value);

};


%extend kDataFrame {
    %template(addColumn_int) addColumn<int>;
    %template(addColumn_bool) addColumn<bool>;
    %template(addColumn_double) addColumn<double>;

    %template(getKmerColumnValue_int) getKmerColumnValue<int>;
    %template(getKmerColumnValue_bool) getKmerColumnValue<bool>;
    %template(getKmerColumnValue_double) getKmerColumnValue<double>;

    %template(setKmerColumnValue_int) setKmerColumnValue<int>;
    %template(setKmerColumnValue_bool) setKmerColumnValue<bool>;
    %template(setKmerColumnValue_double) setKmerColumnValue<double>;
};