class kDataFrameIterator {
private:
  kDataFrame* origin;
public:
  _kDataFrameIterator* iterator;
  using iterator_category = std::forward_iterator_tag;
  kDataFrameIterator() {
    iterator = nullptr;
  }
  kDataFrameIterator(const kDataFrameIterator& other) {
    if (other.iterator != nullptr) {
      iterator = other.iterator->clone();
      origin = other.origin;
    }
    else {
      iterator = nullptr;
      origin = nullptr;
    }
  }
  kDataFrameIterator& operator= (const kDataFrameIterator& other) {
    if (other.iterator != nullptr) {
      iterator = other.iterator->clone();
      origin = other.origin;
    }
    else {
      iterator = nullptr;
      origin = nullptr;
    }
    return *this;
  }
  kDataFrameIterator(_kDataFrameIterator* it, kDataFrame* o) {
    this->origin = o;
    iterator = it;
  }
  /// Increment the iterator to the next kmer
  kDataFrameIterator& operator * () {
    return *this;
  }

  /// Increment the iterator to the next kmer
  kDataFrameIterator& operator ++ (int) {
    (*iterator)++;
    return *this;
  }
  kDataFrameIterator& operator ++ () {
    (*iterator)++;
    return *this;
  }

  /// Increment the iterator to the next kmer (Implemented mainly for python interface)
  kDataFrameIterator& next() {
    (*iterator)++;
    return *this;
  }

  // /// Increment the iterator to the next kmer
  //   kDataFrameIterator operator ++ (int){
  //     kDataFrameIterator temp=*this;
  //     (*iterator)++;
  //     return temp;
  //   }
  // /// Compare the position of each iterator in the underlying datastructure.
  // /*! returns True when other points to kmer points to a further position than the current */
  //   bool operator <(const kDataFrameIterator& other){
  //     return *iterator < *other.iterator;
  //   }
  //   /// Compare the position of each iterator in the underlying datastructure.
  //   /*! returns True when other points to kmer points to a nerarer position than the current */
  //   bool operator >(const kDataFrameIterator& other){
  //     return *iterator > *other.iterator;
  //   }
  //   /// Compare the position of each iterator in the underlying datastructure.
  //   /*! returns True when other points to kmer points to a further or equal position than the current */
  //   bool operator <=(const kDataFrameIterator& other){
  //     return *iterator <= *other.iterator;
  //   }
  //   /// Compare the position of each iterator in the underlying datastructure.
  //   /*! returns True when other points to kmer points to a nearer or equal position than the current */
  //   bool operator >=(const kDataFrameIterator& other){
  //     return *iterator >= *other.iterator;
  //   }
    /// Compare the position of each iterator in the underlying datastructure.
    /*! returns True when current and other points to the same kmer */
  bool operator ==(const kDataFrameIterator& other)
  {
    return *iterator == *other.iterator;
  }
  /// Compare the position of each iterator in the underlying datastructure.
  /*! returns True when current and other points to different kmers */
  bool operator !=(const kDataFrameIterator& other)
  {
    return *iterator != *other.iterator;
  }
  /// Returns the hash value of the current kmer
  std::uint64_t getHashedKmer() {
    return iterator->getHashedKmer();
  };
  /// Returns the current kmer order
  uint64_t getOrder() {
    return iterator->getOrder();
  }
  /// Returns the current kmer
  string getKmer() {
    return iterator->getKmer();
  }
  /// Returns the count of the current kmer
  std::uint64_t getCount();
  /// sets the count of the current kmer
  bool setCount(std::uint64_t count);

  bool setOrder(std::uint64_t count) {
    return iterator->setOrder(count);
  }

  template<typename Container>
  void getColumnValue(const string& colName, typename Container::dataType& res);
  template<typename Container>
  void setColumnValue(const string& colName, typename Container::dataType value);


  void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName);

  ~kDataFrameIterator() {
    delete iterator;
  }
};
