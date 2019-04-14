class kDataFrameIterator{
private:
  kDataFrame* origin;
  _kDataFrameIterator* iterator;
public:
  kDataFrameIterator(){
    iterator=NULL;
  }
  kDataFrameIterator(const kDataFrameIterator& other){
    if(other.iterator!=NULL){
      iterator=other.iterator->clone();
      origin=other.origin;
    }
    else{
      iterator=NULL;
      origin=NULL;
    }
  }
  kDataFrameIterator& operator= (const kDataFrameIterator& other){
    if(other.iterator!=NULL){
      iterator=other.iterator->clone();
      origin=other.origin;
    }
    else{
      iterator=NULL;
      origin=NULL;
    }
    return *this;
  }
  kDataFrameIterator(_kDataFrameIterator* it,kDataFrame* o){
    this->origin=o;
    iterator=it;
  }
/// Increment the iterator to the next kmer
  kDataFrameIterator& operator ++ (int){
    (*iterator)++;
    return *this;
  }

/// IMPLEMENTED ONLY IN SWIG PYTHON INTERFACE
  kDataFrameIterator& next(){
    (*iterator)++;
    return *this;
  }


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
  uint64_t getHashedKmer(){
    return iterator->getHashedKmer();
  };
  /// Returns the current kmer
  string getKmer(){
    return iterator->getKmer();
  }
  /// Returns the count of the current kmer
  uint64_t getKmerCount(){
    return iterator->getKmerCount();
  }
  /// sets the count of the current kmer
  bool setKmerCount(uint64_t count){
    return iterator->setKmerCount(count);
  }
  kmerRow operator*(){
    return kmerRow(iterator->getKmer(),
                   iterator->getHashedKmer(),
                   iterator->getKmerCount()
                  );
  }
  ~kDataFrameIterator(){
    delete iterator;
  }
};