#ifndef _kDataFRAME_H_
#define _kDataFRAME_H_

#include <vector>

#include <iostream>
#include "kmerDecoder.hpp"
#include <any>
#include "defaultColumn.hpp"
#include <cstdint>

using namespace std;

class kDataFrame;


class _kDataFrameIterator{
protected:
    std::uint64_t kSize;
public:
  _kDataFrameIterator()= default;
  explicit _kDataFrameIterator(std::uint64_t k):kSize(k){}
  virtual uint64_t getOrder()=0;
  virtual _kDataFrameIterator& operator ++ (int)=0;
  virtual _kDataFrameIterator* clone()=0;
  virtual std::uint64_t getHashedKmer()=0;
  virtual string getKmer()=0;
  virtual std::uint64_t getCount()=0;
  virtual bool setOrder(std::uint64_t count)=0;
  virtual bool operator ==(const _kDataFrameIterator& other)=0;
  virtual bool operator !=(const _kDataFrameIterator& other)=0;
  virtual ~_kDataFrameIterator()= default;

};

class kDataFrameIterator{
private:
  kDataFrame* origin;
public:
    _kDataFrameIterator* iterator;
    using iterator_category = std::forward_iterator_tag;
  kDataFrameIterator(){
    iterator=nullptr;
  }
  kDataFrameIterator(const kDataFrameIterator& other){
    if(other.iterator!=nullptr){
      iterator=other.iterator->clone();
      origin=other.origin;
    }
    else{
      iterator=nullptr;
      origin=nullptr;
    }
  }
  kDataFrameIterator& operator= (const kDataFrameIterator& other){
    if(other.iterator!=nullptr){
      iterator=other.iterator->clone();
      origin=other.origin;
    }
    else{
      iterator=nullptr;
      origin=nullptr;
    }
    return *this;
  }
  kDataFrameIterator(_kDataFrameIterator* it,kDataFrame* o){
    this->origin=o;
    iterator=it;
  }
/// Increment the iterator to the next kmer
    kDataFrameIterator& operator * (){
        return *this;
    }

/// Increment the iterator to the next kmer
  kDataFrameIterator& operator ++ (int){
    (*iterator)++;
    return *this;
  }
    kDataFrameIterator& operator ++ (){
        (*iterator)++;
        return *this;
    }

  /// Increment the iterator to the next kmer (Implemented mainly for python interface)
  kDataFrameIterator& next(){
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
  std::uint64_t getHashedKmer(){
    return iterator->getHashedKmer();
  };
  /// Returns the current kmer order
    uint64_t getOrder(){
        return iterator->getOrder();
    }
  /// Returns the current kmer
  string getKmer(){
    return iterator->getKmer();
  }
  /// Returns the count of the current kmer
  std::uint64_t getCount();
  /// sets the count of the current kmer
  bool setCount(std::uint64_t count);
  
  bool setOrder(std::uint64_t count){
    return iterator->setOrder(count);
  }

  template<typename Container>
  void getColumnValue(const string& colName, typename Container::dataType & res);
  template<typename Container>
  void setColumnValue(const string& colName, typename Container::dataType value);


  void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName);

  ~kDataFrameIterator(){
    delete iterator;
  }
};



class dbgIterator{
public:
    kDataFrame* frame;
    string currentKmer;
    vector<string> nextFwdKmers;
    vector<string> nextRevKmers;

    dbgIterator();
    dbgIterator(kDataFrame*,string kmer);

    dbgIterator(const dbgIterator& other);
    dbgIterator& operator= (const dbgIterator& other);

    void nextFWD(uint32_t index);
    void nextREV(uint32_t index);
private:
    void generateNextKmers();

};

class kDataFrame{
protected:
  std::uint64_t kSize;
  string class_name; // Default = MQF, change if MAP. Temporary until resolving #17



  unordered_map<uint64_t,uint32_t> orderCheckpoints;
  uint32_t lastCheckpoint;
  vectorColumn<uint32_t>* countColumn;
  kDataFrameIterator* endIterator;
  friend class kDataFrameIterator;
public:
    uint64_t lastKmerOrder;
    unordered_map<string, Column*> columns;
    typedef kDataFrameIterator iterator;
    kmerDecoder * KD;
  virtual string get_class_name(){ return class_name;}  // Temporary until resolving #17
  kDataFrame();
  explicit kDataFrame(uint8_t kSize);


  virtual kDataFrame* clone()=0;
  virtual ~kDataFrame();
/// creates a new kDataframe using the same parameters as the current kDataFrame.
/*! It is like clone but without copying the data */
  virtual kDataFrame* getTwin()=0;

  /// request a capacity change so that the kDataFrame can approximately hold at least n kmers
  void reserve (std::uint64_t n );
/// request a capacity change so that the kDataFrame can approximately hold at least n kmers
  virtual void _reserve (std::uint64_t n )=0;
/// request a capacity change so that the kDataFrame can approximately hold kmers with countHistogram distribution
  virtual void _reserve (vector<std::uint64_t> countHistogram)=0;
/// insert the kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
  virtual uint32_t insert(const string &kmer)=0;
/// insert the hashed kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
  virtual uint32_t insert(std::uint64_t kmer)=0;

 /// set the kmer's count to N time in the kDataFrame
/*! Returns bool value indicating whether the kmer is inserted or not.
The difference between setCount and insert is that setCount set the count to N no matter the previous kmer count was*/
  void addCountColumn();
  virtual bool setCount(const string &kmer, std::uint64_t N);
  virtual bool setCount(std::uint64_t kmer,std::uint64_t N);
  std::uint64_t getCount(const string &kmer);
  std::uint64_t getCount(std::uint64_t kmer);
  void incrementCount(std::uint64_t kmer);
  void incrementCount(const string kmer);
/// returns the count of the kmer in the kDataFrame, i.e. the number of times the kmer is inserted in the kdataFrame.
  virtual std::uint64_t getkmerOrder(const string &kmer)=0;
  virtual std::uint64_t getkmerOrder(std::uint64_t kmer)=0;
  virtual bool setOrder(const string &kmer, std::uint64_t N)=0;
  virtual  bool setOrder(std::uint64_t kmer,std::uint64_t N)=0;
// Removes  a kmer from the kDataFrame
/*! Returns bool value indicating whether the kmer is erased or not*/
  virtual bool erase(const string &kmer)=0;
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
  virtual kDataFrameIterator find(const string &kmer)=0;
  virtual kDataFrameIterator find(uint64_t kmer)=0;


  dbgIterator getDBGIterator(const string &kmer);
  virtual void serialize(string filePath)=0;
/// Returns the kmerDecoder used by kDataframe.
  kmerDecoder* getkmerDecoder() const{
    return KD;
  };



  static kDataFrame* load(string filePath);
  void save(string filePath);

  vector<string> getColumnNames();
  std::uint64_t getkSize() const{return kSize;}

  // duplicate for easier name in python, getkSize won't be wrapped
  std::uint64_t ksize() const{return kSize;}

  void setkSize(std::uint64_t k){
    kSize=k;
    if(KD!= nullptr)
    {
        kmerDecoder* newKD=new Kmers(kSize,KD->hash_mode);
        delete KD;
        KD=newKD;
    }else{
        KD= new Kmers(kSize);
    }
  }


  void addColumn(string columnName, Column*);
  void removeColumn(string columnName);
  template<typename Container>
  typename Container::dataType getKmerColumnValue(const string& columnName,string kmer);

  template<typename Container>
  void setKmerColumnValue(const string& columnName,string kmer, typename Container::dataType value);


  template<typename Container>
  typename Container::dataType getKmerColumnValue(const string& columnName,uint64_t kmer);

  template<typename Container>
  void setKmerColumnValue(const string& columnName,uint64_t kmer, typename Container::dataType value);


  template<typename Container>
  typename Container::dataType getKmerColumnValueByOrder(const string& columnName,uint64_t kmerOrder);

  template<typename Container>
  void setKmerColumnValueByOrder(const string& columnName,uint64_t kmerOrder, typename Container::dataType value);




  virtual bool kmerExist(string kmer)=0;
  virtual bool kmerExist(uint64_t kmer)=0;

  void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName, std::uint64_t kmer);
  void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName, string kmer);

};

template<typename Container>
typename Container::dataType kDataFrame::getKmerColumnValue(const string& columnName,string kmer)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    if(kmerOrder==0)
        throw std::logic_error("kmer not found!");
    return ((Container*)columns[columnName])->get(kmerOrder);
}
template<typename Container>
void kDataFrame::setKmerColumnValue(const string& columnName,string kmer,typename Container::dataType value)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    if(kmerOrder==0)
        throw std::logic_error("kmer not found!");
    ((Container*)columns[columnName])->insert(value,kmerOrder);
}


template<typename Container>
typename Container::dataType kDataFrame::getKmerColumnValue(const string& columnName,uint64_t kmer)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    if(kmerOrder==0)
        throw std::logic_error("kmer not found!");
    return ((Container*)columns[columnName])->get(kmerOrder);
}
template<typename Container>
void kDataFrame::setKmerColumnValue(const string& columnName,uint64_t kmer,typename Container::dataType value)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    if(kmerOrder==0)
        throw std::logic_error("kmer not found!");
    ((Container*)columns[columnName])->insert(value,kmerOrder);
}

template<typename Container>
typename Container::dataType kDataFrame::getKmerColumnValueByOrder(const string& columnName,uint64_t kmerOrder)
{
    return ((Container*)columns[columnName])->get(kmerOrder);
}
template<typename Container>
void kDataFrame::setKmerColumnValueByOrder(const string& columnName,uint64_t kmerOrder,typename Container::dataType value)
{
    ((Container*)columns[columnName])->insert(value,kmerOrder);
}





template<typename Container>
void kDataFrameIterator::getColumnValue(const string& colName, typename Container::dataType& res)
{
    res= origin->getKmerColumnValueByOrder<Container>(colName, iterator->getOrder());
}

template<typename Container>
void kDataFrameIterator::setColumnValue(const string& colName, typename Container::dataType value)
{
    origin->setKmerColumnValueByOrder<Container>(colName, iterator->getOrder(),value);
}



class kDataFrameFactory{
public:
    static kDataFrame* loadMQF(string filePath);
    static kDataFrame* createMQF(uint32_t kSize,uint32_t numKmers=10000);
    static kDataFrame* createMQF(kDataFrame* kframe);

    static kDataFrame* loadPHMAP(string filePath);
    static kDataFrame* createPHMAP(uint32_t kSize,uint32_t numKmers=10000);
    static kDataFrame* createPHMAP(uint64_t ksize, hashingModes hash_mode);


    static kDataFrame* loadMAP(string filePath);
    static kDataFrame* createMAP(uint32_t kSize,uint32_t numKmers=10000);

    static kDataFrame* loadBtree(string filePath);
    static kDataFrame* createBtree(uint32_t kSize,uint32_t numKmers=10000);

    static kDataFrame* loadBMQF(string filePath);
    static kDataFrame* createBMQF(uint32_t kSize,string filePath,uint32_t numKmers=10000);
    static kDataFrame* createBMQF(kDataFrame* kframe,string filePath);


    static kDataFrame* loadBlight(string filePath);
    static kDataFrame* createBlight(uint32_t kSize,string filePath);

    static kDataFrame* loadSSHASH(string filePath);
    static kDataFrame* createSSHASH(uint32_t kSize,string filePath);

};


class kDataFrameUtility {
public:
    static void deleteMemoryBufferBMQF(kDataFrame* frame);
};


#endif
