#ifndef _kDataFRAME_H_
#define _kDataFRAME_H_

#include "blight.h"
#include <vector>

#include "gqf.h"
#include <iostream>
#include <parallel_hashmap/phmap.h>
#include "kmerDecoder.hpp"
#include <any>
#include "bufferedMQF.h"

#include "defaultColumn.hpp"
#include <cstdint>

using phmap::flat_hash_map;
using namespace std;

class kDataFrame;
class kDataFrameMAP;
class kDataFramePHMAP;


class kmerRow{
public:
  string kmer;
  uint64_t hashedKmer;
  uint64_t count;
  uint64_t order;
  kDataFrame * origin;
  kmerRow(){
    kmer="";
    hashedKmer=0;
    count=0;
    order=0;
    origin=nullptr;
  }
  kmerRow(string kmer,std::uint64_t hashedKmer,std::uint64_t count,std::uint64_t order,kDataFrame* o)
  {
    this->kmer=kmer;
    this->hashedKmer=hashedKmer;
    this->count=count;
    this->order=order;
    this->origin = o;
  }
  kmerRow(const kmerRow& other)
  {
    kmer=other.kmer;
    hashedKmer=other.hashedKmer;
    count=other.count;
    origin = other.origin ;
    order=0;
  }

  static kmerRow copy(const kmerRow& other)
  {
    return * (new kmerRow(other));
  }

  template<typename T,typename Container>
  void getColumnValue(const string& colName, T& res);

  template<typename T,typename Container>
  void setColumnValue(const string& colName, T value);

  bool operator==(kmerRow &other) const
  {
    return hashedKmer == other.hashedKmer;
  }

  bool operator < (kmerRow& other) const
  {
    return hashedKmer<other.hashedKmer;
  }
  bool operator > ( kmerRow& other) const
  {
    return hashedKmer>other.hashedKmer;
  }
  bool operator < (const kmerRow& other) const
  {
    return hashedKmer<other.hashedKmer;
  }
  bool operator > (const kmerRow& other) const
  {
    return hashedKmer>other.hashedKmer;
  }

};

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
  kmerRow getKmerRow(){
    return kmerRow(iterator->getKmer(),
                   iterator->getHashedKmer(),
                   getCount(),
                   iterator->getOrder(),
                   origin
                  );
  }
  template<typename T,typename Container>
  void getColumnValue(const string& colName, T& res);
  template<typename T,typename Container>
  void setColumnValue(const string& colName, T value);


  void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName);

  ~kDataFrameIterator(){
    delete iterator;
  }
};


class kDataFrameMQFIterator:public _kDataFrameIterator{
private:
  QFi* qfi;
  kmerDecoder * KD;
public:
  kDataFrameMQFIterator(QF*,std::uint64_t kSize,kmerDecoder* KD);
  kDataFrameMQFIterator(QFi*,std::uint64_t kSize,kmerDecoder* KD);
  kDataFrameMQFIterator(const kDataFrameMQFIterator&);
  kDataFrameMQFIterator& operator ++ (int);
  _kDataFrameIterator* clone();
  std::uint64_t getHashedKmer();
  string getKmer();
  uint64_t getOrder();
  std::uint64_t getCount();
  bool setOrder(std::uint64_t count);
  void endIterator();
  bool operator ==(const _kDataFrameIterator& other);
  bool operator !=(const _kDataFrameIterator& other);
   ~kDataFrameMQFIterator();
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

  uint64_t lastKmerOrder;

  unordered_map<uint64_t,uint32_t> orderCheckpoints;
  uint32_t lastCheckpoint;
  vectorColumn<uint32_t>* countColumn;
  kDataFrameIterator* endIterator;
  friend class kDataFrameIterator;
public:
    unordered_map<string, Column*> columns;
    typedef kDataFrameIterator iterator;
    typedef kmerRow value_type;
    kmerDecoder * KD;
  virtual string get_class_name(){ return class_name;}  // Temporary until resolving #17
  kDataFrame();
  explicit kDataFrame(uint8_t kSize);


  virtual ~kDataFrame();
/// creates a new kDataframe using the same parameters as the current kDataFrame.
/*! It is like clone but without copying the data */
  virtual kDataFrame* getTwin()=0;
/// request a capacity change so that the kDataFrame can approximately hold at least n kmers
  virtual void reserve (std::uint64_t n )=0;
/// request a capacity change so that the kDataFrame can approximately hold kmers with countHistogram distribution
  virtual void reserve (vector<std::uint64_t> countHistogram)=0;
/// insert the kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
  virtual bool insert(const string &kmer)=0;
/// insert the hashed kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
  virtual bool insert(std::uint64_t kmer)=0;
  /// insert the kmer in the kmer row time in the kDataFrame, or increment the kmer count with the count in the row if it is already exists.
  /*! Returns bool value indicating whether the kmer is inserted or not*/
  bool insert(kmerRow k);
  kDataFrame::iterator insert(kDataFrame::iterator& it,kmerRow k);
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
  bool setOrder(const string &kmer, std::uint64_t N);
  bool setOrder(std::uint64_t kmer,std::uint64_t N);
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


  std::uint64_t getkSize() const{return kSize;}

  // duplicate for easier name in python, getkSize won't be wrapped
  std::uint64_t ksize() const{return kSize;}

  void setkSize(std::uint64_t k){

    kSize=k;
    if(KD!= nullptr)
        delete KD;
    KD= new Kmers(kSize);
  }


  void addColumn(string columnName, Column*);
  void removeColumn(string columnName);
  template<typename T,typename Container>
  T getKmerColumnValue(const string& columnName,string kmer);

  template<typename T,typename Container>
  void setKmerColumnValue(const string& columnName,string kmer, T value);


  template<typename T,typename Container>
  T getKmerColumnValue(const string& columnName,uint64_t kmer);

  template<typename T,typename Container>
  void setKmerColumnValue(const string& columnName,uint64_t kmer, T value);


  template<typename T,typename Container>
  T getKmerColumnValueByOrder(const string& columnName,uint64_t kmerOrder);

  template<typename T,typename Container>
  void setKmerColumnValueByOrder(const string& columnName,uint64_t kmerOrder, T value);




  virtual bool kmerExist(string kmer)=0;
  virtual bool kmerExist(uint64_t kmer)=0;

  void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName, std::uint64_t kmer);
  void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName, string kmer);

};

template<typename T,typename Container>
T kDataFrame::getKmerColumnValue(const string& columnName,string kmer)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    return ((Container*)columns[columnName])->get(kmerOrder);
}
template<typename T,typename Container>
void kDataFrame::setKmerColumnValue(const string& columnName,string kmer,T value)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    ((Container*)columns[columnName])->insert(value,kmerOrder);
}


template<typename T,typename Container>
T kDataFrame::getKmerColumnValue(const string& columnName,uint64_t kmer)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    return ((Container*)columns[columnName])->get(kmerOrder);
}
template<typename T,typename Container>
void kDataFrame::setKmerColumnValue(const string& columnName,uint64_t kmer,T value)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    ((Container*)columns[columnName])->insert(value,kmerOrder);
}

template<typename T,typename Container>
T kDataFrame::getKmerColumnValueByOrder(const string& columnName,uint64_t kmerOrder)
{
    return ((Container*)columns[columnName])->get(kmerOrder);
}
template<typename T,typename Container>
void kDataFrame::setKmerColumnValueByOrder(const string& columnName,uint64_t kmerOrder,T value)
{
    ((Container*)columns[columnName])->insert(value,kmerOrder);
}



template<typename T,typename Container>
void kmerRow::getColumnValue(const string& colName, T& res)
{
    res= origin->getKmerColumnValueByOrder<T,Container>(colName, order);
}

template<typename T,typename Container>
void kmerRow::setColumnValue(const string& colName, T value)
{
    origin->setKmerColumnValueByOrder<T,Container>(colName, order,value);
}


template<typename T,typename Container>
void kDataFrameIterator::getColumnValue(const string& colName, T& res)
{
    res= origin->getKmerColumnValueByOrder<T,Container>(colName, iterator->getOrder());
}

template<typename T,typename Container>
void kDataFrameIterator::setColumnValue(const string& colName, T value)
{
    origin->setKmerColumnValueByOrder<T,Container>(colName, iterator->getOrder(),value);
}


class kDataFrameMQF: public kDataFrame{

private:
  QF* mqf;
  double falsePositiveRate;
  std::uint64_t hashbits;
  __uint128_t range;
  static bool isEnough(vector<std::uint64_t> histogram,std::uint64_t noSlots,std::uint64_t fixedSizeCounter,std::uint64_t slotSize);
  friend class kDataframeMQF;
protected:
  
public:
  kDataFrameMQF();
  explicit kDataFrameMQF(std::uint64_t kSize);
  kDataFrameMQF(std::uint64_t kSize, hashingModes hash_mode);
  kDataFrameMQF(QF *mqf, readingModes RM, hashingModes HM, map<string, int> params);
  kDataFrameMQF(std::uint64_t ksize,uint8_t q,uint8_t fixedCounterSize,uint8_t tagSize
    ,double falsePositiveRate);

  kDataFrameMQF(std::uint64_t ksize, uint8_t q, hashingModes hash_mode);

  kDataFrameMQF(QF* mqf,std::uint64_t ksize,double falsePositiveRate);
  //count histogram is array where count of kmers repeated n times is found at index n. index 0 holds number of distinct kmers.
  kDataFrameMQF(std::uint64_t ksize,vector<std::uint64_t> countHistogram,uint8_t tagSize
    ,double falsePositiveRate);
  kDataFrameMQF(std::uint64_t kSize,vector<std::uint64_t> kmersHistogram);
  kDataFrameMQF(std::uint64_t kSize,uint64_t nKmers);

  ~kDataFrameMQF(){
    qf_destroy(mqf);
    delete mqf;
  }
  void reserve (std::uint64_t n);
  void reserve (vector<std::uint64_t> countHistogram);
  kDataFrame* getTwin();

  static std::uint64_t estimateMemory(std::uint64_t nslots,std::uint64_t slotSize,
    std::uint64_t fcounter, std::uint64_t tagSize);


  static void estimateParameters(vector<std::uint64_t> countHistogram,
    std::uint64_t numHashBits,std::uint64_t tagSize,
  std::uint64_t *res_noSlots,std::uint64_t *res_fixedSizeCounter, std::uint64_t *res_memory);


  bool setOrder(const string &kmer, std::uint64_t count);
  bool setOrder(std::uint64_t kmer, std::uint64_t count);
  
  bool insert(const string &kmer);
  bool insert(std::uint64_t kmer);
  std::uint64_t getkmerOrder(const string &kmer);
  std::uint64_t getkmerOrder(std::uint64_t kmer);


  bool erase(const string &kmer);
  bool erase(std::uint64_t kmer);

  std::uint64_t size();
/// max_size function returns the estimated maximum number of kmers that the kDataframeMQF can hold.
/*! The number of kmers is estimated as if all the kmers repeated 2^(fixed counter size)-1 times.*/
  std::uint64_t max_size();
  float load_factor();
  float max_load_factor();


  QF* getMQF(){
    return mqf;
  }

  void serialize(string filePath);
  static kDataFrame* load(string filePath);

  kDataFrameIterator begin();
 // kDataFrameIterator end();
  kDataFrameIterator find(const string &kmer);
    kDataFrameIterator find(uint64_t kmer);

    bool kmerExist(string kmerS);
    bool kmerExist(uint64_t kmer);
};


// kDataFrameMAPIterator

class kDataFrameMAPIterator:public _kDataFrameIterator{
private:
    kDataFrameMAP* origin;
    kmerDecoder * KD;
public:
    std::map<std::uint64_t, std::uint64_t>::iterator iterator;
    kDataFrameMAPIterator(std::map<std::uint64_t, std::uint64_t>::iterator,kDataFrameMAP* origin,std::uint64_t kSize);
    kDataFrameMAPIterator(const kDataFrameMAPIterator&);
    kDataFrameMAPIterator& operator ++ (int);
    _kDataFrameIterator* clone();
    std::uint64_t getHashedKmer();
    string getKmer();
    uint64_t getOrder();
    std::uint64_t getCount();
    bool setOrder(std::uint64_t count);
    void endIterator();
    bool operator ==(const _kDataFrameIterator& other);
    bool operator !=(const _kDataFrameIterator& other);
    ~kDataFrameMAPIterator();
};

// kDataFrameBMQF _____________________________

class kDataFrameBMQF: public kDataFrame{

private:
  bufferedMQF* bufferedmqf;
  double falsePositiveRate;
  std::uint64_t hashbits;
  __uint128_t range;
  static bool isEnough(vector<std::uint64_t> histogram,std::uint64_t noSlots,std::uint64_t fixedSizeCounter,std::uint64_t slotSize);
  friend class kDataframeBMQF;
  string fileName;
public:
  kDataFrameBMQF();
  kDataFrameBMQF(std::uint64_t kSize,string path);
  kDataFrameBMQF(std::uint64_t kSize,uint64_t nKmers,string path);
  kDataFrameBMQF(std::uint64_t ksize,uint8_t q,uint8_t fixedCounterSize,uint8_t tagSize,double falsePositiveRate,string path);
  kDataFrameBMQF(bufferedMQF* bufferedmqf,std::uint64_t ksize,double falsePositiveRate);
  kDataFrameBMQF(bufferedMQF* bufferedmqf, readingModes RM, hashingModes HM, map<string, int> params);
  //count histogram is array where count of kmers repeated n times is found at index n. index 0 holds number of distinct kmers.
  kDataFrameBMQF(std::uint64_t ksize,vector<std::uint64_t> countHistogram,uint8_t tagSize
    ,double falsePositiveRate);
  ~kDataFrameBMQF(){
    delete bufferedmqf;
  }
  void reserve (std::uint64_t n);
  void reserve (vector<std::uint64_t> countHistogram);

  kDataFrame* getTwin();

  static std::uint64_t estimateMemory(std::uint64_t nslots,std::uint64_t slotSize,
    std::uint64_t fcounter, std::uint64_t tagSize);


  static void estimateParameters(vector<std::uint64_t> countHistogram,
    std::uint64_t numHashBits,std::uint64_t tagSize,
  std::uint64_t *res_noSlots,std::uint64_t *res_fixedSizeCounter, std::uint64_t *res_memory);



  bool insert(const string &kmer);
  bool insert(std::uint64_t kmer);
  bool setOrder(const string &kmer, std::uint64_t count);
  bool setOrder(std::uint64_t kmer, std::uint64_t count);


  bool setCount(const string &kmer, std::uint64_t N) override;
  bool setCount(std::uint64_t kmer,std::uint64_t N) override;
  



  std::uint64_t getkmerOrder(const string &kmer);
  std::uint64_t getkmerOrder(std::uint64_t kmer);


  bool erase(const string &kmer);
  bool erase(std::uint64_t kmer);

  std::uint64_t size();
/// max_size function returns the estimated maximum number of kmers that the kDataframeBMQF can hold.
/*! The number of kmers is estimated as if all the kmers repeated 2^(fixed counter size)-1 times.*/
  std::uint64_t max_size();
  float load_factor();
  float max_load_factor();


  bufferedMQF* getBMQF(){
    return bufferedmqf;
  }

  void serialize(string filePath);
  static kDataFrame* load(string filePath);

  kDataFrameIterator begin();
 // kDataFrameIterator end();
  kDataFrameIterator find(const string &kmer);
  kDataFrameIterator find(uint64_t kmer);
  string getFilename();

  void deleteMemoryBuffer();
  bool kmerExist(string kmer);
  bool kmerExist(uint64_t kmer);
};

// kDataFrameMAP _____________________________

class kDataFrameMAP : public kDataFrame
{
private:
  std::map<std::uint64_t, std::uint64_t> MAP;
public:
  kDataFrameMAP();
  kDataFrameMAP(std::uint64_t ksize);
  kDataFrameMAP(std::uint64_t kSize,vector<std::uint64_t> kmersHistogram);
  kDataFrameMAP(std::uint64_t kSize,uint64_t nKmers);
  kDataFrameMAP(readingModes RM, hashingModes HM, map<string, int> params);
  kDataFrame* getTwin();
  void reserve (std::uint64_t n);
  void reserve (vector<std::uint64_t> countHistogram);

  bool kmerExist(string kmer);
  bool kmerExist(uint64_t kmer);

  bool setOrder(const string &kmer, std::uint64_t count);
  bool setOrder(std::uint64_t kmer, std::uint64_t count);
  bool insert(const string &kmer);
  bool insert(std::uint64_t kmer);
  std::uint64_t getkmerOrder(const string &kmer);
  std::uint64_t getkmerOrder(std::uint64_t kmerS);
  bool erase(const string &kmer);
  bool erase(std::uint64_t kmer);

  std::uint64_t size();
  std::uint64_t max_size();
  float load_factor();
  float max_load_factor();
  kDataFrameIterator begin();
  kDataFrameIterator end();
  kDataFrameIterator find(const string &kmer);
  kDataFrameIterator find(uint64_t kmer);
  std::uint64_t bucket(string kmer);
  void serialize(string filePath);
  static kDataFrame *load(string filePath);

    ~kDataFrameMAP() {
        this->MAP.clear();
    }
};


// kDataFramePHMAPIterator

class kDataFramePHMAPIterator : public _kDataFrameIterator {
private:
    kDataFramePHMAP *origin;
    kmerDecoder * KD;
public:
    flat_hash_map<std::uint64_t, std::uint64_t>::iterator iterator;
    kDataFramePHMAPIterator(flat_hash_map<std::uint64_t, std::uint64_t>::iterator, kDataFramePHMAP *origin, std::uint64_t kSize);

    kDataFramePHMAPIterator(const kDataFramePHMAPIterator &);

    kDataFramePHMAPIterator &operator++(int);

    _kDataFrameIterator *clone();

    std::uint64_t getHashedKmer();

    string getKmer();
    uint64_t getOrder();
    std::uint64_t getCount();

    bool setOrder(std::uint64_t count);

    void endIterator();

    bool operator==(const _kDataFrameIterator &other);

    bool operator!=(const _kDataFrameIterator &other);

    ~kDataFramePHMAPIterator();
};


// kDataFramePHMAP _____________________________

class kDataFramePHMAP : public kDataFrame {
private:
    flat_hash_map<std::uint64_t, std::uint64_t> MAP;
public:
    kDataFramePHMAP();

    explicit kDataFramePHMAP(uint64_t ksize);
    kDataFramePHMAP(std::uint64_t kSize,uint64_t nKmers);
    kDataFramePHMAP(readingModes RM, hashingModes hash_mode, map<string, int> params);
    kDataFramePHMAP(uint64_t ksize, hashingModes hash_mode);
    kDataFramePHMAP(uint64_t kSize,vector<uint64_t> kmersHistogram);







    kDataFrame *getTwin();

    void reserve(std::uint64_t n);
    void reserve (vector<std::uint64_t> countHistogram);

    bool kmerExist(string kmer);
    bool kmerExist(uint64_t kmer);

    bool setOrder(const string &kmer, std::uint64_t count);
    bool setOrder(std::uint64_t kmer, std::uint64_t count);

    bool insert(const string &kmer);
    bool insert(std::uint64_t kmer);

    std::uint64_t getkmerOrder(const string &kmer);
    std::uint64_t getkmerOrder(std::uint64_t kmer);

    bool erase(const string &kmer);
    bool erase(std::uint64_t kmer);

    std::uint64_t size();

    std::uint64_t max_size();

    float load_factor();

    float max_load_factor();

    kDataFrameIterator begin();

    kDataFrameIterator end();
    kDataFrameIterator find(const string &kmer);
    kDataFrameIterator find(uint64_t kmer);

    std::uint64_t bucket(string kmer);

    void serialize(string filePath);

    static kDataFrame *load(string filePath);

    ~kDataFramePHMAP() {
        this->MAP.clear();
    }
};

class kDataFrameBMQFIterator:public _kDataFrameIterator{
private:
    bufferedMQFIterator* qfi;
    kmerDecoder * KD;
    bufferedMQF* mqf;
public:
    kDataFrameBMQFIterator(bufferedMQF*,std::uint64_t kSize,kmerDecoder* h);
    kDataFrameBMQFIterator(const kDataFrameBMQFIterator&);
    kDataFrameBMQFIterator(bufferedMQF*,bufferedMQFIterator* qfi,std::uint64_t kSize,kmerDecoder* KD);
    kDataFrameBMQFIterator& operator ++ (int);
    _kDataFrameIterator* clone();
    std::uint64_t getHashedKmer();
    string getKmer();
    kDataFrame *getTwin();
    uint64_t getOrder();
    std::uint64_t getCount();
    bool setOrder(std::uint64_t count);
    void endIterator();
    bool operator ==(const _kDataFrameIterator& other);
    bool operator !=(const _kDataFrameIterator& other);
    ~kDataFrameBMQFIterator();
};


class kDataFrameBlight;


class kDataFrameBlightIterator : public _kDataFrameIterator {
private:
    kmer_Set_Light_iterator iterator;
    kDataFrameBlight *origin;
public:
    kDataFrameBlightIterator(kmer_Set_Light_iterator it, kDataFrameBlight *origin, std::uint64_t kSize);

    kDataFrameBlightIterator(const kDataFrameBlightIterator &);

    kDataFrameBlightIterator &operator++(int);

    _kDataFrameIterator *clone();

    std::uint64_t getHashedKmer();

    string getKmer();

    uint64_t getOrder();
    uint64_t getCount();

    bool setOrder(uint64_t count);

    void endIterator();

    bool operator==(const _kDataFrameIterator &other);

    bool operator!=(const _kDataFrameIterator &other);

    ~kDataFrameBlightIterator();
};


// kDataFrameBlight _____________________________

class kDataFrameBlight : public kDataFrame {
private:
    kmer_Set_Light* blight_index;

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
    bool kmerExist(uint64_t kmer);

    bool setOrder(const string &kmer, std::uint64_t count);
    bool setOrder(std::uint64_t kmer, std::uint64_t count);

    bool insert(const string &kmer);
    bool insert(std::uint64_t kmer);

    std::uint64_t getkmerOrder(const string &kmer);
    std::uint64_t getkmerOrder(std::uint64_t kmer);

    bool erase(const string &kmer);
    bool erase(std::uint64_t kmer);

    std::uint64_t size();

    std::uint64_t max_size();

    float load_factor();

    float max_load_factor();

    kDataFrameIterator begin();

    // kDataFrameIterator end();
    kDataFrameIterator find(const string &kmer);
    kDataFrameIterator find(uint64_t kmer);


    void serialize(string filePath);

    static kDataFrame *load(string filePath);

    ~kDataFrameBlight() = default;


    template<typename T,typename Container>
    T getKmerColumnValue(const string& columnName,string kmer);

    template<typename T,typename Container>
    void setKmerColumnValue(const string& columnName,string kmer, T value);






};
template<typename T,typename Container>
T kDataFrameBlight::getKmerColumnValue(const string& columnName,string kmer)
{
  std::uint64_t kmerOrder=getkmerOrder(kmer);
  return ((Container*)columns[columnName])->get(kmerOrder);
}
template<typename T,typename Container>
void kDataFrameBlight::setKmerColumnValue(const string& columnName,string kmer,T value)
{
  std::uint64_t kmerOrder=getkmerOrder(kmer);
  ((Container*)columns[columnName])->insert(value,kmerOrder);
}






#endif
