#ifndef _kDataFRAME_H_
#define _kDataFRAME_H_

#include <vector>
#include <stdint.h>
#include "gqf.h"
#include <iostream>
#include <parallel_hashmap/phmap.h>
#include "kmerDecoder.hpp"
#include <any>
#include "bufferedMQF.h"
#include "defaultColumn.hpp"

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
  kmerRow(){
    kmer="";
    hashedKmer=0;
    count=0;
  }
  kmerRow(string kmer,uint64_t hashedKmer,uint64_t count)
  {
    this->kmer=kmer;
    this->hashedKmer=hashedKmer;
    this->count=count;
  }
  kmerRow(const kmerRow& other)
  {
    kmer=other.kmer;
    hashedKmer=other.hashedKmer;
    count=other.count;
  }

  kmerRow copy(const kmerRow& other)
  {
    return * (new kmerRow(other));
  }

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

class _kDataFrameIterator{
protected:
  uint64_t kSize;
public:
  _kDataFrameIterator(){}
  _kDataFrameIterator(uint64_t k):kSize(k){}
  virtual _kDataFrameIterator& operator ++ (int)=0;
  virtual _kDataFrameIterator* clone()=0;
  virtual uint64_t getHashedKmer()=0;
  virtual string getKmer()=0;
  virtual uint64_t getCount()=0;
  virtual bool setCount(uint64_t count)=0;
  virtual bool operator ==(const _kDataFrameIterator& other)=0;
  virtual bool operator !=(const _kDataFrameIterator& other)=0;
  virtual ~_kDataFrameIterator(){};

};

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
  uint64_t getHashedKmer(){
    return iterator->getHashedKmer();
  };
  /// Returns the current kmer
  string getKmer(){
    return iterator->getKmer();
  }
  /// Returns the count of the current kmer
  uint64_t getCount(){
    return iterator->getCount();
  }
  /// sets the count of the current kmer
  bool setCount(uint64_t count){
    return iterator->setCount(count);
  }
  kmerRow operator*(){
    return kmerRow(iterator->getKmer(),
                   iterator->getHashedKmer(),
                   iterator->getCount()
                  );
  }
  ~kDataFrameIterator(){
    delete iterator;
  }
};


class kDataFrameMQFIterator:public _kDataFrameIterator{
private:
  QFi* qfi;
  kmerDecoder * KD;
public:
  kDataFrameMQFIterator(QF*,uint64_t kSize,kmerDecoder* KD);
  kDataFrameMQFIterator(QFi*,uint64_t kSize,kmerDecoder* KD);
  kDataFrameMQFIterator(const kDataFrameMQFIterator&);
  kDataFrameMQFIterator& operator ++ (int);
  _kDataFrameIterator* clone();
  uint64_t getHashedKmer();
  string getKmer();
  uint64_t getCount();
  bool setCount(uint64_t count);
  void endIterator();
  bool operator ==(const _kDataFrameIterator& other);
  bool operator !=(const _kDataFrameIterator& other);
   ~kDataFrameMQFIterator();
};


class kDataFrame{
protected:
  uint64_t kSize;
  string class_name; // Default = MQF, change if MAP. Temporary until resolving #17
  bool isStatic;
  bool isKmersOrderComputed;
  unordered_map<string, Column*> columns;

  unordered_map<string,uint32_t> orderCheckpoints;
  Column* defaultColumn;
  virtual void preprocessKmerOrder();
  virtual uint64_t getkmerOrder(string kmer);
public:
    kmerDecoder * KD;
  virtual string get_class_name(){ return class_name;}  // Temporary until resolving #17
  kDataFrame();
  kDataFrame(uint8_t kSize);


  virtual ~kDataFrame(){
    delete KD;
  }
/// creates a new kDataframe using the same parameters as the current kDataFrame.
/*! It is like clone but without copying the data */
  virtual kDataFrame* getTwin()=0;
/// request a capacity change so that the kDataFrame can approximately hold at least n kmers
  virtual void reserve (uint64_t n )=0;
/// request a capacity change so that the kDataFrame can approximately hold kmers with countHistogram distribution
  virtual void reserve (vector<uint64_t> countHistogram)=0;
/// insert the kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
  virtual bool insert(string kmer)=0;
/// insert the kmer N time in the kDataFrame, or increment the kmer count with N if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
  virtual bool insert(string kmer,uint64_t N)=0;
/// insert the hashed kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
    virtual bool insert(uint64_t kmer)=0;
/// insert the hashed kmer N time in the kDataFrame, or increment the kmer count with N if it is already exists.
/*! Returns bool value indicating whether the hashed kmer is inserted or not*/
    virtual bool insert(uint64_t kmer,uint64_t N)=0;
  /// insert the kmer in the kmer row time in the kDataFrame, or increment the kmer count with the count in the row if it is already exists.
  /*! Returns bool value indicating whether the kmer is inserted or not*/
  bool insert(kmerRow k);
/// set the kmer's count to N time in the kDataFrame
/*! Returns bool value indicating whether the kmer is inserted or not.
The difference between setCount and insert is that setCount set the count to N no matter the previous kmer count was*/
  virtual bool setCount(string kmer,uint64_t N)=0;
  virtual bool setCount(uint64_t kmer,uint64_t N)=0;
/// returns the count of the kmer in the kDataFrame, i.e. the number of times the kmer is inserted in the kdataFrame.
  virtual uint64_t getCount(string kmer)=0;
  virtual uint64_t getCount(uint64_t kmer)=0;
// Removes  a kmer from the kDataFrame
/*! Returns bool value indicating whether the kmer is erased or not*/
  virtual bool erase(string kmer)=0;
  virtual bool erase(uint64_t kmer)=0;

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
///Returns an iterator at the specific kmer.
  virtual kDataFrameIterator find(string kmer)=0;

  virtual void serialize(string filePath)=0;
/// Returns the kmerDecoder used by kDataframe.
  kmerDecoder* getkmerDecoder(){
    return KD;
  };



  static kDataFrame* load(string filePath);
  void save(string filePath);


  uint64_t getkSize(){return kSize;}

  // duplicate for easier name in python, getkSize won't be wrapped
  uint64_t ksize(){return kSize;}

  void setkSize(uint64_t k){kSize=k;}


  void addColumn(string columnName, Column*);

  template<typename T,typename Container>
  T getKmerColumnValue(string columnName,string kmer);

  template<typename T,typename Container>
  void setKmerColumnValue(string columnName,string kmer, T value);



  void changeDefaultColumnType(Column*);

  template<typename T,typename Container>
  T getKmerDefaultColumnValue(string kmer);

  template<typename T,typename Container>
  T getKmerDefaultColumnValue(uint64_t kmer);

  template<typename T,typename Container>
  void setKmerDefaultColumnValue(string kmer, T value);

  template<typename T,typename Container>
  void setKmerDefaultColumnValue(uint64_t kmer, T value);


};



class kDataFrameMQF: public kDataFrame{

private:
  QF* mqf;
  double falsePositiveRate;
  uint64_t hashbits;
  __uint128_t range;
  static bool isEnough(vector<uint64_t> histogram,uint64_t noSlots,uint64_t fixedSizeCounter,uint64_t slotSize);
  friend class kDataframeMQF;
protected:
  void preprocessKmerOrder();
  uint64_t getkmerOrder(string kmer);
public:
  kDataFrameMQF();
  kDataFrameMQF(uint64_t kSize);
  kDataFrameMQF(uint64_t kSize, int mode);
  kDataFrameMQF(uint64_t ksize,uint8_t q,uint8_t fixedCounterSize,uint8_t tagSize
    ,double falsePositiveRate);

  kDataFrameMQF(uint64_t ksize, uint8_t q, int mode);

  kDataFrameMQF(QF* mqf,uint64_t ksize,double falsePositiveRate);
  //count histogram is array where count of kmers repeated n times is found at index n. index 0 holds number of distinct kmers.
  kDataFrameMQF(uint64_t ksize,vector<uint64_t> countHistogram,uint8_t tagSize
    ,double falsePositiveRate);
  kDataFrameMQF(uint64_t kSize,vector<uint64_t> kmersHistogram);

  ~kDataFrameMQF(){
    qf_destroy(mqf);
    delete mqf;
  }
  void reserve (uint64_t n);
  void reserve (vector<uint64_t> countHistogram);
  kDataFrame* getTwin();

  static uint64_t estimateMemory(uint64_t nslots,uint64_t slotSize,
    uint64_t fcounter, uint64_t tagSize);


  static void estimateParameters(vector<uint64_t> countHistogram,
    uint64_t numHashBits,uint64_t tagSize,
  uint64_t *res_noSlots,uint64_t *res_fixedSizeCounter, uint64_t *res_memory);


  bool setCount(string kmer,uint64_t count);
  bool setCount(uint64_t kmer, uint64_t count);
  bool insert(string kmer,uint64_t count);
  bool insert(string kmer);
  bool insert(uint64_t kmer, uint64_t count);
  bool insert(uint64_t kmer);
  uint64_t getCount(string kmer);
  uint64_t getCount(uint64_t kmer);


  bool erase(string kmer);
  bool erase(uint64_t kmer);

  uint64_t size();
/// max_size function returns the estimated maximum number of kmers that the kDataframeMQF can hold.
/*! The number of kmers is estimated as if all the kmers repeated 2^(fixed counter size)-1 times.*/
  uint64_t max_size();
  float load_factor();
  float max_load_factor();


  QF* getMQF(){
    return mqf;
  }

  void serialize(string filePath);
  static kDataFrame* load(string filePath);

  kDataFrameIterator begin();
  kDataFrameIterator end();
  kDataFrameIterator find(string kmer);
};


// kDataFrameMAPIterator

class kDataFrameMAPIterator:public _kDataFrameIterator{
private:
    std::map<uint64_t, uint64_t>::iterator iterator;
    kDataFrameMAP* origin;
    kmerDecoder * KD;
public:
    kDataFrameMAPIterator(std::map<uint64_t, uint64_t>::iterator,kDataFrameMAP* origin,uint64_t kSize);
    kDataFrameMAPIterator(const kDataFrameMAPIterator&);
    kDataFrameMAPIterator& operator ++ (int);
    _kDataFrameIterator* clone();
    uint64_t getHashedKmer();
    string getKmer();
    uint64_t getCount();
    bool setCount(uint64_t count);
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
  uint64_t hashbits;
  __uint128_t range;
  static bool isEnough(vector<uint64_t> histogram,uint64_t noSlots,uint64_t fixedSizeCounter,uint64_t slotSize);
  friend class kDataframeBMQF;
  kDataFrameIterator* endIterator;
  string fileName;
public:
  kDataFrameBMQF();
  kDataFrameBMQF(uint64_t kSize,string path);
  kDataFrameBMQF(uint64_t ksize,uint8_t q,uint8_t fixedCounterSize,uint8_t tagSize,double falsePositiveRate,string path);
  kDataFrameBMQF(bufferedMQF* bufferedmqf,uint64_t ksize,double falsePositiveRate);
  //count histogram is array where count of kmers repeated n times is found at index n. index 0 holds number of distinct kmers.
  kDataFrameBMQF(uint64_t ksize,vector<uint64_t> countHistogram,uint8_t tagSize
    ,double falsePositiveRate);
  ~kDataFrameBMQF(){
    delete bufferedmqf;
    delete endIterator;
  }
  void reserve (uint64_t n);
  void reserve (vector<uint64_t> countHistogram);

  kDataFrame* getTwin();

  static uint64_t estimateMemory(uint64_t nslots,uint64_t slotSize,
    uint64_t fcounter, uint64_t tagSize);


  static void estimateParameters(vector<uint64_t> countHistogram,
    uint64_t numHashBits,uint64_t tagSize,
  uint64_t *res_noSlots,uint64_t *res_fixedSizeCounter, uint64_t *res_memory);



  bool insert(string kmer,uint64_t count);
  bool insert(string kmer);
  bool insert(uint64_t kmer,uint64_t count);
  bool insert(uint64_t kmer);
  bool setCount(string kmer,uint64_t count);
  bool setCount(uint64_t kmer, uint64_t count);





  uint64_t getCount(string kmer);
  uint64_t getCount(uint64_t kmer);


  bool erase(string kmer);
  bool erase(uint64_t kmer);

  uint64_t size();
/// max_size function returns the estimated maximum number of kmers that the kDataframeBMQF can hold.
/*! The number of kmers is estimated as if all the kmers repeated 2^(fixed counter size)-1 times.*/
  uint64_t max_size();
  float load_factor();
  float max_load_factor();


  bufferedMQF* getBMQF(){
    return bufferedmqf;
  }

  void serialize(string filePath);
  static kDataFrame* load(string filePath);

  kDataFrameIterator begin();
  kDataFrameIterator end();
  kDataFrameIterator find(string kmer);

  string getFilename();
};

// kDataFrameMAP _____________________________

class kDataFrameMAP : public kDataFrame
{
private:
  std::map<uint64_t, uint64_t> MAP;
public:
  kDataFrameMAP();
  kDataFrameMAP(uint64_t ksize);
  kDataFrameMAP(uint64_t kSize,vector<uint64_t> kmersHistogram);
  kDataFrame* getTwin();
  void reserve (uint64_t n);
  void reserve (vector<uint64_t> countHistogram);

  inline bool kmerExist(string kmer);

  bool setCount(string kmer, uint64_t count);
  bool setCount(uint64_t kmer, uint64_t count);
  bool insert(string kmer);
  bool insert(string kmer, uint64_t count);
  bool insert(uint64_t kmer, uint64_t count);
  bool insert(uint64_t kmer);
  uint64_t getCount(string kmer);
  uint64_t getCount(uint64_t kmerS);
  bool erase(string kmer);
  bool erase(uint64_t kmer);

  uint64_t size();
  uint64_t max_size();
  float load_factor();
  float max_load_factor();
  kDataFrameIterator begin();
  kDataFrameIterator end();
  kDataFrameIterator find(string kmer);

  uint64_t bucket(string kmer);
  void serialize(string filePath);
  static kDataFrame *load(string filePath);

    ~kDataFrameMAP() {
        this->MAP.clear();
    }
};


// kDataFramePHMAPIterator

class kDataFramePHMAPIterator : public _kDataFrameIterator {
private:
    flat_hash_map<uint64_t, uint64_t>::iterator iterator;
    kDataFramePHMAP *origin;
    kmerDecoder * KD;
public:
    kDataFramePHMAPIterator(flat_hash_map<uint64_t, uint64_t>::iterator, kDataFramePHMAP *origin, uint64_t kSize);

    kDataFramePHMAPIterator(const kDataFramePHMAPIterator &);

    kDataFramePHMAPIterator &operator++(int);

    _kDataFrameIterator *clone();

    uint64_t getHashedKmer();

    string getKmer();

    uint64_t getCount();

    bool setCount(uint64_t count);

    void endIterator();

    bool operator==(const _kDataFrameIterator &other);

    bool operator!=(const _kDataFrameIterator &other);

    ~kDataFramePHMAPIterator();
};


// kDataFramePHMAP _____________________________

class kDataFramePHMAP : public kDataFrame {
private:
    flat_hash_map<uint64_t, uint64_t> MAP;
public:
    kDataFramePHMAP();

    kDataFramePHMAP(uint64_t ksize);
    kDataFramePHMAP(uint64_t ksize, int mode);
    kDataFramePHMAP(uint64_t kSize,vector<uint64_t> kmersHistogram);

    kDataFrame *getTwin();

    void reserve(uint64_t n);
    void reserve (vector<uint64_t> countHistogram);

    inline bool kmerExist(string kmer);

    bool setCount(string kmer, uint64_t count);
    bool setCount(uint64_t kmer, uint64_t count);

    bool insert(string kmer);

    bool insert(string kmer, uint64_t count);

    bool insert(uint64_t kmer, uint64_t count);

    bool insert(uint64_t kmer);

    uint64_t getCount(string kmer);
    uint64_t getCount(uint64_t kmer);

    bool erase(string kmer);
    bool erase(uint64_t kmer);

    uint64_t size();

    uint64_t max_size();

    float load_factor();

    float max_load_factor();

    kDataFrameIterator begin();

    kDataFrameIterator end();
    kDataFrameIterator find(string kmer);

    uint64_t bucket(string kmer);

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
    kDataFrameBMQFIterator(bufferedMQF*,uint64_t kSize,kmerDecoder* h);
    kDataFrameBMQFIterator(const kDataFrameBMQFIterator&);
    kDataFrameBMQFIterator(bufferedMQF*,bufferedMQFIterator* qfi,uint64_t kSize,kmerDecoder* KD);
    kDataFrameBMQFIterator& operator ++ (int);
    _kDataFrameIterator* clone();
    uint64_t getHashedKmer();
    string getKmer();
    kDataFrame *getTwin();
    uint64_t getCount();
    bool setCount(uint64_t count);
    void endIterator();
    bool operator ==(const _kDataFrameIterator& other);
    bool operator !=(const _kDataFrameIterator& other);
    ~kDataFrameBMQFIterator();
};


#endif
