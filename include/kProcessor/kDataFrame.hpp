#include <HashUtils/hashutil.h>
#include <vector>
#include <stdint.h>
#include "gqf.hpp"
#include "KmerCounter/kmer.h"
#include <map>
#include <unordered_map>

using namespace std;

class kmerRow{
public:
  string kmer;
  uint64_t kmerHash;
  uint64_t count;
  uint64_t tag;
};


class _kDataFrameIterator{
protected:
  bool end;
  kmerRow current;
public:
  _kDataFrameIterator(){
    end=false;
  }
  virtual void operator ++ (int)=0;
  virtual kmerRow operator * ()=0;
  virtual bool isEnd(){return end;}
};

class kDataFrameIterator{
private:
  _kDataFrameIterator* iterator;
public:
  kDataFrameIterator(_kDataFrameIterator* it){
    iterator=it;
  }
  void operator ++ (int){(*iterator)++;}
  kmerRow operator * (){
    return *(*iterator);
  }
  bool isEnd(){
    return iterator->isEnd();
  }
};


class kDataFrameMQFIterator:public _kDataFrameIterator{
private:
  QFi* qfi;
public:
  kDataFrameMQFIterator(QF* mqf);
  void operator ++ (int);
  kmerRow operator * ();
};

class kDataFrame{
protected:
  vector<Hasher*> hashFunctions;
  uint64_t kSize;
  bool autoResize;
public:
  kDataFrame();
  kDataFrame(double falsePositiveRate,uint8_t kSize);

  virtual ~kDataFrame(){

  }
  uint64_t hashKmer(string kmer);
/// insert the kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
  virtual bool insert(string kmer)=0;
/// insert the kmer N time in the kDataFrame, or increment the kmer count with N if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
  virtual bool insert(string kmer,uint64_t N)=0;
/// set the kmer's count to N time in the kDataFrame
/*! Returns bool value indicating whether the kmer is inserted or not.
The difference between setCount and insert is that setCount set the count to N no matter the previous kmer count was*/
  virtual bool setCount(string kmer,uint64_t N)=0;
/// returns the count of the kmer in the kDataFrame, i.e. the number of times the kmer is inserted in the kdataFrame.
  virtual uint64_t count(string kmer)=0;
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

  virtual bool isFull()=0;

  virtual kDataFrameIterator begin()=0;

  virtual void save(string filePath)=0;
  static kDataFrame* load(string filePath, string method);


  uint64_t getkSize(){return kSize;}
  uint64_t setkSize(uint64_t k){kSize=k;}

  inline string getCanonicalKmer(string kmer){
    uint64_t kmerI = kmer::str_to_int(kmer);
    uint64_t kmerIR = kmer::reverse_complement(kmerI, kSize);
    uint64_t item;
    if (kmer::compare_kmers(kmerI, kmerIR))
      return kmer;
    else
      return kmer::int_to_str(kmerIR,kSize);
  }

};

class kDataFrameMQF: public kDataFrame{
private:
  QF* mqf;
  static bool isEnough(vector<uint64_t> histogram,uint64_t noSlots,uint64_t fixedSizeCounter,uint64_t slotSize);
public:
  kDataFrameMQF();
  kDataFrameMQF(uint64_t ksize,uint8_t q,uint8_t fixedCounterSize,uint8_t tagSize
    ,double falsePositiveRate);
  kDataFrameMQF(QF* mqf,uint64_t ksize,double falsePositiveRate);
  //count histogram is array where count of kmers repeated n times is found at index n. index 0 holds number of distinct kmers.
  kDataFrameMQF(uint64_t ksize,vector<uint64_t> countHistogram,uint8_t tagSize
    ,double falsePositiveRate);
  ~kDataFrameMQF(){
    qf_destroy(mqf);
    delete mqf;
  }
  static uint64_t estimateMemory(uint64_t nslots,uint64_t slotSize,
    uint64_t fcounter, uint64_t tagSize);


  static void estimateParameters(vector<uint64_t> countHistogram,
    uint64_t numHashBits,uint64_t tagSize,
  uint64_t *res_noSlots,uint64_t *res_fixedSizeCounter, uint64_t *res_memory);


  bool setCount(string kmer,uint64_t count);
  bool insert(string kmer,uint64_t count);
  bool insert(string kmer);
  uint64_t count(string kmer);


  bool erase(string kmer);

  uint64_t size();
/// max_size function returns the estimated maximum number of kmers that the kDataframeMQF can hold.
/*! The number of kmers is estimated as if all the kmers repeated 2^(fixed counter size)-1 times.*/
  uint64_t max_size();
  float load_factor();
  float max_load_factor();
  bool isFull();

  void save(string filePath);
  static kDataFrame* load(string filePath);

  kDataFrameIterator begin();







};

// kDataFrameMAP _____________________________

class kDataFrameMAP : public kDataFrame
{
private:
  unordered_map<string, uint64_t> MAP;

public:
  kDataFrameMAP();
  kDataFrameMAP(uint64_t ksize);


  inline bool kmerExist(string kmer);

  bool setCount(string kmer, uint64_t count);
  bool insert(string kmer);
  bool insert(string kmer, uint64_t count);
  uint64_t count(string kmer);


  bool erase(string kmer);

  uint64_t size();
  uint64_t max_size();
  float load_factor();
  float max_load_factor();
  bool isFull();

  kDataFrameIterator begin();

  void save(string filePath);
  static kDataFrame *load(string filePath);


    ~kDataFrameMAP() {
        this->MAP.clear();
    }
};
