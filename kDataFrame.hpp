#include "HashUtils/hashutil.h"
#include <vector>
#include <stdint.h>
#include "MQF/gqf.h"
#include "KmerCounter/kmer.h"
#include <map>

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
  map<uint64_t,string> tagsLegend;
public:
  kDataFrame();
  kDataFrame(double falsePositiveRate,uint8_t kSize);

  uint64_t hashKmer(string kmer);
  virtual bool setCounter(string kmer,uint64_t count)=0;
  virtual bool incrementCounter(string kmer,uint64_t count)=0;
  virtual uint64_t getCounter(string kmer)=0;

  virtual bool setTag(string kmer,uint64_t tag)=0;
  virtual uint64_t getTag(string kmer)=0;

  virtual bool removeKmer(string kmer)=0;

  virtual uint64_t size()=0;
  virtual uint64_t filled_space()=0;
  virtual bool isFull()=0;

  virtual kDataFrameIterator begin()=0;

  virtual void save(string filePath)=0;
  static kDataFrame* load(string filePath);

};

class kDataFrameMQF: public kDataFrame{
private:
  QF* mqf;

public:
  kDataFrameMQF();
  kDataFrameMQF(uint64_t ksize,uint8_t q,uint8_t fixedCounterSize,uint8_t tagSize,double falsePositiveRate);
  kDataFrameMQF(QF* mqf,uint64_t ksize,double falsePositiveRate);
  bool setCounter(string kmer,uint64_t count);
  bool incrementCounter(string kmer,uint64_t count);
  uint64_t getCounter(string kmer);

  bool setTag(string kmer,uint64_t tag);
  uint64_t getTag(string kmer);

  bool removeKmer(string kmer);

  uint64_t size();
  uint64_t filled_space();
  bool isFull();

  void save(string filePath);
  static kDataFrame* load(string filePath);

  kDataFrameIterator begin();
};
