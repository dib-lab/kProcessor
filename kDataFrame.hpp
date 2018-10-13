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
  map<uint64_t,vector<int> >* colorsMap;
public:
  kDataFrame();
  kDataFrame(double falsePositiveRate,uint8_t kSize);

  virtual ~kDataFrame(){

  }
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

  vector<int> getColors(string kmer);
  uint64_t getkSize(){return kSize;}
  uint64_t setkSize(uint64_t k){kSize=k;}

  string getCanonicalKmer(string kmer){
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

  static kDataFrameMQF* index(vector<kDataFrameMQF*> kframes);

  void loadIntoFastq(string sequenceFilename,int noThreads);


  std::map<uint64_t, std::vector<int> > * get_legend(){
    return mqf->metadata->tags_map;
  }

  void set_legend(std::map<uint64_t, std::vector<int> > * t){
      mqf->metadata->tags_map=t;
      this->colorsMap=t;
    }
};

class kDataFrameMAP : public kDataFrame
{
private:
  map<string, uint64_t> MAP;

public:
  kDataFrameMAP(uint64_t ksize);


  inline bool kmerExist(string kmer);

  bool setCounter(string kmer, uint64_t count);
  bool incrementCounter(string kmer, uint64_t count);
  uint64_t getCounter(string kmer);

  bool setTag(string kmer, uint64_t tag);
  uint64_t getTag(string kmer);

  bool removeKmer(string kmer);

  uint64_t size();
  uint64_t filled_space();
  bool isFull();

  kDataFrameIterator begin();

  void save(string filePath);
  static kDataFrame *load(string filePath);

  void set_legend(std::map<uint64_t, std::vector<int>> *t)
  {
    this->colorsMap = t;
  }

  std::map<uint64_t, std::vector<int>> *get_legend()
  {
    return this->colorsMap;
  }

    ~kDataFrameMAP() {
        this->MAP.clear();
        this->colorsMap->clear();
    }
};
