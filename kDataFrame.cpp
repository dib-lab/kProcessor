#include "kDataFrame.hpp"
#include "KmerCounter/kmer.h"
#include <iostream>
#include <fstream>
using namespace std;

kDataFrameMQFIterator::kDataFrameMQFIterator(QF* mqf)
:_kDataFrameIterator()
{
  qfi= new QFi();
  qf_iterator(mqf,qfi,0);
}
void kDataFrameMQFIterator::operator ++ (int){
  qfi_next(qfi);
  end=qfi_end(qfi)==1;
}
kmerRow kDataFrameMQFIterator::operator * (){
  kmerRow res;
  qfi_get(qfi,&res.kmerHash,&res.tag,&res.count);
  return res;
}

kDataFrame::kDataFrame(){
  kSize=31;
  autoResize=false;
  hashFunctions.push_back(new IntegerHasher(BITMASK(2*kSize)));
}
kDataFrame::kDataFrame(double falsePositiveRate,uint8_t k_size){
  kSize=k_size;
  autoResize=false;
  Hasher* hasher;
  if(falsePositiveRate==0){
    hasher=new IntegerHasher(BITMASK(2*kSize));
  }
  else if(falsePositiveRate<1){
  hasher=new MumurHasher(2038074761);
  }
  hashFunctions.push_back(hasher);
}

uint64_t kDataFrame::hashKmer(string kmer){
  return hashFunctions[0]->hash(kmer::str_to_int(kmer));
}
kDataFrame* kDataFrame::load(string filePath){
  return kDataFrameMQF::load(filePath);
}

kDataFrameMQF::kDataFrameMQF():kDataFrame(){
  mqf=new QF();
  qf_init(mqf, (1ULL<<16), 2*kSize, 0,2, true, "", 2038074761);
}
kDataFrameMQF::kDataFrameMQF(uint64_t ksize,uint8_t q,uint8_t fixedCounterSize,uint8_t tagSize,double falsePositiveRate):
kDataFrame(falsePositiveRate,ksize){
  mqf=new QF();
  qf_init(mqf, (1ULL<<q), 2*ksize,tagSize,fixedCounterSize, true, "", 2038074761);
}
kDataFrameMQF::kDataFrameMQF(QF* mqf,uint64_t ksize,double falsePositiveRate):
kDataFrame(falsePositiveRate,ksize)
{
  this->mqf=mqf;
}

bool kDataFrameMQF::setCounter(string kmer,uint64_t count){
  uint64_t hash=hashKmer(kmer)%mqf->metadata->range;
  uint64_t currentCount=qf_count_key(mqf,hash);
  if(currentCount>count){
    qf_remove(mqf,hash,currentCount-count,false,false);
  }
  else{
    qf_insert(mqf,hash,count-currentCount,false,false);
  }
  return true;
}
bool kDataFrameMQF::incrementCounter(string kmer,uint64_t count){
  uint64_t hash=hashKmer(kmer)%mqf->metadata->range;
  qf_insert(mqf,hash,count,true,true);
  return true;
}
uint64_t kDataFrameMQF::getCounter(string kmer){
  uint64_t hash=hashKmer(kmer)%mqf->metadata->range;
  return qf_count_key(mqf,hash);
}

bool kDataFrameMQF::setTag(string kmer,uint64_t tag){
  uint64_t hash=hashKmer(kmer)%mqf->metadata->range;
  return qf_add_tag(mqf,hash,tag,false,false)==1;
}
uint64_t kDataFrameMQF::getTag(string kmer){
  uint64_t hash=hashKmer(kmer)%mqf->metadata->range;
  return qf_get_tag(mqf,hash);
}

bool kDataFrameMQF::removeKmer(string kmer){
  uint64_t hash=hashKmer(kmer)%mqf->metadata->range;
  uint64_t currentCount=qf_count_key(mqf,hash);
  qf_remove(mqf,hash,currentCount,true,true);
  return true;
}

uint64_t kDataFrameMQF::size(){
  return mqf->metadata->xnslots;
}
uint64_t kDataFrameMQF::filled_space(){
  return (uint64_t)qf_space(mqf);
}
bool kDataFrameMQF::isFull(){
  return mqf->metadata->noccupied_slots>=mqf->metadata->maximum_occupied_slots;
}

void kDataFrameMQF::save(string filePath){
  ofstream file(filePath+".extra");
  file<<kSize<<endl;
  // uint64_t legendSize=tagsLegend.size();
  // file<<legendSize<<endl;
  // auto it = tagsLegend.begin();
  // while(it==tagsLegend.end())
  // {
  //   file<<it->first<<" "<<it->second<<endl;
  //   it++;
  // }
  // file.close();
  qf_serialize(mqf,(filePath+".mqf").c_str());
}
kDataFrame* kDataFrameMQF::load(string filePath){
  ifstream file(filePath+".extra");
  uint64_t filekSize;
  file>>filekSize;
  QF* mqf=new QF();
  qf_deserialize(mqf,(filePath+".mqf").c_str());
  return new kDataFrameMQF(mqf,filekSize,0);
}

kDataFrameIterator kDataFrameMQF::begin(){
    return *(new kDataFrameIterator(new kDataFrameMQFIterator(mqf)));
}
