#include "kDataFrame.hpp"
#include "Utils/kmer.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <limits>
#include <sstream>


#include "algorithms.hpp"
using namespace std;



inline bool fileExists(const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

kDataFrameMQFIterator::kDataFrameMQFIterator(QF* mqf,uint64_t kSize,Hasher* h)
:_kDataFrameIterator(kSize)
{
  qfi= new QFi();
  qf_iterator(mqf,qfi,0);
  this->hasher=h;
}

kDataFrameMQFIterator::kDataFrameMQFIterator(const kDataFrameMQFIterator& other):
_kDataFrameIterator(other.kSize)
{
  qfi=new QFi();
  qfi->qf=other.qfi->qf;
  qfi->run=other.qfi->run;
  qfi->current=other.qfi->current;
  qfi->cur_start_index=other.qfi->cur_start_index;
  qfi->cur_length=other.qfi->cur_length;
  qfi->num_clusters=other.qfi->num_clusters;
  qfi->c_info=other.qfi->c_info;
  hasher=other.hasher;
}
_kDataFrameIterator* kDataFrameMQFIterator::clone(){
  return new kDataFrameMQFIterator(*this);
}

kDataFrameMQFIterator& kDataFrameMQFIterator::operator ++ (int){
  qfi_next(qfi);
  return *this;
}
uint64_t kDataFrameMQFIterator::getHashedKmer(){
  uint64_t key,value,count;
  qfi_get(qfi,&key,&value,&count);
  return key;

}
string kDataFrameMQFIterator::getKmer(){
  return hasher->Ihash(getHashedKmer());
}
uint64_t kDataFrameMQFIterator::getKmerCount(){
  uint64_t key,value,count;
  qfi_get(qfi,&key,&value,&count);
  return count;
}
bool kDataFrameMQFIterator::setKmerCount(uint64_t count){
  uint64_t key,value,currentCount;
  qfi_get(qfi,&key,&value,&currentCount);
  if(currentCount>count){
    qf_remove(qfi->qf,key,currentCount-count,false,false);
  }
  else{
      qf_insert(qfi->qf,key,count-currentCount,false,false);
  }
  return true;
}
void kDataFrameMQFIterator::endIterator(){
  qfi->current=qfi->qf->metadata->xnslots;
}

bool kDataFrameMQFIterator::operator ==(const _kDataFrameIterator& other){
  if(qfi->current >= qfi->qf->metadata->xnslots &&
    ((kDataFrameMQFIterator*)&other)->qfi->current >=
     ((kDataFrameMQFIterator*)&other)->qfi->qf->metadata->xnslots )
     return true;

  return qfi->current == ((kDataFrameMQFIterator*)&other)->qfi->current;
}
bool kDataFrameMQFIterator::operator !=(const _kDataFrameIterator& other){
  if(qfi->current >= qfi->qf->metadata->xnslots &&
    ((kDataFrameMQFIterator*)&other)->qfi->current >=
     ((kDataFrameMQFIterator*)&other)->qfi->qf->metadata->xnslots )
     return false;
  if(qfi->current >= qfi->qf->metadata->xnslots &&
    ((kDataFrameMQFIterator*)&other)->qfi->current <
     ((kDataFrameMQFIterator*)&other)->qfi->qf->metadata->xnslots )
     return true;
  if(qfi->current < qfi->qf->metadata->xnslots &&
       ((kDataFrameMQFIterator*)&other)->qfi->current >=
        ((kDataFrameMQFIterator*)&other)->qfi->qf->metadata->xnslots )
        return true;
  return qfi->current != ((kDataFrameMQFIterator*)&other)->qfi->current;
}
kDataFrameMQFIterator::~kDataFrameMQFIterator(){
   delete qfi;
}


kDataFrameMAPIterator::kDataFrameMAPIterator(unordered_map<string,uint64_t>::iterator it,kDataFrameMAP* origin,uint64_t kSize)
:_kDataFrameIterator(kSize)
{
  iterator=it;
  this->origin=origin;
}

kDataFrameMAPIterator::kDataFrameMAPIterator(const kDataFrameMAPIterator& other):
_kDataFrameIterator(other.kSize)
{
  iterator=other.iterator;
  this->origin=other.origin;
}
_kDataFrameIterator* kDataFrameMAPIterator::clone(){
  return new kDataFrameMAPIterator(*this);
}

kDataFrameMAPIterator& kDataFrameMAPIterator::operator ++ (int){
  iterator++;
  return *this;
}
uint64_t kDataFrameMAPIterator::getHashedKmer(){
  return origin->getHasher()->hash(iterator->first);

}
string kDataFrameMAPIterator::getKmer(){
  return iterator->first;
}
uint64_t kDataFrameMAPIterator::getKmerCount(){
  return iterator->second;
}
bool kDataFrameMAPIterator::setKmerCount(uint64_t count){
  iterator->second=count;
}


bool kDataFrameMAPIterator::operator ==(const _kDataFrameIterator& other){
  return iterator == ((kDataFrameMAPIterator*)&other)->iterator;
}
bool kDataFrameMAPIterator::operator !=(const _kDataFrameIterator& other){
  return iterator != ((kDataFrameMAPIterator*)&other)->iterator;
}
kDataFrameMAPIterator::~kDataFrameMAPIterator(){

}



kDataFrame::kDataFrame(){
  kSize=31;
}
kDataFrame::kDataFrame(uint8_t k_size){
  kSize=k_size;

}

bool kDataFrame::empty(){
  return this->size()==0;
}

bool kDataFrame::insert(kmerRow k){
  return this->insert(k.kmer,k.count);
}


kDataFrame *kDataFrame::load(string filePath) {
  if(fileExists(filePath+".mqf"))
    return kDataFrameMQF::load(filePath);
  else if (fileExists(filePath+".map"))
    return kDataFrameMAP::load(filePath);
  else
    throw std::runtime_error("Could not open kDataFrame file");
}



kDataFrameMQF::kDataFrameMQF():kDataFrame(){
  mqf=new QF();
  qf_init(mqf, (1ULL<<16), 2*kSize, 0,2,0, true, "", 2038074761);
  hasher=(new IntegerHasher(kSize));
  falsePositiveRate=0;
  hashbits=2*kSize;
  range=(1ULL<<hashbits);
}
kDataFrameMQF::kDataFrameMQF(uint64_t ksize,uint8_t q,uint8_t fixedCounterSize,uint8_t tagSize,double falsePositiveRate):
kDataFrame(ksize){
  mqf=new QF();
  qf_init(mqf, (1ULL<<q), 2*ksize,tagSize,fixedCounterSize, 0,true, "", 2038074761);
  this->falsePositiveRate=falsePositiveRate;
  if(falsePositiveRate==0){
    hasher=(new IntegerHasher(kSize));
  }
  else if(falsePositiveRate<1){
  hasher=(new MumurHasher(2038074761));
  }
  hashbits=2*kSize;
  range=(1ULL<<hashbits);

}
kDataFrameMQF::kDataFrameMQF(uint64_t ksize):
kDataFrame(ksize){
  this->falsePositiveRate=0.0;
  hasher=(new IntegerHasher(kSize));
  hashbits=2*kSize;
  range=(1ULL<<hashbits);
  mqf=NULL;
  reserve(10000);
}

kDataFrameMQF::kDataFrameMQF(QF* mqf,uint64_t ksize,double falsePositiveRate):
kDataFrame(ksize)
{
  this->mqf=mqf;
this->falsePositiveRate=falsePositiveRate;
  if(falsePositiveRate==0){
    hasher=(new IntegerHasher(kSize));
  }
  else if(falsePositiveRate<1){
  hasher=(new MumurHasher(2038074761));
  }
  hashbits=this->mqf->metadata->key_bits;
  hashbits=2*kSize;
  range=(1ULL<<hashbits);
}
kDataFrame* kDataFrameMQF::getTwin(){
  uint64_t q=log2(mqf->metadata->nslots);
  return ((kDataFrame*)new kDataFrameMQF(kSize,q,mqf->metadata->fixed_counter_size,
    mqf->metadata->tag_bits,falsePositiveRate));
}

void kDataFrameMQF::reserve(uint64_t n)
{
  QF* old=mqf;
  mqf=new QF();
  uint64_t q=(uint64_t)ceil(log2((double)n*1.4));
  qf_init(mqf, (1ULL<<q), hashbits,0,2, 0,true, "", 2038074761);
  if(old!=NULL)
  {
    qf_migrate(old,mqf);
    qf_destroy(old);
    delete old;
  }
}

kDataFrameMQF::kDataFrameMQF(uint64_t ksize,vector<uint64_t> countHistogram,uint8_t tagSize
  ,double falsePositiveRate):
  kDataFrame(ksize)
  {
    uint64_t nSlots;
    uint64_t fixedCounterSize;
    uint64_t memory;
    kDataFrameMQF::estimateParameters(countHistogram,2*ksize,tagSize,
    &nSlots,&fixedCounterSize,&memory);
    qf_init(mqf, nSlots, 2*ksize,tagSize,fixedCounterSize, 0,true,"", 2038074761);
  }

uint64_t kDataFrameMQF::estimateMemory(uint64_t nslots,uint64_t slotSize, uint64_t fcounter, uint64_t tagSize)
     {
       uint64_t SLOTS_PER_BLOCK=64;
       uint64_t xnslots = nslots + 10*sqrt((double)nslots);
     	uint64_t nblocks = (xnslots + SLOTS_PER_BLOCK - 1) / SLOTS_PER_BLOCK;
       uint64_t blocksize=17;

       return ((nblocks)*(blocksize+8*(slotSize+fcounter+tagSize)))/1024;

     }

void kDataFrameMQF::estimateParameters(vector<uint64_t> countHistogram,
  uint64_t numHashBits,uint64_t tagSize,
uint64_t *res_noSlots,uint64_t *res_fixedSizeCounter, uint64_t *res_memory){

  uint64_t noDistinctKmers=countHistogram[0];
  *res_memory=numeric_limits<uint64_t>::max();
  for(int i=8;i<64;i++)
  {
    uint64_t noSlots=(1ULL)<<i;
    if(noSlots<noDistinctKmers)
    continue;
    bool moreWork=false;
    uint64_t slotSize=numHashBits-log2((double)noSlots);
    for(uint64_t fixedSizeCounter=1;fixedSizeCounter<slotSize;fixedSizeCounter++)
    {
      if(isEnough(countHistogram,noSlots,fixedSizeCounter,slotSize))
      {
        uint64_t tmpMem=estimateMemory(noSlots,slotSize,fixedSizeCounter,tagSize);
        if(*res_memory>tmpMem)
        {
          *res_memory=tmpMem;
          *res_fixedSizeCounter=fixedSizeCounter;
          *res_noSlots=noSlots;
          moreWork=true;
        }
        else{
          break;
        }
      }

    }
    if(!moreWork && *res_memory!=numeric_limits<uint64_t>::max())
    break;
  }
  if(*res_memory==numeric_limits<uint64_t>::max())
  {
    throw std::overflow_error("Data limits exceeds MQF capabilities(> uint64). Check if ntcard file is corrupted");
  }

}

bool kDataFrameMQF::isEnough(vector<uint64_t> histogram,uint64_t noSlots,uint64_t fixedSizeCounter,uint64_t slotSize)
{
  // cout<<"noSlots= "<<noSlots<<endl
  //     <<"fcounter= "<<fixedSizeCounter<<endl
  //     <<"slot size= "<<numHashBits<<endl;

  noSlots=(uint64_t)((double)noSlots*0.90);
  for(uint64_t i=1;i<1000;i++)
  {
    uint64_t usedSlots=1;

    if(i>((1ULL)<<fixedSizeCounter)-1)
    {
      uint64_t nSlots2=0;
      __uint128_t capacity;
      do{
        nSlots2++;
        capacity=((__uint128_t)(1ULL)<<(nSlots2*slotSize+fixedSizeCounter))-1;
      //  cout<<"slots num "<<nSlots2<<" "<<capacity<<endl;
    }while((__uint128_t)i>capacity);
      usedSlots+=nSlots2;
    }
    //cout<<"i= "<<i<<"->"<<usedSlots<<" * "<<histogram[i]<<endl;
    if(noSlots>=(usedSlots*histogram[i]))
    {
      noSlots-=(usedSlots*histogram[i]);
    }
    else
    {
    //  cout<<"failed"<<endl<<endl;
      return false;
    }

  }
  //cout<<"success"<<endl<<endl;
  return true;
}



bool kDataFrameMQF::setCount(string kmer,uint64_t count){
  uint64_t hash=hasher->hash(kmer)% mqf->metadata->range;
  uint64_t currentCount=qf_count_key(mqf,hash);
  if(currentCount>count){
    qf_remove(mqf,hash,currentCount-count,false,false);
  }
  else{
    try{
      qf_insert(mqf,hash,count-currentCount,false,false);
    }
    catch(overflow_error & e)
    {
      reserve(mqf->metadata->nslots);
      return setCount(kmer,count);
    }
  }
  return true;
}
bool kDataFrameMQF::insert(string kmer,uint64_t count){
  uint64_t hash=hasher->hash(kmer)% mqf->metadata->range;
  try{
  qf_insert(mqf,hash,count,true,true);
  }
  catch(overflow_error & e)
  {
    reserve(mqf->metadata->nslots);
    return insert(kmer,count);
  }
  return true;
}
bool kDataFrameMQF::insert(string kmer){
  if(load_factor()>0.9)
    reserve(mqf->metadata->nslots);
  uint64_t hash=hasher->hash(kmer)% mqf->metadata->range;
  try{
    qf_insert(mqf,hash,1,false,false);
  }
  catch(overflow_error & e)
  {
    reserve(mqf->metadata->nslots);
    return insert(kmer);
  }
  return true;
}
uint64_t kDataFrameMQF::count(string kmer){
  uint64_t hash=hasher->hash(kmer)% mqf->metadata->range;
  return qf_count_key(mqf,hash);
}



bool kDataFrameMQF::erase(string kmer){
  uint64_t hash=hasher->hash(kmer)% mqf->metadata->range;
  uint64_t currentCount=qf_count_key(mqf,hash);

  qf_remove(mqf,hash,currentCount,false,false);
  return true;
}

uint64_t kDataFrameMQF::size(){
  return mqf->metadata->ndistinct_elts;
}
uint64_t kDataFrameMQF::max_size(){
  return mqf->metadata->xnslots;
}
float kDataFrameMQF::load_factor(){
  return (float)qf_space(mqf)/100.0;
}
float kDataFrameMQF::max_load_factor(){
  return 0.9;
}


void kDataFrameMQF::save(string filePath){
  //filePath += ".mqf";
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
    return (kDataFrameIterator(
      (_kDataFrameIterator*)new kDataFrameMQFIterator(mqf,kSize,hasher),
      (kDataFrame*)this));
}
kDataFrameIterator kDataFrameMQF::end(){
  kDataFrameMQFIterator* it=new kDataFrameMQFIterator(mqf,kSize,hasher);
  it->endIterator();
  return (kDataFrameIterator(it,(kDataFrame*)this));
}




// kDataFrameMAP _____________________________

kDataFrameMAP::kDataFrameMAP(uint64_t ksize) {
    this->kSize = ksize;
    this->MAP=unordered_map<string,uint64_t>(1000);
    hasher=new wrapperHasher<unordered_map<string,uint64_t>::hasher >(MAP.hash_function(),ksize);
}
kDataFrameMAP::kDataFrameMAP() {
    this->kSize = 23;
    this->MAP=unordered_map<string,uint64_t>(1000);
    hasher=new wrapperHasher<unordered_map<string,uint64_t>::hasher >(MAP.hash_function(),kSize);
}
inline bool kDataFrameMAP::kmerExist(string kmerS) {
    return (this->MAP.find(kmer::canonicalKmer(kmerS)) == this->MAP.end()) ? 0 : 1;
}



bool kDataFrameMAP::insert(string kmerS, uint64_t count) {
    this->MAP[kmer::canonicalKmer(kmerS)]+=count;
    return true;
}
bool kDataFrameMAP::insert(string kmerS) {
    this->MAP[kmer::canonicalKmer(kmerS)]++;
    return true;
}


bool kDataFrameMAP::setCount(string kmerS, uint64_t tag) {
    this->MAP[kmer::canonicalKmer(kmerS)]=tag;
    return true;
}

uint64_t kDataFrameMAP::count(string kmerS)
{
  return this->MAP[kmer::canonicalKmer(kmerS)];
}

uint64_t kDataFrameMAP::bucket(string kmerS)
{
  return this->MAP.bucket(kmer::canonicalKmer(kmerS));
}


bool kDataFrameMAP::erase(string kmerS)
{
  return this->MAP.erase(kmer::canonicalKmer(kmerS));
}

uint64_t kDataFrameMAP::size() {
    return (uint64_t)this->MAP.size();
}

uint64_t kDataFrameMAP::max_size() {
    return (uint64_t)this->MAP.max_size();
}

float kDataFrameMAP::load_factor() {
    return this->MAP.load_factor();
}

float kDataFrameMAP::max_load_factor() {
    return this->MAP.max_load_factor();
}



void kDataFrameMAP::save(string filePath) {
    ofstream myfile;
    myfile.open(filePath + ".map", ios::out);
    unordered_map<string, uint64_t>::iterator it;
    myfile << this->kSize << endl;
    for (auto const &it : this->MAP) {
        myfile << it.first << ":" << it.second << endl;
    }



}

kDataFrame *kDataFrameMAP::load(string filePath) {
    filePath += ".map";
    ifstream myfile(filePath);
    string key, value;

    getline(myfile, key, '\n'); // Get Kmer_size from first line
    kDataFrameMAP *KMAP = new kDataFrameMAP(std::stoull(key));

    while (!myfile.eof()) {
        getline(myfile, key, ':');
        if (getline(myfile, value, '\n')) {
            // cout << "Key:" << key << "| Value: " << value << endl;
            KMAP->setCount(key, std::stoull(value));
        } else break;
    }
    myfile.close();
    return KMAP;
}
kDataFrame* kDataFrameMAP::getTwin(){
  return ((kDataFrame*)new kDataFrameMAP(kSize));
}
void kDataFrameMAP::reserve(uint64_t n)
{
  this->MAP.reserve(n);
}
kDataFrameIterator kDataFrameMAP::begin() {
  return *(new kDataFrameIterator(
    (_kDataFrameIterator*)new kDataFrameMAPIterator(MAP.begin(),this,kSize),
    (kDataFrame*)this));
}
kDataFrameIterator kDataFrameMAP::end() {
  return *(new kDataFrameIterator(
    (_kDataFrameIterator*)new kDataFrameMAPIterator(MAP.end(),this,kSize),
    (kDataFrame*)this));
}
