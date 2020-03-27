#include "kDataFrame.hpp"
#include "Utils/kmer.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <limits>
#include <sstream>

using namespace std;

inline bool fileExists(const std::string &name) {
    ifstream f(name.c_str());
    return f.good();
}

kDataFrame::kDataFrame() {
    kSize = 31;
    isStatic=false;
    isKmersOrderComputed=false;
}

kDataFrame::kDataFrame(uint8_t k_size) {
    kSize = k_size;
    isStatic=false;
    isKmersOrderComputed=false;

}

bool kDataFrame::empty() {
    return this->size() == 0;
}

bool kDataFrame::insert(kmerRow k) {
    return this->insert(k.kmer, k.count);
}


kDataFrame * kDataFrame::load(string filePath) {
    if (fileExists(filePath + ".mqf"))
        return kDataFrameMQF::load(filePath);
    else if (fileExists(filePath + ".map"))
        return kDataFrameMAP::load(filePath);
    else if (fileExists(filePath + ".phmap"))
        return kDataFramePHMAP::load(filePath);
    else if (fileExists(filePath))
        return kDataFrameBMQF::load(filePath);
    else
        throw std::runtime_error("Could not open kDataFrame file");
}

void kDataFrame::preprocessKmerOrder()
{
  int prevOrder=0;
  int checkpointsDistance=64;
  int index=0;
  kDataFrameIterator it=this->begin();
  while(it!=this->end())
  {
    string kmer=it.getKmer();
    if(index%checkpointsDistance==0)
    {
      orderCheckpoints[kmer]=prevOrder;
      prevOrder=index;
    }
    index++;
    it++;
  }
  orderCheckpoints["THEEND"]=prevOrder;
}
uint64_t kDataFrame::getkmerOrder(string kmer)
{
  kDataFrameIterator it=this->find(kmer);
  kmer=it.getKmer();
  uint32_t offset=0;
  while(it!=this->end()&&
  orderCheckpoints.find(kmer) == orderCheckpoints.end())
  {
    offset++;
    it++;
    kmer=it.getKmer();
  }
  if(it==this->end())
  {
    kmer="THEEND";
  }
  return orderCheckpoints[kmer]+offset;
}


template void kDataFrame::addColumn<int>(string columnName);
template void kDataFrame::addColumn<double>(string columnName);
template void kDataFrame::addColumn<bool>(string columnName);

template int kDataFrame::getKmerColumnValue<int>(string columnName,string kmer);
template double kDataFrame::getKmerColumnValue<double>(string columnName,string kmer);
template bool kDataFrame::getKmerColumnValue<bool>(string columnName,string kmer);

template void kDataFrame::setKmerColumnValue<int>(string columnName,string kmer, int value);
template void kDataFrame::setKmerColumnValue<double>(string columnName,string kmer, double value);
template void kDataFrame::setKmerColumnValue<bool>(string columnName,string kmer, bool value);

template<typename T>
void kDataFrame::addColumn(string columnName)
{
  if(!isKmersOrderComputed)
  {
    this->preprocessKmerOrder();
    isKmersOrderComputed=true;
    isStatic=true;
  }
  columns[columnName]=new vector<T>(this->size());

}



template<typename T>
T kDataFrame::getKmerColumnValue(string columnName,string kmer)
{
  uint64_t kmerOrder=getkmerOrder(kmer);
  vector<T>* col=any_cast<vector<T>* >(columns[columnName]);
  return (*col)[kmerOrder];
}
template<typename T>
void kDataFrame::setKmerColumnValue(string columnName,string kmer,T value)
{
  uint64_t kmerOrder=getkmerOrder(kmer);
  vector<T>* col=any_cast<vector<T>* >(columns[columnName]);
  (*col)[kmerOrder]=value;
}
