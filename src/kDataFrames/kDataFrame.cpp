#include "kDataFrame.hpp"
#include "Utils/kmer.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <limits>
#include <sstream>
#include "defaultColumn.hpp"
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


template void kDataFrame::addColumn<vector<int>  >(string columnName,vector<int>* ptr);
template void kDataFrame::addColumn<vector<double> >(string columnName,vector<double>* ptr);
template void kDataFrame::addColumn<vector<bool> >(string columnName,vector<bool>* ptr);

template int kDataFrame::getKmerColumnValue<int,vector<int> >(string columnName,string kmer);
template double kDataFrame::getKmerColumnValue<double,vector<double> >(string columnName,string kmer);
template bool kDataFrame::getKmerColumnValue<bool,vector<bool> >(string columnName,string kmer);

template void kDataFrame::setKmerColumnValue<int,vector<int> >(string columnName,string kmer, int value);
template void kDataFrame::setKmerColumnValue<double, vector<double> >(string columnName,string kmer, double value);
template void kDataFrame::setKmerColumnValue<bool, vector<bool> >(string columnName,string kmer, bool value);

template<typename T>
void kDataFrame::addColumn(string columnName,T* ptr)
{
  if(!isKmersOrderComputed)
  {
    this->preprocessKmerOrder();
    isKmersOrderComputed=true;
    isStatic=true;
  }
  columns[columnName]=ptr;

}



template<typename T, typename Container>
T kDataFrame::getKmerColumnValue(string columnName,string kmer)
{
  uint64_t kmerOrder=getkmerOrder(kmer);
  Container* col=any_cast<Container* >(columns[columnName]);
  return (*col)[kmerOrder];
}
template<typename T, typename Container>
void kDataFrame::setKmerColumnValue(string columnName,string kmer,T value)
{
  uint64_t kmerOrder=getkmerOrder(kmer);
  Container* col=any_cast<Container* >(columns[columnName]);
  (*col)[kmerOrder]=value;
}

template void kDataFrame::changeDefaultColumnType<vectorColumn<double> >(vectorColumn<double>* ptr);
template double kDataFrame::getKmerDefaultColumnValue<double,vectorColumn<double> >(string kmer);
template void kDataFrame::setKmerDefaultColumnValue<double,vectorColumn<double> >(string kmer, double value);


template<typename T>
void kDataFrame::changeDefaultColumnType(T* ptr)
{
    defaultColumn=ptr;
}

template<typename T, typename Container>
T kDataFrame::getKmerDefaultColumnValue(string kmer)
{
    Container* col=any_cast<Container* >(defaultColumn);
    return col->getWithIndex(getCount(kmer));
}

template<typename T, typename Container>
void kDataFrame::setKmerDefaultColumnValue(string kmer, T value)
{
    Container* col=any_cast<Container* >(defaultColumn);
    uint32_t i=col->insertAndGetIndex(value);
    setCount(kmer,i);
}
