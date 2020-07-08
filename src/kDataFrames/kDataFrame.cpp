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
    defaultColumn=NULL;
}

kDataFrame::kDataFrame(uint8_t k_size) {
    kSize = k_size;
    isStatic=false;
    isKmersOrderComputed=false;
    defaultColumn=NULL;
}

kDataFrame::~kDataFrame(){
    delete KD;
    if(defaultColumn!=NULL)
        delete defaultColumn;
    for(auto c:columns)
        delete c.second;
    delete endIterator;
}

bool kDataFrame::empty() {
    return this->size() == 0;
}

bool kDataFrame::insert(kmerRow k) {
    return this->insert(k.kmer, k.count);
}

void kDataFrame::save(string filePath)
{
    ofstream out(filePath+".multiColumn");
    out<<"isStatic=\t"<<isStatic<<endl;
    if(isStatic)
    {
        out<<"numColumns=\t"<<columns.size()<<endl;
        for(auto c:columns)
        {
	    string suffix=".multiColumn."+c.first;
            string filename=filePath+ suffix;
            size_t columnType=typeid(*(c.second)).hash_code();
            out<<c.first<<"\t"<<columnType<<"\t"<<suffix<<endl;
            c.second->serialize(filename);
        }

    }
    if(defaultColumn!=NULL)
    {
        string suffix=".defaultColumn";
        string filename=filePath+suffix;
        size_t columnType=typeid(*(defaultColumn)).hash_code();
        out<<"default\t"<<columnType<<"\t"<<suffix<<endl;
        defaultColumn->serialize(filename);
        cout<<"saving default finsihed"<<endl;
    }
    else{
        out<<"default\t"<<0<<"\tNULL"<<endl;
    }
    out.close();
    this->serialize(filePath);
}

kDataFrame * kDataFrame::load(string filePath) {
    kDataFrame* res;
    if (fileExists(filePath + ".mqf"))
        res=kDataFrameMQF::load(filePath);
    else if (fileExists(filePath + ".map"))
        res=kDataFrameMAP::load(filePath);
    else if (fileExists(filePath + ".phmap"))
        res=kDataFramePHMAP::load(filePath);
    else if (fileExists(filePath+ ".bmqf"))
        res=kDataFrameBMQF::load(filePath);
    else
        throw std::runtime_error("Could not open kDataFrame file");

    if (fileExists(filePath +".multiColumn") ){
        ifstream inp(filePath + ".multiColumn");
        string tmp;
        string name, path;
        uint64_t type;
        inp >> tmp >> res->isStatic;
        if (res->isStatic) {
            uint32_t numColumns;
            inp >> tmp >> numColumns;
            for (uint32_t i = 0; i < numColumns; i++) {
                inp >> name >> type >> path;
                Column *c = Column::getContainerByName(type);
                c->deserialize(filePath+path);
                res->columns[name] = c;
            }
            res->preprocessKmerOrder();
        }
        inp >> name >> type >> path;
        if (type != 0) {
            Column *c = Column::getContainerByName(type);
            c->deserialize(filePath+path);
            res->defaultColumn = c;
        }
    }

    return res;

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
  isKmersOrderComputed=true;
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



//template int kDataFrame::getKmerColumnValue<int, vectorColumn<int> >(string columnName,string kmer);
//template double kDataFrame::getKmerColumnValue<double, vectorColumn<double> >(string columnName,string kmer);
//template bool kDataFrame::getKmerColumnValue<bool, vectorColumn<bool> >(string columnName,string kmer);
//
//template void kDataFrame::setKmerColumnValue<int, vectorColumn<int>  >(string columnName,string kmer, int value);
//template void kDataFrame::setKmerColumnValue<double, vectorColumn<double>  >(string columnName,string kmer, double value);
//template void kDataFrame::setKmerColumnValue<bool, vectorColumn<bool>  >(string columnName,string kmer, bool value);
//
//



void kDataFrame::addColumn(string columnName,Column* ptr)
{
  if(!isKmersOrderComputed)
  {
    this->preprocessKmerOrder();
    isStatic=true;
  }
  columns[columnName]=ptr;

}




//template double kDataFrame::getKmerDefaultColumnValue<double, vectorColumn<double>  >(string kmer);
//template void kDataFrame::setKmerDefaultColumnValue<double, vectorColumn<double>>(string kmer, double value);
//
//
//
//template vector<uint32_t > kDataFrame::getKmerDefaultColumnValue<vector<uint32_t >, colorColumn>(string kmer);
//template void kDataFrame::setKmerDefaultColumnValue<vector<uint32_t >, colorColumn>(string kmer, vector<uint32_t > value);
//
//template double kDataFrame::getKmerDefaultColumnValue<double, vectorColumn<double>  >(uint64_t kmer);
//template void kDataFrame::setKmerDefaultColumnValue<double, vectorColumn<double>>(uint64_t kmer, double value);
//
//template vector<uint32_t > kDataFrame::getKmerDefaultColumnValue<vector<uint32_t >, colorColumn>(uint64_t kmer);
//template void kDataFrame::setKmerDefaultColumnValue<vector<uint32_t >, colorColumn>(uint64_t kmer, vector<uint32_t > value);
//
//
//template vector<string> kDataFrame::getKmerDefaultColumnValue<vector<string>, StringColorColumn>(string kmer);
//template vector<string> kDataFrame::getKmerDefaultColumnValue<vector<string>, StringColorColumn>(uint64_t kmer);


void kDataFrame::changeDefaultColumnType(Column* ptr)
{
    defaultColumn=ptr;
}

//template<typename T,typename Container>
//T kDataFrame::getKmerDefaultColumnValue(string kmer)
//{
//    return ((Container*)defaultColumn)->getWithIndex(getCount(kmer));
//}
//
//template<typename T,typename Container>
//void kDataFrame::setKmerDefaultColumnValue(string kmer, T value)
//{
//    uint32_t i=((Container*)defaultColumn)->insertAndGetIndex(value);
//    setCount(kmer,i);
//}
//
//template<typename T,typename Container>
//T kDataFrame::getKmerDefaultColumnValue(uint64_t kmer)
//{
//    return ((Container*)defaultColumn)->getWithIndex(getCount(kmer));
//}
//
//template<typename T,typename Container>
//void kDataFrame::setKmerDefaultColumnValue(uint64_t kmer, T value)
//{
//    uint32_t i=((Container*)defaultColumn)->insertAndGetIndex(value);
//    setCount(kmer,i);
//}
//

kDataFrameIterator kDataFrame::end(){
//    kDataFrameBMQFIterator* it=new kDataFrameBMQFIterator(bufferedmqf,kSize,KD);
//    it->endIterator();
//    return (kDataFrameIterator(it,(kDataFrame*)this));
    return *endIterator;
}