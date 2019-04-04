#include "gtest/gtest.h"
#include "kDataFrame.hpp"
#include <unordered_map>
#include "colorTable.hpp"
using namespace std;



class kDataFrameTest : public ::testing::TestWithParam<kDataFrame* >{
};

class algorithmsTest : public ::testing::TestWithParam<tuple<kDataFrame*,string> >{
};

//string is the name of color table class
//and the first integer is the number of samples and the second one is the number of colors
class colorsTableTest : public ::testing::TestWithParam<tuple<string, uint64_t,uint64_t> >{
  public:
    unordered_map<uint64_t,vector<uint32_t> > simColors;
    virtual void SetUp();
};


class setFunctionsTest : public ::testing::TestWithParam<vector<kDataFrame*>  >{
};

class kmersGenerator{
private:
  unordered_map<int,unordered_map<string,int>* > db;
public:
  kmersGenerator(){}
  unordered_map<string,int>* getKmers(int kSize)
  {
    auto it=db.find(kSize);
    if(it==db.end())
    {
      srand (time(NULL));
      db[kSize]=new unordered_map<string,int>();
      size_t nKmers=100000;
      uint64_t range=(1ULL<<kSize);
      while(db[kSize]->size() < nKmers)
      {
        uint64_t kmerInt=rand()%range;
        string kmerStr=kmer::int_to_str(kmerInt,kSize);
        uint64_t kmerCount=(rand()%1000)+1;
        db[kSize]->insert(make_pair(kmerStr,kmerInt));
      }
      it=db.find(kSize);
    }
    return it->second;
  }
};
