#include "gtest/gtest.h"
#include "kDataFrame.hpp"
#include <unordered_map>

using namespace std;



class kDataFrameTest : public ::testing::TestWithParam<kDataFrame* >{


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
