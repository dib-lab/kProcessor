#include "testkDataframe.h"
#include "Utils/kmer.h"
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <iostream>
#include <vector>
using namespace std;
//add the new kDataframes here to be tested
// using MyTypes = ::testing::Types<kDataFrameMQF,kDataFrameMAP>;
// TEST_P_SUITE(kDataFrameTest, MyTypes);
//static vector<kDataFrame*> framesToBeTested;

vector<kDataFrame*> BuildTestFrames()
{
  vector<kDataFrame*> framesToBeTested;
  vector<int> kSizes={21,31};
  for(auto k:kSizes)
  {
    framesToBeTested.push_back(new kDataFrameMQF(k));
  }
  for(auto k:kSizes)
  {
    framesToBeTested.push_back(new kDataFrameMAP(k));
  }
  return framesToBeTested;
}

INSTANTIATE_TEST_SUITE_P(testFrames,
                        kDataFrameTest,
                        ::testing::ValuesIn(BuildTestFrames()));

static  kmersGenerator kmersGen;


TEST_P(kDataFrameTest,emptykDataFrame)
{

    kDataFrame* kframe=GetParam();
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
    kframe->insert(kmers->begin()->first);
    EXPECT_EQ(kframe->empty(),false);
}

TEST_P(kDataFrameTest,insertOneTime)
{

    kDataFrame* kframe=GetParam();
    EXPECT_EQ(kframe->empty(), true);
    int insertedKmers=0;
    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
    for(auto k:*kmers)
    {
        kframe->insert(k.first);
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8)
        {
          break;
        }
        insertedKmers++;
    }
    int checkedKmers=0;
    for(auto k:*kmers)
    {
        int c=kframe->count(k.first);
        EXPECT_GE(c,1);
        if(checkedKmers==insertedKmers)
        {
          break;
        }
        checkedKmers++;
    }

}

TEST_P(kDataFrameTest,insertNTimes)
{

    kDataFrame* kframe=GetParam();
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
    int insertedKmers=0;
    for(auto k:*kmers)
    {
        kframe->insert(k.first,k.second);
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8)
        {
          break;
        }
        insertedKmers++;
    }
    int checkedKmers=0;
    for(auto k:*kmers)
    {
        int c=kframe->count(k.first);
        EXPECT_GE(c,k.second);
        if(checkedKmers==insertedKmers)
        {
          break;
        }
        checkedKmers++;
    }
}
TEST_P(kDataFrameTest,eraseKmers)
{

    kDataFrame* kframe=GetParam();
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
    int insertedKmers=0;
    for(auto k:*kmers)
    {
        kframe->insert(k.first);
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8)
        {
          break;
        }
        insertedKmers++;
    }
    int checkedKmers=0;
    for(auto k:*kmers)
    {
        if(checkedKmers%2==0)
          kframe->erase(k.first);

        if(checkedKmers==insertedKmers)
        {
          break;
        }
        checkedKmers++;

    }

    checkedKmers=0;
    for(auto k:*kmers)
    {
        int c=kframe->count(k.first);
        if(checkedKmers%2==0)
        {
          ASSERT_EQ(c,0);
        }
        else{
          ASSERT_GE(c,1);
        }

        if(checkedKmers==insertedKmers)
        {
          break;
        }
        checkedKmers++;
    }

}


TEST_P(kDataFrameTest,autoResize)
{

    kDataFrame* kframe=GetParam();
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
    int insertedKmers=0;
    for(auto k:*kmers)
    {
        kframe->insert(k.first,k.second);
        insertedKmers++;
    }
    int checkedKmers=0;
    for(auto k:*kmers)
    {
        int c=kframe->count(k.first);
        EXPECT_GE(c,k.second);
        if(checkedKmers==insertedKmers)
        {
          break;
        }
        checkedKmers++;
    }
}

TEST_P(kDataFrameTest,iterateOverAllKmers)
{

    kDataFrame* kframe=GetParam();
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
    unordered_map<string,int> insertedKmers;
    for(auto k:*kmers)
    {
      kframe->insert(k.first,k.second);
      insertedKmers[k.first]+=k.second;
      if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
        break;
      }
    }
    int checkedKmers=0;
    kDataFrameIterator it=kframe->begin();
    while(it!=kframe->end())
    {
      string kmer=it.getKmer();
      uint64_t count=it.getKmerCount();
      ASSERT_EQ(count,insertedKmers[kmer]);
      insertedKmers.erase(kmer);
      it++;
    }
    EXPECT_EQ(insertedKmers.size(),0);

}
