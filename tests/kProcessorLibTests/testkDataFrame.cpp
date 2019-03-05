#include "testkDataframe.h"
#include "KmerCounter/kmer.h"
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <iostream>
using namespace std;
//add the new kDataframes here to be tested
using MyTypes = ::testing::Types<kDataFrameMQF,kDataFrameMAP>;
TYPED_TEST_SUITE(kDataFrameTest, MyTypes);

template<class T_kDataFrame>
unordered_map<string,int> kDataFrameTest<T_kDataFrame>::kmers;

template<class T_kDataFrame>
void kDataFrameTest<T_kDataFrame>::SetUp() {
  if(kmers.size()==0)
  {
    srand (time(NULL));
    int kSize=kFrame.getkSize();
    int nKmers=100000;
    uint64_t range=(1ULL<<kSize);
    while(kmers.size() < nKmers)
    {
      uint64_t kmerInt=rand()%range;
      string kmerStr=kmer::int_to_str(kmerInt,kSize);
      uint64_t kmerCount=(rand()%1000)+1;
      kmers[kmerStr]=kmerInt;
    }
  }
}


TYPED_TEST(kDataFrameTest,emptykDataFrame)
{

    kDataFrame* kframe=&this->kFrame;
    EXPECT_EQ(kframe->empty(), true);
    kframe->insert(this->kmers.begin()->first);
    EXPECT_EQ(kframe->empty(),false);
}

TYPED_TEST(kDataFrameTest,insertOneTime)
{

    kDataFrame* kframe=&this->kFrame;
    EXPECT_EQ(kframe->empty(), true);
    int insertedKmers=0;
    for(auto k:this->kmers)
    {
        kframe->insert(k.first);
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8)
        {
          break;
        }
        insertedKmers++;
    }
    int checkedKmers=0;
    for(auto k:this->kmers)
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

TYPED_TEST(kDataFrameTest,insertNTimes)
{

    kDataFrame* kframe=&this->kFrame;
    EXPECT_EQ(kframe->empty(), true);
    int insertedKmers=0;
    for(auto k:this->kmers)
    {
        kframe->insert(k.first,k.second);
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8)
        {
          break;
        }
        insertedKmers++;
    }
    int checkedKmers=0;
    for(auto k:this->kmers)
    {
        int c=kframe->count(k.first);
        EXPECT_GE(c,k.second);
        if(checkedKmers==insertedKmers)
        {
          break;
        }
        checkedKmers++;
    }
    cout<<kframe->size()<<endl;
    cout<<kframe->load_factor()<<" "<<kframe->max_load_factor()<<endl;

}
TYPED_TEST(kDataFrameTest,eraseKmers)
{

    kDataFrame* kframe=&this->kFrame;
    EXPECT_EQ(kframe->empty(), true);
    int insertedKmers=0;
    for(auto k:this->kmers)
    {
        kframe->insert(k.first);
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8)
        {
          break;
        }
        insertedKmers++;
    }
    int checkedKmers=0;
    for(auto k:this->kmers)
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
    for(auto k:this->kmers)
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
