#include "testkDataframe.h"
#include "Utils/kmer.h"
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <seqan/seq_io.h>
#include "algorithms.hpp"
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


vector<string> fastqFiles={"test.noN.fastq","test2.noN.fastq","test2.noN.fastq"};
INSTANTIATE_TEST_SUITE_P(testcounting,
                         algorithmsTest,
                        ::testing::Combine(
                        ::testing::ValuesIn(BuildTestFrames()),
                        ::testing::ValuesIn(fastqFiles)
                      ));

vector<vector<string> > setFunctionsTestInput={{"test.noN.fastq","test2.noN.fastq","test2.noN.fastq"}};

vector<vector<kDataFrame*> > BuildTestFramesForSetFunctions()
{
  vector<vector<kDataFrame*> > framesToBeTested(2);
  int k=31;
  for(auto file:setFunctionsTestInput[0])
  {
    framesToBeTested[0].push_back(new kDataFrameMQF(k));
    kProcessor::parseSequences(file,1,framesToBeTested[0].back());
    framesToBeTested[1].push_back(new kDataFrameMAP(k));
    kProcessor::parseSequences(file,1,framesToBeTested[1].back());
  }
  return framesToBeTested;
}
INSTANTIATE_TEST_SUITE_P(testSetFunctions,
                         setFunctionsTest,
                        ::testing::ValuesIn(BuildTestFramesForSetFunctions())
                      );


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

string gen_random(const int len) {
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    string s="";
    for (int i = 0; i < len; ++i) {
        s+= alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    return s;
}

TEST_P(kDataFrameTest,saveAndIterateOverAllKmers)
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
    string fileName="tmp.kdataframe."+gen_random(4);
    kframe->save(fileName);
    kDataFrame* kframeLoaded=kDataFrame::load(fileName);
    int checkedKmers=0;
    kDataFrameIterator it=kframeLoaded->begin();
    while(it!=kframeLoaded->end())
    {
      string kmer=it.getKmer();
      uint64_t count=it.getKmerCount();
      ASSERT_EQ(count,insertedKmers[kmer]);
      insertedKmers.erase(kmer);
      it++;
    }
    EXPECT_EQ(insertedKmers.size(),0);

}


TEST_P(kDataFrameTest,transformPlus10)
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
    kDataFrame* kframe2=kProcessor::transform(kframe,[](kmerRow k)
    {
      k.count+=10;
      return k;
    });
    kDataFrameIterator it=kframe2->begin();
    while(it!=kframe2->end())
    {
      string kmer=it.getKmer();
      uint64_t count=it.getKmerCount();
      ASSERT_EQ(count,insertedKmers[kmer]+10);
      insertedKmers.erase(kmer);
      it++;
    }
    EXPECT_EQ(insertedKmers.size(),0);

}


TEST_P(algorithmsTest,parsingTest)
{
  kDataFrame* kframe=get<0>(GetParam());
  string fileName=get<1>(GetParam());
  int kSize=kframe->getkSize();
  kProcessor::parseSequences(fileName,1,kframe);
  seqan::SeqFileIn seqIn(fileName.c_str());
  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::CharString> reads;
  int chunkSize=1000;
  while(!atEnd(seqIn)){
    clear(reads);
    clear(ids);

    seqan::readRecords(ids, reads, seqIn,chunkSize);
    for(int j=0;j<length(reads);j++)
    {
      string seq=string((char*)seqan::toCString(reads[j]));
      for(int i=0;i<seq.size()-kSize+1;i++)
      {
          string kmer=seq.substr(i,kSize);
          ASSERT_GE(kframe->count(kmer),1);
      }
    }

  }
  seqan::close(seqIn);
}

TEST_P(setFunctionsTest,unioinTest)
{
  vector<kDataFrame*> input=GetParam();
  kDataFrame* unioinResult=kProcessor::kFrameUnion(input);
  for(auto kframe:input)
  {
    auto it=kframe->begin();
    while(it!=kframe->end())
    {
      int count=unioinResult->count((*it).kmer);
      ASSERT_GE(count,((*it).count));
      it++;
    }
  }
  delete unioinResult;

}
TEST_P(setFunctionsTest,intersectTest)
{
  vector<kDataFrame*> input=GetParam();

  kDataFrame* intersectResult=kProcessor::kFrameIntersect(input);
  auto it=intersectResult->begin();
  while(it!=intersectResult->end())
  {
    for(auto kframe:input)
    {
      int count=kframe->count((*it).kmer);
      ASSERT_GE(count,((*it).count));
    }
    it++;
  }
  delete intersectResult;

}

TEST_P(setFunctionsTest,differenceTest)
{
  vector<kDataFrame*> input=GetParam();

  kDataFrame* diffResult=kProcessor::kFrameDiff(input);
  auto it=diffResult->begin();
  while(it!=diffResult->end())
  {
      int count=input[0]->count((*it).kmer);
      ASSERT_EQ(count,((*it).count));
      for(int i=1;i<input.size();i++)
      {
        int count=input[i]->count((*it).kmer);
        ASSERT_EQ(count,0);
      }
      it++;
  }
  delete diffResult;


}
