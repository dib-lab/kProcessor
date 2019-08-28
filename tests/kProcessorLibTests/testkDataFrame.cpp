#include "testkDataframe.h"
#include "Utils/kmer.h"
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <seqan/seq_io.h>
#include "algorithms.hpp"
#include <iterator>
#include <algorithm>


#include <set>
using namespace std;


INSTANTIATE_TEST_SUITE_P(testcolorsTable,
                         colorsTableTest,
                        ::testing::Combine(
                        ::testing::Values("bitVector","intVector"),
                        ::testing::Values(10,20,100),
                        ::testing::Values(10,100,1000)
                      ));
void colorsTableTest::SetUp(){
  uint64_t numSamples=get<1>(GetParam());
  uint64_t numColors=get<2>(GetParam());
  for(int i=1;i<=numColors;i++)
  {
    int n=rand()%(numSamples*2);
    set<uint32_t> colors;
    colors.clear();
    for(int j=0;j<n;j++)
    {
      colors.insert(rand()%numSamples);
    }
    vector<uint32_t> colorsVec;
    colorsVec.clear();
    copy(colors.begin(),colors.end(),back_inserter(colorsVec));
    simColors[i]=colorsVec;
  }
}
colorTable* createColorTables(string name,uint64_t numSamples,uint64_t numColors)
{
  if(name=="bitVector")
    return new BitVectorsTable(numSamples);
  if(name=="intVector")
      return new intVectorsTable();
  return NULL;
}


TEST_P(colorsTableTest,insertAndQuery)
{
  string colorTableName=get<0>(GetParam());
  uint64_t numSamples=get<1>(GetParam());
  uint64_t numColors=get<2>(GetParam());
  colorTable* table=createColorTables(colorTableName,numSamples,numColors);
  for(int i=1;i<=numColors;i++)
  {
    table->setColor(i,simColors[i]);
  }
  for(auto it:simColors)
  {
    vector<uint32_t> res;
    res.clear();
    table->getSamples(it.first,res);
    EXPECT_EQ(res,it.second);
  }
}

TEST_P(colorsTableTest,saveAndLoad)
{
  string colorTableName=get<0>(GetParam());
  uint64_t numSamples=get<1>(GetParam());
  uint64_t numColors=get<2>(GetParam());
  colorTable* table=createColorTables(colorTableName,numSamples,numColors);
  for(int i=1;i<=numColors;i++)
  {
    table->setColor(i,simColors[i]);
  }
  string fileName="colorTable.test";
  table->save(fileName);
  delete table;
  colorTable* loadedTable=colorTable::load(fileName);
  for(auto it:simColors)
  {
    vector<uint32_t> res;
    res.clear();
    loadedTable->getSamples(it.first,res);
    EXPECT_EQ(res,it.second);
  }
  delete loadedTable;
}


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
  for(auto k:kSizes)
  {
    framesToBeTested.push_back(new kDataFrameBMQF(k));
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

INSTANTIATE_TEST_SUITE_P(testntCard,
                         estimateTest,
                         ::testing::ValuesIn(fastqFiles)
                                            );

vector<vector<string> > setFunctionsTestInput={{"test.noN.fastq","test2.noN.fastq","test2.noN.fastq"}};

vector<vector<kDataFrame*> > BuildTestFramesForSetFunctions()
{
  vector<vector<kDataFrame*> > framesToBeTested(3);
  int k=31;
  for(auto file:setFunctionsTestInput[0])
  {
    framesToBeTested[0].push_back(new kDataFrameMQF(k));
//    kProcessor::parseSequences(file,1,framesToBeTested[0].back());
    framesToBeTested[1].push_back(new kDataFrameMQF(k)); // Temporary until resolving #17
    framesToBeTested[1].push_back(new kDataFrameMAP(k)); // Set functions should now work fine on sorted map
///    kProcessor::parseSequences(file,1,framesToBeTested[1].back());
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

TEST_P(kDataFrameTest,transformFilterLessThan5)
{

    kDataFrame* kframe=GetParam();
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
    for(auto k:*kmers)
    {
      kframe->insert(k.first,k.second);
      if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
        break;
      }
    }
    int checkedKmers=0;
    kDataFrame* kframe2=kProcessor::filter(kframe,[](kmerRow k)
    {
      return k.count>=5;
    });
    kDataFrameIterator it=kframe2->begin();
    while(it!=kframe2->end())
    {
      string kmer=it.getKmer();
      uint64_t count=it.getKmerCount();
      ASSERT_GE(count,5);
      it++;
    }

}

TEST_P(kDataFrameTest,aggregateSum)
{

    kDataFrame* kframe=GetParam();
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
    uint64_t goldSum=0;
    for(auto k:*kmers)
    {
      kframe->insert(k.first,k.second);
      goldSum+=k.second;
      if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
        break;
      }
    }
    int checkedKmers=0;
    any initial=(uint64_t)0;
    any sum=kProcessor::aggregate(kframe,initial,[](kmerRow k,any v)
    {
      uint64_t tmp=any_cast<uint64_t>(v);
      any result=tmp+k.count;
      return result;
    });
    ASSERT_EQ(any_cast<uint64_t>(sum),goldSum);


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

TEST_P(estimateTest,estimateTestTest)
{
  string fileName=GetParam();
  int kSize=31;
  vector<uint64_t> res= kProcessor::estimateKmersHistogram(fileName,kSize,4);
  ifstream gold((fileName+".ntCardRes").c_str());
  uint64_t g;
  int i=0;
  while(gold>>g)
  {
    ASSERT_EQ(g,res[i++]);
  }
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

vector<colorTableInv*> BuildColorTableInv()
{
  vector<colorTableInv*> tablesToBeTested;
  tablesToBeTested.push_back(new stringColorTableInv());
  return tablesToBeTested;
}

void colorsTableInvTest::SetUp(){
  for(int i=1;i<=numColors;i++)
  {
    int n=rand()%(numSamples*2)+1;
    set<uint32_t> colors;
    colors.clear();
    for(int j=0;j<n;j++)
    {
      colors.insert(rand()%numSamples);
    }
    vector<uint32_t> colorsVec;
    colorsVec.clear();
    copy(colors.begin(),colors.end(),back_inserter(colorsVec));
    simColors[i]=colorsVec;
  }
}

INSTANTIATE_TEST_SUITE_P(testcolorsTableInvTest,
                        colorsTableInvTest,
                        ::testing::ValuesIn(BuildColorTableInv()));


TEST_P(colorsTableInvTest,setAndGet)
{
    colorTableInv* table=GetParam();
    for(int i=1;i<=numColors;i++)
    {
      if(table->getColorId(simColors[i])==0)
        table->setColorId(i,simColors[i]);
      else{
        simColors.erase(i);
      }
    }
    for(auto it:simColors)
    {
      vector<uint32_t> res;
      res.clear();
      uint64_t colorId =table->getColorId(it.second);
      EXPECT_EQ(colorId,it.first);
    }

}

INSTANTIATE_TEST_SUITE_P(testIndexing,
                        indexingTest,
                        ::testing::Values("test1.fa"));


TEST_P(indexingTest,index)
{
  string filename=GetParam();
  int chunkSize = 1000;

  kDataFrame *KF = new kDataFrameMQF(25, 28, 1);
  kmerDecoder *KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});
  colored_kDataFrame* res= kProcessor::index(KMERS, filename+".names", KF);

  seqan::SeqFileIn seqIn2(filename.c_str());
  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::CharString> reads;
  uint64_t kSize=res->getkSize();
  int readCount=0;
  int correct=0,wrong=0;
  vector<uint32_t> colors;
  string kmer;
  while(!atEnd(seqIn2)){
    clear(reads);
    clear(ids);
    seqan::readRecords(ids, reads, seqIn2,chunkSize);
    for(int j=0;j<length(reads);j++){
      string readName=string((char*)seqan::toCString(ids[j]));
      uint32_t sampleID=res->namesMapInv[readName];
      ASSERT_NE(sampleID,0);


      string seq=string((char*)seqan::toCString(reads[j]));
      if(seq.size()<kSize)
        continue;
      for (int i = 0; i < seq.size() - kSize + 1; i++)
      {
        kmer=seq.substr(i,kSize);
        colors.clear();
        res->getSamplesIDForKmer(kmer,colors);
        ASSERT_NE(colors.size(),0);
        auto colorIt=colors.end();
        colorIt=find(colors.begin(),colors.end(),sampleID);
        ASSERT_NE(colorIt,colors.end());

      }
    }

  }

}
TEST_P(indexingTest,saveAndLoad)
{
  string filename=GetParam();
  int chunkSize = 1000;
  kDataFrame *KF = new kDataFrameMQF(25, 28, 1);
  kmerDecoder *KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});
  colored_kDataFrame* res1= kProcessor::index(KMERS,filename+".names", KF);
  res1->save("tmp.coloredKdataFrame");
  colored_kDataFrame* res=colored_kDataFrame::load("tmp.coloredKdataFrame");
  seqan::SeqFileIn seqIn2(filename.c_str());
  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::CharString> reads;
  uint64_t kSize=res->getkSize();
  int readCount=0;
  int correct=0,wrong=0;
  vector<uint32_t> colors;
  string kmer;
  while(!atEnd(seqIn2)){
    clear(reads);
    clear(ids);
    seqan::readRecords(ids, reads, seqIn2,chunkSize);
    for(int j=0;j<length(reads);j++){
      string readName=string((char*)seqan::toCString(ids[j]));
      uint32_t sampleID=res->namesMapInv[readName];
      ASSERT_NE(sampleID,0);


      string seq=string((char*)seqan::toCString(reads[j]));
      if(seq.size()<kSize)
        continue;
      for (int i = 0; i < seq.size() - kSize + 1; i++)
      {
        kmer=seq.substr(i,kSize);
        colors.clear();
        res->getSamplesIDForKmer(kmer,colors);
        ASSERT_NE(colors.size(),0);
        auto colorIt=colors.end();
        colorIt=find(colors.begin(),colors.end(),sampleID);
        ASSERT_NE(colorIt,colors.end());

      }
    }

  }

}
