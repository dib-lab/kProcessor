#include "testkDataframe.h"
#include "Utils/kmer.h"
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include "algorithms.hpp"
#include <iterator>
#include <algorithm>
#include <tuple>


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
  for(uint32_t i=1;i<=numColors;i++)
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
  return framesToBeTested;
}

vector<kDataFrameBMQF*> BuildTestBufferedFrames()
{
    vector<kDataFrameBMQF*> framesToBeTested;
    vector<int> kSizes={21,31};
    for(auto k:kSizes)
    {
        framesToBeTested.push_back(new kDataFrameBMQF((uint64_t)k));
    }
    return framesToBeTested;
}

INSTANTIATE_TEST_SUITE_P(testFrames,
                        kDataFrameTest,
                        ::testing::ValuesIn(BuildTestFrames()));


INSTANTIATE_TEST_SUITE_P(testFrames,
                         kDataFrameBufferedTest,
                         ::testing::ValuesIn(BuildTestBufferedFrames()));

vector<string> fastqFiles={"test.noN.fastq","test2.noN.fastq"};
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

vector<vector<string> > setFunctionsTestInput={{"test.noN.fastq","test2.noN.fastq"}};

vector<vector<kDataFrame*> > BuildTestFramesForSetFunctions()
{
  vector<vector<kDataFrame*> > framesToBeTested(2);
  int k=31;
  for(auto file:setFunctionsTestInput[0])
  {
    framesToBeTested[0].push_back(new kDataFrameMQF(k));
//    kProcessor::countKmersFromFile(framesToBeTested[0].back(), {{"mode", 1}}, file, 1000); // Mode 1 : kmers, KmerSize will be cloned from the kFrame

    framesToBeTested[1].push_back(new kDataFrameMQF(k)); // Temporary until resolving #17
//    kProcessor::countKmersFromFile(framesToBeTested[1].back(), {{"mode", 1}}, file, 1000); // Mode 1 : kmers, KmerSize will be cloned from the kFrame
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

    kDataFrame* kframe=GetParam()->getTwin();
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
    kframe->insert(kmers->begin()->first);
    EXPECT_EQ(kframe->empty(),false);
    delete kframe;
}

TEST_P(kDataFrameTest,insertOneTime)
{

    kDataFrame* kframe=GetParam()->getTwin();
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
        int c=kframe->getCount(k.first);
        EXPECT_GE(c,1);
        if(checkedKmers==insertedKmers)
        {
          break;
        }
        checkedKmers++;
    }
    delete kframe;

}

TEST_P(kDataFrameTest,insertNTimes)
{

    kDataFrame* kframe=GetParam()->getTwin();
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
        int c=kframe->getCount(k.first);
        EXPECT_GE(c,k.second);
        if(checkedKmers==insertedKmers)
        {
          break;
        }
        checkedKmers++;
    }
    delete kframe;
}
TEST_P(kDataFrameTest,eraseKmers)
{

    kDataFrame* kframe=GetParam()->getTwin();
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
        int c=kframe->getCount(k.first);
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
    delete kframe;

}


TEST_P(kDataFrameTest,autoResize)
{

    kDataFrame* kframe=GetParam()->getTwin();
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
        int c=kframe->getCount(k.first);
        EXPECT_GE(c,k.second);
        if(checkedKmers==insertedKmers)
        {
          break;
        }
        checkedKmers++;
    }
    delete kframe;
}

TEST_P(kDataFrameTest,iterateOverAllKmers)
{

    kDataFrame* kframe=GetParam()->getTwin();
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
    unordered_map<string,int> insertedKmers;
    insertedKmers.clear();
    bool first=true;
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
      uint64_t count=it.getCount();
      ASSERT_EQ(count,insertedKmers[kmer]);
      insertedKmers.erase(kmer);
      it++;
    }
    EXPECT_EQ(insertedKmers.size(),0);
    delete kframe;

}

TEST_P(kDataFrameTest,multiColumns)
{

    kDataFrame* kframe=GetParam()->getTwin();
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
    kframe->addColumn<int>("intColumn");
    kframe->addColumn<double>("doubleColumn");
    kframe->addColumn<bool>("boolColumn");

    map<string,tuple<int,double,bool> > simColumns;
    kDataFrameIterator it=kframe->begin();
    while(it!=kframe->end())
    {
      string kmer=it.getKmer();

      int randInt=rand()%1000000;
      double randDouble=(double)(rand()%1000000);
      bool randBool=rand()%2==0;

      simColumns[kmer]=make_tuple(randInt,randDouble,randBool);

      kframe->setKmerColumnValue<int>("intColumn",kmer,randInt);
      kframe->setKmerColumnValue<double>("doubleColumn",kmer,randDouble);
      kframe->setKmerColumnValue<bool>("boolColumn",kmer,randBool);
      it++;
    }
    for(auto simRow:simColumns)
    {
      string kmer=simRow.first;
      int randInt=get<0>(simRow.second);
      double randDouble=get<1>(simRow.second);
      bool randBool=get<2>(simRow.second);

      int retInt=kframe->getKmerColumnValue<int>("intColumn",kmer);
      double retDouble=kframe->getKmerColumnValue<double>("doubleColumn",kmer);
      bool retBool=kframe->getKmerColumnValue<bool>("boolColumn",kmer);

      ASSERT_EQ(randInt,retInt);
      ASSERT_EQ(randDouble,retDouble);
      ASSERT_EQ(randBool,retBool);

    }
    delete kframe;

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

    kDataFrame* kframe=GetParam()->getTwin();
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
//    int checkedKmers=0;
    kDataFrameIterator it=kframeLoaded->begin();
    while(it!=kframeLoaded->end())
    {
      string kmer=it.getKmer();
      uint64_t count=it.getCount();
      ASSERT_EQ(count,insertedKmers[kmer]);
      insertedKmers.erase(kmer);
      it++;
    }
    EXPECT_EQ(insertedKmers.size(),0);
    delete kframe;

}


TEST_P(kDataFrameTest,transformPlus10)
{

    kDataFrame* kframe=GetParam()->getTwin();
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
      uint64_t count=it.getCount();
      ASSERT_EQ(count,insertedKmers[kmer]+10);
      insertedKmers.erase(kmer);
      it++;
    }
    EXPECT_EQ(insertedKmers.size(),0);
    delete kframe;

}

TEST_P(kDataFrameTest,transformFilterLessThan5)
{

    kDataFrame* kframe=GetParam()->getTwin();
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
      uint64_t count=it.getCount();
      ASSERT_GE(count,5);
      it++;
    }
    delete kframe;

}

TEST_P(kDataFrameTest,aggregateSum)
{

    kDataFrame* kframe=GetParam()->getTwin();
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
    delete kframe;

}

TEST_P(algorithmsTest,parsingTest)
{
  kDataFrame* kframe=get<0>(GetParam())->getTwin();
  string fileName=get<1>(GetParam());
  int kSize=kframe->getkSize();
  int chunkSize=1000;

  kProcessor::countKmersFromFile(kframe, {{"mode", 1}}, fileName, 1000); // Mode 1 : kmers, KmerSize will be cloned from the kFrame
  kmerDecoder *KMERS = kProcessor::initialize_kmerDecoder(fileName, 1000, "kmers", {{"k_size", kSize}});

    while (!KMERS->end()) {
        KMERS->next_chunk();
        for (const auto &seq : *KMERS->getKmers()) {
            for (const auto &kmer : seq.second) {
                ASSERT_GE(kframe->getCount(kmer.str),1);
            }
        }
    }

    delete kframe;
    
    
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
    kDataFrame *unioinResult;

  try {
      unioinResult = kProcessor::kFrameUnion(input);
  }
  catch (const logic_error& expected)
  {
      SUCCEED();
      //FAIL();
  }
  for(auto kframe:input)
  {
    auto it=kframe->begin();
    while(it!=kframe->end())
    {
      int count=unioinResult->getCount((*it).kmer);
      ASSERT_GE(count,((*it).count));
      it++;
    }
  }
  delete unioinResult;


}
TEST_P(setFunctionsTest,intersectTest)
{
  vector<kDataFrame*> input=GetParam();

  kDataFrame* intersectResult;
    try {
        intersectResult=kProcessor::kFrameIntersect(input);
    }
    catch (const logic_error& expected)
    {
        SUCCEED();
        //FAIL();
    }

  auto it=intersectResult->begin();
  while(it!=intersectResult->end())
  {
    for(auto kframe:input)
    {
      int count=kframe->getCount((*it).kmer);
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
  try {
      kDataFrame* diffResult=kProcessor::kFrameDiff(input);
  }
  catch (const logic_error& expected)
  {
        SUCCEED();

  }
  auto it=diffResult->begin();
  while(it!=diffResult->end())
  {
      int count=input[0]->getCount((*it).kmer);
      ASSERT_EQ(count,((*it).count));
      for(int i=1;i<input.size();i++)
      {
        int count=input[i]->getCount((*it).kmer);
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

    uint64_t kSize=res->getkSize();
    int readCount=0;
    int correct=0,wrong=0;
    vector<uint32_t> colors;


    while (!KMERS->end()) {
            KMERS->next_chunk();
            for (const auto &seq : *KMERS->getKmers()) {
                string readName = seq.first;
                uint32_t sampleID = res->namesMapInv[readName];
                ASSERT_NE(sampleID, 0);
                for (const auto &kmer : seq.second) {
                    colors.clear();
                    res->getKmerSource(kmer.str, colors);
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

  uint64_t kSize=res->getkSize();
  int readCount=0;
  int correct=0,wrong=0;
  vector<uint32_t> colors;

    while (!KMERS->end()) {
        KMERS->next_chunk();
        for (const auto &seq : *KMERS->getKmers()) {
            string readName = seq.first;
            uint32_t sampleID = res->namesMapInv[readName];
            ASSERT_NE(sampleID, 0);
            for (const auto &kmer : seq.second) {
                colors.clear();
                res->getKmerSource(kmer.str, colors);
                ASSERT_NE(colors.size(),0);
                auto colorIt=colors.end();
                colorIt=find(colors.begin(),colors.end(),sampleID);
                ASSERT_NE(colorIt,colors.end());
            }
        }
    }

}


TEST_P(kDataFrameBufferedTest,iterateOverAllKmers)
{

    kDataFrameBMQF* kframe=(kDataFrameBMQF*)GetParam();
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
    unordered_map<string,int> insertedKmers;
    insertedKmers.clear();
    bool first=true;
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
        uint64_t count=it.getCount();
        ASSERT_EQ(count,insertedKmers[kmer]);
        insertedKmers.erase(kmer);
        it++;
    }
    EXPECT_EQ(insertedKmers.size(),0);
    //delete kframe;

}

TEST_P(kDataFrameBufferedTest,autoResize)
{

    kDataFrame* kframe=(kDataFrameBMQF*)GetParam();
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
    unordered_map<string,int> insertedKmers;

    for(auto k:*kmers)
    {
        kframe->insert(k.first,k.second);
        insertedKmers[k.first]+=k.second;


    }
    kDataFrameIterator it=kframe->begin();
    while(it!=kframe->end())
    {
        string kmer=it.getKmer();
        uint64_t count=it.getCount();
        ASSERT_EQ(count,insertedKmers[kmer]);
        insertedKmers.erase(kmer);
        it++;
    }
    EXPECT_EQ(insertedKmers.size(),0);
    //delete kframe;
}


//TEST_P(kDataFrameBufferedTest,transformPlus10)
//{
//
//    kDataFrame* kframe=GetParam()->getTwin();
//    EXPECT_EQ(kframe->empty(), true);
//    unordered_map<string,int>* kmers=kmersGen.getKmers((int)kframe->getkSize());
//    unordered_map<string,int> insertedKmers;
//    for(auto k:*kmers)
//    {
//        kframe->insert(k.first,k.second);
//        insertedKmers[k.first]+=k.second;
//        if(kframe->load_factor()>=kframe->max_load_factor()*0.3){
//            break;
//        }
//    }
//    int checkedKmers=0;
//    kDataFrame* kframe2=kProcessor::transform(kframe,[](kmerRow k)
//    {
//        k.count+=10;
//        return k;
//    });
//    kDataFrameIterator it=kframe2->begin();
//    while(it!=kframe2->end())
//    {
//        string kmer=it.getKmer();
//        uint64_t count=it.getCount();
//        ASSERT_EQ(count,insertedKmers[kmer]+10);
//        insertedKmers.erase(kmer);
//        it++;
//    }
//    EXPECT_EQ(insertedKmers.size(),0);
//
//}

TEST_P(kDataFrameBufferedTest,parsingTest)
{
    kDataFrame* kframe=GetParam();
    string fileName="test2.noN.fastq";
//kDataFrame* kframe=get<0>(GetParam())->getTwin();
//string fileName=get<1>(GetParam());
    int kSize=kframe->getkSize();
    kProcessor::countKmersFromFile(kframe, {{"mode", 1}}, fileName, 1000); // Mode 1 : kmers, KmerSize will be cloned from the kFrame

    seqan::SeqFileIn seqIn(fileName.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::CharString> reads;
    int chunkSize=1000;
    unordered_map<string,uint64_t > insertedKmers;
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
                kmer=kmer::canonicalKmer(kmer);
                insertedKmers[kmer]++;
            }
        }

    }
    seqan::close(seqIn);
    kDataFrameIterator it=kframe->begin();
    while(it!=kframe->end())
    {
        string kmer=it.getKmer();
        uint64_t count=it.getCount();
        if(count != insertedKmers[kmer])
        {
            cout<<kmer<<endl;

        }
        EXPECT_EQ(count,insertedKmers[kmer]);
        insertedKmers.erase(kmer);
        it++;
    }
    EXPECT_EQ(insertedKmers.size(),0);
   // delete kframe;
}

TEST_P(algorithmsTest,parsingTest2)
{

    kDataFrame* kframe=get<0>(GetParam())->getTwin();
    string fileName=get<1>(GetParam());
    int kSize=kframe->getkSize();
    RecordProperty("kdataFrame Type", kframe->get_class_name());
    int chunkSize=1000;
    unordered_map<string,uint64_t > insertedKmers;
    kProcessor::countKmersFromFile(kframe, {{"mode", 1}}, fileName, 1000); // Mode 1 : kmers, KmerSize will be cloned from the kFrame

    seqan::SeqFileIn seqIn(fileName.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::CharString> reads;


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
                kmer=kmer::canonicalKmer(kmer);
                insertedKmers[kmer]++;
            }
        }

    }
    seqan::close(seqIn);
    kDataFrameIterator it=kframe->begin();
    while(it!=kframe->end())
    {
        string kmer=it.getKmer();
        uint64_t count=it.getCount();
        if(count != insertedKmers[kmer])
        {
            cout<<kmer<<endl;
        }
        ASSERT_EQ(count,insertedKmers[kmer]);
        insertedKmers.erase(kmer);
        it++;
    }
    EXPECT_EQ(insertedKmers.size(),0);
    //delete kframe;
}