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
#include "defaultColumn.hpp"


#include <set>
using namespace std;

string gen_random(const int len) {
    srand (time(NULL));
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
map< pair<uint64_t,uint64_t>, insertColorColumn*> insertColumns;
INSTANTIATE_TEST_SUITE_P(testcolorsTable,
                         queryColumnTest,
                        ::testing::Combine(
                        ::testing::Values("mixVectors","prefixTrie"),
                        ::testing::Values(10,20,100),
                        ::testing::Values(10,100)
                      ));

void queryColumnTest::SetUp(){
  uint64_t numSamples=get<1>(GetParam());
  uint64_t numColors=get<2>(GetParam());
  auto searchPair=make_pair(numSamples,numColors);
  auto it=insertColumns.find(searchPair);
  insertColorColumn* insertColumn;
  if(it == insertColumns.end()){
      string tmpFolder = "tmp.colorColumn."+gen_random(8);
      insertColumn =new insertColorColumn(numSamples,tmpFolder);
      while(simColors.size() < numColors)
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
          uint64_t index=insertColumn->insertAndGetIndex(colorsVec);
          simColors[index]=colorsVec;

      }
      insertColumns[searchPair]=insertColumn;

  }

  testColumn=nullptr;
  testColumnLoaded=nullptr;




}
void queryColumnTest::TearDownTestSuite(){
    for(auto it: insertColumns)
    {
        it.second->cleanFiles();
        delete it.second;
    }
}
void queryColumnTest::TearDown(){
    if(testColumn!=nullptr)
        delete testColumn;
    if(testColumnLoaded!=nullptr)
        delete testColumnLoaded;
}
queryColorColumn* createIndexingColumn(string name, insertColorColumn* col)
{
  if(name=="mixVectors")
    return new mixVectors(col);
  if(name=="prefixTrie")
      return new prefixTrie(col);
  return NULL;
}


TEST_P(queryColumnTest, insertAndQuery)
{
  string colorTableName=get<0>(GetParam());
  uint64_t numSamples=get<1>(GetParam());
  uint64_t numColors=get<2>(GetParam());
  auto inputPair= make_pair(numSamples,numColors);
  testColumn= createIndexingColumn(colorTableName,insertColumns[inputPair]);
  EXPECT_EQ(numSamples,testColumn->noSamples);
  EXPECT_EQ(numColors,testColumn->numColors);

  for(auto it:simColors)
  {
    vector<uint32_t> res=testColumn->getWithIndex(it.first);
    EXPECT_EQ(res,it.second);
  }

}

TEST_P(queryColumnTest, saveAndLoad)
{
  string colorTableName=get<0>(GetParam());
  uint64_t numSamples=get<1>(GetParam());
  uint64_t numColors=get<2>(GetParam());
  auto inputPair= make_pair(numSamples,numColors);
  testColumn= createIndexingColumn(colorTableName,insertColumns[inputPair]);
  EXPECT_EQ(numSamples,testColumn->noSamples);
  EXPECT_EQ(numColors,testColumn->numColors);


  string fileName="colorTable.test."+ gen_random(8);
  testColumn->serialize(fileName);
  testColumnLoaded=(queryColorColumn*)Column::getContainerByName(typeid(*(testColumn)).hash_code());
  delete testColumn;
  testColumn=nullptr;
  testColumnLoaded->deserialize(fileName);

  EXPECT_EQ(numSamples,testColumnLoaded->noSamples);
  EXPECT_EQ(numColors,testColumnLoaded->numColors);

  for(auto it:simColors)
  {
      vector<uint32_t> res=testColumnLoaded->getWithIndex(it.first);
      EXPECT_EQ(res,it.second);
  }

}


vector<kDataFrame*> BuildTestFrames()
{
  vector<kDataFrame*> framesToBeTested;
  vector<int> kSizes={21};
  for(auto k:kSizes)
  {
   framesToBeTested.push_back(new kDataFrameMQF(k));
  }
  for(auto k:kSizes)
  {
   // framesToBeTested.push_back(new kDataFrameMAP(k));
  }
  return framesToBeTested;
}
kDataFrame* getFrame(tuple<string,int> input)
{
    string type=get<0>(input);
    uint64_t  kSize=get<1>(input);


    if(type=="MQF")
    {
        return new kDataFrameMQF(kSize);
    }
    else if(type=="MAP")
    {
        return new kDataFrameMAP(kSize);
    }
    else if(type=="PHMAP")
    {
        return new kDataFramePHMAP(kSize);
    }
    else if(type=="BMQF")
    {
        string fileName="tmp.kdataframeMQF."+gen_random(8);
        return new kDataFrameBMQF((uint64_t)kSize,fileName);
    }
    else{
        throw std::logic_error("Unknown kdataframe type");
    }
    return NULL;
}

vector<kDataFrameBMQF*> BuildTestBufferedFrames()
{
    vector<kDataFrameBMQF*> framesToBeTested;
    vector<int> kSizes={21};
    for(auto k:kSizes)
    {
        int randNum=rand();
        string fileName="tmp"+to_string(randNum);
        framesToBeTested.push_back(new kDataFrameBMQF((uint64_t)k,fileName));
    }
    return framesToBeTested;
}

INSTANTIATE_TEST_SUITE_P(testFrames,
                        kDataFrameTest,
                         ::testing::Combine(
                                 ::testing::Values("MQF","MAP","PHMAP"),
                                 ::testing::Values(21,31))
);


INSTANTIATE_TEST_SUITE_P(testFrames,
                         kDataFrameBufferedTest,
                         ::testing::Values(31));

vector<string> fastqFiles={"test.noN.fastq"};
INSTANTIATE_TEST_SUITE_P(testcounting,
                         algorithmsTest,
                        ::testing::Combine(
                                ::testing::Values("MQF","MAP","PHMAP"),
                                ::testing::Values(21,31),
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
    kProcessor::countKmersFromFile(framesToBeTested[0].back(), {{"mode", 1}}, file, 1000); // Mode 1 : kmers, KmerSize will be cloned from the kFrame

    framesToBeTested[1].push_back(new kDataFrameMQF(k)); // Temporary until resolving #17
    kProcessor::countKmersFromFile(framesToBeTested[1].back(), {{"mode", 1}}, file, 1000); // Mode 1 : kmers, KmerSize will be cloned from the kFrame
  }
  return framesToBeTested;
}
//INSTANTIATE_TEST_SUITE_P(testSetFunctions,
//                         setFunctionsTest,
//                        ::testing::ValuesIn(BuildTestFramesForSetFunctions())
//                      );





TEST_P(kDataFrameTest,emptykDataFrame)
{


    kDataFrame* kframe=getFrame(GetParam());
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
    kframe->insert(kmers->begin()->first);
    EXPECT_EQ(kframe->empty(),false);
    delete kframe;
    kframe=nullptr;
}

TEST_P(kDataFrameTest,insertOneTime)
{
    EXPECT_EQ(kframe->empty(), true);
    int insertedKmers=0;
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
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
    kframe= nullptr;

}

TEST_P(kDataFrameTest,insertNTimes)
{
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
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
    kframe= nullptr;
}
TEST_P(kDataFrameTest,eraseKmers)
{
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
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
    kframe= nullptr;
}


TEST_P(kDataFrameTest,autoResize)
{
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
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
    kframe= nullptr;
}

TEST_P(kDataFrameTest,iterateOverAllKmers)
{
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
   // unordered_set<string> insertedKmers;
//    insertedKmers.clear();
    bool first=true;
    int numInsertedKmers=0;
    for(auto k:*kmers)
    {
        numInsertedKmers++;
      kframe->insert(k.first,k.second);

      if(kframe->load_factor()>=kframe->max_load_factor()){
        break;
      }
    }
    int checkedKmers=0;

    for(auto it:*kframe)
    {
      string kmer=it.getKmer();
      uint64_t count=it.getCount();
      ASSERT_EQ(count,(*kmers)[kmer]);
      //insertedKmers.erase(kmer);
      checkedKmers++;
    }
    EXPECT_EQ(checkedKmers,numInsertedKmers);
    delete kframe;
    kframe= nullptr;

}

TEST_P(kDataFrameTest,multiColumns)
{
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
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
    kframe->addColumn("intColumn",new vectorColumn<int>(kframe->size()));
    kframe->addColumn("doubleColumn",new vectorColumn<double>(kframe->size()));
    kframe->addColumn("boolColumn",new vectorColumn<bool>(kframe->size()));

    map<string,tuple<int,double,bool> > simColumns;
    kDataFrameIterator it=kframe->begin();
    while(it!=kframe->end())
    {
      string kmer=it.getKmer();

      int randInt=rand()%1000000;
      double randDouble=(double)(rand()%1000000);
      bool randBool=rand()%2==0;

      simColumns[kmer]=make_tuple(randInt,randDouble,randBool);

      kframe->setKmerColumnValue<int, vectorColumn<int> >("intColumn",kmer,randInt);
      kframe->setKmerColumnValue<double, vectorColumn<double> >("doubleColumn",kmer,randDouble);
      kframe->setKmerColumnValue<bool, vectorColumn<bool> >("boolColumn",kmer,randBool);
      it++;
    }
    for(auto simRow:simColumns)
    {
      string kmer=simRow.first;
      int randInt=get<0>(simRow.second);
      double randDouble=get<1>(simRow.second);
      bool randBool=get<2>(simRow.second);

      int retInt=kframe->getKmerColumnValue<int, vectorColumn<int> >("intColumn",kmer);
      double retDouble=kframe->getKmerColumnValue<double, vectorColumn<double> >("doubleColumn",kmer);
      bool retBool=kframe->getKmerColumnValue<bool, vectorColumn<bool>>("boolColumn",kmer);

      ASSERT_EQ(randInt,retInt);
      ASSERT_EQ(randDouble,retDouble);
      ASSERT_EQ(randBool,retBool);

    }
    delete kframe;
    kframe= nullptr;

}


TEST_P(kDataFrameTest,saveAndLoadMultiColumns)
{
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
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
    kframe->addColumn("boolColumn",new vectorColumn<bool>(kframe->size()));
    kframe->addColumn("intColumn",new vectorColumn<int>(kframe->size()));
    kframe->addColumn("doubleColumn",new vectorColumn<double>(kframe->size()));


    map<string,tuple<int,double,bool> > simColumns;
    kDataFrameIterator it=kframe->begin();
    while(it!=kframe->end())
    {
        string kmer=it.getKmer();

        int randInt=rand()%1000000;
        double randDouble=(double)(rand()%1000000);
        bool randBool=rand()%2==0;

        simColumns[kmer]=make_tuple(randInt,randDouble,randBool);

        kframe->setKmerColumnValue<int, vectorColumn<int> >("intColumn",kmer,randInt);
        kframe->setKmerColumnValue<double, vectorColumn<double> >("doubleColumn",kmer,randDouble);
        kframe->setKmerColumnValue<bool, vectorColumn<bool> >("boolColumn",kmer,randBool);
        it++;
    }
    string fileName="tmp.kdataframe."+gen_random(4);
    kframe->save(fileName);
    delete kframe;
    kframe=nullptr;
    kframeLoaded=kDataFrame::load(fileName);
    for(auto simRow:simColumns)
    {
        string kmer=simRow.first;
        int randInt=get<0>(simRow.second);
        double randDouble=get<1>(simRow.second);
        bool randBool=get<2>(simRow.second);

        int retInt=kframeLoaded->getKmerColumnValue<int, vectorColumn<int> >("intColumn",kmer);
        double retDouble=kframeLoaded->getKmerColumnValue<double, vectorColumn<double> >("doubleColumn",kmer);
        bool retBool=kframeLoaded->getKmerColumnValue<bool, vectorColumn<bool> >("boolColumn",kmer);

        ASSERT_EQ(randInt,retInt);
        ASSERT_EQ(randDouble,retDouble);
        EXPECT_EQ(randBool,retBool);

    }
    delete kframeLoaded;
    kframeLoaded=nullptr;

}



TEST_P(kDataFrameTest,changeDefaultColumn)
{
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
    kframe->changeDefaultColumnType(new vectorColumn<double>());
    map<string,tuple<int,double,bool> > simColumns;

    int insertedKmers=0;
    for(auto k:*kmers)
    {

        int randInt=rand()%1000000;
        double randDouble=(double)(rand()%1000000);
        bool randBool=rand()%2==0;

        simColumns[k.first]=make_tuple(randInt,randDouble,randBool);

        kframe->setKmerDefaultColumnValue<double, vectorColumn<double> >(k.first,randDouble);

        if(kframe->load_factor()>=kframe->max_load_factor()*0.8)
        {
            break;
        }
        insertedKmers++;
    }
    int checkedKmers=0;



    for(auto simRow:simColumns)
    {
        string kmer=simRow.first;
        int randInt=get<0>(simRow.second);
        double randDouble=get<1>(simRow.second);
        bool randBool=get<2>(simRow.second);

        //int retInt=kframe->getKmerColumnValue<int,vector<int> >("intColumn",kmer);
        double retDouble=kframe->getKmerDefaultColumnValue<double, vectorColumn<double> >(kmer);
        //bool retBool=kframe->getKmerColumnValue<bool,vector<bool> >("boolColumn",kmer);

       // ASSERT_EQ(randInt,retInt);
        ASSERT_EQ(randDouble,retDouble);
       // ASSERT_EQ(randBool,retBool);

    }
    delete kframe;
    kframe=nullptr;
}


TEST_P(kDataFrameTest,saveAndLoadChangeDefaultColumn)
{

    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
    kframe->changeDefaultColumnType(new vectorColumn<double>());
    map<string,tuple<int,double,bool> > simColumns;

    int insertedKmers=0;
    for(auto k:*kmers)
    {

        int randInt=rand()%1000000;
        double randDouble=(double)(rand()%1000000);
        bool randBool=rand()%2==0;

        simColumns[k.first]=make_tuple(randInt,randDouble,randBool);

        kframe->setKmerDefaultColumnValue<double, vectorColumn<double> >(k.first,randDouble);

        if(kframe->load_factor()>=kframe->max_load_factor()*0.8)
        {
            break;
        }
        insertedKmers++;
    }
    int checkedKmers=0;
    string fileName="tmp.kdataframe."+gen_random(4);
    kframe->save(fileName);
    delete kframe;
    kframe=nullptr;
    kframeLoaded=kDataFrame::load(fileName);

    for(auto simRow:simColumns)
    {
        string kmer=simRow.first;
        double randDouble=get<1>(simRow.second);
        double retDouble=kframeLoaded->getKmerDefaultColumnValue<double, vectorColumn<double> >(kmer);
        ASSERT_EQ(randDouble,retDouble);
    }
    delete kframeLoaded;
    kframeLoaded= nullptr;

}




TEST_P(kDataFrameTest,saveAndIterateOverAllKmers)
{

    kDataFrame* kframe=getFrame(GetParam());
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
    int numInsertedKmers=0;
  //  unordered_map<string,int> insertedKmers;
    for(auto k:*kmers)
    {
      numInsertedKmers++;
      kframe->insert(k.first,k.second);
      if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
        break;
      }
    }
    string fileName="tmp.kdataframe."+gen_random(8);
    kframe->save(fileName);
    delete kframe;
    kframe=nullptr;
    kDataFrame* kframeLoaded=kDataFrame::load(fileName);
    int checkedKmers=0;
    kDataFrameIterator it=kframeLoaded->begin();
    while(it!=kframeLoaded->end())
    {
      string kmer=it.getKmer();
      uint64_t count=it.getCount();
      //cout<<kmer<<endl;
      ASSERT_EQ(count,(*kmers)[kmer]);
      checkedKmers++;
      it++;
    }
    EXPECT_EQ(checkedKmers,numInsertedKmers);
    delete kframeLoaded;
    kframeLoaded=nullptr;


}


TEST_P(kDataFrameTest,transformPlus10)
{
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
    //unordered_map<string,int> insertedKmers;
    int numInsertedkmers=0;
    for(auto k:*kmers)
    {
      numInsertedkmers++;
      kframe->insert(k.first,k.second);
      if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
        break;
      }

    }
    kframe2=kProcessor::transform(kframe,[](kmerRow k)
    {
      k.count+=10;
      return k;
    });
    kDataFrameIterator it=kframe2->begin();
    int numCheckedKmers=0;
    while(it!=kframe2->end())
    {
      numCheckedKmers++;
      string kmer=it.getKmer();
      uint64_t count=it.getCount();
      ASSERT_EQ(count,(*kmers)[kmer]+10);
      it++;
    }
    EXPECT_EQ(numCheckedKmers,numInsertedkmers);
    delete kframe2;
    kframe2= nullptr;

}

TEST_P(kDataFrameTest,transformFilterLessThan5)
{
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
    for(auto k:*kmers)
    {
      kframe->insert(k.first,k.second);
      if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
        break;
      }
    }
    int checkedKmers=0;
    kframe2=kProcessor::filter(kframe,[](kmerRow k)
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

    delete kframe2;
    kframe2=nullptr;



}

TEST_P(kDataFrameTest,transformFilterLessThan5MultipleColumns)
{
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
    for(auto k:*kmers)
    {
        kframe->insert(k.first,k.second);
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
            break;
        }
    }
    int checkedKmers=0;
    kProcessor::createCountColumn(kframe);
    kframe2=kProcessor::filter(kframe,[](kDataFrameIterator& k) -> bool
    {
        uint32_t count;
        k.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count",count);
        return count>=5;
    });
    kDataFrameIterator it=kframe2->begin();
    while(it!=kframe2->end())
    {
        string kmer=it.getKmer();
        uint64_t count=it.getCount();
        ASSERT_GE(count,5);
        it++;
    }
    delete kframe2;
    kframe2=nullptr;



}

TEST_P(kDataFrameTest,aggregateSum)
{
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
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
  string kframeType=get<0>(GetParam());
  int kSize=get<1>(GetParam());
  kDataFrame* kframe=getFrame(make_tuple(kframeType,kSize));
  string fileName=get<2>(GetParam());
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

    delete KMERS;
    delete kframe;
    kframe=nullptr;
    
    
}

TEST_P(algorithmsTest,loadingKMCTest)
{
    string kframeType=get<0>(GetParam());
    int kSize=get<1>(GetParam());
    kDataFrame* kframe=getFrame(make_tuple(kframeType,kSize));
    string fileName=get<2>(GetParam());
    string db=fileName+"."+std::to_string(kSize);


    kProcessor::loadFromKMC(kframe,db);

    kmerDecoder *KMERS = kProcessor::initialize_kmerDecoder(fileName, 1000, "kmers", {{"k_size", kSize}});

    while (!KMERS->end()) {
        KMERS->next_chunk();
        for (const auto &seq : *KMERS->getKmers()) {
            for (const auto &kmer : seq.second) {
                ASSERT_GE(kframe->getCount(kmer.str),1);
            }
        }
    }

    delete KMERS;
    delete kframe;
    kframe=nullptr;


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

//TEST_P(setFunctionsTest,unioinTest)
//{
//  vector<kDataFrame*> input=GetParam();
//    kDataFrame *unioinResult;
//
//  try {
//      unioinResult = kProcessor::kFrameUnion(input);
//  }
//  catch (const logic_error& expected)
//  {
//      SUCCEED();
//      //FAIL();
//  }
//  for(auto kframe:input)
//  {
//    auto it=kframe->begin();
//    while(it!=kframe->end())
//    {
//      int count=unioinResult->getCount((*it).getKmer());
//      ASSERT_GE(count,((*it).getCount()));
//      it++;
//    }
//  }
//  delete unioinResult;
//
//
//}
//TEST_P(setFunctionsTest,intersectTest)
//{
//  vector<kDataFrame*> input=GetParam();
//
//  kDataFrame* intersectResult;
//    try {
//        intersectResult=kProcessor::kFrameIntersect(input);
//    }
//    catch (const logic_error& expected)
//    {
//        SUCCEED();
//        //FAIL();
//    }
//
//  auto it=intersectResult->begin();
//  while(it!=intersectResult->end())
//  {
//    for(auto kframe:input)
//    {
//      int count=kframe->getCount((*it).getKmer());
//      ASSERT_GE(count,((*it).getCount()));
//    }
//    it++;
//  }
//  delete intersectResult;
//
//
//}
//
//TEST_P(setFunctionsTest,differenceTest)
//{
//  vector<kDataFrame*> input=GetParam();
//
//  kDataFrame* diffResult=kProcessor::kFrameDiff(input);
//  try {
//      kDataFrame* diffResult=kProcessor::kFrameDiff(input);
//  }
//  catch (const logic_error& expected)
//  {
//        SUCCEED();
//
//  }
//  auto it=diffResult->begin();
//  while(it!=diffResult->end())
//  {
//      int count=input[0]->getCount((*it).getKmer());
//      ASSERT_EQ(count,((*it).getCount()));
//      for(int i=1;i<input.size();i++)
//      {
//        int count=input[i]->getCount((*it).getKmer());
//        ASSERT_EQ(count,0);
//      }
//      it++;
//  }
//  delete diffResult;
//
//}

//vector<colorTableInv*> BuildColorTableInv()
//{
//  vector<colorTableInv*> tablesToBeTested;
//  tablesToBeTested.push_back(new stringColorTableInv());
//  return tablesToBeTested;
//}
//
//void colorsTableInvTest::SetUp(){
//  for(int i=1;i<=numColors;i++)
//  {
//    int n=rand()%(numSamples*2)+1;
//    set<uint32_t> colors;
//    colors.clear();
//    for(int j=0;j<n;j++)
//    {
//      colors.insert(rand()%numSamples);
//    }
//    vector<uint32_t> colorsVec;
//    colorsVec.clear();
//    copy(colors.begin(),colors.end(),back_inserter(colorsVec));
//    simColors[i]=colorsVec;
//  }
//}
//
//INSTANTIATE_TEST_SUITE_P(testcolorsTableInvTest,
//                        colorsTableInvTest,
//                        ::testing::ValuesIn(BuildColorTableInv()));
//
//
//TEST_P(colorsTableInvTest,setAndGet)
//{
//    colorTableInv* table=GetParam();
//    for(int i=1;i<=numColors;i++)
//    {
//      if(table->getColorId(simColors[i])==0)
//        table->setColorId(i,simColors[i]);
//      else{
//        simColors.erase(i);
//      }
//    }
//    for(auto it:simColors)
//    {
//      vector<uint32_t> res;
//      res.clear();
//      uint64_t colorId =table->getColorId(it.second);
//      EXPECT_EQ(colorId,it.first);
//    }
//
//}

INSTANTIATE_TEST_SUITE_P(testIndexing,
                        indexingTest,
                        ::testing::Values("test1.fa"));


TEST_P(indexingTest,index)
{
  string filename=GetParam();
  int chunkSize = 1000;

  kDataFrame *KF = new kDataFrameMQF(25, 25, 1);
  kmerDecoder *KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});
  kProcessor::index(KMERS, filename+".names", KF);

  uint64_t kSize=KF->getkSize();

  delete KMERS;
  string names_fileName=filename+".names";
  ifstream namesFile(names_fileName.c_str());
  string seqName, groupName,line;
    flat_hash_map<string, string> namesMap;
    while (std::getline(namesFile, line)) {
        std::vector<string> tokens;
        std::istringstream iss(line);
        std::string token;
        while (std::getline(iss, token, '\t'))   // but we can specify a different one
            tokens.push_back(token);
        seqName = tokens[0];
        groupName = tokens[1];
        namesMap.insert(make_pair(seqName, groupName));
    }
    namesFile.close();

    KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});

    while (!KMERS->end()) {
            KMERS->next_chunk();
            for (const auto &seq : *KMERS->getKmers()) {
                string readName = seq.first;
                string groupName=namesMap[readName];
                for (const auto &kmer : seq.second) {
                    vector<string> colors=KF->getKmerDefaultColumnValue<vector<string> ,StringColorColumn>(kmer.hash);
                    ASSERT_NE(colors.size(),0);
                    auto colorIt=colors.end();
                    colorIt=find(colors.begin(),colors.end(),groupName);
                    ASSERT_NE(colorIt,colors.end());
                }
            }
    }
    delete KF;
    delete KMERS;

}

TEST_P(indexingTest,indexPriorityQSaveAndLoad)
{
    string filename=GetParam();
    int chunkSize = 1000;

    kDataFrame *KF = new kDataFrameMQF(25, 25, 1);
    kmerDecoder *KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});


    vector<kDataFrame*> inputFrames;
    while (!KMERS->end()) {
        KMERS->next_chunk();
        for (const auto &seq : *KMERS->getKmers()) {
            kDataFrame* curr=new kDataFrameMAP(KMERS->get_kSize());
            for (const auto &kmer : seq.second) {
                curr->insert(kmer.hash);
            }
            inputFrames.push_back(curr);
        }
    }

    kProcessor::indexPriorityQueue(inputFrames,"", KF);
    string fileName="tmp.kdataframe."+gen_random(4);
    KF->save(fileName);
    delete KF;
    kDataFrame* kframeLoaded=kDataFrame::load(fileName);
    for(int i=0;i<inputFrames.size();i++)
    {
        kDataFrameIterator it=inputFrames[i]->begin();
        while(it!=inputFrames[i]->end())
        {
            vector<uint32_t> colors=kframeLoaded->getKmerDefaultColumnValue<vector<uint32_t >, mixVectors>(it.getHashedKmer());
            ASSERT_NE(colors.size(),0);
            auto colorIt=colors.end();
            colorIt=find(colors.begin(),colors.end(),i);
            ASSERT_NE(colorIt,colors.end());
            it.next();
        }
        delete inputFrames[i];
    }

    delete kframeLoaded;
    kframeLoaded=nullptr;
    delete KMERS;

}

TEST_P(indexingTest,indexPriorityQ)
{
    string filename=GetParam();
    int chunkSize = 1000;

    kDataFrame *KF = new kDataFrameMQF(25, 25, 1);
    kmerDecoder *KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});


    vector<kDataFrame*> inputFrames;
    while (!KMERS->end()) {
        KMERS->next_chunk();
        for (const auto &seq : *KMERS->getKmers()) {
            kDataFrame* curr=new kDataFrameMAP(KMERS->get_kSize());
            for (const auto &kmer : seq.second) {
                curr->insert(kmer.hash);
            }
            inputFrames.push_back(curr);
        }
    }

    kProcessor::indexPriorityQueue(inputFrames,"", KF);
    string fileName="tmp.kdataframe."+gen_random(4);


    for(int i=0;i<inputFrames.size();i++)
    {
        kDataFrameIterator it=inputFrames[i]->begin();
        while(it!=inputFrames[i]->end())
        {
            vector<uint32_t> colors=KF->getKmerDefaultColumnValue<vector<uint32_t >, mixVectors>(it.getHashedKmer());
            ASSERT_NE(colors.size(),0);
            auto colorIt=colors.end();
            colorIt=find(colors.begin(),colors.end(),i);
            ASSERT_NE(colorIt,colors.end());
            it.next();
        }
        delete inputFrames[i];
    }

    delete KF;
    delete KMERS;

}

TEST_P(indexingTest,mergeIndexes)
{
    string filename=GetParam();
    int chunkSize = 1000;

    kDataFrame *KF = new kDataFrameMQF(25, 25, 1);
    kmerDecoder *KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});


    vector<kDataFrame*> inputFrames;
    while (!KMERS->end()) {
        KMERS->next_chunk();
        for (const auto &seq : *KMERS->getKmers()) {
            kDataFrame* curr=new kDataFrameMAP(KMERS->get_kSize());
            for (const auto &kmer : seq.second) {
                curr->insert(kmer.hash);
            }
            inputFrames.push_back(curr);
        }
    }
    int numIndexes=10;
    int sizeOfIndexes=inputFrames.size()/numIndexes;
    vector<kDataFrame*> indexes(numIndexes);

    for(int i=0;i<numIndexes;i++)
    {
        vector<kDataFrame*> input;
        input.clear();
        for(int j=i*sizeOfIndexes; j<(i+1)*sizeOfIndexes;j++) {
            input.push_back(inputFrames[j]);
        }
        if(i==numIndexes-1)
        {
            for(int j=(i+1)*sizeOfIndexes; j<inputFrames.size();j++)
                input.push_back(inputFrames[j]);
        }
        kDataFrame *KF2 = new kDataFrameMQF(25,25,2);
        kProcessor::indexPriorityQueue(input,"", KF2);
        indexes[i]=KF2;
    }


    kProcessor::mergeIndexes(indexes, "",KF);

    for(int i=0;i<numIndexes;i++)
        delete indexes[i];

    for(int i=0;i<inputFrames.size();i++)
    {
        kDataFrameIterator it=inputFrames[i]->begin();
        while(it!=inputFrames[i]->end())
        {
            vector<uint32_t> colors=KF->getKmerDefaultColumnValue<vector<uint32_t >, mixVectors>(it.getHashedKmer());
            ASSERT_NE(colors.size(),0);
            auto colorIt=colors.end();
            colorIt=find(colors.begin(),colors.end(),i);
            ASSERT_NE(colorIt,colors.end());
            it.next();
        }
    }

    delete KF;
    delete KMERS;

}

TEST_P(indexingTest,saveAndLoad)
{
    string filename=GetParam();
    int chunkSize = 1000;

    kDataFrame *KF = new kDataFrameMQF(25, 25, 1);
    kmerDecoder *KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});
    kProcessor::index(KMERS, filename+".names", KF);
    string fileName="tmp.kdataframe."+gen_random(4);
    KF->save(fileName);
    delete KF;
    kDataFrame* kframeLoaded=kDataFrame::load(fileName);


    uint64_t kSize=kframeLoaded->getkSize();
    //vector<uint32_t> colors;
    delete KMERS;
    string names_fileName=filename+".names";
    ifstream namesFile(names_fileName.c_str());
    string seqName, groupName,line;
    flat_hash_map<string, string> namesMap;
    while (std::getline(namesFile, line)) {
        std::vector<string> tokens;
        std::istringstream iss(line);
        std::string token;
        while (std::getline(iss, token, '\t'))   // but we can specify a different one
            tokens.push_back(token);
        seqName = tokens[0];
        groupName = tokens[1];
        namesMap.insert(make_pair(seqName, groupName));
    }
    namesFile.close();

    KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});

    while (!KMERS->end()) {
        KMERS->next_chunk();
        for (const auto &seq : *KMERS->getKmers()) {
            string readName = seq.first;
            string groupName=namesMap[readName];
            for (const auto &kmer : seq.second) {
                vector<string> colors=kframeLoaded->getKmerDefaultColumnValue<vector<string> ,StringColorColumn>(kmer.hash);
                ASSERT_NE(colors.size(),0);
                auto colorIt=colors.end();
                colorIt=find(colors.begin(),colors.end(),groupName);
                ASSERT_NE(colorIt,colors.end());
            }
        }
    }
    delete kframeLoaded;
    kframeLoaded=nullptr;
    delete KMERS;

}


TEST_P(kDataFrameBufferedTest,iterateOverAllKmers)
{
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
    int numInsertedKmers=0;
   // unordered_map<string,int> insertedKmers;
    //insertedKmers.clear();

    for(auto k:*kmers)
    {
        numInsertedKmers++;
        kframe->insert(k.first,k.second);
      //  insertedKmers[k.first]+=k.second;
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
            break;
        }
    }


    int testKmers=0;
    for(auto it:*kframe)
    {
        string kmer=it.getKmer();
        uint64_t count=it.getCount();
        ASSERT_EQ(count,(*kmers)[kmer]);
        //kmers->erase(kmer);
        testKmers++;
        //if(testKmers%1000==0)
          //  cout<<testKmers<<endl;
    }
    EXPECT_EQ(numInsertedKmers,testKmers);
//    delete kframe;
//    kframe=nullptr;

}

TEST_P(kDataFrameBufferedTest,autoResize)
{

    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
    //unordered_map<string,int> insertedKmers;

    int numInsertedKmers=0;
    for(auto k:*kmers)
    {
        numInsertedKmers++;
        kframe->insert(k.first,k.second);
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
            break;
        }
      //  insertedKmers[k.first]+=k.second;


    }

    kDataFrameIterator it=kframe->begin();
    int testedKmers=0;
    while(it!=kframe->end())
    {
        string kmer=it.getKmer();
        uint64_t count=it.getCount();
        ASSERT_EQ(count,(*kmers)[kmer]);
        //insertedKmers.erase(kmer);
        it++;
        testedKmers++;
    }
    EXPECT_EQ(numInsertedKmers,testedKmers);
//    delete kframe;
//    kframe=nullptr;
}
TEST_P(kDataFrameBufferedTest,saveAndIterateOverAllKmers)
{
    string filename=kframe->getFilename();;
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
    int numInsertedKmers=0;
    //  unordered_map<string,int> insertedKmers;
    for(auto k:*kmers)
    {
        numInsertedKmers++;
        kframe->insert(k.first,k.second);
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
            break;
        }
    }

    kframe->save(filename);
    delete kframe;
    kframe=nullptr;
    kDataFrame* kframeLoaded=kDataFrame::load(filename);
    int checkedKmers=0;
    kDataFrameIterator it=kframeLoaded->begin();
    while(it!=kframeLoaded->end())
    {
        string kmer=it.getKmer();
        uint64_t count=it.getCount();
        ASSERT_EQ(count,(*kmers)[kmer]);
        checkedKmers++;
        it++;
    }
    EXPECT_EQ(checkedKmers,numInsertedKmers);
//    delete kframeLoaded;
//    delete kframe;
//    kframe=nullptr;
}

TEST_P(kDataFrameBufferedTest,saveAndIterateOverAllKmersNoMemory)
{
    string filename=kframe->getFilename();
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
    int numInsertedKmers=0;
    //  unordered_map<string,int> insertedKmers;
    for(auto k:*kmers)
    {
        numInsertedKmers++;
        kframe->insert(k.first,k.second);
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
            break;
        }
    }

    kframe->save(filename);
    delete kframe;
    kframe=nullptr;
    kframeLoaded=(kDataFrameBMQF*)kDataFrame::load(filename);
    kframeLoaded->deleteMemoryBuffer();
    int checkedKmers=0;
    kDataFrameIterator it=kframeLoaded->begin();
    while(it!=kframeLoaded->end())
    {
        string kmer=it.getKmer();
        uint64_t count=it.getCount();
        ASSERT_EQ(count,(*kmers)[kmer]);
        checkedKmers++;
        it++;
    }
    EXPECT_EQ(checkedKmers,numInsertedKmers);
//    delete kframeLoaded;
//    delete kframe;
//    kframe=nullptr;
//    kframeLoaded= nullptr;
}

TEST_P(kDataFrameBufferedTest,transformPlus10)
{

    kDataFrame* kframe=(kDataFrameBMQF*)getFrame(make_tuple("BMQF",GetParam()));
    EXPECT_EQ(kframe->empty(), true);
    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
    unordered_map<string,int> insertedKmers;
    for(auto k:*kmers)
    {
        kframe->insert(k.first,k.second);
        insertedKmers[k.first]+=k.second;
        if(kframe->load_factor()>=kframe->max_load_factor()*0.3){
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

}

//TEST_P(kDataFrameBufferedTest,parsingTest)
//{
//    kDataFrameBMQF* kframe=(kDataFrameBMQF*)getFrame(make_tuple("BMQF",GetParam()));
//    string fileName="test2.noN.fastq";
////kDataFrame* kframe=get<0>(GetParam())->getTwin();
////string fileName=get<1>(GetParam());
//    int kSize=kframe->getkSize();
//    kProcessor::countKmersFromFile(kframe, {{"mode", 1}}, fileName, 1000); // Mode 1 : kmers, KmerSize will be cloned from the kFrame
//
//    ifstream kmerCountGoldFile("test.noN.dsk.txt");
//    string kmer;
//    uint64_t  count;
//    while(kmerCountGoldFile>>kmer>>count)
//    {
//        uint64_t kDatframe_count=kframe->getCount(kmer);
//        ASSERT_EQ(count,kDatframe_count);
//    }
//    kmerCountGoldFile.close();
//
//    seqan::SeqFileIn seqIn(fileName.c_str());
//    seqan::StringSet<seqan::CharString> ids;
//    seqan::StringSet<seqan::CharString> reads;
//    int chunkSize=1000;
//    unordered_map<string,uint64_t > insertedKmers;
//    while(!atEnd(seqIn)){
//        clear(reads);
//        clear(ids);
//
//        seqan::readRecords(ids, reads, seqIn,chunkSize);
//        for(int j=0;j<length(reads);j++)
//        {
//            string seq=string((char*)seqan::toCString(reads[j]));
//            for(int i=0;i<seq.size()-kSize+1;i++)
//            {
//                string kmer=seq.substr(i,kSize);
//                kmer=kmer::canonicalKmer(kmer);
//                insertedKmers[kmer]++;
//            }
//        }
//
//    }
//    seqan::close(seqIn);
//    kDataFrameIterator it=kframe->begin();
//    while(it!=kframe->end())
//    {
//        string kmer=it.getKmer();
//        uint64_t count=it.getCount();
//        if(count != insertedKmers[kmer])
//        {
//            cout<<kmer<<endl;
//
//        }
//        EXPECT_EQ(count,insertedKmers[kmer]);
//        insertedKmers.erase(kmer);
//        it++;
//    }
//    EXPECT_EQ(insertedKmers.size(),0);
// //   delete kframe;
//}

TEST_P(algorithmsTest,parsingTest2)
{
    string kframeType=get<0>(GetParam());
    int kSize=get<1>(GetParam());
    kDataFrame* kframe=getFrame(make_tuple(kframeType,kSize));
    string fileName=get<2>(GetParam());
    RecordProperty("kdataFrame Type", kframe->get_class_name());
    kProcessor::countKmersFromFile(kframe, {{"mode", 1}}, fileName, 1000); // Mode 1 : kmers, KmerSize will be cloned from the kFrame
    string goldFileName=fileName.substr(0,fileName.size()-6)+"."+to_string(kSize)+".dsk.txt";
    ifstream kmerCountGoldFile(goldFileName);
    string kmer;
    uint64_t  count;
    while(kmerCountGoldFile>>kmer>>count)
    {
        uint64_t kDatframe_count=kframe->getCount(kmer);
        ASSERT_EQ(count,kDatframe_count);
    }
    kmerCountGoldFile.close();
//    seqan::SeqFileIn seqIn(fileName.c_str());
//    seqan::StringSet<seqan::CharString> ids;
//    seqan::StringSet<seqan::CharString> reads;
//
//
//    while(!atEnd(seqIn)){
//        clear(reads);
//        clear(ids);
//
//        seqan::readRecords(ids, reads, seqIn,chunkSize);
//        for(int j=0;j<length(reads);j++)
//        {
//            string seq=string((char*)seqan::toCString(reads[j]));
//            for(int i=0;i<seq.size()-kSize+1;i++)
//            {
//                string kmer=seq.substr(i,kSize);
//                kmer=kmer::canonicalKmer(kmer);
//                insertedKmers[kmer]++;
//            }
//        }
//
//    }
//    seqan::close(seqIn);
//    kDataFrameIterator it=kframe->begin();
//    while(it!=kframe->end())
//    {
//        string kmer=it.getKmer();
//        uint64_t count=it.getCount();
//        if(count != insertedKmers[kmer])
//        {
//            cout<<kmer<<endl;
//        }
//        ASSERT_EQ(count,insertedKmers[kmer]);
//        insertedKmers.erase(kmer);
//        it++;
//    }
//    EXPECT_EQ(insertedKmers.size(),0);
    delete kframe;
    kframe=nullptr;
}