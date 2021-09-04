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
void deleteFiles(string prefix)
{
    if(prefix=="")
        return;
    string command ="rm -f "+prefix+"*";
    system(command.c_str());
}
map< pair<uint64_t,uint64_t>, insertColorColumn*> insertColumns;





vector<tuple<string,uint64_t,uint64_t> > createQueryInputs() {
    vector<tuple<string, uint64_t, uint64_t> > queryInputs;
    vector<string> queryClasses = {"mixVectors", "prefixTrie"};
    for (auto classType:queryClasses)
        for (auto numColors: {10, 100, 1000})
            for (auto numSamples:{10, 20, 100}) {
                if (numColors == 1000 && (numSamples == 20 || numSamples == 100))
                    continue;
                queryInputs.push_back(make_tuple(classType, numSamples,numColors));
            }
    //queryInputs.push_back(make_tuple("prefixTrie", 10,5000));
    return queryInputs;
}
INSTANTIATE_TEST_SUITE_P(testcolorsTable,
                         queryColumnTest,
                        ::testing::ValuesIn(createQueryInputs())
                      );

void queryColumnTest::SetUp(){
  uint64_t numSamples=get<1>(GetParam());
  uint64_t numColors=get<2>(GetParam());
  auto searchPair=make_pair(numSamples,numColors);
  auto it=insertColumns.find(searchPair);
  insertColorColumn* insertColumn;
  fileName="colorTable.test."+ gen_random(8);
  if(it == insertColumns.end()){
      tmpFolder = "tmp.colorColumn."+gen_random(8);
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

    deleteFiles(fileName);
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

TEST_P(queryColumnTest, clone)
{
    string colorTableName=get<0>(GetParam());
    uint64_t numSamples=get<1>(GetParam());
    uint64_t numColors=get<2>(GetParam());
    auto inputPair= make_pair(numSamples,numColors);
    testColumn= createIndexingColumn(colorTableName,insertColumns[inputPair]);
    EXPECT_EQ(numSamples,testColumn->noSamples);
    EXPECT_EQ(numColors,testColumn->numColors);




    testColumnLoaded=(queryColorColumn*)testColumn->clone();
    delete testColumn;
    testColumn=nullptr;

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
        return new kDataFrameMQF(kSize,NKmersTEST);
    }
    else if(type=="MAP")
    {
        return new kDataFrameMAP(kSize,NKmersTEST);
    }
    else if(type=="PHMAP")
    {
        return new kDataFramePHMAP(kSize,NKmersTEST);
    }
    else if(type=="BMQF")
    {
        string fileName="tmp.kdataframeBMQF."+gen_random(8);
        return new kDataFrameBMQF((uint64_t)kSize,NKmersTEST,fileName);
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
                                 ::testing::Values("MAP","PHMAP"),
                                 ::testing::Values(21,31))
);


//INSTANTIATE_TEST_SUITE_P(testFrames,
//                         kDataFrameBufferedTest,
//                         ::testing::Values(31));

vector<string> fastqFiles={"test.noN.fastq"};
INSTANTIATE_TEST_SUITE_P(testcounting,
                         algorithmsTest,
                        ::testing::Combine(
                                ::testing::Values("MAP","PHMAP"),
                                ::testing::Values(21,31),
                             ::testing::ValuesIn(fastqFiles)
                      ));

INSTANTIATE_TEST_SUITE_P(testntCard,
                         estimateTest,
                         ::testing::ValuesIn(fastqFiles)
                                            );

INSTANTIATE_TEST_SUITE_P(testcounting,
                         kDataFrameBlightTest,
                         ::testing::Combine(
                                 ::testing::Values(31),
                                 ::testing::ValuesIn(fastqFiles)
                                 ));

void setFunctionsTest::SetUp()
{
    tuple<string, int> specification=GetParam();
    auto it=input_SET.find(specification);
    if(it==input_SET.end())
    {
        string kFrameType=get<0>(specification);
        uint32_t k=get<1>(specification);
        vector<kDataFrame*> frames(setFunctionsTestInput[0].size());
        for(int i=0;i<setFunctionsTestInput[0].size();i++)
        {
            frames[i]=(getFrame(specification));
            // Mode 1 : kmers, KmerSize will be cloned from the kFrame
            kProcessor::countKmersFromFile(frames[i], {{"mode", 1}}, setFunctionsTestInput[0][i], 1000);
        }
        kDataFrameMAP* indexRes=new kDataFrameMAP(k);
        kProcessor::indexPriorityQueue(frames,"",indexRes);
        frames.push_back(indexRes);
        input_SET[specification]=frames;
        it=input_SET.find(specification);
    }
    input=it->second;
}
INSTANTIATE_TEST_SUITE_P(testSetFunctions,
                         setFunctionsTest,
                        ::testing::ValuesIn({make_tuple(string("MAP"),21),make_tuple(string("MAP"),31)})
                      );





TEST_P(kDataFrameTest,emptykDataFrame)
{


    kframe=getFrame(GetParam());
    EXPECT_EQ(kframe->empty(), true);

    kframe->insert(kmers->begin()->first);
    EXPECT_EQ(kframe->empty(),false);
    delete kframe;
    kframe=nullptr;
}

TEST_P(kDataFrameTest,insertOneTime)
{
    EXPECT_EQ(kframe->empty(), true);
    int insertedKmers=0;

    for(auto k:*kmers)
    {
        kframe->incrementCount(k.first);
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

    int insertedKmers=0;
    for(auto k:*kmers)
    {
        kframe->setCount(k.first,k.second);
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

    int insertedKmers=0;
    for(auto k:*kmers)
    {
        kframe->incrementCount(k.first);
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
    int insertedKmers=0;
    for(auto k:*kmers)
    {
        kframe->setCount(k.first,k.second);
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
   // unordered_set<string> insertedKmers;
//    insertedKmers.clear();
    bool first=true;
    int numInsertedKmers=0;
    for(auto k:*kmers)
    {
        numInsertedKmers++;
      kframe->setCount(k.first,k.second);

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

    int insertedKmers=0;
    for(auto k:*kmers)
    {
        kframe->setCount(k.first,k.second);
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

    vector<string> correctColumNames={"intColumn","doubleColumn","boolColumn","count"};
    sort(correctColumNames.begin(),correctColumNames.end());
    vector<string> columnNames=kframe->getColumnNames();
    sort(columnNames.begin(),columnNames.end());
    ASSERT_EQ(columnNames,correctColumNames);
    delete kframe;
    kframe= nullptr;

}

TEST_P(kDataFrameTest,multiColumnsDefaultValue)
{
    EXPECT_EQ(kframe->empty(), true);

    int insertedKmers=0;
    for(auto k:*kmers)
    {

        kframe->setCount(k.first,k.second);
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

        int randInt=abs(rand()%1000000);
        double randDouble=(double)(rand()%1000000);
        bool randBool=rand()%2==0;

        if(kframe->getKmerColumnValue<int, vectorColumn<int> >("intColumn",kmer) !=0){
            cout<<kframe->getKmerColumnValue<int, vectorColumn<int> >("intColumn",kmer)<<endl;
            cout<<(simColumns.find(kmer)==simColumns.end())<<endl;
        }
        simColumns[kmer]=make_tuple(randInt,randDouble,randBool);

        if(randInt%3==0){
            kframe->setKmerColumnValue<int, vectorColumn<int> >("intColumn",kmer,randInt);
            kframe->setKmerColumnValue<double, vectorColumn<double> >("doubleColumn",kmer,randDouble);
            kframe->setKmerColumnValue<bool, vectorColumn<bool> >("boolColumn",kmer,randBool);
        }
        it++;
    }
    for(auto simRow:simColumns)
    {
        string kmer=simRow.first;
        int randInt=get<0>(simRow.second);
        double randDouble=get<1>(simRow.second);
        bool randBool=get<2>(simRow.second);

        if(randInt%3!=0)
        {
            randInt=0;
            randDouble=0.0;
            randBool=false;
        }
        else
        {
           // cout<<kmer<<endl;
        }

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


TEST_P(kDataFrameTest,convertTOMQF)
{
    EXPECT_EQ(kframe->empty(), true);

    int insertedKmers=0;
    for(auto k:*kmers)
    {
        kframe->setCount(k.first,k.second);
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


    unordered_map<string,tuple<int,double,bool> > simColumns;
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

    kframeMQF=new kDataFrameMQF(kframe);
    delete kframe;
    kframe= nullptr;

    for(auto simRow:simColumns)
    {
        string kmer=simRow.first;
        int randInt=get<0>(simRow.second);
        double randDouble=get<1>(simRow.second);
        bool randBool=get<2>(simRow.second);

        int retInt=kframeMQF->getKmerColumnValue<int, vectorColumn<int> >("intColumn",kmer);
        double retDouble=kframeMQF->getKmerColumnValue<double, vectorColumn<double> >("doubleColumn",kmer);
        bool retBool=kframeMQF->getKmerColumnValue<bool, vectorColumn<bool> >("boolColumn",kmer);

        ASSERT_EQ(randInt,retInt);
        ASSERT_EQ(randDouble,retDouble);
        EXPECT_EQ(randBool,retBool);

    }
}
TEST_P(kDataFrameTest,convertTOBMQF)
{
    EXPECT_EQ(kframe->empty(), true);

    int insertedKmers=0;
    for(auto k:*kmers)
    {
        kframe->setCount(k.first,k.second);
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8)
        {
            break;
        }
        insertedKmers++;
    }
   // int checkedKmers=0;
    kframe->addColumn("boolColumn",new vectorColumn<bool>(kframe->size()));
    kframe->addColumn("intColumn",new vectorColumn<int>(kframe->size()));
    kframe->addColumn("doubleColumn",new vectorColumn<double>(kframe->size()));


    unordered_map<string,tuple<int,double,bool> > simColumns;
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

    kDataFrameBMQF* kframeMQF=new kDataFrameBMQF(kframe,"tmp.bmqf."+ gen_random(10));
    delete kframe;
    kframe= nullptr;
    unsigned checkedKmers=0;
    for(auto k:*kframeMQF)
    {
        bool kmerExists=simColumns.find(k.getKmer()) != simColumns.end();
        ASSERT_EQ(kmerExists, true);
        checkedKmers++;
    }
    ASSERT_EQ(checkedKmers,simColumns.size());
//    for(auto simRow:simColumns)
//    {
//        string kmer=simRow.first;
//        int randInt=get<0>(simRow.second);
//        double randDouble=get<1>(simRow.second);
//        bool randBool=get<2>(simRow.second);
//
//        bool kmerExists=kframeMQF->kmerExist(kmer);
//        ASSERT_EQ(kmerExists, true);
////        int retInt=kframeMQF->getKmerColumnValue<int, vectorColumn<int> >("intColumn",kmer);
////        double retDouble=kframeMQF->getKmerColumnValue<double, vectorColumn<double> >("doubleColumn",kmer);
////        bool retBool=kframeMQF->getKmerColumnValue<bool, vectorColumn<bool> >("boolColumn",kmer);
////
////        ASSERT_EQ(randInt,retInt);
////        ASSERT_EQ(randDouble,retDouble);
////        EXPECT_EQ(randBool,retBool);
//
//    }
}



TEST_P(kDataFrameTest,saveAndLoadMultiColumns)
{
    EXPECT_EQ(kframe->empty(), true);

    int insertedKmers=0;
    for(auto k:*kmers)
    {
        kframe->setCount(k.first,k.second);
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

    vector<string> correctColumNames={"intColumn","doubleColumn","boolColumn","count"};
    sort(correctColumNames.begin(),correctColumNames.end());
    vector<string> columnNames=kframeLoaded->getColumnNames();
    sort(columnNames.begin(),columnNames.end());
    ASSERT_EQ(columnNames,correctColumNames);

    delete kframeLoaded;
    kframeLoaded=nullptr;

}










TEST_P(kDataFrameTest,saveAndIterateOverAllKmers)
{

    EXPECT_EQ(kframe->empty(), true);

    int numInsertedKmers=0;
  //  unordered_map<string,int> insertedKmers;
    for(auto k:*kmers)
    {
      numInsertedKmers++;
      kframe->setCount(k.first,k.second);
      if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
        break;
      }
    }

    kframe->save(fileName);
    delete kframe;
    kframe=nullptr;
    kframeLoaded=kDataFrame::load(fileName);
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

    //unordered_map<string,int> insertedKmers;
    int numInsertedkmers=0;
    for(auto k:*kmers)
    {
      numInsertedkmers++;
      kframe->setCount(k.first,k.second);
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

    for(auto k:*kmers)
    {
      kframe->setCount(k.first,k.second);
      if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
        break;
      }
    }
    int checkedKmers=0;
    kframe2=kProcessor::filter(kframe,[](kDataFrameIterator& k)
    {
      return k.getCount()>=5;
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

TEST_P(kDataFrameTest,FilterLessThan5MultipleColumns)
{
    EXPECT_EQ(kframe->empty(), true);
    kframe->addColumn("boolColumn",new vectorColumn<bool>());
    for(auto k:*kmers)
    {
        uint32_t count=(k.second%9)+1;
        kframe->setCount(k.first,count);
        kframe->setKmerColumnValue<bool, vectorColumn<bool> >("boolColumn",k.first,count%5==0);
        if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
            break;
        }
    }
    int checkedKmers=0;
    kframe2=kProcessor::filter(kframe,[](kDataFrameIterator& k) -> bool
    {
        uint32_t count=0;
        k.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count",count);
        return count>=5;
    });
    kDataFrameIterator it=kframe2->begin();
    while(it!=kframe2->end())
    {
        string kmer=it.getKmer();
        uint32_t count;
        it.getColumnValue<uint32_t, vectorColumn<uint32_t> >("count",count);
        bool boolvalue;
        it.getColumnValue<bool, vectorColumn<bool> >("boolColumn",boolvalue);

        //ASSERT_EQ(count,((*kmers)[kmer]%9)+1);
        ASSERT_GE(count,5);
        ASSERT_EQ(count%5==0,boolvalue);
        it++;
    }
    delete kframe2;
    kframe2=nullptr;



}

TEST_P(kDataFrameTest,aggregateSum)
{
    EXPECT_EQ(kframe->empty(), true);

    uint64_t goldSum=0;
    for(auto k:*kmers)
    {
      kframe->setCount(k.first,k.second);
      goldSum+=k.second;
      if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
        break;
      }
    }
    int checkedKmers=0;
    any initial=(uint64_t)0;
    any sum=kProcessor::aggregate(kframe,initial,[](kDataFrameIterator& k,any v)
    {
      uint64_t tmp=any_cast<uint64_t>(v);
      any result=tmp+k.getCount();
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
  kmerDecoder *KD_KMERS = kmerDecoder::getInstance(fileName, chunkSize, KMERS, kframe->KD->hash_mode, {{"kSize", kSize}});

    while (!KD_KMERS->end()) {
        KD_KMERS->next_chunk();
        for (const auto &seq : *KD_KMERS->getKmers()) {
            for (const auto &kmer : seq.second) {
                ASSERT_GE(kframe->getCount(kmer.str),1);
            }
        }
    }

    delete KD_KMERS;
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

    kmerDecoder *KD_KMERS = kProcessor::initialize_kmerDecoder(fileName, 1000, "kmers", {{"k_size", kSize}});

    while (!KD_KMERS->end()) {
        KD_KMERS->next_chunk();
        for (const auto &seq : *KD_KMERS->getKmers()) {
            for (const auto &kmer : seq.second) {
                ASSERT_GE(kframe->getCount(kmer.str),1);
            }
        }
    }

    delete KD_KMERS;
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

TEST_P(setFunctionsTest,unioinTest)
{



  try {
      result = kProcessor::kFrameUnion(input);
  }
  catch (const logic_error& expected)
  {
      SUCCEED();
      //FAIL();
  }
  for(int i =0 ;i<input.size()-1;i++)
  {
    kDataFrame* kframe=input[i];
    auto it=kframe->begin();
    while(it!=kframe->end())
    {
      uint32_t countGold;
      it.getColumnValue<uint32_t,vectorColumn<uint32_t>>("count",countGold);
      uint32_t count;
      count = result->getKmerColumnValue<uint32_t,vectorColumn<uint32_t>>("count"+ to_string(i),it.getHashedKmer());
      ASSERT_GE(count,countGold);
      it++;
    }
  }


}
TEST_P(setFunctionsTest,intersectTest)
{
    try {
        result=kProcessor::kFrameIntersect(input);
    }
    catch (const logic_error& expected)
    {
        SUCCEED();
        //FAIL();
    }

  auto it=result->begin();
  while(it!=result->end())
  {
    for(int i=0;i<input.size()-1;i++)
    {
        auto kframe=input[i];
        uint32_t countGold=kframe->getKmerColumnValue<uint32_t,vectorColumn<uint32_t>>("count",it.getHashedKmer());
        uint32_t count;
        it.getColumnValue<uint32_t,vectorColumn<uint32_t>>("count"+ to_string(i),count);
        ASSERT_GE(count,countGold);
    }
    it++;
  }



}
//
TEST_P(setFunctionsTest,differenceTest)
{
  try {
      result=kProcessor::kFrameDiff(input);
  }
  catch (const logic_error& expected)
  {
        SUCCEED();

  }
  auto it=result->begin();
  while(it!=result->end())
  {

      for(int i=1;i<input.size()-1;i++)
      {
        int count=input[i]->getCount((*it).getKmer());
        ASSERT_EQ(count,0);
      }
      it++;
  }


}
TEST_P(setFunctionsTest,innerJoinTest)
{
    try {
        result=kProcessor::innerJoin(input,{1});
    }
    catch (const logic_error& expected)
    {
        SUCCEED();

    }
    auto it=result->begin();
    while(it!=result->end())
    {
        int count=input[1]->getCount((*it).getKmer());
        ASSERT_GE(count,1);
        it++;
    }


}

TEST_P(setFunctionsTest,innerJoinTest2)
{
    try {
        result=kProcessor::innerJoin(input,{0,1,2});
    }
    catch (const logic_error& expected)
    {
        SUCCEED();

    }
    auto it=result->begin();
    while(it!=result->end())
    {
        for(unsigned i=0 ; i<2; i++){
            uint32_t countRes;
            it.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count"+to_string(i),countRes);
            uint32_t countGold=input[i]->getCount(it.getHashedKmer());
            ASSERT_EQ(countRes,countGold);
        }
        vector<uint32_t> colorsCorrect=input[2]->getKmerColumnValue<vector<uint32_t >, deduplicatedColumn<vector<uint32_t>, mixVectors> >("color",it.getHashedKmer());
        vector<uint32_t> colors;
        it.getColumnValue<vector<uint32_t >, deduplicatedColumn<vector<uint32_t>, mixVectors> >("color2",colors);
        ASSERT_EQ(colors,colorsCorrect);
        it++;

    }


}


TEST_P(setFunctionsTest,parallelinnerJoinTest)
{
    try {
        vector<string> inputFileNames(input.size());
        for(unsigned i=0;i<input.size();i++)
        {
            inputFileNames[i]="tmp.parrallelinerJoin."+to_string(i)+".";
            input[i]->save(inputFileNames[i]);
        }
        result=kProcessor::parallelJoin(inputFileNames,{0,1,2});
    }
    catch (const logic_error& expected)
    {
        SUCCEED();

    }
    auto it=result->begin();
    while(it!=result->end())
    {
        for(unsigned i=0 ; i<2; i++){
            uint32_t countRes;
            it.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count"+to_string(i),countRes);
            uint32_t countGold=input[i]->getCount(it.getHashedKmer());
            if(countRes!=countGold)
            {
                cout<<it.getHashedKmer()<<endl;
                for(int k=0;k<2;k++)
                {
                    cout<<k<< "->";
                    cout<<input[k]->kmerExist(it.getHashedKmer())<<" ";
                    cout<<input[k]->getCount(it.getHashedKmer())<<endl;
                }

                cout<<"order in output is "<<result->getkmerOrder(it.getHashedKmer())<<endl;
                for(auto c: result->columns)
                    cout<<c.second->size()<<endl;
            }
            ASSERT_EQ(countRes,countGold);
        }
        vector<uint32_t> colorsCorrect=input[2]->getKmerColumnValue<vector<uint32_t >, deduplicatedColumn<vector<uint32_t>, mixVectors> >("color",it.getHashedKmer());
        vector<uint32_t> colors;
        it.getColumnValue<vector<uint32_t >, deduplicatedColumn<vector<uint32_t>, mixVectors> >("color2",colors);
        ASSERT_EQ(colors,colorsCorrect);
        it++;

    }


}



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
  int q = 25;
 // kDataFrame *KF = new kDataFrameMQF(25, q, integer_hasher);
  KF = new kDataFramePHMAP(25,integer_hasher);
  kmerDecoder *KD_KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});
  kProcessor::index(KD_KMERS, filename+".names", KF);

  uint64_t kSize=KF->getkSize();

  delete KD_KMERS;
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

    KD_KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});

    while (!KD_KMERS->end()) {
            KD_KMERS->next_chunk();
            for (const auto &seq : *KD_KMERS->getKmers()) {
                string readName = seq.first;
                string groupName=namesMap[readName];
                for (const auto &kmer : seq.second) {
                    vector<string> colors=KF->getKmerColumnValue<vector<string> ,deduplicatedColumn<vector<string>, StringColorColumn>>("color",kmer.hash);
                    ASSERT_NE(colors.size(),0);
                    auto colorIt=find(colors.begin(),colors.end(),groupName);
                    ASSERT_NE(colorIt,colors.end());
                }
            }
    }
    delete KF;
    KF= nullptr;
    delete KD_KMERS;

}

TEST_P(indexingTest,indexPriorityQSaveAndLoad)
{
    string filename=GetParam();
    int chunkSize = 1000;

    int q = 25;
    //kDataFrame *KF = new kDataFrameMQF(25, q, integer_hasher);
    KF = new kDataFramePHMAP(25,integer_hasher);
    kmerDecoder *KD_KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});


    vector<kDataFrame*> inputFrames;
    while (!KD_KMERS->end()) {
        KD_KMERS->next_chunk();
        for (const auto &seq : *KD_KMERS->getKmers()) {
            kDataFrame* curr=new kDataFrameMAP(KD_KMERS->get_kSize());
            for (const auto &kmer : seq.second) {
                curr->insert(kmer.hash);
            }
            inputFrames.push_back(curr);
        }
    }

    kProcessor::indexPriorityQueue(inputFrames,"", KF);
    KF->save(fileName);
    delete KF;
    KF= nullptr;
    kframeLoaded=kDataFrame::load(fileName);
    for(int i=0;i<inputFrames.size();i++)
    {
        kDataFrameIterator it=inputFrames[i]->begin();
        while(it!=inputFrames[i]->end())
        {
            vector<uint32_t> colors=kframeLoaded->getKmerColumnValue<vector<uint32_t >, deduplicatedColumn<vector<uint32_t>, mixVectors> >("color",it.getHashedKmer());
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
    delete KD_KMERS;

}

TEST_P(indexingTest,indexPriorityQ)
{
    string filename=GetParam();
    int chunkSize = 1000;
    int q = 25;
//    kDataFrame *KF = new kDataFrameMQF(25, q, integer_hasher);
    KF = new kDataFramePHMAP(25,integer_hasher);

    kmerDecoder *KD_KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});


    vector<kDataFrame*> inputFrames;
    while (!KD_KMERS->end()) {
        KD_KMERS->next_chunk();
        for (const auto &seq : *KD_KMERS->getKmers()) {
            kDataFrame* curr=new kDataFrameMAP(KD_KMERS->get_kSize());
            for (const auto &kmer : seq.second) {
                curr->insert(kmer.hash);
            }
            inputFrames.push_back(curr);
        }
    }

    kProcessor::indexPriorityQueue(inputFrames,"", KF);


    for(int i=0;i<inputFrames.size();i++)
    {
        kDataFrameIterator it=inputFrames[i]->begin();
        while(it!=inputFrames[i]->end())
        {
            vector<uint32_t> colors=KF->getKmerColumnValue<vector<uint32_t >, deduplicatedColumn<vector<uint32_t>, mixVectors> >("color",it.getHashedKmer());
            ASSERT_NE(colors.size(),0);
            auto colorIt=colors.end();
            colorIt=find(colors.begin(),colors.end(),i);
            ASSERT_NE(colorIt,colors.end());
            it.next();
        }
        delete inputFrames[i];
    }

    delete KF;
    KF=nullptr;
    delete KD_KMERS;

}

TEST_P(indexingTest,indexPriorityQAndOptimize)
{
    string filename=GetParam();
    int chunkSize = 1000;
    int q = 25;
    //    kDataFrame *KF = new kDataFrameMQF(25, q, integer_hasher);
    KF = new kDataFramePHMAP(25,integer_hasher);

    kmerDecoder *KD_KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});


    vector<kDataFrame*> inputFrames;
    while (!KD_KMERS->end()) {
        KD_KMERS->next_chunk();
        for (const auto &seq : *KD_KMERS->getKmers()) {
            kDataFrame* curr=new kDataFrameMAP(KD_KMERS->get_kSize());
            for (const auto &kmer : seq.second) {
                curr->insert(kmer.hash);
            }
            inputFrames.push_back(curr);
        }
    }

    kProcessor::indexPriorityQueue(inputFrames,"", KF);

    auto prevColor=(deduplicatedColumn<vector<uint32_t>, mixVectors>*)KF->columns["color"];
    auto newColor=new deduplicatedColumn<vector<uint32_t>, prefixTrie>();
    newColor->index=prevColor->index;
    newColor->values=new prefixTrie(prevColor->values);
    KF->columns["color"]=newColor;
    delete prevColor;

    KF->save(fileName);
    delete KF;
    KF= nullptr;
    kframeLoaded=kDataFrame::load(fileName);

    for(int i=0;i<inputFrames.size();i++)
    {
        kDataFrameIterator it=inputFrames[i]->begin();
        while(it!=inputFrames[i]->end())
        {
            vector<uint32_t> colors=kframeLoaded->getKmerColumnValue<vector<uint32_t >, deduplicatedColumn<vector<uint32_t>, prefixTrie> >("color",it.getHashedKmer());
            ASSERT_NE(colors.size(),0);
            auto colorIt=colors.end();
            colorIt=find(colors.begin(),colors.end(),i);
            ASSERT_NE(colorIt,colors.end());
            it.next();
        }
        delete inputFrames[i];
    }


}
//
//TEST_P(indexingTest,mergeIndexes)
//{
//    string filename=GetParam();
//    int chunkSize = 1000;
//
//    int q = 25;
//    kDataFrame *KF = new kDataFramePHMAP(25, integer_hasher);
//    kmerDecoder *KD_KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});
//
//
//    vector<kDataFrame*> inputFrames;
//    while (!KD_KMERS->end()) {
//        KD_KMERS->next_chunk();
//        for (const auto &seq : *KD_KMERS->getKmers()) {
//            kDataFrame* curr=new kDataFrameMAP(KD_KMERS->get_kSize());
//            for (const auto &kmer : seq.second) {
//                curr->insert(kmer.hash);
//            }
//            inputFrames.push_back(curr);
//        }
//    }
//    int numIndexes=10;
//    int sizeOfIndexes=inputFrames.size()/numIndexes;
//    vector<kDataFrame*> indexes(numIndexes);
//
//    for(int i=0;i<numIndexes;i++)
//    {
//        vector<kDataFrame*> input;
//        input.clear();
//        for(int j=i*sizeOfIndexes; j<(i+1)*sizeOfIndexes;j++) {
//            input.push_back(inputFrames[j]);
//        }
//        if(i==numIndexes-1)
//        {
//            for(int j=(i+1)*sizeOfIndexes; j<inputFrames.size();j++)
//                input.push_back(inputFrames[j]);
//        }
//        int q = 25;
//        kDataFrame *KF2 = new kDataFramePHMAP(25, TwoBits_hasher);
//        kProcessor::indexPriorityQueue(input,"", KF2);
//        indexes[i]=KF2;
//    }
//
//
//    kProcessor::mergeIndexes(indexes, "",KF);
//
//    for(int i=0;i<numIndexes;i++)
//        delete indexes[i];
//
//    for(int i=0;i<inputFrames.size();i++)
//    {
//        kDataFrameIterator it=inputFrames[i]->begin();
//        while(it!=inputFrames[i]->end())
//        {
//            vector<uint32_t> colors=KF->getKmerColumnValue<vector<uint32_t >, deduplicatedColumn<vector<uint32_t>, mixVectors> >("color",it.getHashedKmer());
//            ASSERT_NE(colors.size(),0);
//            auto colorIt=colors.end();
//            colorIt=find(colors.begin(),colors.end(),i);
//            ASSERT_NE(colorIt,colors.end());
//            it.next();
//        }
//    }
//
//    delete KF;
//    delete KD_KMERS;
//
//}

TEST_P(indexingTest,saveAndLoad)
{
    string filename=GetParam();
    int chunkSize = 1000;


    int q = 25;
    KF = new kDataFramePHMAP(25, integer_hasher);
    kmerDecoder *KD_KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});
    kProcessor::index(KD_KMERS, filename+".names", KF);

    KF->save(fileName);
    delete KF;
    KF= nullptr;
    kframeLoaded=kDataFrame::load(fileName);


    uint64_t kSize=kframeLoaded->getkSize();
    //vector<uint32_t> colors;
    delete KD_KMERS;
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

    KD_KMERS = kProcessor::initialize_kmerDecoder(filename, chunkSize, "kmers", {{"k_size", 25}});

    while (!KD_KMERS->end()) {
        KD_KMERS->next_chunk();
        for (const auto &seq : *KD_KMERS->getKmers()) {
            string readName = seq.first;
            string groupName=namesMap[readName];
            for (const auto &kmer : seq.second) {
                vector<string> colors=kframeLoaded->getKmerColumnValue<vector<string> ,deduplicatedColumn<vector<string>, StringColorColumn>>("color",kmer.hash);
                ASSERT_NE(colors.size(),0);
                auto colorIt=colors.end();
                colorIt=find(colors.begin(),colors.end(),groupName);
                ASSERT_NE(colorIt,colors.end());
            }
        }
    }
    delete kframeLoaded;
    kframeLoaded=nullptr;
    delete KD_KMERS;

}


//TEST_P(kDataFrameBufferedTest,iterateOverAllKmers)
//{
//    EXPECT_EQ(kframe->empty(), true);
//    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
//    int numInsertedKmers=0;
//   // unordered_map<string,int> insertedKmers;
//    //insertedKmers.clear();
//
//    for(auto k:*kmers)
//    {
//        numInsertedKmers++;
//        kframe->setCount(k.first,k.second);
//      //  insertedKmers[k.first]+=k.second;
//        if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
//            break;
//        }
//    }
//
//
//    int testKmers=0;
//    for(auto it:*kframe)
//    {
//        string kmer=it.getKmer();
//        uint64_t count=it.getCount();
//        ASSERT_EQ(count,(*kmers)[kmer]);
//        //kmers->erase(kmer);
//        testKmers++;
//        //if(testKmers%1000==0)
//          //  cout<<testKmers<<endl;
//    }
//    EXPECT_EQ(numInsertedKmers,testKmers);
////    delete kframe;
////    kframe=nullptr;
//
//}
//
//TEST_P(kDataFrameBufferedTest,autoResize)
//{
//
//    EXPECT_EQ(kframe->empty(), true);
//    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
//    //unordered_map<string,int> insertedKmers;
//
//    int numInsertedKmers=0;
//    for(auto k:*kmers)
//    {
//        numInsertedKmers++;
//        kframe->setCount(k.first,k.second);
//        if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
//            break;
//        }
//      //  insertedKmers[k.first]+=k.second;
//
//
//    }
//
//    kDataFrameIterator it=kframe->begin();
//    int testedKmers=0;
//    while(it!=kframe->end())
//    {
//        string kmer=it.getKmer();
//        uint64_t count=it.getCount();
//        ASSERT_EQ(count,(*kmers)[kmer]);
//        //insertedKmers.erase(kmer);
//        it++;
//        testedKmers++;
//    }
//    EXPECT_EQ(numInsertedKmers,testedKmers);
////    delete kframe;
////    kframe=nullptr;
//}
//TEST_P(kDataFrameBufferedTest,saveAndIterateOverAllKmers)
//{
//    string filename=kframe->getFilename();;
//    EXPECT_EQ(kframe->empty(), true);
//    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
//    int numInsertedKmers=0;
//    //  unordered_map<string,int> insertedKmers;
//    for(auto k:*kmers)
//    {
//        numInsertedKmers++;
//        kframe->setCount(k.first,k.second);
//        if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
//            break;
//        }
//    }
//
//    kframe->save(filename);
//    delete kframe;
//    kframe=nullptr;
//    kDataFrame* kframeLoaded=kDataFrame::load(filename);
//    int checkedKmers=0;
//    kDataFrameIterator it=kframeLoaded->begin();
//    while(it!=kframeLoaded->end())
//    {
//        string kmer=it.getKmer();
//        uint64_t count=it.getCount();
//        ASSERT_EQ(count,(*kmers)[kmer]);
//        checkedKmers++;
//        it++;
//    }
//    EXPECT_EQ(checkedKmers,numInsertedKmers);
////    delete kframeLoaded;
////    delete kframe;
////    kframe=nullptr;
//}
//
//TEST_P(kDataFrameBufferedTest,saveAndIterateOverAllKmersNoMemory)
//{
//    string filename=kframe->getFilename();
//    EXPECT_EQ(kframe->empty(), true);
//    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
//    int numInsertedKmers=0;
//    //  unordered_map<string,int> insertedKmers;
//    for(auto k:*kmers)
//    {
//        numInsertedKmers++;
//        kframe->setCount(k.first,k.second);
//        if(kframe->load_factor()>=kframe->max_load_factor()*0.8){
//            break;
//        }
//    }
//
//    kframe->save(filename);
//    delete kframe;
//    kframe=nullptr;
//    kframeLoaded=(kDataFrameBMQF*)kDataFrame::load(filename);
//    kframeLoaded->deleteMemoryBuffer();
//    int checkedKmers=0;
//    kDataFrameIterator it=kframeLoaded->begin();
//    while(it!=kframeLoaded->end())
//    {
//        string kmer=it.getKmer();
//        uint64_t count=it.getCount();
//        ASSERT_EQ(count,(*kmers)[kmer]);
//        checkedKmers++;
//        it++;
//    }
//    EXPECT_EQ(checkedKmers,numInsertedKmers);
////    delete kframeLoaded;
////    delete kframe;
////    kframe=nullptr;
////    kframeLoaded= nullptr;
//}
//
//TEST_P(kDataFrameBufferedTest,transformPlus10)
//{
//
//    EXPECT_EQ(kframe->empty(), true);
//    unordered_map<string,int>* kmers=kmersGen->getKmers((int)kframe->getkSize());
//    unordered_map<string,int> insertedKmers;
//    for(auto k:*kmers)
//    {
//        kframe->setCount(k.first,k.second);
//        insertedKmers[k.first]+=k.second;
//        if(kframe->load_factor()>=kframe->max_load_factor()*0.3){
//            break;
//        }
//    }
//    int checkedKmers=0;
//    kframeLoaded=(kDataFrameBMQF*)kProcessor::transform(kframe,[](kmerRow k)
//    {
//        k.count+=10;
//        return k;
//    });
//    kDataFrameIterator it=kframeLoaded->begin();
//    while(it!=kframeLoaded->end())
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

INSTANTIATE_TEST_SUITE_P(testprefixColumn,
        prefixColumnTest,
        ::testing::Values("colorColumn.mixVector"));

TEST_P(prefixColumnTest,optimizeAndCheck)
{
    string filename=GetParam();
    mixVectors *qColumn= new mixVectors();
    qColumn->deserialize(filename);
    prefixTrie *pColumn = new prefixTrie(qColumn);

    for(uint32_t i=1; i<qColumn->numColors; i++)
    {
        auto correctColor=qColumn->getWithIndex(i);
        auto queryColor=pColumn->getWithIndex(i);
        if(correctColor != queryColor)
        {
            cout<<"Not Matching Colors "<<endl;
            cout<<i<<endl;
        }
        EXPECT_EQ(correctColor,queryColor);
        EXPECT_EQ(correctColor.size(),queryColor.size());


    }

    delete pColumn;
    delete qColumn;


}

TEST_P(kDataFrameBlightTest,parsingTest)
{
    int kSize=get<0>(GetParam());
    string fileName=get<1>(GetParam());
    kDataFrame* kframe =new kDataFrameBlight(kSize,fileName);
    int chunkSize=1000;
 // Mode 1 : kmers, KmerSize will be cloned from the kFrame
    kmerDecoder *KD_KMERS = kmerDecoder::getInstance(fileName, chunkSize, KMERS, TwoBits_hasher, {{"kSize", kSize}});

    while (!KD_KMERS->end()) {
        KD_KMERS->next_chunk();
        for (const auto &seq : *KD_KMERS->getKmers()) {
            for (const auto &kmer : seq.second) {
                ASSERT_TRUE(kframe->kmerExist(kmer.str));
            }
        }
    }

    delete KD_KMERS;
    delete kframe;
    kframe=nullptr;


}

