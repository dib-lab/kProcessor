#include "gtest/gtest.h"
#include "kDataFrame.hpp"
#include <unordered_map>
#include "colorTable.hpp"
#include <unistd.h>

using namespace std;
kDataFrame* getFrame(tuple<string,int> input);

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
//        uint64_t kmerCount=(rand()%1000)+1;
                db[kSize]->insert(make_pair(kmerStr,kmerInt));
            }
            it=db.find(kSize);
        }
        return it->second;
    }
    ~kmersGenerator(){
        for(auto i: db)
            delete i.second;
    }
};

static kmersGenerator* kmersGen;

class kDataFrameTest : public ::testing::TestWithParam<tuple<string,int> >{
public:
    kDataFrame* kframe;
    kDataFrame* kframeLoaded;
    kDataFrame* kframe2;
    virtual void SetUp()
    {
        kframe=getFrame(GetParam());
        kframeLoaded=nullptr;
        kframe2=nullptr;
    }
    static void SetUpTestSuite()
    {
        kmersGen=new kmersGenerator();
    }
    virtual void TearDown()
    {

        if(kframe!= nullptr)
            delete kframe;
        if(kframe2!= nullptr)
            delete kframe2;
        if(kframeLoaded!= nullptr)
            delete kframeLoaded;

    }
    static void TearDownTestSuite()
    {
        delete kmersGen;
    }

};

class kDataFrameBufferedTest : public ::testing::TestWithParam<int >{
public:
    kDataFrameBMQF* kframe;
    kDataFrameBMQF* kframeLoaded;
    virtual void SetUp()
    {
        kframe=(kDataFrameBMQF*)getFrame(make_tuple("BMQF",GetParam()));
        kframeLoaded=nullptr;
    }
    static void SetUpTestSuite()
    {
        kmersGen=new kmersGenerator();
    }
    void deleteFiles(string prefix)
    {
        string tmp=prefix+".bmqf";
        unlink(tmp.c_str());
        tmp=prefix+".bmqf.bufferedMem.metadata";
        unlink(tmp.c_str());
        tmp=prefix+".extra";
        unlink(tmp.c_str());
        tmp=prefix+".bmqf.ondisk.metadata";
        unlink(tmp.c_str());
        tmp=prefix+".multiColumn";
        unlink(tmp.c_str());

    }
    virtual void TearDown()
    {

        if(kframe!= nullptr) {
            deleteFiles(kframe->getFilename());
            delete kframe;
        }
        if(kframeLoaded!= nullptr) {
            deleteFiles(kframeLoaded->getFilename());
            delete kframeLoaded;
        }
    }
    static void TearDownTestSuite()
    {
        delete kmersGen;
    }
};

class algorithmsTest : public ::testing::TestWithParam<tuple<string,int,string> >{
};
class estimateTest : public ::testing::TestWithParam<string >{
};
//string is the name of color table class
//and the first integer is the number of samples and the second one is the number of colors
class colorsTableTest : public ::testing::TestWithParam<tuple<string, uint64_t,uint64_t> >{
  public:
    unordered_map<uint64_t,vector<uint32_t> > simColors;
    virtual void SetUp();
};

class colorsTableInvTest : public ::testing::TestWithParam<colorTableInv* >{
public:
  const uint64_t numColors=10000;
  const uint64_t numSamples=1000;
  unordered_map<uint64_t,vector<uint32_t> > simColors;
  virtual void SetUp();
 // virtual void TearDown();
};

class setFunctionsTest : public ::testing::TestWithParam<vector<kDataFrame*>  >{

    virtual void TearDown()
    {
        cout<<"here"<<endl;
    }
};

class indexingTest : public ::testing::TestWithParam<string>{

};

