#ifndef TESTING_HPP
#define TESTING_HPP
#include "gtest/gtest.h"
#include "kDataFrame.hpp"
#include <unordered_map>
#include <unistd.h>

using namespace std;
const uint64_t NKmersTEST=100000;
kDataFrame* getFrame(tuple<string,int> input);
string gen_random(const int len);
void deleteFiles(string prefix);
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
            size_t nKmers=NKmersTEST;
            uint64_t range=(1ULL<<(2*kSize));
            while(db[kSize]->size() < nKmers)
            {
                uint64_t kmerInt=rand()%range;
                string kmerStr=kmer::int_to_str(kmerInt,kSize);
                kmerStr=kmer::canonicalKmer(kmerStr);
                uint64_t count=kmerInt%10000+1;

                if(db[kSize]->find(kmerStr) == db[kSize]->end())
                    db[kSize]->insert(make_pair(kmerStr,count));
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
    unordered_map<string,int>* kmers;
    string fileName;
    virtual void SetUp()
    {
        kframe=getFrame(GetParam());
	kframe->addCountColumn();
        kframeLoaded=nullptr;
        kframe2=nullptr;
        kmers=kmersGen->getKmers((int)kframe->getkSize());
        fileName="tmp.kdataframe."+gen_random(8);
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
        deleteFiles(fileName);
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
	kframe->addCountColumn();
	kframeLoaded=nullptr;
    }
    static void SetUpTestSuite()
    {
        kmersGen=new kmersGenerator();
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

class queryColumnTest : public ::testing::TestWithParam<tuple<string, uint64_t,uint64_t> >{
  public:
    unordered_map<uint64_t,vector<uint32_t> > simColors;
    queryColorColumn* testColumn;
    queryColorColumn* testColumnLoaded;
    string tmpFolder;
    string fileName;


    static void TearDownTestSuite();
    virtual void SetUp();
    virtual void TearDown();

};

//class colorsTableInvTest : public ::testing::TestWithParam<colorTableInv* >{
//public:
//  const uint64_t numColors=10000;
//  const uint64_t numSamples=1000;
//  unordered_map<uint64_t,vector<uint32_t> > simColors;
//  virtual void SetUp();
// // virtual void TearDown();
//};
static vector<vector<string> > setFunctionsTestInput={{"test.noN.fastq","test2.noN.fastq"}};
static map<tuple<string,int> ,vector<kDataFrame*>> input_SET;

class setFunctionsTest : public ::testing::TestWithParam<tuple<string, int>>{

public:
    kDataFrame* result;
    vector<kDataFrame*> input;
    virtual void SetUp();

    virtual void TearDown()
    {

        if(result!= nullptr) {
            delete result;
        }
    }
    static void TearDownTestSuite()
    {
        for(auto it:input_SET)
        {
            for(auto frame:it.second)
                delete frame;
        }
    }
};

class indexingTest : public ::testing::TestWithParam<string>{
public:
    kDataFrame* KF;
    kDataFrame* kframeLoaded;
    string fileName;
    virtual void SetUp()
    {
        KF= nullptr;
        kframeLoaded=nullptr;
        fileName="tmp.kdataframe."+gen_random(8);
    }

    virtual void TearDown()
    {

        if(KF!= nullptr)
            delete KF;
        if(kframeLoaded!= nullptr)
            delete kframeLoaded;
        deleteFiles(fileName);
    }
};

class prefixColumnTest : public ::testing::TestWithParam<string>{

};

#endif
