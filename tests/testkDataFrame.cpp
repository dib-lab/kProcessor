#include "../ThirdParty/catch.hpp"
#include "../algorithms.hpp"
#include <gqf.hpp>
#include <stdint.h>
#include <fstream>
#include <map>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <iostream>
#include "../kDataFrame.hpp"
#include "../KmerDecoder/FastqReader.hpp"
#include <algorithm>
using namespace std;

TEST_CASE("MAP load fasta and query")
{

  map<uint64_t, int> gold;
  int kSize = 31;
  vector<kDataFrame *> kframes;
  kframes.push_back(new kDataFrameMAP(kSize));

  for (auto kframe : kframes)
  {
    FastqReader reader("tests/testData/test.fastq");
    deque<pair<string, string>> sequences;

    while (!reader.isEOF())
    {
      reader.readNSeq(&sequences);
    }
    string kmer;

    for (auto seqPair : sequences)
    {

      string seq = seqPair.first;

      for (int i = 0; i < seq.size() - kSize; i++)
      {

        kmer = seq.substr(i, kSize);
        kframe->incrementCounter(kmer, 1);

        uint64_t kmerHash = kframe->hashKmer(kmer);

        auto goldIT = gold.find(kmerHash);
        if (goldIT == gold.end())
          gold.insert(make_pair(kmerHash, 1));
        // else
        // {
        //   goldIT->second++; // Comment this
        // }
      }
    }

    for (auto seqPair : sequences)
    {

      string seq = seqPair.first;
      for (int i = 0; i < seq.size() - kSize; i++)
      {
        kmer = seq.substr(i, kSize);
        uint64_t kmerHash = kframe->hashKmer(kmer);
        auto goldIT = gold.find(kmerHash);
        // REQUIRE(goldIT->second == kframe->getCounter(kmer));
      }
    }

    kframe->removeKmer(kmer);

    REQUIRE(kframe->getCounter(kmer) == 0);

    kframe->setCounter(kmer, 10);
    REQUIRE(kframe->getCounter(kmer));

    kframe->setCounter(kmer, 4);
    REQUIRE(kframe->getCounter(kmer));
    kframe->setCounter(kmer, 15);
    REQUIRE(kframe->getCounter(kmer));

    kframe->setTag(kmer, 1);
    REQUIRE(kframe->getTag(kmer) == 1);
  }
}

TEST_CASE( "load fasta and query" ) {

    map<uint64_t,int> gold;
    int kSize=31;
    vector<kDataFrame*> kframes;
    kframes.push_back(new kDataFrameMQF(kSize,20,2,2,0));

    for(auto kframe: kframes){
      FastqReader reader("tests/testData/test.fastq");
      deque<pair<string,string> > sequences;
      while(!reader.isEOF())
      {
        reader.readNSeq(&sequences);
      }
      string kmer;
      for(auto seqPair:sequences){

          string seq=seqPair.first;
          for(int i=0;i<seq.size()-kSize;i++)
          {
            kmer=seq.substr(i,kSize);
            kframe->incrementCounter(kmer,1);
            uint64_t kmerHash=kframe->hashKmer(kmer);
            auto goldIT=gold.find(kmerHash);
            if(goldIT==gold.end())
              gold.insert(make_pair(kmerHash,1));
            else{
              goldIT->second++; // Comment this
            }
          }
      }

      for(auto seqPair:sequences){

          string seq=seqPair.first;
          for(int i=0;i<seq.size()-kSize;i++)
          {
            kmer=seq.substr(i,kSize);
            uint64_t kmerHash=kframe->hashKmer(kmer);
            auto goldIT=gold.find(kmerHash);
            REQUIRE(goldIT->second==kframe->getCounter(kmer));
          }
      }





      kframe->removeKmer(kmer);
      REQUIRE(kframe->getCounter(kmer)==0);


      kframe->setCounter(kmer,10);
      REQUIRE(kframe->getCounter(kmer)==10);

      kframe->setCounter(kmer,4);
      REQUIRE(kframe->getCounter(kmer)==4);
      kframe->setCounter(kmer,15);
      REQUIRE(kframe->getCounter(kmer)==15);



      kframe->setTag(kmer,1);
      REQUIRE(kframe->getTag(kmer)==1);
    }

}

TEST_CASE( "bufferedMQF load fasta and query" ) {

    map<uint64_t,int> gold;
    int kSize=31;
    vector<kDataFrame*> kframes;
    kframes.push_back(new kDataFrameBMQF(kSize,20,2,2,0));

    for(auto kframe: kframes){
      FastqReader reader("tests/testData/test.fastq");
      deque<pair<string,string> > sequences;
      while(!reader.isEOF())
      {
        reader.readNSeq(&sequences);
      }
      string kmer;
      for(auto seqPair:sequences){

          string seq=seqPair.first;
          for(int i=0;i<seq.size()-kSize;i++)
          {
            kmer=seq.substr(i,kSize);
            kframe->incrementCounter(kmer,1);
            uint64_t kmerHash=kframe->hashKmer(kmer);
            auto goldIT=gold.find(kmerHash);
            if(goldIT==gold.end())
              gold.insert(make_pair(kmerHash,1));
            else{
              goldIT->second++; // Comment this
            }
          }
      }

      for(auto seqPair:sequences){

          string seq=seqPair.first;
          for(int i=0;i<seq.size()-kSize;i++)
          {
            kmer=seq.substr(i,kSize);
            uint64_t kmerHash=kframe->hashKmer(kmer);
            auto goldIT=gold.find(kmerHash);
            REQUIRE(goldIT->second==kframe->getCounter(kmer));
          }
      }





      kframe->removeKmer(kmer);
      REQUIRE(kframe->getCounter(kmer)==0);


      kframe->setCounter(kmer,10);
      REQUIRE(kframe->getCounter(kmer)==10);

      kframe->setCounter(kmer,4);
      REQUIRE(kframe->getCounter(kmer)==4);
      kframe->setCounter(kmer,15);
      REQUIRE(kframe->getCounter(kmer)==15);



      kframe->setTag(kmer,1);
      REQUIRE(kframe->getTag(kmer)==1);
    }

}

TEST_CASE( "test iterator" ) {

    map<uint64_t,int> gold;
    int kSize=31;
    kDataFrameMQF kframe(kSize,20,2,2,0);
    FastqReader reader("tests/testData/test.fastq");
    deque<pair<string,string> > sequences;
    while(!reader.isEOF())
    {
      reader.readNSeq(&sequences);
    }
    string kmer;
    for(auto seqPair:sequences){

        string seq=seqPair.first;
        for(int i=0;i<seq.size()-kSize;i++)
        {
          kmer=seq.substr(i,kSize);
          kframe.incrementCounter(kmer,1);
          uint64_t kmerhash=kframe.hashKmer(kmer);
          auto goldIT=gold.find(kmerhash);
          if(goldIT==gold.end())
            gold.insert(make_pair(kmerhash,1));
          else{
            goldIT->second++;
          }
        }
    }


    auto kframeIt=kframe.begin();
    uint64_t numberOfKmers=0;
    while(!kframeIt.isEnd())
    {
      kmerRow res=*kframeIt;
      auto goldIT=gold.find(res.kmerHash);
      REQUIRE(res.count==goldIT->second);
      kframeIt++;
      numberOfKmers++;
    }
    REQUIRE(numberOfKmers==gold.size());



    kframe.removeKmer(kmer);
    REQUIRE(kframe.getCounter(kmer)==0);


    kframe.setCounter(kmer,10);
    REQUIRE(kframe.getCounter(kmer)==10);

    kframe.setCounter(kmer,4);
    REQUIRE(kframe.getCounter(kmer)==4);
    kframe.setCounter(kmer,15);
    REQUIRE(kframe.getCounter(kmer)==15);



    kframe.setTag(kmer,1);
    REQUIRE(kframe.getTag(kmer)==1);


}


TEST_CASE( "save and load" ) {

  map<uint64_t,int> gold;
  int kSize=31;
  vector<kDataFrameMQF*> kframes;
  kframes.push_back(new kDataFrameMQF(kSize,20,2,2,0));
  for(auto kframe: kframes){
    FastqReader reader("tests/testData/test.fastq");
    deque<pair<string,string> > sequences;
    while(!reader.isEOF())
    {
      reader.readNSeq(&sequences);
    }
    string kmer;
    for(auto seqPair:sequences){

        string seq=seqPair.first;
        for(int i=0;i<seq.size()-kSize;i++)
        {
          kmer=seq.substr(i,kSize);
          kframe->incrementCounter(kmer,1);
          uint64_t kmerHash=kframe->hashKmer(kmer);
          auto goldIT=gold.find(kmerHash);
          if(goldIT==gold.end())
            gold.insert(make_pair(kmerHash,1));
          else{
            goldIT->second++;
          }
        }
    }

    string filePath="tests/testData/tmp.kDataFrame";
    kframe->save(filePath);
    kDataFrame* kframe2=kDataFrame::load(filePath, "MQF");
    for(auto seqPair:sequences){

        string seq=seqPair.first;
        for(int i=0;i<seq.size()-kSize;i++)
        {
          kmer=seq.substr(i,kSize);
          uint64_t kmerHash=kframe->hashKmer(kmer);
          auto goldIT=gold.find(kmerHash);
          REQUIRE(goldIT->second==kframe->getCounter(kmer));
        }
    }





    kframe->removeKmer(kmer);
    REQUIRE(kframe->getCounter(kmer)==0);


    kframe->setCounter(kmer,10);
    REQUIRE(kframe->getCounter(kmer)==10);

    kframe->setCounter(kmer,4);
    REQUIRE(kframe->getCounter(kmer)==4);
    kframe->setCounter(kmer,15);
    REQUIRE(kframe->getCounter(kmer)==15);



    kframe->setTag(kmer,1);
    REQUIRE(kframe->getTag(kmer)==1);
  }



}


// Needs to be modified for the new change in indexMain.cpp
TEST_CASE( "indexing" ) {

    int kSize=31;

    FastqReader reader("tests/testData/test.fastq");
    deque<pair<string,string> > sequences;
    while(!reader.isEOF())
    {
      reader.readNSeq(&sequences);
    }
    string kmer;
    vector<kDataFrameMQF*> frames;
    for(auto seqPair:sequences){
        string seq=seqPair.first;
        std::vector<uint64_t> countHistogram(2);
        countHistogram[0]=seq.size()*3;
        countHistogram[1]=seq.size()*3;
        kDataFrameMQF* curr=new kDataFrameMQF(kSize,20,2,0,0);
        frames.push_back(curr);
        for(int i=0;i<seq.size()-kSize;i++)
        {
          kmer=seq.substr(i,kSize);
          curr->incrementCounter(kmer,1);
        }
    }

    kDataFrameMQF* indexFrame=kDataFrameMQF::index(frames);
    std::map<uint64_t, std::vector<int> > * legend=indexFrame->get_legend();


    // auto it2=indexFrame->begin();
    // while(!it2.isEnd()){
    //   cout<<(*it2).kmerHash<<" "<<(*it2).tag<<endl;
    //   it2++;
    // }
    for(int j=0;j<sequences.size();j++){
      string seq=sequences[j].first;
      for(int i=0;i<seq.size()-kSize;i++)
      {
        kmer=seq.substr(i,kSize);
        uint64_t tag=indexFrame->getTag(kmer);
        auto colors=legend->find(tag)->second;
        auto colorIt=find(colors.begin(),colors.end(),j);
        // REQUIRE(colorIt!=colors.end());
      }
    }

}
