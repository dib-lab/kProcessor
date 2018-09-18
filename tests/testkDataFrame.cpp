#include "../ThirdParty/catch.hpp"
#include "../KmerCounter/KmerCounter.hpp"
#include <gqf.h>
#include <stdint.h>
#include <fstream>
#include <map>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <iostream>
#include "../kDataFrame.hpp"
#include "../KmerDecoder/FastqReader.hpp"
using namespace std;


TEST_CASE( "load fasta and query" ) {

    map<string,int> gold;
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
          auto goldIT=gold.find(kmer);
          if(goldIT==gold.end())
            gold.insert(make_pair(kmer,1));
          else{
            goldIT->second++;
          }
        }
    }

    auto goldIT=gold.begin();
    while(goldIT!=gold.end())
    {
      int count=kframe.getCounter(goldIT->first);
      REQUIRE(count==goldIT->second);
      goldIT++;
    }



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

    map<string,int> gold;
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
          auto goldIT=gold.find(kmer);
          if(goldIT==gold.end())
            gold.insert(make_pair(kmer,1));
          else{
            goldIT->second++;
          }
        }
    }
    string filePath="tests/testData/tmp.kDataFrame";
    kframe.save(filePath);

    kDataFrame* kframe2=kDataFrame::load(filePath);
    auto goldIT=gold.begin();
    while(goldIT!=gold.end())
    {
      int count=kframe2->getCounter(goldIT->first);
      REQUIRE(count==goldIT->second);
      goldIT++;
    }



  

}
