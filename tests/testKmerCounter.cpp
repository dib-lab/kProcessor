#include "../catch.hpp"
#include "../KmerCounter/KmerCounter.hpp"
#include <gqf.h>
#include <stdint.h>
#include <fstream>
#include <map>
#include <seqan/modifier.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <iostream>
using namespace std;
TEST_CASE( "loadIntoMQF", "[KmerCounter]" ) {
    string testfilename="tests/testData/test.noN.fastq";
    string tmpOutput="tests/testData/tmp.test.noN.count";
    int k=15;
    int noThreads=1;
    QF qf;
    uint64_t nslots=32768;
    uint64_t fixed_size_counter=1;
    qf_init(&qf, nslots, 2*k+15, 0,fixed_size_counter, true, "", 2038074761);
    loadIntoMQF(testfilename, k, noThreads,&qf);
    dumpMQF(&qf,k,tmpOutput);

    map<string,int> gold;
    ifstream goldFile("tests/testData/test.noN.dsk.txt");
    string kmer;
    int count;
    int tmp=0;
    while(goldFile>>kmer>>count){
      gold.insert(make_pair(kmer,count));
      tmp++;
    }
    ifstream result(tmpOutput.c_str());
    int resultKmerNumber=0;
    while(result>>kmer>>count)
    {
      resultKmerNumber+=1;
      auto it=gold.find(kmer);
      if(it==gold.end())
      {
        seqan::Dna5String revKmer(kmer);
        seqan::reverseComplement(revKmer);
        string revKmer2;
        assign(revKmer2, revKmer);
        it=gold.find(revKmer2);
        INFO("Reverse Complement")
      }
      INFO("Checking count of"<<it->first);
      REQUIRE(count==it->second);
    }
    REQUIRE(resultKmerNumber==gold.size());


}
