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
using namespace std;
void KmerCounter_main(int argc, char *argv[]);


TEST_CASE( "loadIntoMQF", "[KmerCounter]" ) {
    char* argv[]={
      (char*)"count",
      (char*)"-i", (char*)"tests/testData/test.noN.fastq",
      (char*)"-k", (char*)"15",
      (char*)"-o", (char*)"tests/testData/tmp.test.noN.mqf",
      (char*)"-s", (char*)"32768",
      (char*)"-t", (char*)"4",
      (char*)"-u", (char*)"tests/testData/tmp.test.noN.count"
    };
    int argc=13;


    KmerCounter_main(argc, argv);

    map<string,int> gold;
    ifstream goldFile("tests/testData/test.noN.dsk.txt");
    string kmer;
    int count;
    int tmp=0;
    while(goldFile>>kmer>>count){
      gold.insert(make_pair(kmer,count));
      tmp++;
    }
    ifstream result("tests/testData/tmp.test.noN.count");
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
