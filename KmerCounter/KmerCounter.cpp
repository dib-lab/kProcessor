#include "KmerCounter.hpp"
#include <iostream>
#include "kmer.h"
#include <fstream>

#include <seqan/seq_io.h>
#include "../HashUtils/hashutil.h"

#include <gqf.h>
using namespace std;
using namespace seqan;


void loadIntoMQF(string sequenceFilename,int ksize,int noThreads,QF * memoryMQF){
  SeqFileIn seqFileIn(sequenceFilename.c_str());
  CharString id;
  Dna5String read;
  while(!atEnd(seqFileIn))
  {
    seqan::readRecord(id, read, seqFileIn);
    if(length(read)<ksize)
    {
      continue;
    }

    uint64_t first = 0;
    uint64_t first_rev = 0;
    uint64_t item = 0;
    //cout << "K " << read.substr(0,K) << endl;
    for(int i=0; i<ksize; i++) {
      //First kmer
      uint8_t curr = kmer::map_base(read[i]);
      if (curr > DNA_MAP::G) {
        // 'N' is encountered
        //  read = read.substr(i+1, length(read));
        continue;
        //goto start_read;
      }
      first = first | curr;
      first = first << 2;
    }
    first = first >> 2;
    first_rev = kmer::reverse_complement(first, ksize);

    //cout << "kmer: "; cout << int_to_str(first);
    //cout << " reverse-comp: "; cout << int_to_str(first_rev) << endl;

    if (kmer::compare_kmers(first, first_rev))
    item = first;
    else
    item = first_rev;
    item = HashUtil::hash_64(item, BITMASK(2*ksize));

    qf_insert(memoryMQF,item,1);

    uint64_t next = (first << 2) & BITMASK(2*ksize);
    uint64_t next_rev = first_rev >> 2;

    for(uint32_t i=ksize; i<length(read); i++) {
      //next kmers
      //cout << "K: " << read.substr(i-K+1,K) << endl;
      uint8_t curr = kmer::map_base(read[i]);
      if (curr > DNA_MAP::G) {
        // 'N' is encountered
        continue;
        //read = read.substr(i+1, length(read));
      }
      next |= curr;
      uint64_t tmp = kmer::reverse_complement_base(curr);
      tmp <<= (ksize*2-2);
      next_rev = next_rev | tmp;
      if (kmer::compare_kmers(next, next_rev))
      item = next;
      else
      item = next_rev;


      item = HashUtil::hash_64(item, BITMASK(2*ksize));
      qf_insert(memoryMQF,item,1);

      next = (next << 2) & BITMASK(2*ksize);
      next_rev = next_rev >> 2;
    }


  }

}

void dumpMQF(QF * MQF,int ksize,std::string outputFilename){
  ofstream output(outputFilename.c_str());
  QFi qfi;
  qf_iterator(MQF, &qfi, 0);
  do {
    uint64_t key, value, count;
    qfi_get(&qfi, &key, &value, &count);
    string kmer=int_to_str(HashUtil::hash_64i(key,BITMASK(2*ksize)),ksize);
    output<<kmer<<"\t"<<count<<endl;
  } while(!qfi_next(&qfi));
}
