#include "KmerCounter.hpp"
#include <iostream>
#include "kmer.h"
#include <fstream>

#include <seqan/seq_io.h>
#include "../HashUtils/hashutil.h"
#include <omp.h>

#include <gqf.h>
using namespace std;
using namespace seqan;
#define QBITS_LOCAL_QF 16

/* dump the contents of a local QF into the main QF */
static inline void dump_local_qf_to_main(QF* local, QF* main )
{
  #pragma omp critical
  {
    QFi local_cfi;
    if (qf_iterator(local, &local_cfi, 0)) {
      do {
        uint64_t key = 0, value = 0, count = 0;
        qfi_get(&local_cfi, &key, &value, &count);
        //qf_spin_lock((int*)&main_qf_lock,true);
        //main_qf_lock=false;
        qf_insert(main, key, count, true, true);
      } while (!qfi_next(&local_cfi));
      qf_reset(local);
    }
  }
}

static inline void insertToLevels(uint64_t item,QF* local,QF* main)
{
  if(!qf_insert(main, item%main->metadata->range, 1,
										 true, false)) {
				qf_insert(local, item%local->metadata->range, 1,
									false, false);
				// check of the load factor of the local QF is more than 50%
				if (qf_space(local)>95) {
					dump_local_qf_to_main(local,main);
				}
			}


}


void loadIntoMQF(string sequenceFilename,int ksize,int noThreads, Hasher *hasher,QF * memoryMQF){
  cout<<noThreads<<endl;
  SeqFileIn seqFileIn(sequenceFilename.c_str());
  StringSet<CharString> ids;
  StringSet<Dna5String> reads;
  omp_set_num_threads(noThreads);
  QF* localMQF;
  #pragma omp parallel private(ids,reads,localMQF) shared(seqFileIn)
  {
    localMQF= new QF();
    qf_init(localMQF, (1ULL << QBITS_LOCAL_QF), memoryMQF->metadata->key_bits,
    0,memoryMQF->metadata->fixed_counter_size, true,"", 2038074761);


    while(!atEnd(seqFileIn))
    {
      #pragma omp critical
      {
        seqan::readRecords(ids, reads, seqFileIn,10000);
      }

      for(auto read:reads){
        if(length(read)<ksize)
        {
          continue;
        }

        uint64_t first = 0;
        uint64_t first_rev = 0;
        uint64_t item = 0;

        for(int i=0; i<ksize; i++) {
          //First kmer
          uint8_t curr = kmer::map_base(read[i]);
          if (curr > DNA_MAP::G) {
            // 'N' is encountered
            //    read = read.substr(i+1, length(read));
            continue;
            //goto start_read;
          }
          first = first | curr;
          first = first << 2;
        }
        first = first >> 2;
        first_rev = kmer::reverse_complement(first, ksize);


        if (kmer::compare_kmers(first, first_rev))
        item = first;
        else
        item = first_rev;
        item = hasher->hash(item)%memoryMQF->metadata->range;
        insertToLevels(item,localMQF,memoryMQF);

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


          item = hasher->hash(item)%memoryMQF->metadata->range;
          insertToLevels(item,localMQF,memoryMQF);

          next = (next << 2) & BITMASK(2*ksize);
          next_rev = next_rev >> 2;
        }
      }

    }
    dump_local_qf_to_main(localMQF,memoryMQF);
    qf_destroy(localMQF);
  }

}

void dumpMQF(QF * MQF,int ksize,std::string outputFilename){
  IntegerHasher Ihasher(BITMASK(2*ksize));
  ofstream output(outputFilename.c_str());
  QFi qfi;
  qf_iterator(MQF, &qfi, 0);
  do {
    uint64_t key, value, count;
    qfi_get(&qfi, &key, &value, &count);
    string kmer=kmer::int_to_str(Ihasher.Ihash(key),ksize);
    output<<kmer<<"\t"<<count<<endl;
  } while(!qfi_next(&qfi));
}
