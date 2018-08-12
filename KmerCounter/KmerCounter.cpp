#include "KmerCounter.hpp"
#include <iostream>
#include "kmer.h"
#include <fstream>

#include <seqan/seq_io.h>
#include "../HashUtils/hashutil.h"
#include <seqan/parallel.h>
#include "../KmerDecoder/FastqReader.hpp"
#include <limits>
#include <omp.h>
#include <stdexcept>
#include <math.h>
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <chrono>
#include <thread>
#include"../ThirdParty/spsc.hpp"

#include <gqf.h>
using namespace std;
using namespace seqan;
#define numElementsInBuffer 1000000
#define QBITS_LOCAL_QF 16
#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))
uint64_t bitmaskLookup2[]={0,1, 3, 7, 15, 31, 63, 127, 255, 511, 1023,
  2047, 4095, 8191, 16383, 32767, 65535, 131071, 262143, 524287, 1048575,
  2097151, 4194303, 8388607, 16777215, 33554431, 67108863, 134217727, 268435455, 536870911, 1073741823,
  2147483647, 4294967295, 8589934591, 17179869183, 34359738367, 68719476735, 137438953471, 274877906943, 549755813887, 1099511627775,
  2199023255551, 4398046511103, 8796093022207, 17592186044415, 35184372088831, 70368744177663, 140737488355327, 281474976710655, 562949953421311, 1125899906842623,
  2251799813685247, 4503599627370495, 9007199254740991, 18014398509481983, 36028797018963967, 72057594037927935, 144115188075855871, 288230376151711743, 576460752303423487, 1152921504606846975,
  2305843009213693951, 4611686018427387903, 9223372036854775807, 18446744073709551615
};

//#define BITMASK(nbits) ((nbits) == 64 ? 0xffffffffffffffff : (1ULL << (nbits)) \
- 1ULL)
#define BITMASK(nbits) bitmaskLookup2[nbits]


static inline void insertToLevels(uint64_t item,QF* local,QF* main,QF * diskMQF=NULL)
{
  if(!qf_insert(main, item%main->metadata->range, 1,true, false)) {
    qf_insert(local, item%local->metadata->range, 1,false, false);
    // check of the load factor of the local QF is more than 50%
    if (qf_space(local)>90) {
      if(main->metadata->noccupied_slots+local->metadata->noccupied_slots
        < main->metadata->maximum_occupied_slots){
          qf_migrate(local,main);
        }
        else if(diskMQF!=NULL){
          SEQAN_OMP_PRAGMA(critical){
            qf_general_lock(main,true);
            qf_migrate(main,diskMQF);
            qf_reset(main);
            qf_general_unlock(main);
            qf_migrate(local,main);
          }
        }
        else{
          throw overflow_error("memory MQF doesn't have enough space");
        }
        qf_reset(local);
      }
    }
    else{
      if (qf_space(main)>90) {
        SEQAN_OMP_PRAGMA(critical)
        {
          if (qf_space(main)>90) {
            if(diskMQF!=NULL){
              qf_general_lock(main,true);
              qf_migrate(main,diskMQF);
              qf_reset(main);
              qf_general_unlock(main);
            }else{
              throw overflow_error("memory MQF doesn't have enough space");
            }
          }
        }
      }
    }
  }


  void loadIntoMQF(string sequenceFilename,int ksize,int noThreads, Hasher *hasher,QF * memoryMQF,QF * diskMQF){
    FastqReader reader(sequenceFilename);
    omp_set_num_threads(noThreads);
    const int numBuffers=500;

    boost::lockfree::queue<uint64_t*> freeBuffers(numBuffers);
    for(int i=0;i<numBuffers;i++)
      freeBuffers.push(new uint64_t[numElementsInBuffer]);


    SpScLockFreeQueue<pair<uint64_t*,int>,1000>** kmersBuffer;
    kmersBuffer=new SpScLockFreeQueue<pair<uint64_t*,int>,1000>*[noThreads];
    for(int i=0;i<noThreads;i++)
      kmersBuffer[i]=new SpScLockFreeQueue<pair<uint64_t*,int>,1000>();

    QF* localMQF;
    bool moreWork=true;
    uint64_t numReads=0;
    vector<pair<string,string> > reads;
    string read,tag;
    #pragma omp parallel private(reads,localMQF,read,tag) shared(reader,moreWork,numReads,kmersBuffer)
    {
      uint64_t threadBits=LOG2(noThreads);
      uint64_t totalBits=LOG2(memoryMQF->metadata->range);
      uint64_t threadId= omp_get_thread_num();
      uint64_t item;
      localMQF= new QF();
      qf_init(localMQF, (1ULL<<QBITS_LOCAL_QF) , memoryMQF->metadata->key_bits,
        0,memoryMQF->metadata->fixed_counter_size, true,"", 2038074761);
      if(threadId==0)//parser thread
      {
        int bufferTops[noThreads];
        uint64_t* buffers[noThreads];
        for(int i=0;i<noThreads;i++){
          bufferTops[i]=0;
          freeBuffers.pop(buffers[i]);
        }

        reads=vector<pair<string,string> >(100000);
        while(moreWork)
        {
          reader.readNSeq(&reads,10000);
          numReads+=10000;
          //std::cout << numReads << endl;
          for(int j=0;j<reads.size();j++){
            read=reads[j].first;
            start_read:
            if(read.size()<ksize)
            {
              continue;
            }
            uint64_t first = 0;
            uint64_t first_rev = 0;
            item = 0;

            for(int i=0; i<ksize; i++) {
              //First kmer
              uint8_t curr = kmer::map_base(read[i]);
              if (curr > DNA_MAP::G) {
                // 'N' is encountered

                read=read.substr(i+1, read.length());

                //continue;
                goto start_read;
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
            // 1) split the bits of the items
            // 2) push to the queue
            // 3)
            uint64_t workerID = item >> (totalBits-threadBits);
            //uint64_t newItem = item & BITMASK(totalBits-threadBits);
            //cout<<"item = "<<item<<" workerID= "<<workerID<<" new item= "<<newItem<<endl;
            if(workerID==0){
              insertToLevels(item,localMQF,memoryMQF,diskMQF);
            }else{
              buffers[workerID][bufferTops[workerID]++]=item;
              if(bufferTops[workerID]==numElementsInBuffer)
              {
                while(!kmersBuffer[workerID]->push(make_pair(buffers[workerID],bufferTops[workerID])));
                bufferTops[workerID]=0;
                while(!freeBuffers.pop(buffers[workerID]));
              }
            }

            //insertToLevels(item,localMQF,memoryMQF,diskMQF);

            uint64_t next = (first << 2) & BITMASK(2*ksize);
            uint64_t next_rev = first_rev >> 2;

            for(uint32_t i=ksize; i<length(read); i++) {
              //next kmers
              //cout << "K: " << read.substr(i-K+1,K) << endl;
              uint8_t curr = kmer::map_base(read[i]);
              if (curr > DNA_MAP::G) {
                // 'N' is encountered
                //continue;
                //read = read.substr(i+1, length(read));
                read=read.substr(i+1, read.length());
                //erase(read,0,i+1);

                goto start_read;
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
              //insertToLevels(item,localMQF,memoryMQF,diskMQF);
              workerID = item >> (totalBits-threadBits);
              //newItem = item & BITMASK(totalBits-threadBits);
              if(workerID==0){
                insertToLevels(item,localMQF,memoryMQF,diskMQF);
              }
              else{
                buffers[workerID][bufferTops[workerID]++]=item;
                if(bufferTops[workerID]==numElementsInBuffer)
                {
                  while(!kmersBuffer[workerID]->push(make_pair(buffers[workerID],bufferTops[workerID])));
                  bufferTops[workerID]=0;
                  while(!freeBuffers.pop(buffers[workerID]));
                }
              }
              next = (next << 2) & BITMASK(2*ksize);
              next_rev = next_rev >> 2;
            }
          }

          bool tmp=!reader.isEOF();
          if(!tmp){
            for(int i=1;i<noThreads;i++)
            {
                while(!kmersBuffer[i]->push(make_pair(buffers[i],bufferTops[i])));
            }

          }
          moreWork=tmp;

        }
        qf_migrate(localMQF,memoryMQF);
        qf_destroy(localMQF);

      }
      else{
      // worker thread
      uint64_t *local_buffer;
      int top=0;
      pair<uint64_t*,int> buffer_pair;
      while(moreWork)
      {
        if(kmersBuffer[threadId]->pop(buffer_pair))
        {
          local_buffer=buffer_pair.first;
          for(int i=0;i<buffer_pair.second;i++){
          insertToLevels(local_buffer[i],localMQF,memoryMQF,diskMQF);
          }
          while(!freeBuffers.push(local_buffer));
        }
        else
        {
          std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
      }
      while(kmersBuffer[threadId]->pop(buffer_pair))
        {
          local_buffer=buffer_pair.first;
          for(int i=0;i<buffer_pair.second;i++){
            insertToLevels(local_buffer[i],localMQF,memoryMQF,diskMQF);
          }
        }
        qf_migrate(localMQF,memoryMQF);
        qf_destroy(localMQF);
    }

    }
      if(diskMQF!=NULL){
        qf_migrate(memoryMQF,diskMQF);
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
        output<<kmer<<" "<<count<<"\n";
      } while(!qfi_next(&qfi));
    }

    bool isEnough(vector<uint64_t> histogram,uint64_t noSlots,uint64_t fixedSizeCounter,uint64_t slotSize)
    {
      // cout<<"noSlots= "<<noSlots<<endl
      //     <<"fcounter= "<<fixedSizeCounter<<endl
      //     <<"slot size= "<<numHashBits<<endl;

      noSlots=(uint64_t)((double)noSlots*0.90);
      for(uint64_t i=1;i<1000;i++)
      {
        uint64_t usedSlots=1;

        if(i>((1ULL)<<fixedSizeCounter)-1)
        {
          uint64_t nSlots2=0;
          __uint128_t capacity;
          do{
            nSlots2++;
            capacity=((__uint128_t)(1ULL)<<(nSlots2*slotSize+fixedSizeCounter))-1;
            //  cout<<"slots num "<<nSlots2<<" "<<capacity<<endl;
          }while((__uint128_t)i>capacity);
          usedSlots+=nSlots2;
        }
        //cout<<"i= "<<i<<"->"<<usedSlots<<" * "<<histogram[i]<<endl;
        if(noSlots>=(usedSlots*histogram[i]))
        {
          noSlots-=(usedSlots*histogram[i]);
        }
        else
        {
          //  cout<<"failed"<<endl<<endl;
          return false;
        }

      }
      //cout<<"success"<<endl<<endl;
      return true;
    }


    void estimateMemRequirement(std::string ntcardFilename,
      uint64_t numHashBits,uint64_t tagSize,
      uint64_t *res_noSlots,uint64_t *res_fixedSizeCounter, uint64_t *res_memory)
      {
        uint64_t noDistinctKmers=0,totalNumKmers=0;
        vector<uint64_t> histogram(1000,0);
        ifstream ntcardFile(ntcardFilename);
        string f;
        uint64_t count;
        while(ntcardFile>>f>>count)
        {
          if(count==numeric_limits<uint64_t>::max())
          continue;
          if(f=="F0")
          noDistinctKmers=count;
          else if(f=="F1")
          totalNumKmers=count;
          else{
            f=f.substr(1,f.size());
            int n=atoi(f.c_str());
            histogram[n]=count;
          }
        }
        *res_memory=numeric_limits<uint64_t>::max();
        for(int i=8;i<64;i++)
        {
          uint64_t noSlots=(1ULL)<<i;
          if(noSlots<noDistinctKmers)
          continue;
          bool moreWork=false;
          uint64_t slotSize=numHashBits-log2((double)noSlots);
          for(uint64_t fixedSizeCounter=1;fixedSizeCounter<slotSize;fixedSizeCounter++)
          {
            if(isEnough(histogram,noSlots,fixedSizeCounter,slotSize))
            {
              uint64_t tmpMem=estimateMemory(noSlots,slotSize,fixedSizeCounter,tagSize);
              if(*res_memory>tmpMem)
              {
                *res_memory=tmpMem;
                *res_fixedSizeCounter=fixedSizeCounter;
                *res_noSlots=noSlots;
                moreWork=true;
              }
              else{
                break;
              }
            }

          }
          if(!moreWork && *res_memory!=numeric_limits<uint64_t>::max())
          break;
        }
        if(*res_memory==numeric_limits<uint64_t>::max())
        {
          throw std::overflow_error("Data limits exceeds MQF capabilities(> uint64). Check if ntcard file is corrupted");
        }


      }
