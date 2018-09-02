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
#include <deque>
#include <gqf.h>

#include <boost/fusion/support/pair.hpp>
#include <boost/fusion/include/pair.hpp>
using namespace std;
using namespace seqan;
#define QBITS_LOCAL_QF 16
#define numElementsInBuffer 1000000
#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

static inline void insertToLevels(uint64_t item,QF* local,QF* main,QF * diskMQF=NULL)
{
  if(!qf_insert(main, item, 1,
										 true, false)) {
				qf_insert(local, item, 1,
									false, false);
				// check of the load factor of the local QF is more than 50%
				if (qf_space(local)>50) {

          {
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

class bufferPair
{
public:
  uint64_t* buff;
  int length;

    bufferPair(){
      buff=0;
      length=0;
    }
    bufferPair(const bufferPair& rhs) { buff=rhs.buff;
      length=rhs.length;
     }
     bufferPair(uint64_t* inbuff,int inlength) {
       buff=inbuff;
       length=inlength;
      }
    //Implicitly defined destructor for itself and all member variables
    //Implicitly defined operator= for itself and all member variables

};

void loadIntoMQF(string sequenceFilename,int ksize,int noThreads, Hasher *hasher,QF * memoryMQF,QF * diskMQF){
  FastqReaderSqueker reader(sequenceFilename);
  omp_set_num_threads(noThreads);
  QF* localMQF;
  bool moreWork=true;
  uint64_t numReads=0;
  deque<pair<string,string> > reads;
  string read,tag;
  const int numBuffers=noThreads*noThreads*3;
  int threadsToJoin=0;

  boost::lockfree::queue<uint64_t*> freeBuffers(numBuffers);
  for(int i=0;i<numBuffers;i++)
    freeBuffers.push(new uint64_t[numElementsInBuffer]);

  boost::lockfree::queue<bufferPair >** kmersBuffer;
  kmersBuffer=new boost::lockfree::queue<bufferPair >*[noThreads];
  for(int i=0;i<noThreads;i++)
    kmersBuffer[i]=new boost::lockfree::queue<bufferPair >(noThreads*4);

#pragma omp parallel private(reads,localMQF,read,tag) shared(reader,moreWork,numReads)  firstprivate(ksize,noThreads,memoryMQF,diskMQF)
  {
    int threadNum=omp_get_thread_num();
    uint64_t threadBits=LOG2(noThreads);
    uint64_t totalBits=LOG2(memoryMQF->metadata->range);

    uint64_t workerID;
    auto localHasher=hasher->clone();
    localMQF= new QF();

    reads=deque<pair<string,string> >(15000);
    qf_init(localMQF, (1ULL << QBITS_LOCAL_QF), memoryMQF->metadata->key_bits,
    0,memoryMQF->metadata->fixed_counter_size, true,"", 2038074761);
    bufferPair tmp_buffer;
    int bufferTops[noThreads];
    uint64_t* buffers[noThreads];
    for(int i=0;i<noThreads;i++){
      bufferTops[i]=0;
      freeBuffers.pop(buffers[i]);
    }



    while(moreWork)
    {
      SEQAN_OMP_PRAGMA(critical)
      {
        reader.readNSeq(&reads,15000);
        numReads+=15000;
        bool tmp=!reader.isEOF();
        moreWork=tmp;
        //cout<<threadNum<<" was here"<<endl;
      }

      for(int j=0;j<reads.size();j++){

        while(kmersBuffer[threadNum]->pop(tmp_buffer))
        {
          for(int i=0;i<tmp_buffer.length;i++)
            insertToLevels(tmp_buffer.buff[i],localMQF,memoryMQF,diskMQF);
        }
        read=reads[j].first;
start_read:
        if(read.size()<ksize)
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

        item = localHasher->hash(item)%memoryMQF->metadata->range;
        workerID = item >> (totalBits-threadBits);
        if(workerID==threadNum){
             insertToLevels(item,localMQF,memoryMQF,diskMQF);
        }else{
            buffers[workerID][bufferTops[workerID]++]=item;
            if(bufferTops[workerID]==numElementsInBuffer)
             {
               while(!kmersBuffer[workerID]->push(bufferPair(buffers[workerID],bufferTops[workerID])));
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


          item = localHasher->hash(item)%memoryMQF->metadata->range;
          workerID = item >> (totalBits-threadBits);
          if(workerID==threadNum){
               insertToLevels(item,localMQF,memoryMQF,diskMQF);
          }else{
              buffers[workerID][bufferTops[workerID]++]=item;
              if(bufferTops[workerID]==numElementsInBuffer)
               {
                 while(!kmersBuffer[workerID]->push(bufferPair(buffers[workerID],bufferTops[workerID])));
                 bufferTops[workerID]=0;
                 while(!freeBuffers.pop(buffers[workerID]));
               }
          }
          //insertToLevels(item,localMQF,memoryMQF,diskMQF);
          next = (next << 2) & BITMASK(2*ksize);
          next_rev = next_rev >> 2;
        }
      }

    }
    cout<<threadNum<<"about to finish"<<endl;
    #pragma omp atomic
      threadsToJoin++;
    while(threadsToJoin<noThreads)
      {
        while(kmersBuffer[threadNum]->pop(tmp_buffer))
        {
          for(int i=0;i<tmp_buffer.length;i++)
            insertToLevels(tmp_buffer.buff[i],localMQF,memoryMQF,diskMQF);
        }
      }
    #pragma omp barier
    for(int i=0;i<noThreads;i++){
      kmersBuffer[i]->push(bufferPair(buffers[i],bufferTops[i]));
    }
    #pragma omp barier
    //    #pragma omp critical
    {
      while(kmersBuffer[threadNum]->pop(tmp_buffer))
      {
        for(int i=0;i<tmp_buffer.length;i++)
          insertToLevels(tmp_buffer.buff[i],localMQF,memoryMQF,diskMQF);
      }

      qf_migrate(localMQF,memoryMQF);
    }
    qf_destroy(localMQF);
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
