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
#include <deque>
#include <gqf.h>
#include <stdio.h>
#include <string.h>
using namespace std;
//using namespace seqan;
#define QBITS_LOCAL_QF 16

const uint64_t maxLocal=(uint64_t)((double)(1ULL<<16)*0.5);


static inline void insertToLevels2(uint64_t item,QF* local,QF* main,uint64_t *local_capacity)
{
  if(!qf_insert(main, item, 1,true, false)) {
    qf_insert(local, item, 1,false, false);
    *local_capacity++;
    if (*local_capacity>maxLocal) {
      qf_migrate(local,main);
      *local_capacity=0;
      qf_reset(local);
    }
  }
}
//
static inline void insertToLevels(uint64_t item,QF* local,QF* main,QF * diskMQF=NULL)
{
  if(!qf_insert(main, item, 1,true, false)) {
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
              throw overflow_error("memory MQF doesn't have enough space to migrate");
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

std::vector<std::string> string_split(std::string s, const char delimiter)
{
    size_t start=0;
    size_t end=s.find_first_of(delimiter);

    std::vector<std::string> output;

    while (end <= std::string::npos)
    {
	    output.emplace_back(s.substr(start, end-start));

	    if (end == std::string::npos)
	    	break;

    	start=end+1;
    	end = s.find_first_of(delimiter, start);
    }

    return output;
}


void loadIntoMQF(string sequenceFilename,int ksize,int noThreads, Hasher *hasher,QF * memoryMQF,QF * diskMQF){

  FastqReaderSqueker reader(sequenceFilename);
  omp_set_num_threads(noThreads);
  QF* localMQF;
  bool moreWork=true;
  uint64_t numReads=0;
  deque<pair<string,string> > reads;


  string read,tag;
#pragma omp parallel private(reads,localMQF,read,tag) shared(reader,moreWork,numReads)  firstprivate(ksize,noThreads,memoryMQF,diskMQF)
  {
    uint32_t OVERHEAD_SIZE = 65535;
    uint64_t part_size = 1ULL << 23;
    char* chunk= (char *)malloc((part_size + OVERHEAD_SIZE)*sizeof(char));
    uint64_t chunkSize;
    auto localHasher=hasher->clone();
    localMQF= new QF();
    uint64_t local_capacity=0;
    reads=deque<pair<string,string> >(15000);
    qf_init(localMQF, (1ULL << QBITS_LOCAL_QF), memoryMQF->metadata->key_bits,
    0,memoryMQF->metadata->fixed_counter_size, true,"", 2038074761);

    uint64_t  *LocalkmersBuffer,LocalkmersBufferSize=150,LocalkmersBufferTop=0;
    LocalkmersBuffer=new uint64_t[LocalkmersBufferSize];

    while(moreWork)
    {
      #pragma omp critical
      {
        //reader.readNSeq(&reads,15000);
        //numReads+=15000;
        reader.readChunk(chunk,&chunkSize);
        bool tmp=!reader.isEOF();
        moreWork=tmp;
      }
      auto fs = chunk;
      auto fe = chunk;
      auto end = chunk + chunkSize;
      while (fs && fs!=end){
        fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore the first line
        fs++; // increment the pointer

        fe = static_cast<char*>(memchr(fs, '\n', end-fs)); // read the read
        read=string(fs, fe-fs);
start_read:
        if(read.size()<ksize)
        {
	  fs = ++fe;		// increment the pointer
	  fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one line
	  fs++; // increment the pointer
	  fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one more line
	  fs++; // increment the pointer
          continue;
        }
        LocalkmersBufferTop=0;
        if(read.size()>LocalkmersBufferSize)
        {
          LocalkmersBufferSize=read.size();
          delete [] LocalkmersBuffer;
          LocalkmersBuffer=new uint64_t[LocalkmersBufferSize];
        }
<<<<<<< HEAD

        vector<string> readParts=string_split(read,,'N');
        for(auto readPart:readParts){
          if(readPart.size()<k)
            continue;
          uint64_t first = 0;
          uint64_t first_rev = 0;
          uint64_t item = 0;

          for(int i=0; i<ksize; i++) {
            //First kmer
            uint8_t curr = kmer::map_base_vectorized(readPart[i]);
            first = first | curr;
            first = first << 2;
=======
        first = first >> 2;
        first_rev = kmer::reverse_complement(first, ksize);


        if (kmer::compare_kmers(first, first_rev))
        item = first;
        else
        item = first_rev;

        item = localHasher->hash(item)%memoryMQF->metadata->range;
        insertToLevels2(item,localMQF,memoryMQF,&local_capacity);
	//  insertToLevels(item,localMQF,memoryMQF,diskMQF);
        uint64_t next = (first << 2) & BITMASK(2*ksize);
        uint64_t next_rev = first_rev >> 2;

        for(uint32_t i=ksize; i<(read.size()); i++) {
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
>>>>>>> b6dfe42318ca828b0183a66290abba37924469b8
          }
          first = first >> 2;
          first_rev = kmer::reverse_complement(first, ksize);
          LocalkmersBuffer[LocalkmersBufferTop++]=kmer::smallerKmer(first,first_rev);


<<<<<<< HEAD
        //  item = localHasher->hash(item)%memoryMQF->metadata->range;
          //insertToLevels2(item,localMQF,memoryMQF,&local_capacity);
          insertToLevels(item,localMQF,memoryMQF,diskMQF);
          uint64_t next = (first << 2) & BITMASK(2*ksize);
          uint64_t next_rev = first_rev >> 2;

          for(uint32_t i=ksize; i<(read.size()); i++) {
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
            insertToLevels(item,localMQF,memoryMQF,diskMQF);
            //insertToLevels2(item,localMQF,memoryMQF,&local_capacity);
            next = (next << 2) & BITMASK(2*ksize);
            next_rev = next_rev >> 2;
          }
=======
          item = localHasher->hash(item)%memoryMQF->metadata->range;
	  // insertToLevels(item,localMQF,memoryMQF,diskMQF);
          insertToLevels2(item,localMQF,memoryMQF,&local_capacity);
          next = (next << 2) & BITMASK(2*ksize);
          next_rev = next_rev >> 2;
>>>>>>> b6dfe42318ca828b0183a66290abba37924469b8
        }
        fs = ++fe;		// increment the pointer
        fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one line
        fs++; // increment the pointer
        fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one more line
        fs++; // increment the pointer
      }

    }
    //    #pragma omp critical
    {
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


void estimateMemRequirement_2Structures(std::string ntcardFilename,
  uint64_t numHashBits,uint64_t tagSize,
  uint64_t *res_SingleQ,uint64_t *res_noSlots,uint64_t *res_fixedSizeCounter, uint64_t *res_memory)
{
  uint64_t noDistinctKmers=0,totalNumKmers=0;
  vector<uint64_t> histogram(1000,0);
  ifstream ntcardFile(ntcardFilename);
  string f;
  uint64_t count;
  uint64_t singletons;
  while(ntcardFile>>f>>count)
  {
    if(count==numeric_limits<uint64_t>::max())
      continue;
    if(f=="F0")
      noDistinctKmers=count;
    else if(f=="F1")
      totalNumKmers=count;
    else if(f=="f1"){
      singletons=count;
    }
      //else if(f=="f2")
      //cout<<"two "<<count<<endl;
    else{
      f=f.substr(1,f.size());
      int n=atoi(f.c_str());
      histogram[n]=count;
    }
  }

  noDistinctKmers-=singletons;

  uint64_t singletonsQbits=(uint64_t)log2((double)singletons)+1ULL;
  *res_SingleQ=singletonsQbits;
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
  *res_memory+=estimateMemory((1ULL<<singletonsQbits),numHashBits-singletonsQbits,1,tagSize);
  if(*res_memory==numeric_limits<uint64_t>::max())
  {
    throw std::overflow_error("Data limits exceeds MQF capabilities(> uint64). Check if ntcard file is corrupted");
  }


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
