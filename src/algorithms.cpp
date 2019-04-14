#include "algorithms.hpp"
#include <iostream>
#include "Utils/kmer.h"
#include <fstream>

#include <seqan/seq_io.h>
#include "HashUtils/hashutil.h"
#include <seqan/parallel.h>
#include "KmerDecoder/FastqReader.hpp"
#include <limits>
#include <omp.h>
#include <stdexcept>
#include <math.h>
#include <deque>
#include <gqf.hpp>
#include <queue>
#include <functional>
#include <limits>
using namespace std;
using namespace seqan;
#define QBITS_LOCAL_QF 16

namespace kProcessor{
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
    auto localHasher=hasher->clone();
    localMQF= new QF();
    reads=deque<pair<string,string> >(15000);
    qf_init(localMQF, (1ULL << QBITS_LOCAL_QF), memoryMQF->metadata->key_bits,
    0,memoryMQF->metadata->fixed_counter_size,0 ,true,"", 2038074761);
    while(moreWork)
    {
      SEQAN_OMP_PRAGMA(critical)
      {
        reader.readNSeq(&reads,15000);
        numReads+=15000;
        bool tmp=!reader.isEOF();
        moreWork=tmp;
      }

      for(int j=0;j<reads.size();j++){
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
        insertToLevels(item,localMQF,memoryMQF,diskMQF);

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
          insertToLevels(item,localMQF,memoryMQF,diskMQF);
          next = (next << 2) & BITMASK(2*ksize);
          next_rev = next_rev >> 2;
        }
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
    string kmer=Ihasher.Ihash(key);
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

kDataFrame* transform(kDataFrame* input,kmerRow (*fn)(kmerRow it)){
  kDataFrame* res=input->getTwin();
  kDataFrameIterator it=input->begin();
  while(it!=input->end())
  {
    kmerRow newkmer=fn(*it);
    res->insert(newkmer);
    it++;
  }
  return res;

}


void parseSequences(string seqFileName,int nThreads,kDataFrame* output){
//  if(dynamic_cast<kDataFrameMQF*>(output))
//  {
//    loadIntoMQF(seqFileName,output->getkSize(),nThreads, output->getHasher(),((kDataFrameMQF*)output)->getMQF(),NULL);
//    return;
//  }
  FastqReaderSqueker reader(seqFileName);
  deque<pair<string,string> > reads;
  int k=output->getkSize();
  while(!reader.isEOF())
  {
    reader.readNSeq(&reads,5000);
    for(auto seqPair:reads)
    {
      string seq=seqPair.first;
      for(int i=0;i<seq.size()-k+1;i++){
        string kmer=seq.substr(i,k);
        output->insert(kmer);
      }
    }
    reads.clear();
  }
}
void  parseSequencesFromString(string sequence,kDataFrame* output){
  int k=output->getkSize();
  for(int i=0;i<sequence.size()-k+1;i++){
    string kmer=sequence.substr(i,k);
    output->insert(kmer);
  }
}


struct CustomKmerRow
{
    bool operator()(const pair<kmerRow,int>& lhs, const pair<kmerRow,int>& rhs)
    {
        return lhs.first.hashedKmer > rhs.first.hashedKmer;
    }
};
void merge(const vector<kDataFrame*>& input,kDataFrame* res,kmerRow (*fn)(vector<kmerRow>& i)){
  priority_queue <pair<kmerRow,int>,vector<pair<kmerRow,int> > , CustomKmerRow  > Q;
  vector<kDataFrameIterator> iterators(input.size());
  for(int i=0;i<input.size();i++)
  {
    iterators[i]=input[i]->begin();
    if(iterators[i]!=input[i]->end()){
    //  cout<<i<<" "<<(*iterators[i]).hashedKmer<<endl;
      Q.push(make_pair(*iterators[i],i));
    }
  }
  vector<kmerRow> current(input.size());
  while(Q.size()>0)
  {
    for(int i=0;i<current.size();i++)
      current[i]=kmerRow();
    pair<kmerRow,int> top=Q.top();
    Q.pop();
    iterators[top.second]++;
    current[top.second]=top.first;
  //  cout<<top.first.hashedKmer<<" "<<top.second<<endl;
    if(iterators[top.second]!=input[top.second]->end())
      Q.push(make_pair(*iterators[top.second],top.second));
    while(Q.size()>0 && top.first.hashedKmer==Q.top().first.hashedKmer)
    {
      top=Q.top();
      Q.pop();
      current[top.second]=top.first;
      iterators[top.second]++;
      if(iterators[top.second]!=input[top.second]->end())
        Q.push(make_pair(*iterators[top.second],top.second));
    }

    kmerRow newRow=fn(current);
    res->insert(newRow);
  }

}
kDataFrame* kFrameUnion(const vector<kDataFrame*>& input){
  kDataFrame* res=input[0]->getTwin();
  uint64_t numKmers=0;
  for(auto kframe:input)
  {
    numKmers+=kframe->size();
  }
  res->reserve((uint64_t)((double)numKmers*0.75));
  merge(input,res,[](vector<kmerRow>& input)->kmerRow {
    kmerRow res;
    for(int i=0;i<input.size();i++)
    {
      if(input[i].count!=0)
      {
        res.kmer=input[i].kmer;
        res.hashedKmer=input[i].hashedKmer;
      }
      res.count+=input[i].count;
    }
    return res;
  });
  return res;
}

kDataFrame* kFrameIntersect(const vector<kDataFrame*>& input){
  kDataFrame* res=input[0]->getTwin();
  uint64_t numKmers=numeric_limits<uint64_t>::max();
  for(auto kframe:input)
  {
    numKmers=min(numKmers,(uint64_t)kframe->size());
  }
  res->reserve((uint64_t)((double)numKmers*1.2));
  merge(input,res,[](vector<kmerRow>& input)->kmerRow {
    kmerRow res=input[0];
  //  cout<<"Start"<<endl;
  //  cout<<input[0].kmer<<endl;
    for(int i=1;i<input.size();i++)
    {
      //cout<<input[i].kmer<<endl;
      res.count=min(res.count,input[i].count);
    }
    if(res.count==0)
      return kmerRow();
    //cout<<res.count<<endl;
    return res;
  });
  //cout<<"Size "<<res->size()<<endl;
  return res;
}

kDataFrame* kFrameDiff(const vector<kDataFrame*>& input){
  kDataFrame* res=input[0]->getTwin();
  merge(input,res,[](vector<kmerRow>& input)->kmerRow {
    kmerRow res=input[0];
    bool found=false;
    for(int i=1;i<input.size();i++)
    {
      if(input[i].count>0)
      {
        found=true;
        break;
      }
    }
    if(found)
      return kmerRow();
    return res;
  });
  //cout<<"Size "<<res->size()<<endl;
  return res;
}
colored_kDataFrame* index(string input_file,string names_fileName,uint64_t kSize)
{
  map<string,string> namesMap;
  map<string,uint64_t> tagsMap;
  map<string,uint64_t> groupNameMap;
  std::map<uint64_t, std::vector<uint32_t> >  *legend=new std::map<uint64_t, std::vector<uint32_t> >();
  map<uint64_t,uint64_t> colorsCount;
  uint64_t readID=0,groupID=1;
  ifstream namesFile(names_fileName.c_str());
  string seqName,groupName;
  priority_queue <uint64_t, vector<uint64_t>, std::greater<uint64_t> > freeColors;
  map<string,uint64_t> groupCounter;
  while(namesFile>>seqName>>groupName){
    namesMap.insert(make_pair(seqName,groupName));
    auto it=groupNameMap.find(groupName);
    groupCounter[groupName]++;
    if(it==groupNameMap.end())
    {
      groupNameMap.insert(make_pair(groupName,groupID));
      tagsMap.insert(make_pair(to_string(groupID),groupID));
      vector<uint32_t> tmp;
      tmp.clear();
      tmp.push_back(groupID);
      legend->insert(make_pair(groupID,tmp));
      colorsCount.insert(make_pair(groupID,0));
      groupID++;
    }
  }



  seqan::SeqFileIn seqIn(input_file.c_str());

  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::CharString> reads;
  int chunkSize=1000;
  vector<kDataFrameMQF*> frames;
  int currIndex=0;
  string kmer;
  uint64_t tagBits=0;
  uint64_t maxTagValue=(1ULL<<tagBits)-1;

  kDataFrame *frame;
  frame = new kDataFrameMQF(kSize, 27, 2, tagBits, 0);
//frame = new kDataFrameMAP(kSize);


  uint64_t lastTag=0;
  readID=0;

  while(!atEnd(seqIn)){
    clear(reads);
    clear(ids);

    seqan::readRecords(ids, reads, seqIn,chunkSize);

    map<uint64_t,uint64_t> convertMap;


    for(int j=0;j<length(reads);j++)
    {
      string readName=string((char*)seqan::toCString(ids[j]));

      auto it=namesMap.find(readName);
      if(it==namesMap.end())
      {
        cout<<"read "<<readName<<"dont have group. Please check the group names file."<<endl;
      }
      string groupName=it->second;

      uint64_t readTag=groupNameMap.find(groupName)->second;

      string seq=string((char*)seqan::toCString(reads[j]));
      if(seq.size()<kSize)
        continue;
      convertMap.clear();
      convertMap.insert(make_pair(0,readTag));
      convertMap.insert(make_pair(readTag,readTag));
  //    cout<<readName<<"   "<<seq.size()<<endl;
      for(int i=0;i<seq.size()-kSize+1;i++)
      {
        kmer=seq.substr(i,kSize);
      //  cout<<i<<" "<<kmer<<" "<<frame->count(kmer)<<endl;
        //frame->incrementCounter(kmer,1);
        uint64_t currentTag=frame->count(kmer);

        auto itc=convertMap.find(currentTag);
        if(itc==convertMap.end())
        {
          vector<uint32_t> colors=legend->find(currentTag)->second;
          auto tmpiT=find(colors.begin(),colors.end(),readTag);
          if(tmpiT==colors.end()){
            colors.push_back(readTag);
            sort(colors.begin(),colors.end());
          }

          string colorsString=to_string(colors[0]);
          for(int k=1;k<colors.size();k++)
          {
            colorsString+=";"+to_string(colors[k]);
          }

          auto itTag=tagsMap.find(colorsString);
          if(itTag==tagsMap.end())
          {
            uint64_t newColor;
            if(freeColors.size()==0){
              newColor=groupID++;
            }
            else{
              newColor=freeColors.top();
              freeColors.pop();
            }

            tagsMap.insert(make_pair(colorsString,newColor));
            legend->insert(make_pair(newColor,colors));
            itTag=tagsMap.find(colorsString);
            colorsCount[newColor]=0;
            // if(groupID>=maxTagValue){
            //   cerr<<"Tag size is not enough. ids reached "<<groupID<<endl;
            //   return -1;
            // }
          }
          uint64_t newColor=itTag->second;

          convertMap.insert(make_pair(currentTag,newColor));
          itc=convertMap.find(currentTag);
        }

        if(itc->second!=currentTag)
        {

          colorsCount[currentTag]--;
          if(colorsCount[currentTag]==0  && currentTag!=0 ){
            freeColors.push(currentTag);
            legend->erase(currentTag);
            if(convertMap.find(currentTag)!=convertMap.end())
              convertMap.erase(currentTag);
          }
          colorsCount[itc->second]++;
        }



        frame->setCount(kmer,itc->second);
        if(frame->count(kmer)!=itc->second)
        {
              //frame->setC(kmer,itc->second);
          cout<<"Error Founded "<<kmer<<" from sequence "<<readName<<" expected "
          <<itc->second<<" found "<<frame->count(kmer)<<endl;
          return NULL;
        }

      }
      readID+=1;
      if(colorsCount[readTag]==0)
      {
        groupCounter[groupName]--;
        if(groupCounter[groupName]==0){
          freeColors.push(readTag);
          legend->erase(readTag);
        }
      }

    }

  }
  colorTable* colors=new intVectorsTable();
  for(auto it:*legend)
  {
    colors->setColor(it.first,it.second);
  }

  colored_kDataFrame* res=new colored_kDataFrame();
  res->setColorTable(colors);
  res->setkDataFrame(frame);
  for(auto iit=namesMap.begin(); iit!=namesMap.end();iit++){
    uint32_t sampleID=groupNameMap[iit->second];
    res->namesMap[sampleID]=iit->first;
    res->namesMapInv[iit->first]=sampleID;
  }
  return res;

  // frame->set_legend(legend);
  //
  // cerr << "[!] Done Indexing, writing files ... " << endl;
  //
  // frame->save(outDB);
  // ofstream outNames(outDB+".namesMap");



}

//kDataFrame* intersect(vector<kDataFrame*> input);
//kDataFrame* diff(vector<kDataFrame*> input);
}