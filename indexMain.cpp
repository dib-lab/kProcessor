#include <iostream>
#include <string>
#include "ThirdParty/CLI11.hpp"
#include <vector>
#include <stdint.h>
#include <gqf.h>
#include "KmerCounter/KmerCounter.hpp"
#include "KmerCounter/kmer.h"
#include "fstream"
#include <cmath>
#include <map>
#include <seqan/seq_io.h>
#include "kDataFrame.hpp"
using namespace std;


int index_main(int argc, char *argv[]){
  CLI::App app;
  string input_file;
  string names_fileName="";
  int kSize;
  int numThreads=1;

  app.add_option("-i,--input", input_file,
   "Fasta file containing the sequences to create the cDBG.")->required()
  ->check(CLI::ExistingFile);

  app.add_option("-n,--names", input_file,
   "TSV file of two columns: fasta sequences header and group name. If not supplied the sequence header is used as the group name. ")
  ->check(CLI::ExistingFile);

  app.add_option("-k,--kmer-length",kSize,"kmer length")->required();


  app.add_option("-t,--threads", numThreads,
   "Number of threads used in kmer counting. Default 1");

CLI11_PARSE(app, argc, argv);
  map<string,string> namesMap;

  if(names_fileName!=""){
    ifstream namesFile(names_fileName.c_str());
    string seqName,groupName;
    while(namesFile>>seqName>>groupName){
      namesMap.insert(make_pair(seqName,groupName));
    }
  }

  omp_set_num_threads(numThreads);
  seqan::SeqFileIn seqIn(input_file.c_str());

  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::CharString> reads;
  int chunkSize=1000;
  vector<kDataFrameMQF*> frames;
  int currIndex=0;
  string kmer;
  while(!atEnd(seqIn)){
    clear(reads);
    clear(ids);
    seqan::readRecords(ids, reads, seqIn,chunkSize);
    if(names_fileName==""){
      for(auto id:ids){
        string tmp=string((char*)seqan::toCString(id));
        namesMap.insert(make_pair(tmp,tmp));
      }
    }
    vector<kDataFrameMQF*> framesBuffer(chunkSize);

    #pragma omp parallel for private(kmer) firstprivate(kSize)
    for(int j=0;j<length(reads);j++)
    {
      kDataFrameMQF* frame=new kDataFrameMQF(kSize,15,2,0,0);
      string seq=string((char*)seqan::toCString(reads[j]));
      for(int i=0;i<seq.size()-kSize;i++)
      {
        kmer=seq.substr(i,kSize);
        frame->incrementCounter(kmer,1);
      }
      framesBuffer[j]=frame;
    }
    for(int j=0;j<length(reads);j++)
    {
      frames.push_back(framesBuffer[j]);
    }
  }
  cout<<"Loaded Sequences= "<<frames.size()<<endl;
  vector<kDataFrameMQF*> tmpFrames;
  // while(frames.size()>1)
  // {
  //   tmpFrames.clear();
  //   cout<<"Curr size= "<<frames.size()<<endl;
  //   #pragma omp parallel for ordered
  //   for(int i=0;i<frames.size()-1;i+=2)
  //   {
  //     vector<kDataFrameMQF*> localFrames;
  //     localFrames.clear();
  //     localFrames.push_back(frames[i]);
  //     localFrames.push_back(frames[i+1]);
  //     if(i==frames.size()-3&&frames.size()%2==1)
  //       localFrames.push_back(frames[i+2]);
  //
  //     kDataFrameMQF* subRes=kDataFrameMQF::index(localFrames);
  //
  //     for(auto a: localFrames){
  //       delete a;
  //     }
  //     #pragma omp ordered
  //     {
  //       tmpFrames.push_back(subRes);
  //     }
  //   }
  //   frames=tmpFrames;
  //   cout<<frames.size()<<endl;
  // }
  kDataFrameMQF* indexFrame2=kDataFrameMQF::index(frames);
//  kDataFrameMQF* indexFrame2=frames[0];
  string filePath="tests/testData/tmp.kDataFrame";
  indexFrame2->save(filePath);
  kDataFrameMQF* indexFrame=(kDataFrameMQF*)kDataFrame::load(filePath);
  std::map<uint64_t, std::vector<int> > * legend=indexFrame->get_legend();



  cout<<"Number of Groups= "<<legend->size()<<endl;
  auto it=legend->begin();
  while(it!=legend->end())
  {
    cout<<"map "<<it->first<<" ->  ";
    for(auto a:it->second)
      cout<<a<<" ";
    cout<<endl;
    it++;
  }
  // auto it2=indexFrame->begin();
  // while(!it2.isEnd()){
  //   cout<<(*it2).kmerHash<<" "<<(*it2).count<<endl;
  //   it2++;
  // }

  seqan::SeqFileIn seqIn2(input_file.c_str());
  int readCount=0;
  while(!atEnd(seqIn2)){
    clear(reads);
    clear(ids);
    seqan::readRecords(ids, reads, seqIn2,chunkSize);
    for(int j=0;j<length(reads);j++){
      string seq=string((char*)seqan::toCString(reads[j]));
      for(int i=0;i<seq.size()-kSize;i++)
      {
        kmer=seq.substr(i,kSize);
        uint64_t tag=indexFrame->getTag(kmer);
        //cout<<tag<<endl;
        auto colors=legend->find(tag)->second;
        auto colorIt=find(colors.begin(),colors.end(),j+readCount);
        if(colorIt==colors.end()){
          cout<<"Failed"<<endl;
          return -1;
        }
      }
    }
    readCount+=length(reads);
  }
  cout<<"Tested "<<readCount<<endl;



  return 0;
}
