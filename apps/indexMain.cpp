#include <iostream>
#include <string>
#include "CLI11.hpp"
#include <vector>
#include <stdint.h>
#include <gqf.hpp>
#include "KmerCounter/KmerCounter.hpp"
#include "KmerCounter/kmer.h"
#include "fstream"
#include <cmath>
#include <map>
#include <seqan/seq_io.h>
#include "kDataFrame.hpp"
#include <algorithm>
#include <queue>
using namespace std;

int index_main(int argc, char *argv[]){
  CLI::App app;
  string input_file;
  string names_fileName="";
  string outDB;
  string method;
  bool _DEBUG = 0;
  int kSize;
  int numThreads=1;

  app.add_option("-i,--input", input_file,
   "Fasta file containing the sequences to create the cDBG.")->required()
  ->check(CLI::ExistingFile);

  app.add_option("-o,--output", outDB,
   "Output kDataFrame filename.")->required();

  app.add_option("-n,--names", names_fileName,
   "TSV file of two columns: fasta sequences header and group name. If not supplied the sequence header is used as the group name. ")
  ->required()->check(CLI::ExistingFile);

  app.add_option("-k,--kmer-length",kSize,"kmer length")->required();



  app.add_option("-t,--threads", numThreads,
   "Number of threads used in kmer counting. Default 1");

  app.add_option("-m,--method", method,
                 "Data strucure used in indexing <MAP> or <MQF>, default: ksize < 32 = MQF, Ksize > 31 = MAP.");

  app.add_option("-d,--debug", _DEBUG,
                 "to activate debugging, set -d 1.");


  CLI11_PARSE(app, argc, argv);
  map<string,string> namesMap;

  if (method.empty())
  {
    method = (kSize <= 31) ? "MQF" : "MAP";
  }
  else
  {
    if (method.compare("MAP") && method.compare("MQF"))
    {
      cerr << "method can be only set to <MAP> or <MQF>" << endl;
      return 0;
    }
    else if (kSize > 31 && method == "MQF")
    {
      cerr << "can't proceed using MQF with ksize > 31" << endl;
      return 0;
    }
  }

  map<string,uint64_t> tagsMap;
  map<string,uint64_t> groupNameMap;
  std::map<uint64_t, std::vector<int> >  *legend=new std::map<uint64_t, std::vector<int> >();
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
      vector<int> tmp;
      tmp.clear();
      tmp.push_back(groupID);
      legend->insert(make_pair(groupID,tmp));
      colorsCount.insert(make_pair(groupID,0));
      groupID++;
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
  uint64_t tagBits=0;
  uint64_t maxTagValue=(1ULL<<tagBits)-1;

  kDataFrame *frame;
  if (!method.compare("MQF"))
    frame = new kDataFrameMQF(kSize, 29, 2, tagBits, 0);

  else
    frame = new kDataFrameMAP(kSize);

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
      //  cout<<i<<" "<<kmer<<" "<<frame->getCounter(kmer)<<endl;
        //frame->incrementCounter(kmer,1);
        uint64_t currentTag=frame->getCounter(kmer);

        auto itc=convertMap.find(currentTag);
        if(itc==convertMap.end())
        {
          vector<int> colors=legend->find(currentTag)->second;
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
          if(colorsCount[currentTag]==0  && currentTag!=0 && currentTag>groupNameMap.size()){
            freeColors.push(currentTag);
            vector<int> colors = legend->find(currentTag)->second;
            string colorsString = to_string(colors[0]);
            for (int k = 1; k < colors.size(); k++)
            {
              colorsString += ";" + to_string(colors[k]);
            }
            tagsMap.erase(colorsString);
            legend->erase(currentTag);
            if(convertMap.find(currentTag)!=convertMap.end())
              convertMap.erase(currentTag);
          }
          colorsCount[itc->second]++;
        }



        frame->setCounter(kmer,itc->second);
        if(frame->getCounter(kmer)!=itc->second)
        {
              //frame->setC(kmer,itc->second);
          cout<<"Error Founded "<<kmer<<" from sequence "<<readName<<" expected "
          <<itc->second<<" found "<<frame->getTag(kmer)<<" "<<frame->getTag(kmer)<<endl;
          return -1;
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
  cerr<<"Loaded Sequences= "<<readID<<endl;
  //kDataFrameMQF* indexFrame2=frame;
  //string filePath="tests/testData/tmp.kDataFrame";
  //indexFrame2->save(filePath);
  //kDataFrameMQF* indexFrame=(kDataFrameMQF*)kDataFrame::load(filePath, method);
  //kDataFrameMQF* indexFrame=frame;
  cerr<<"Number of Groups= "<<legend->size()<<endl;

  frame->set_legend(legend);

  cerr << "[!] Done Indexing, writing files ... " << endl;

  frame->save(outDB);
  ofstream outNames(outDB+".namesMap");
  for(auto iit=groupNameMap.begin(); iit!=groupNameMap.end();iit++){
    outNames<<iit->first<<"\t"<<iit->second<<endl;
  }

 cerr << "[!] Done writing the index on disk ..." << endl;

  if (!_DEBUG)
    return 0;

  cerr << "[!] Double checking the index quality ..." << endl;
  // auto it=legend->begin();
  // while(it!=legend->end())
  // {
  //   cout<<"map "<<it->first<<" ->  ";
  //   for(auto a:it->second){
  //     auto git=groupNameMap.begin();
  //     string res;
  //     while(git!=groupNameMap.end())
  //     {
  //       if(git->second==a)
  //         {
  //           res=git->first;
  //           break;
  //         }
  //       git++;
  //     }
  //     cout<<res<<" ";
  //   }
  //   cout<<endl;
  //   it++;
  // }
  //
  // kDataFrameMQF* indexFrame=(kDataFrameMQF*)kDataFrame::load(outDB);
  // legend=indexFrame->get_legend();
  kDataFrame* indexFrame;

  if (!method.compare("MQF"))
    indexFrame = (kDataFrameMQF*) frame;

  else
    indexFrame = (kDataFrameMAP *)frame;


  //
  // auto it2=indexFrame->begin();
  // while(!it2.isEnd()){
  //   cout<<(*it2).kmerHash<<" "<<(*it2).count<<endl;
  //   it2++;
  // }
  //
  seqan::SeqFileIn seqIn2(input_file.c_str());
  int readCount=0;
  int correct=0,wrong=0;
  while(!atEnd(seqIn2)){
    clear(reads);
    clear(ids);
    seqan::readRecords(ids, reads, seqIn2,chunkSize);
    for(int j=0;j<length(reads);j++){
      string readName=string((char*)seqan::toCString(ids[j]));

      auto it=namesMap.find(readName);
      if(it==namesMap.end())
      {
        cout<<"read "<<readName<<"dont have group. Please check the group names file."<<endl;
      }
      string groupName=it->second;

      uint64_t readTag=groupNameMap.find(groupName)->second;
      //cout<<groupName<<"-> "<<readTag<<endl;
      string seq=string((char*)seqan::toCString(reads[j]));
      if(seq.size()<kSize)
        continue;
      for (int i = 0; i < seq.size() - kSize + 1; i++)
      {
        kmer=seq.substr(i,kSize);

        uint64_t tag=indexFrame->getCounter(kmer);
      //  cout<<">>"<<tag<<endl;
        auto colors=legend->find(tag);
        if(colors==legend->end())
        {
          wrong++;
          cout<<"Colors not found "<<kmer<<" "<<tag<<endl;
          return -1;
        }
        else{
          auto colorIt=colors->second.end();
          colorIt=find(colors->second.begin(),colors->second.end(),readTag);
          if(colorIt==colors->second.end()){
      //      cout<<"Failed"<<endl;
            cout<<"Found colors dont include the read target "<<kmer<<" readTag = "<<readTag<<endl;
            cout<<"Tag= "<<tag<<endl;
            for(auto a:colors->second){
              cout<<a<<" ";
            }
            cout<<endl;
            wrong++;
          return -1;
          }
          else{
            correct++;
          }
        }
      }
      readCount+=1;
    }

  }
  cerr<<"Tested "<<readCount<<endl;
  cerr<<"Correct "<<correct<<endl;
  cerr<<"Wrong "<<wrong<<endl;
  delete frame;
  //delete legend;


  return 0;
}
