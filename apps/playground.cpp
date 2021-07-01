#include <iostream>
#include <string>
#include "CLI11.hpp"
#include "kDataFrame.hpp"
#include <cstdlib>
#include "set"
#include <fstream>
#include <algorithm>
#include "algorithms.hpp"
using namespace std;


int main(int argc, char *argv[]){

    string inputColumn=argv[1];
    prefixTrie *pColumn = new prefixTrie();
    pColumn->deserialize(inputColumn);
    pColumn->explainSize();
    unordered_map<uint32_t,uint32_t> itemCount;
    for(auto vec:pColumn->edges)
      {
	for(auto i :*vec)
        itemCount[i]++;
      }
    deque<pair<uint32_t,uint32_t> > itemCountVec;
    for(auto i:itemCount)
      itemCountVec.push_back(make_pair(i.second,i.first));
    sort(itemCountVec.begin(),itemCountVec.end());

    sdsl::int_vector<> translate(itemCountVec.size()+1);
    unordered_map<uint32_t,uint32_t> reverse(itemCountVec.size()+1);
    cout<<"Unique Items = "<<itemCountVec.size()<<endl;
    uint32_t ii=0;
    for(;ii<pColumn->noSamples;ii++)
      {
	translate[ii]=ii;
	reverse[ii]=ii;
      }
    
    for(auto i:itemCountVec)
      {
	if(i.second < pColumn->noSamples)
	  continue;
	translate[ii]=i.second;
	reverse[i.second]=ii;
	ii++;
      }
    double eSize=0.0;
    for(auto vec:pColumn->edges)
      {
	sdsl::int_vector<> tmp(vec->size());
	for(uint32_t i=0;i<vec->size();i++)
	  {
	    tmp[i]=reverse[(*vec)[i]];
	  }
	//sdsl::enc_vector<> tmpCompressed(tmp);
	sdsl::util::bit_compress(tmp);
	//sdsl::enc_vector<> tmpCompressed(tmp);
	eSize += sdsl::size_in_mega_bytes(tmp);
	
	  
        
      }
    cout<<"new size = "<<eSize<<"MB"<<endl;
    
    delete pColumn;
    return 0;

    
  //    CLI::App app;
//    string input_file;
//
//
//    app.add_option("-i,--input", input_file,
//                   "KdataFrame")->required();
//
//
//
//    CLI11_PARSE(app, argc, argv);
//    //kDataFrame* frame=kDataFrame::load(input_file);
//    kDataFrame* frame =new kDataFrameMQF();
//    kProcessor::loadFromKMC(frame, "input.unitigs");
//
//    cout<<"size "<<frame->size()<<endl;
//    frame->addColumn("visited",new vectorColumn<bool>(frame->size()));
//
//    for(auto kmer:*frame)
//    {
//        cout<<kmer.getKmer()<<endl;
//        if(!frame->getKmerColumnValue<bool,vectorColumn<bool> >("visited",kmer.kmer)) {
//            dbgIterator it = frame->getDBGIterator(kmer.kmer);
//            bool moreKmers=true;
//            while(moreKmers)
//            {
//                cout<<" "<<it.currentKmer<<endl;
//                frame->setKmerColumnValue<bool,vectorColumn<bool> >("visited",it.currentKmer,true);
//                moreKmers=false;
//                for(unsigned int i=0;i<it.nextFwdKmers.size();i++)
//                {
//                    if(!frame->getKmerColumnValue<bool,vectorColumn<bool> >("visited",it.nextFwdKmers[i]))
//                    {
//                        moreKmers=true;
//                        it.nextFWD(i);
//                        break;
//                    }
//                }
//            }
//        }
//
//    }
//    cout<<"Navigate dbg is completed"<<endl;






    return 0;
}
