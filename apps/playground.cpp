#include <iostream>
#include <string>
#include "CLI11.hpp"
#include "kDataFrame.hpp"
#include <cstdlib>
#include "set"
#include <fstream>
#include <algorithm>
#include "algorithms.hpp"
#include <parallel_hashmap/btree.h>

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/archives/binary.hpp>
using namespace std;


int main(int argc, char *argv[]){


    string columnPrefix=argv[1];
    auto newColor=new deduplicatedColumn<prefixTrie,phmap::btree_map<uint32_t,uint32_t> >();
    newColor->deserialize(columnPrefix);
    cout<<"IDS map size = "<<sdsl::size_in_mega_bytes(newColor->values->idsMap)<<" MB"<<endl;
    sdsl::util::bit_compress(newColor->values->idsMap);
    cout<<"Bit Compressed IDS map size = "<<sdsl::size_in_mega_bytes(newColor->values->idsMap)<<" MB"<<endl;
//    uint64_t nkmers=2722829058;
    uint64_t nkmers=2722829058;
    uint32_t maximum=0;
    sdsl::int_vector<32> intVector(nkmers+1);
    for(auto i:newColor->index)
    {
        maximum=max(maximum,i.second);
        intVector[i.first]=i.second;
    }
    deque<pair<double,string> > results;
    double tmpSize=newColor->index.size()*4*2/(1024*1024);
    results.push_back(make_pair(tmpSize,"btree_map<uint32_t,uint32_t>"));

    cout<<"Maximum value = "<<maximum<<endl;

    tmpSize=sdsl::size_in_mega_bytes(intVector);
    results.push_back(make_pair(tmpSize,"Plain Int Vector"));
    sdsl::util::bit_compress(intVector);
    tmpSize=sdsl::size_in_mega_bytes(intVector);
    results.push_back(make_pair(tmpSize,"bit compress Int Vector"));
    auto  vlcVector=new sdsl::vlc_vector<>(intVector);
    tmpSize=sdsl::size_in_mega_bytes(*vlcVector);
    results.push_back(make_pair(tmpSize,"vlc Vector"));
    delete vlcVector;
    auto encVector = new sdsl::enc_vector<>(intVector);
    tmpSize=sdsl::size_in_mega_bytes(*encVector);
    results.push_back(make_pair(tmpSize,"enc Vector"));
    delete encVector;
    auto dacVector=new sdsl::dac_vector<>(intVector);
    tmpSize=sdsl::size_in_mega_bytes(*dacVector);
    results.push_back(make_pair(tmpSize,"dac Vector"));
    delete dacVector;

    sort(results.begin(),results.end());
    cout<<"Winner: "<<results[0].second<< " size = "<<results[0].first<<" MB"<<endl;

    for(auto r:results)
        cout<<r.second<<" : "<<r.first<<" MB"<<endl;

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
