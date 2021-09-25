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

    string indexFileName=argv[1];
    string tmpFolder=argv[2];
    kDataFrame* index=kDataFrame::load(indexFileName);
    kProcessor::createPrefixForest(index,tmpFolder);

    vector<deduplicatedColumn<prefixTrie,phmap::btree_map<uint32_t,uint32_t>>*> columns;
    uint32_t i=0;
    auto it=index->columns.find("color"+to_string(i));
    while(it!=index->columns.end())
    {
        columns.push_back((deduplicatedColumn<prefixTrie,phmap::btree_map<uint32_t,uint32_t>>*) it->second);
        i++;
        it=index->columns.find("color"+to_string(i));
    }
    uint32_t noColumns=i;
    uint32_t noErrors=0;

    prefixForest* forest=(prefixForest*)index->columns["color"];
    for(uint64_t o=0;o<index->lastKmerOrder;o++)
    {
        uint64_t samplesCount=0;
        vector<uint32_t> correct;
        for(uint32_t i=0;i<noColumns;i++)
        {
            vector<uint32_t> tmp=columns[i]->get(o);
            for(unsigned j=0;j<tmp.size();j++)
                correct.push_back(tmp[j]+samplesCount);
            samplesCount+= columns[i]->values->noSamples;

        }
        vector<uint32_t> query=forest->get(o);

        if(query!=correct){
            query=forest->get(o);
            cout<<"query: ";
            for(auto c: query)
                cout<<c<<" ";
            cout<<endl;
            cout<<"correct: ";
            for(auto c: correct)
                cout<<c<<" ";
            cout<<endl;
            noErrors++;
            if(noErrors>10)
             return -1;
        }

    }


//    string columnPrefix=argv[1];
//    auto newColor=new deduplicatedColumn<prefixTrie,phmap::btree_map<uint32_t,uint32_t> >();
//    newColor->deserialize(columnPrefix);
//    cout<<"IDS map size = "<<sdsl::size_in_mega_bytes(newColor->values->idsMap)<<" MB"<<endl;
//    sdsl::util::bit_compress(newColor->values->idsMap);
//    cout<<"Bit Compressed IDS map size = "<<sdsl::size_in_mega_bytes(newColor->values->idsMap)<<" MB"<<endl;
////    uint64_t nkmers=2722829058;
//    uint64_t nkmers=2722829058;
//    uint32_t maximum=0;
//   // sdsl::int_vector<32> intVector(nkmers+1);
//    for(auto i:newColor->index)
//    {
//        maximum=max(maximum,i.second);
//     //   intVector[i.first]=i.second;
//    }
//    deque<pair<double,string> > results;
//    double tmpSize=newColor->index.size()*4*2/(1024*1024);
//    results.push_back(make_pair(tmpSize,"btree_map<uint32_t,uint32_t>"));
//
//    cout<<"Maximum value = "<<maximum<<endl;
//
//    cout<<"Maximum log = "<<log2(maximum)<<endl;
//    if(log2(maximum)<16.0)
//    {
//        tmpSize=(newColor->index.size()*4+newColor->index.size()*2)/(1024.0*1024.0);
//        results.push_back(make_pair(tmpSize,"btree_map<uint32_t,uint16_t>"));
//    }
//
//
//
////    tmpSize=sdsl::size_in_mega_bytes(intVector);
////    results.push_back(make_pair(tmpSize,"Plain Int Vector"));
//    sdsl::vlc_vector<>* vlcVector;
//    sdsl::enc_vector<>* encVector;
////    sdsl::util::bit_compress(intVector);
////    tmpSize=sdsl::size_in_mega_bytes(intVector);
////   results.push_back(make_pair(tmpSize,"bit compress Int Vector"));
////    vlcVector=new sdsl::vlc_vector<>(intVector);
////    tmpSize=sdsl::size_in_mega_bytes(*vlcVector);
////    results.push_back(make_pair(tmpSize,"vlc Vector"));
////    delete vlcVector;
////    auto encVector = new sdsl::enc_vector<>(intVector);
////    tmpSize=sdsl::size_in_mega_bytes(*encVector);
////    results.push_back(make_pair(tmpSize,"enc Vector"));
////    delete encVector;
////    auto dacVector=new sdsl::dac_vector<>(intVector);
////    tmpSize=sdsl::size_in_mega_bytes(*dacVector);
////    results.push_back(make_pair(tmpSize,"dac Vector"));
////    delete dacVector;
//
//    vector<uint32_t> orders(newColor->index.size());
//    uint32_t ii=0;
//
//
//    for(auto i:newColor->index)
//    {
//       orders[ii++]=i.first;
//    };
//    sort(orders.begin(),orders.end());
//    unordered_map<uint32_t,uint32_t> deltas;
//
//
//    for(ii=1;ii<orders.size();ii++)
//    {
//        deltas[orders[ii]-orders[ii-1]]++;
//    }
//    uint32_t delta=0;
//    uint32_t deltaSize=0;
//    for(auto d: deltas)
//    {
//        if(d.second>deltaSize)
//        {
//            delta=d.first;
//            deltaSize=d.second;
//        }
//    }
//    deltas.clear();
//    cout<<"Delta = "<<delta<<" size= "<<deltaSize<<"total = "<<orders.size()<<endl;
//    string mixResult="Mix:";
//    sdsl::int_vector<> deltaValues(deltaSize+1);
//    uint32_t dIt=0;
//    uint32_t maxOld=0;
//    for(ii=1;ii<orders.size();ii++)
//    {
//        uint32_t currDelta=orders[ii]-orders[ii-1];
//        if(currDelta==delta)
//        {
//            deltaValues[dIt++]=newColor->index[orders[ii]];
//        }
//        else{
//            maxOld=max(maxOld,newColor->index[orders[ii]]);
//        }
//    }
//    if(log2(maxOld)<16.0)
//    {
//        uint32_t numItems=(orders.size()-deltaSize);
//        tmpSize=(numItems*2+numItems*4)/(1024.0*1024.0);
//        mixResult+="btree_map<uint32_t,uint16_t> - ";
//    }
//    else{
//        uint32_t numItems=(orders.size()-deltaSize);
//        tmpSize=(numItems*4+numItems*4)/(1024.0*1024.0);
//        mixResult+="btree_map<uint32_t,uint32_t> - ";
//    }
//    sdsl::util::bit_compress(deltaValues);
//    results.push_back(make_pair(tmpSize+sdsl::size_in_mega_bytes(deltaValues),mixResult+" bit compressed values"));
//
//    encVector=new sdsl::enc_vector<>(deltaValues);
//    results.push_back(make_pair(tmpSize+sdsl::size_in_mega_bytes(*encVector),mixResult+" enc values"));
//    delete vlcVector;
//
//
//
//    sdsl::int_vector intVector2(orders.size());
//    for(ii=0;ii<orders.size();ii++)
//        intVector2[ii]=newColor->index[orders[ii]];
//
//    vlcVector=new sdsl::vlc_vector<>(intVector2);
//    encVector= new sdsl::enc_vector<>(orders);
//    cout<<"Order enc = "<<sdsl::size_in_mega_bytes(*encVector)<<endl;
//    cout<<"values vlc = "<<sdsl::size_in_mega_bytes(*vlcVector)<<endl;
//
//    tmpSize=sdsl::size_in_mega_bytes(*vlcVector)+sdsl::size_in_mega_bytes(*encVector);
//    results.push_back(make_pair(tmpSize,"order(enc) values(vlc)"));
//    delete encVector;
//    delete vlcVector;
//    vlcVector=new sdsl::vlc_vector<>(orders);
//    encVector= new sdsl::enc_vector<>(intVector2);
//    cout<<"Order vlc = "<<sdsl::size_in_mega_bytes(*vlcVector)<<endl;
//    cout<<"values enc = "<<sdsl::size_in_mega_bytes(*encVector)<<endl;
//
//
//    tmpSize=sdsl::size_in_mega_bytes(*vlcVector)+sdsl::size_in_mega_bytes(*encVector);
//    results.push_back(make_pair(tmpSize,"order(vlc) values(enc)"));
//
//
//    sort(results.begin(),results.end());
//    cout<<"Winner: "<<results[0].second<< " size = "<<results[0].first<<" MB"<<endl;
//
//
//
//
//    for(auto r:results)
//        cout<<r.second<<" : "<<r.first<<" MB"<<endl;
//
//    return 0;






    return 0;
}
