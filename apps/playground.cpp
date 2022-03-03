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
#include <omp.h>
#include <headers/codecfactory.h>

using namespace std;


bool cmp(pair<uint32_t,uint32_t> &a,
         pair<uint32_t,uint32_t> &b) {
    return a.second > b.second;
}

// Function to sort the map according
// to value in a (key-value) pairs
vector<pair<uint32_t,uint32_t>> sortMap(unordered_map<uint32_t, uint32_t> &M) {

    // Declare vector of pairs
    vector<pair<uint32_t,uint32_t> > A;

    // Copy key-value pair from Map
    // to vector of pairs
    for (auto &it : M) {
        A.push_back(it);
    }

    // Sort using comparator function
    sort(A.begin(), A.end(), cmp);



    return A;
}

int main(int argc, char *argv[]) {
    using namespace FastPForLib;

    // We pick a CODEC
    IntegerCODEC &codec = *CODECFactory::getFromName("simdfastpfor128");
    // could use others, e.g., "simdbinarypacking", "varintg8iu"
    ////////////

    string indexFileName = argv[1];
    string outputPrefix=argv[2];
    unsigned chunkSize = 1;
    auto newColor = new prefixTrie();
    newColor->deserialize(indexFileName);
    cout << "Number of colors " << newColor->size() << endl;
    newColor->explainSize();

    deque<sdsl::int_vector<> *> unCompressedEdges(newColor->noSamples);
 //   deque<sdsl::int_vector<> *> treeIndexes(newColor->noSamples);
    for (unsigned i = 0; i < newColor->edges.size(); i++) {
        unCompressedEdges[i] = new sdsl::int_vector<>(newColor->edges[i]->size());

//        treeIndexes[i] = new sdsl::int_vector<>(newColor->edges[i]->size());
        for (unsigned j = 0; j < newColor->edges[i]->size(); j++) {
            unsigned newNode = newColor->translateEdges[(*(newColor->edges[i]))[j]];
//            prefixTrieIterator it(newColor,newNode);
            (*(unCompressedEdges[i]))[j] = newNode;
//            (*(treeIndexes[i]))[j] = it.treeIndex;
//
        }
    }

    double encSize = 0;
    double vlcSize = 0;
    double dacSize = 0;
    double bitCompress = 0;
    double translateSize = 0;


//
//    for(auto t:treeIndexes)
//    {
//        sdsl::enc_vector<> encVector(*t);
//        sdsl::vlc_vector<> vlcVector(*t);
//        sdsl::dac_vector<> dacVector(*t);
//        encSize += sdsl::size_in_mega_bytes(encVector);
//        vlcSize += sdsl::size_in_mega_bytes(vlcVector);
//        dacSize += sdsl::size_in_mega_bytes(dacVector);
//
//        sdsl::util::bit_compress(*t);
//        bitCompress += sdsl::size_in_mega_bytes(*t);
//    }
//
//    cout << "Total Tree enc = " << encSize << endl;
//    cout << "Total Tree vlc = " << vlcSize << endl;
//    cout << "Total Tree dac = " << dacSize << endl;
//    cout << "Total Tree bit compress = " << bitCompress << endl;


    encSize = 0;
    vlcSize = 0;
    dacSize = 0;
    bitCompress = 0;
    translateSize = 0;


    unsigned start = 0;
    uint64_t pFORSize=0;
    while (start < newColor->noSamples) {
        uint32_t end = min((uint32_t) (start + chunkSize), (uint32_t) newColor->noSamples);


        unordered_map<uint32_t, uint32_t> nodesCount;
        unordered_map<uint32_t, uint32_t> logNodeCount;
        for (unsigned j = start; j < end; j++) {
            for (auto i: *unCompressedEdges[j]){
                nodesCount[i]++;
                logNodeCount[(uint32_t)log((uint32_t)i)]++;
            }
        }
        string outputFilename=outputPrefix+"."+to_string(start)+".csv";
        ofstream out(outputFilename);
        for(auto i:logNodeCount)
            out<<i.first<<","<<i.second<<"\n";
        out.close();
        sdsl::int_vector<> translateEdges(nodesCount.size());
        cout<<"Tree "<<start<<" Number of Distinct Items = "<<nodesCount.size()<<" Number of Items = "<<unCompressedEdges[start]->size()<<endl;
        auto sortedNodesCount=sortMap(nodesCount);
        uint32_t uniqueNodeID = 0;
        unordered_map<uint32_t, uint32_t> reverse;
        for (auto n:sortedNodesCount) {
            translateEdges[uniqueNodeID] = n.first;
            reverse[n.first] = uniqueNodeID;
            uniqueNodeID++;
        }

        translateSize += sdsl::size_in_mega_bytes(translateEdges);
        double unCompressedSize = 0.0;

        for (unsigned j = start; j < end; j++) {

            uint32_t index = 0;
            for(unsigned k=0;k<unCompressedEdges[j]->size();k++)
                (*(unCompressedEdges[j]))[k]= reverse[(*(unCompressedEdges[j]))[k]];



//            sdsl::util::bit_compress(newVec);
//            bitCompress += sdsl::size_in_mega_bytes(newVec);
        }
        start += chunkSize;
    }
    for (unsigned i = 0; i < newColor->edges.size(); i++) {
        std::vector<uint32_t> newVec(unCompressedEdges[i]->size());
        std::copy(unCompressedEdges[i]->begin(),unCompressedEdges[i]->end(),newVec.begin());
        std::vector<uint32_t> compressedVec(unCompressedEdges[i]->size()+1024);
        size_t compressedsize = compressedVec.size();
        codec.encodeArray(newVec.data(), newVec.size(), compressedVec.data(),
                          compressedsize);
        compressedVec.resize(compressedsize);
        std::cout << std::setprecision(3);
        std::cout << "Tree "<<i<< " are using "
        << 32.0 * static_cast<double>(compressedVec.size()) /
        static_cast<double>(newVec.size())
        << " bits per integer. " << std::endl;

        pFORSize+=(compressedsize*4);

        uint32_t startIndex=0;
        uint32_t VECTOR_SIZE=1000000;
        uint32_t EndIndex=startIndex+VECTOR_SIZE;
        EndIndex=min(EndIndex,(uint32_t)unCompressedEdges[i]->size());

        while(startIndex< unCompressedEdges[i]->size())
        {
            EndIndex=startIndex+VECTOR_SIZE;
            EndIndex=min(EndIndex,(uint32_t)unCompressedEdges[i]->size());
            sdsl::int_vector<> tmpVec(EndIndex-startIndex);
            std::copy(unCompressedEdges[i]->begin()+startIndex,unCompressedEdges[i]->begin()+EndIndex,tmpVec.begin());
            sdsl::enc_vector<> encVector(tmpVec);
            sdsl::vlc_vector<> vlcVector(tmpVec);
            sdsl::dac_vector<> dacVector(tmpVec);
            encSize += sdsl::size_in_mega_bytes(encVector);
            vlcSize += sdsl::size_in_mega_bytes(vlcVector);
            dacSize += sdsl::size_in_mega_bytes(dacVector);
            startIndex+=VECTOR_SIZE;
        }

//        sdsl::enc_vector<> encVector(*(unCompressedEdges[i]));
//        sdsl::vlc_vector<> vlcVector(*(unCompressedEdges[i]));
//        sdsl::dac_vector<> dacVector(*(unCompressedEdges[i]));
//        encSize += sdsl::size_in_mega_bytes(encVector);
//        vlcSize += sdsl::size_in_mega_bytes(vlcVector);
//        dacSize += sdsl::size_in_mega_bytes(dacVector);

    }

//

//
//
//        cout<<"tree = "<<i<<
//        " ,enc size="<<sdsl::size_in_mega_bytes(encVector)<<
//        " ,vlc size="<<sdsl::size_in_mega_bytes(vlcVector)<<
//        " ,dac size="<<sdsl::size_in_mega_bytes(dacVector)<<
//        endl;
    //    cout<<"tree = "<<i<<" ,# integers = "<<newColor->edges[i]->size()<<" ,size = "<<sdsl::size_in_mega_bytes(*(newColor->edges[i]))<<" mb"<<endl;

    cout << "Total enc = " << (uint64_t)encSize << " MB"<<endl;
    cout << "Total vlc = " << (uint64_t)vlcSize << " MB" << endl;
    cout << "Total dac = " << (uint64_t)dacSize << " MB" << endl;
    cout << "Total bit compress = " << bitCompress << endl;
    cout << "Total translate = " << translateSize << endl;
    cout << "Total pFor = " << (double)pFORSize/(1024.0*1024.0) << endl;
   // cout<<pFORSize<<endl;
//    string indexFileName=argv[1];
//    string tmpFolder=argv[2];
//    int nThreads=atoi(argv[3]);
//    omp_set_num_threads(nThreads);
//    kDataFrame* index=kDataFrame::load(indexFileName);
//    kProcessor::createPrefixForest(index,tmpFolder);
//
//    vector<deduplicatedColumn<prefixTrie,phmap::btree_map<uint32_t,uint32_t>>*> columns;
//    uint32_t i=0;
//    auto it=index->columns.find("color"+to_string(i));
//    while(it!=index->columns.end())
//    {
//        columns.push_back((deduplicatedColumn<prefixTrie,phmap::btree_map<uint32_t,uint32_t>>*) it->second);
//        i++;
//        it=index->columns.find("color"+to_string(i));
//    }
//    uint32_t noColumns=i;
//    uint32_t noErrors=0;
//
//    prefixForest* forest=(prefixForest*)index->columns["color"];
//    forest->serialize("forest");
//    cout<<"Kmers to check = "<<index->lastKmerOrder<<endl;
//    uint64_t chunk=index->lastKmerOrder/20;
////#pragma omp parallel for firstprivate (noErrors,noColumns,columns)
//    for(uint64_t o=1;o<index->lastKmerOrder;o++)
//    {
//        if(o%chunk==0)
//            cout<<((double)o/(double)index->lastKmerOrder)*100 <<"% of the kmers are tested"<<endl;
//        uint64_t samplesCount=0;
//        vector<uint32_t> correct;
//        for(uint32_t i=0;i<noColumns;i++)
//        {
//            vector<uint32_t> tmp=columns[i]->get(o);
//            for(unsigned j=0;j<tmp.size();j++)
//                correct.push_back(tmp[j]+samplesCount);
//            samplesCount+= columns[i]->values->noSamples;
//
//        }
//        vector<uint32_t> query=forest->get(o);
//
//        if(query!=correct && noErrors<=10){
//            query=forest->get(o);
//            cout<<"query: ";
//            for(auto c: query)
//                cout<<c<<" ";
//            cout<<endl;
//            cout<<"correct: ";
//            for(auto c: correct)
//                cout<<c<<" ";
//            cout<<endl;
//            noErrors++;
//        }
//
//    }
//

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
