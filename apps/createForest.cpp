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
using namespace std;


int main(int argc, char *argv[]){

    string indexFileName=argv[1];
    string tmpFolder=argv[2];
    int nThreads=atoi(argv[3]);
    int numVectos=atoi(argv[4]);
    int vectorLen=atoi(argv[5]);
    omp_set_num_threads(nThreads);
    kDataFrame* index=kDataFrame::load(indexFileName);
    kProcessor::createPrefixForest(index,tmpFolder,numVectos,vectorLen);

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
    forest->serialize("forest");
    cout<<"Kmers to check = "<<index->lastKmerOrder<<endl;
    uint64_t chunk=index->lastKmerOrder/20;
//#pragma omp parallel for firstprivate (noErrors,noColumns,columns)
    for(uint64_t o=1;o<index->lastKmerOrder;o++)
    {
        if(o%chunk==0)
            cout<<((double)o/(double)index->lastKmerOrder)*100 <<"% of the kmers are tested"<<endl;
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

        if(query!=correct && noErrors<=10){
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
        }

    }






    return 0;
}
