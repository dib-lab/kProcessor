//
// Created by mostafa on 4/6/20.
//

#include <iostream>
#include <string>
#include "kDataFrame.hpp"
#include <vector>
#include "algorithms.hpp"
#include <typeinfo>
using namespace std;

int main(int argc, char *argv[])
{
    string inputPath=argv[1];
    uint64_t  q=atoi(argv[2]);
    string outPath=argv[3];

    vector<string> filenames;
    vector<kDataFrame*> frames;
    string sample;
    ifstream input(inputPath);


    while(input>>sample)
    {
        filenames.push_back(sample);
        frames.push_back(kDataFrame::load(sample));
        if(dynamic_cast<kDataFrameBMQF*>(frames.back()))
        {
            ((kDataFrameBMQF*)frames.back())->deleteMemoryBuffer();
        }
	cerr<<"sample "<<sample<<" loaded"<<endl; 
    }
    uint64_t  kSize=frames[0]->getkSize();
    for(auto f:frames)
    {
        if(f->getkSize()!=kSize)
        {
            cerr<<"All Kdataframes should have the same kSize "<<endl;
            return -1;
        }
    }

    kDataFrame* output= new kDataFrameMQF(kSize,q,1);
    kProcessor::indexPriorityQueue(frames,"",output);
    cout<<"Indexing Finished"<<endl;
//    uint64_t testedKmers=0;
//    uint64_t failedKmers=0;
//    uint64_t notFoundKmers =0;
//
//    for(int i=0;i<filenames.size();i++)
//    {
//        ifstream inp(filenames[i]+".testkmers");
//        string kmer;
//        while(inp>>kmer)
//        {
//            testedKmers++;
//            vector<uint32_t> colors=output->getKmerDefaultColumnValue<vector<uint32_t >, queryColorColumn >(kmer);
//            if(colors.size()==0)
//            {
//                cout<<filenames[i]<<" KMER "<<kmer<<endl;
//                notFoundKmers++;
//                continue;
//            }
//            auto colorIt=find(colors.begin(),colors.end(),i);
//            if(colorIt==colors.end())
//            {
//                cerr<<"Error detected in sample #"<<i <<" "<<
//                    filenames[i]<<" at kmer "<<kmer<<" Combination got "<<endl;
//                for(auto c: colors)
//                    cerr<<c <<" ";
//                cerr<<endl;
//
//                failedKmers++;
//            }
//
//        }
//        inp.close();
//    }
//    cout<<"Numbers of tested Kmers = "<<testedKmers<<endl;
//    cout<<"Numbers of non found kmers = "<<notFoundKmers<<endl;
//    cout<<"Numbers of wrong combination = "<<failedKmers<<endl;



    output->save(outPath);
    cout<<"Saving Finished"<<endl;
    delete output;
    for(auto f:frames)
        delete f;

    return 0;
}
