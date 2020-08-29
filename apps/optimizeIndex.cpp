//
// Created by mostafa on 4/6/20.
//

#include <iostream>
#include <string>
#include "kDataFrame.hpp"
#include <vector>
#include "algorithms.hpp"
using namespace std;

int main(int argc, char *argv[])
{
    string inputList=argv[1];
    string framePath=argv[2];

    vector<string> filenames;
    vector<kDataFrame*> frames;

    string sample;
    ifstream input(inputList);
    while(input>>sample)
    {
        filenames.push_back(sample);
    }

    kDataFrame* indexFrame=kDataFrame::load(framePath);
    ((queryColorColumn*)indexFrame->getDefaultColumn())->optimizeRLE();

    uint64_t testedKmers=0;
    uint64_t failedKmers=0;
    uint64_t notFoundKmers =0;

    for(int i=0;i<filenames.size();i++)
    {
        cerr<<"Testing "<<filenames[i]+".testkmers"<<endl;
        ifstream inp(filenames[i]+".testkmers");
        string kmer;
        uint64_t count;
        while(inp>>kmer>>count)
        {
            testedKmers++;
            vector<uint32_t> colors=indexFrame->getKmerDefaultColumnValue<vector<uint32_t >, queryColorColumn >(kmer);
            if(colors.size()==0)
            {
                notFoundKmers++;
                continue;
            }
            auto colorIt=find(colors.begin(),colors.end(),i);
            if(colorIt==colors.end())
            {
                cerr<<"Error detected in sample #"<<i <<" "<<
                    filenames[i]<<" at kmer "<<kmer<<" Combination got "<<endl;
                for(auto c: colors)
                    cerr<<c <<" ";
                cerr<<endl;

                failedKmers++;
            }

        }
        inp.close();
    }
    cout<<"Numbers of tested Kmers = "<<testedKmers<<endl;
    cout<<"Numbers of non found kmers = "<<notFoundKmers<<endl;
    cout<<"Numbers of wrong combination = "<<failedKmers<<endl;
    return 0;
}
