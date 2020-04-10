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

    kDataFrame* indexFrame=kDataFrame::load(sample);

    uint64_t testedKmers=0;
    uint64_t failedKmers=0;

    for(int i=0;i<filenames.size();i++)
    {
        ifstream inp(filenames[i]+".testkmers");
        string kmer;
        while(inp>>kmer)
        {
            testedKmers++;
            vector<uint32_t> colors=indexFrame->getKmerDefaultColumnValue<vector<uint32_t >, colorColumn>(kmer);
            auto colorIt=find(colors.begin(),colors.end(),i);
            if(colorIt==colors.end())
            {
                cerr<<"Error detected in sample "<<i <<endl;
                failedKmers++;
            }

        }
        inp.close();
    }
    cout<<"Numbers of tested Kmers = "<<testedKmers<<endl;
    cout<<"Numbers of failed kmers = "<<failedKmers<<endl;

    return 0;
}
