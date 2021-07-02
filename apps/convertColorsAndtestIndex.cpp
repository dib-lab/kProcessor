//
// Created by mostafa on 4/6/20.
//

#include <iostream>
#include <string>
#include "kDataFrame.hpp"
#include <vector>
#include "algorithms.hpp"
# include <omp.h>
#include <stdlib.h>

using namespace std;

int main(int argc, char *argv[])
{
    string inputList=argv[1];
    string framePath=argv[2];
    int numThreads=atoi(argv[3]);

    omp_set_num_threads(numThreads);

    vector<string> filenames;
    vector<kDataFrame*> frames;

    string sample;
    ifstream input(inputList);
    while(input>>sample)
    {
        filenames.push_back(sample);
    }

    kDataFrame* indexFrame=kDataFrame::load(framePath);
    mixVectors* col=new mixVectors((insertColorColumn*)indexFrame->getDefaultColumn());
    double s=(double)col->sizeInBytes();
    cout<<"Size = "<<s/(1024.0*1024.0)<<"MB"<<endl;
    indexFrame->changeDefaultColumnType(col);

    uint64_t testedKmers=0;
    uint64_t failedKmers=0;
    uint64_t notFoundKmers =0;

    for(uint32_t i=0;i<filenames.size();i++)
    {
        ifstream inp(filenames[i]+".bmqf.testkmers");
        string kmer;
        while(inp>>kmer)
        {
            testedKmers++;
            vector<uint32_t> colors=indexFrame->getKmerColumnValue<vector<uint32_t >, mixVectors>("color",kmer);

	    if(colors.size()==0)
	      {
		notFoundKmers++;
		continue;
	      }
            auto colorIt=find(colors.begin(),colors.end(),i);
           // if(kmer=="GTGCGTGGTCACCCGAGATT")
            if(colorIt==colors.end())
            {
	            cerr<<"Error detected in sample #"<<i <<" "<<
		filenames[i]<<" at kmer "<<kmer<<" Combination "<<indexFrame->getkmerOrder(kmer)<<" got "<<endl;
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
