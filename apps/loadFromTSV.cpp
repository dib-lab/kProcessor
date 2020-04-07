//
// Created by mostafa on 4/6/20.
//

#include <iostream>
#include <string>
#include "kDataFrame.hpp"
using namespace std;

int main(int argc, char *argv[])
{
    int k=atoi(argv[1]);
    uint64_t  numSlots=atoi(argv[2]);
    string outPath=argv[3];

    kDataFrameBMQF frame(k,outPath);
    //kDataFrameMQF frame(k);
    frame.reserve(numSlots);

    string kmer;
    uint64_t  count;
    uint64_t nKmers=0;
    while (cin>>kmer>>count)
    {
        frame.insert(kmer,count);
	nKmers++;
	if(nKmers%100000==0)
	  cout<<nKmers<<endl;
    }
    cout<<"Finished"<<endl;
    frame.serialize(outPath);
    return 0;

}
