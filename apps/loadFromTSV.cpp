//
// Created by mostafa on 4/6/20.
//

#include <iostream>
#include <string>
#include "kDataFrame.hpp"
#include <math.h>
using namespace std;

int main(int argc, char *argv[])
{
    int k=atoi(argv[1]);
    uint64_t  numSlots=atoi(argv[2]);
    string outPath=argv[3];
    //    numSlots=(double)numSlots*1.4;
    cout<<"Expected number of slots "<<numSlots<<endl;
    uint64_t s=log2(numSlots);
    cout<<"S = "<<s<<endl;
    kDataFrameBMQF frame(k,s,2,0,0,outPath);
    //kDataFrameMQF frame(k);
    cout<<"reserve completed"<<endl;
    string kmer;
    uint64_t  count;
    uint64_t nKmers=0;
    while (cin>>kmer>>count)
    {
        frame.insert(kmer,count);
	nKmers++;
	if(nKmers%10000==0)
	 cout<<nKmers<<endl;
    }
    cout<<"finished "<<nKmers<<endl;
    cout<<"serialize to  "<<outPath<<endl; 
    frame.serialize(outPath);
    return 0;

}
