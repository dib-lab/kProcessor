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
    }
    uint64_t  kSize=frames[0]->getkSize();
    for(auto f:frames)
    {
        if(f->getkSize()!=kSize)
        {
            cerr<<"All Kdataframes should ahve tthe same kSize "<<endl;
            return -1;
        }
    }

    kDataFrame* output= new kDataFrameMQF(kSize,q,1);
    kProcessor::indexPriorityQueue(frames,output);

    output->serialize(outPath);

    delete output;
    for(auto f:frames)
        delete f;

    return 0;
}
