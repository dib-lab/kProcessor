//
// Created by mostafa on 4/6/20.
//

#include <iostream>
#include <string>
#include "kDataFrame.hpp"
#include "algorithms.hpp"
#include <math.h>
using namespace std;

int main(int argc, char *argv[])
{
    string inputPath=argv[1];
    string kdataframeType=argv[2];
    string outPath=argv[3];

    kDataFrame* kframe;
    if(kdataframeType=="PHMAP" || kdataframeType=="BMQF")
        kframe=new kDataFramePHMAP();
    else if(kdataframeType=="MAP")
        kframe=new kDataFrameMAP();

    kProcessor::loadFromKMC(kframe,inputPath);

    cout<<"Number of kmers ="<<kframe->size()<<endl;

    if(kdataframeType=="BMQF")
    {
        kDataFrame* kframe2=new kDataFrameBMQF(kframe,outPath);
        kframe2->save(outPath);
        return 0;
    }
    kframe->save(outPath);

    return 0;

}
