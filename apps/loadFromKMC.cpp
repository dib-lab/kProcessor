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
        kframe=kDataFrameFactory::createPHMAP(21);
    else if(kdataframeType=="MAP")
        kframe=kDataFrameFactory::createMAP(21);
    else if(kdataframeType=="BTREE")
        kframe=kDataFrameFactory::createBtree(21);
    else if(kdataframeType=="MQF")
        kframe=kDataFrameFactory::createMQF(21);

    kProcessor::loadFromKMC(kframe,inputPath);

    cout<<"Number of kmers ="<<kframe->size()<<endl;

    if(kdataframeType=="BMQF")
    {
        kDataFrame* kframe2=kDataFrameFactory::createBMQF(kframe,outPath);
        kframe2->save(outPath);
        return 0;
    }
    kframe->save(outPath);

    return 0;

}
