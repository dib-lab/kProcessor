//
// Created by mostafa on 4/6/20.
//

#include <iostream>
#include <string>
#include "kDataFrame.hpp"
#include "algorithms.hpp"
#include <math.h>
#include "CLI11.hpp"
using namespace std;

int main(int argc, char *argv[])
{
    string inputPath;
    string kdataframeType;
    string outPath;
    CLI::App app;

    app.add_option("-i,--input", inputPath,
                   "Prefix of source kdataframe")->required();

    app.add_option("-o,--output", outPath,
                   "Output Path prefix for destination kdataframe")->required();
    app.add_option("-t,--type", kdataframeType,
                 "Type fo destination kdataframe")->required();


    CLI11_PARSE(app, argc, argv);


    kDataFrame* src=kDataFrame::load(inputPath);

    if(kdataframeType=="BMQF")
    {
        kDataFrame* kframe2=kDataFrameFactory::createBMQF(src,outPath);
        kframe2->save(outPath);
        return 0;
    }

    if(kdataframeType=="MQF")
    {
        kDataFrame* kframe2=kDataFrameFactory::createMQF(src);
        kframe2->save(outPath);
        return 0;
    }


    kDataFrame* kframe;
    if(kdataframeType=="PHMAP")
        kframe=kDataFrameFactory::createPHMAP(src->ksize(),src->size());
    else if(kdataframeType=="MAP")
        kframe=kDataFrameFactory::createMAP(src->ksize(),src->size());

    for(auto k:*src)
    {
        kframe->insert(k.getKmer());
        kframe->setOrder(k.getKmer(),k.getOrder());
    }
    for(auto c:src->columns)
    {
        kframe->columns[c.first]=c.second->clone();
    }


    kframe->save(outPath);

    return 0;

}
