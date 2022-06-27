#include <iostream>
#include <string>
#include "CLI11.hpp"
#include "kDataFrame.hpp"
#include <vector>
#include <algorithm>
#include "algorithms.hpp"
#include <any>
using namespace std;




void upsetPlot(vector<string>& genomeFileNames)
{
    vector<kDataFrame*> genomes_kframes;
    for(auto fileName: genomeFileNames)
    {
        genomes_kframes.push_back(kDataFrame::load(fileName));
    }
    uint32_t kSize=genomes_kframes[0]->getkSize();
    kDataFrame* kframe=kDataFrameFactory::createPHMAP(kSize);
    kProcessor::indexPriorityQueue(genomes_kframes, "", kframe);

    auto colorsCounts=new unordered_map<uint32_t ,uint32_t >();
    kProcessor::aggregate(kframe, colorsCounts, [=](kDataFrameIterator& it, any v) -> any {
        auto dict=any_cast<unordered_map<uint32_t,uint32_t >*>(v);
        uint32_t colorID=it.getColorID();
        (*dict)[colorID]++;
    });


    for(auto currColor: *colorsCounts)
    {
        uint32_t colorID=currColor.first;
        uint32_t colorCount=currColor.second;
        vector<uint32_t> color=kframe->getColorByColorID(colorID);

        cout<<genomeFileNames[color[0]];
        for(unsigned i=1; i<color.size(); i++ )
            cout<<"&"<<genomeFileNames[color[0]];
        cout<<"="<<colorCount<<",";
    }
    cout<<endl;
}

int main(int argc, char *argv[]){
    CLI::App app;
    vector<string> genomesFileNames;
    string outputFilename;


    app.add_option("-g,--genome", genomesFileNames,
                   "Genome File fasta")->required();


    CLI11_PARSE(app, argc, argv);


    //differntialExpression(genes_file,samples,control,outputFilename);
    return 0;
}
