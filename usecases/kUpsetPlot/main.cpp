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
	cerr<<fileName<<" loaded"<<endl; 
    }
    uint32_t kSize=genomes_kframes[0]->getkSize();
    kDataFrame* kframe=kDataFrameFactory::createPHMAP(kSize);
    kProcessor::indexPriorityQueue(genomes_kframes, "", kframe);
    cerr<<"Indexing done"<<endl;
    auto colorsCounts=kframe->getColorHistogram();


    for(auto currColor: colorsCounts)
    {
        uint32_t colorID=currColor.first;
        uint32_t colorCount=currColor.second;
        vector<uint32_t> color=kframe->getColorByColorID(colorID);

        cout<<genomeFileNames[color[0]];
        for(unsigned i=1; i<color.size(); i++ )
            cout<<"&"<<genomeFileNames[color[i]];
        cout<<"="<<colorCount/1000<<", ";
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
    upsetPlot(genomesFileNames);

    //differntialExpression(genes_file,samples,control,outputFilename);
    return 0;
}
