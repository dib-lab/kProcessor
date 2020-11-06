#include <iostream>
#include <string>
#include "CLI11.hpp"
#include "kDataFrame.hpp"
#include <cstdlib>
#include "set"
#include <fstream>
#include <vector>
#include <algorithm>
#include "algorithms.hpp"
using namespace std;


void differntialExpression(string genes_file,
    vector<string> samplesInput,
    vector<string> controlInput,
    string outputFilename)
{
    vector<kDataFrame*> samples,control,kFrames;
    vector<uint32_t> requiredIndices;
    unsigned index=0;
//load kmers,count from kmc DB for sample and control
    for(auto filename:samplesInput) {
        samples.push_back(new kDataFrameMQF());
        kProcessor::loadFromKMC(samples.back(),filename);
        kFrames.push_back(samples.back());
        kProcessor::createCountColumn(kFrames.back());
        cout<<"Load "<<filename<< " kmers: "<<kFrames.back()->size()<<endl;
    }
    for(auto filename:controlInput) {
        control.push_back(new kDataFrameMQF());
        kProcessor::loadFromKMC(control.back(),filename);
        kFrames.push_back(control.back());
        kProcessor::createCountColumn(kFrames.back());
        cout<<"Load "<<filename<< " kmers: "<<kFrames.back()->size()<<endl;
    }
    int kSize=samples.back()->getkSize();
    int chunkSize = 1000;
    kDataFrame * genesFrame = new kDataFrameMQF(kSize);
    kmerDecoder * KMERS = kProcessor::initialize_kmerDecoder(genes_file, chunkSize, "kmers", {{"k_size", kSize}});
    kProcessor::index(KMERS, genes_file+".names", genesFrame);
    kProcessor::createColorColumn(genesFrame);
    kFrames.push_back(genesFrame);
    requiredIndices.push_back(kFrames.size()-1);
    cout<<"Load "<<genes_file<< " kmers: "<<kFrames.back()->size()<<endl;


    for(unsigned int j=0; j<kFrames.size(); j++)
    {
        string outputFilename="out."+ to_string(j);
        ofstream outTmp(outputFilename);
        for( auto k :*kFrames[j])
            outTmp<<k.kmer<<"\n";
        outTmp.close();
    }

    kDataFrame* res= kProcessor::innerJoin(kFrames, requiredIndices);

    kDataFrame* res2=kProcessor::filter(res,[](kmerRow r) -> bool {
        uint32_t count0,count1;
        r.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count.0",count0);
        r.getColumnValue<uint32_t,vectorColumn<uint32_t>>("count.1",count1);
        return (count0>0 || count1>0);

    });
    for(auto k:*res2)
    {
        uint32_t count0,count1;
        vector<uint32_t> color;
        k.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count.0",count0);
        k.getColumnValue<uint32_t,vectorColumn<uint32_t>>("count.1",count1);
        k.getColumnValue<vector<uint32_t>,deduplicatedColumn<vector<uint32_t>,StringColorColumn> >("color.2",color);

        cout<<k.kmer<<"\t"<<count0<<"\t"<<count1;
        for(auto c:color)
        {
            cout<<"\t"<<c;
        }
        cout<<"\n";

    }
    cout<< res->size() <<endl;

    return;
}

int main(int argc, char *argv[]){
    CLI::App app;
    string genes_file;
    vector<string> samples;
    vector<string> control;
    string outputFilename;


    app.add_option("-g,--genes_file", genes_file,
                   "Genes File fasta")->required();
    app.add_option("-s,--sample_file", samples,
                   "Sample KMC db ")->required();
    app.add_option("-c,--control", control,
                   "Control KMC db")->required();

    app.add_option("-o,--output", outputFilename,
                   "Output file of median change per genes")->required();


    CLI11_PARSE(app, argc, argv);
    differntialExpression(genes_file,samples,control,outputFilename);
    return 0;
}
