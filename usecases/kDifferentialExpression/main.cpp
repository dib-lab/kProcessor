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
#include <any>
using namespace std;


void differntialExpression(string genes_file,
    vector<string> samplesInput,
    vector<string> controlInput,
    string outputFilename)
{
    uint32_t nSamples = samplesInput.size();
    uint32_t nControl = controlInput.size();
    uint32_t allDatasets = nSamples+nControl;
    vector<kDataFrame*> kFrames;
    vector<string>  allSamples;
    allSamples.insert(allSamples.end(),samplesInput.begin(),samplesInput.end());
    allSamples.insert(allSamples.end(),controlInput.begin(),controlInput.end());
    vector<uint32_t> requiredIndices;
    unsigned index=0;
//load kmers,count from kmc DB for sample and control
    for(auto filename:allSamples) {
        kDataFrame* currentFrame=new kDataFrameMQF();
        kProcessor::loadFromKMC(currentFrame,filename);
        kProcessor::createCountColumn(currentFrame);
        cout<<"Load "<<filename<< " kmers: "<<currentFrame->size()<<endl;
        any totalCountAny=kProcessor::aggregate(currentFrame,(uint64_t)0,  [](kmerRow it, any v) -> any {
                uint32_t count0;
                it.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count",count0);
                return (any)(any_cast<uint64_t>(v) + (uint64_t)count0);
        });
        double totalCount=  (double)any_cast<uint64_t>(totalCountAny);
        cout <<"Total count = "<<totalCount<<endl;
        currentFrame= kProcessor::transform(currentFrame,  [=](kmerRow it) -> kmerRow {
            uint32_t count0;
            it.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count",count0);
            double normalized = (double)count0*(100000000.0) / totalCount;
            it.setColumnValue<uint32_t,vectorColumn<uint32_t> >("count",(uint32_t)normalized);
            return it;
        });
        kFrames.push_back(currentFrame);

    }

    int kSize=kFrames.back()->getkSize();
    int chunkSize = 1000;
    kDataFrame * genesFrame = new kDataFrameMQF(kSize);
    kmerDecoder * KMERS = kProcessor::initialize_kmerDecoder(genes_file, chunkSize, "kmers", {{"k_size", kSize}});
    kProcessor::index(KMERS, genes_file+".names", genesFrame);
    kProcessor::createColorColumn(genesFrame);
    kFrames.push_back(genesFrame);
    requiredIndices.push_back(kFrames.size()-1);
    string colorColumn="color."+to_string(kFrames.size()-1);
    cout<<"Load "<<genes_file<< " kmers: "<<kFrames.back()->size()<<endl;

    kDataFrame* res= kProcessor::innerJoin(kFrames, requiredIndices);

    res=kProcessor::filter(res,[=](kmerRow r) -> bool {
        for(unsigned i=0; i < allDatasets ;i++ ){
            uint32_t count;
            r.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count.0",count);
            if(count>0)
                return true;
        }
        return false;
    });
    res->addColumn("foldChange",new vectorColumn<double >(res->size()));

    res= kProcessor::transform(res,  [=](kmerRow it) -> kmerRow {
        unsigned i=0;
        uint32_t sampleSum=0;
        for(;i < nSamples ; i++)
        {
            uint32_t count;
            string colName = "count."+to_string(i);
            it.getColumnValue<uint32_t,vectorColumn<uint32_t> >(colName,count);
            sampleSum+=count;
        }
        uint32_t controlSum=0;
        for(;i < allDatasets ; i++)
        {
            uint32_t count;
            string colName = "count."+to_string(i);
            it.getColumnValue<uint32_t,vectorColumn<uint32_t> >(colName,count);
            controlSum+=count;
        }
        double sampleAVG= (double)sampleSum / (double)nSamples;
        double controlAVG= (double)controlSum / (double)nControl;
        double foldChange= sampleAVG / controlAVG;
        it.setColumnValue<double,vectorColumn<double> >("foldChange",foldChange);
        return it;
    });

    auto foldChangeByGene=new unordered_map<uint32_t ,vector<double> >();
    any genesGatherAny=kProcessor::aggregate(res,foldChangeByGene,  [=](kmerRow it, any v) -> any {
        auto dict=any_cast<unordered_map<uint32_t ,vector<double>>*>(v);
        double foldChange;
        it.getColumnValue<double,vectorColumn<double> >("foldChange",foldChange);
        vector<uint32_t> color;
        it.getColumnValue<vector<uint32_t> ,deduplicatedColumn<vector<uint32_t>,StringColorColumn> >(colorColumn,color);
        for(auto c: color)
        {
            (*dict)[c].push_back(foldChange);
        }
        return (any)(dict);
    });
    ofstream output(outputFilename.c_str());

    for(auto k:*foldChangeByGene)
    {
        if(!k.second.empty()) {
            sort(k.second.begin(), k.second.end());
            double median = k.second[k.second.size() / 2];
            output<<k.first<<"\t"<<median<<endl;
        }

    }
    output.close();


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
