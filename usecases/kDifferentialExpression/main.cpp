#include <iostream>
#include <string>
#include "CLI11.hpp"
#include "kDataFrame.hpp"
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
    
// //load kmers,count from kmc DB for sample and control
     for(const auto& filename:allSamples) {
         kDataFrame* currentFrame=new kDataFrameMAP();
         kProcessor::loadFromKMC(currentFrame,filename);
 	//        kProcessor::createCountColumn(currentFrame);
         cout<<"Load "<<filename<< " kmers: "<<currentFrame->size()<<endl;

         const string countColName="count";
         any totalCountAny=kProcessor::aggregate(currentFrame,(uint64_t)0,  [countColName](kDataFrameIterator& it, any v) -> any {
                 uint32_t count0;

                 it.getColumnValue<uint32_t,vectorColumn<uint32_t> >(countColName,count0);
                 return (any)(any_cast<uint64_t>(v) + (uint64_t)count0);
         });

         uint64_t totalCount=  any_cast<uint64_t>(totalCountAny);
         cout <<"Total count = "<<totalCount<<endl;
         kProcessor::transformInPlace(currentFrame,  [=](kDataFrameIterator& it) -> void {
             uint32_t count0;
             it.getColumnValue<uint32_t,vectorColumn<uint32_t> >(countColName,count0);
             double normalized = (double)count0*(100000000.0) / totalCount;
             it.setColumnValue<uint32_t,vectorColumn<uint32_t> >(countColName,(uint32_t)normalized);
         });
         cout<<currentFrame->size()<<endl;
         kFrames.push_back(currentFrame);
    }



     uint64_t kSize=kFrames.back()->getkSize();
     int chunkSize = 1000;
     kDataFrame * genesFrame = new kDataFrameFactory::createMAP(kSize);
     kmerDecoder * KMERS = kProcessor::initialize_kmerDecoder(genes_file, chunkSize, "kmers", {{"k_size", kSize}});
     kProcessor::index(KMERS, genes_file+".names", genesFrame);
     //   kProcessor::createColorColumn(genesFrame);
     kFrames.push_back(genesFrame);
     requiredIndices={kFrames.size()-1};
     string colorColumn="color"+to_string(kFrames.size()-1);
     cout<<"Load "<<genes_file<< " kmers: "<<kFrames.back()->size()<<endl;
     cout<<"Load "<<genes_file<< " kmers: "<<genesFrame->size()<<endl;


     cout<<"Load "<<genes_file<< " kmers: "<<genesFrame->size()<<endl;
     kDataFrame* res= kProcessor::innerJoin(kFrames, requiredIndices);

     cout<<"Joined "<<res->size()<<" kmers"<<endl;
     res=kProcessor::filter(res,[=](kDataFrameIterator& r) -> bool {
         for(unsigned i=0; i < allDatasets ;i++ ){
             uint32_t count;
             r.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count"+ to_string(i),count);
             if(count>0)
                 return true;
         }
         return false;
     });
     const string foldChangeColName="foldChange";
     res->addColumn(foldChangeColName,new vectorColumn<double >(res->size()));

     kProcessor::transformInPlace(res,  [=](kDataFrameIterator& it) -> void {
         unsigned i=0;
         uint32_t sampleSum=0;
         for(;i < nSamples ; i++)
         {
             uint32_t count;
             string colName = "count"+to_string(i);
             it.getColumnValue<uint32_t,vectorColumn<uint32_t> >(colName,count);
             sampleSum+=count;
         }
         uint32_t controlSum=0;
         for(;i < allDatasets ; i++)
         {
             uint32_t count;
             string colName = "count"+to_string(i);
             it.getColumnValue<uint32_t,vectorColumn<uint32_t> >(colName,count);
             controlSum+=count;
         }
         double sampleAVG= (double)sampleSum / (double)nSamples;
         double controlAVG= (double)controlSum / (double)nControl;
         double foldChange= sampleAVG / controlAVG;
         it.setColumnValue<double,vectorColumn<double> >(foldChangeColName,foldChange);

     });

     // if we calculate average we dont need to save all the kmers
     auto foldChangeByGene=new unordered_map<string ,vector<double> >();
     any genesGatherAny=kProcessor::aggregate(res,foldChangeByGene,  [=](kDataFrameIterator& it, any v) -> any {
         auto dict=any_cast<unordered_map<string ,vector<double>>*>(v);
         double foldChange;
         it.getColumnValue<double,vectorColumn<double> >(foldChangeColName,foldChange);
         vector<string> color;
         it.getColumnValue<vector<string> ,deduplicatedColumn<StringColorColumn>>(colorColumn,color);
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
             output<<k.first<<"\t"<<median<<"\n";
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
