#include <iostream>
#include <string>
#include "CLI11.hpp"
#include "kDataFrame.hpp"
#include <vector>
#include <algorithm>
#include "algorithms.hpp"
#include <any>
using namespace std;


// Function to find mean.
float Mean(const vector<uint32_t>& arr)
{
    float sum = 0;
    for (auto a:arr)
        sum = sum + a;
    return sum / (float)arr.size();
}

// Function to find standard
// deviation of given array.
float standardDeviation(const vector<uint32_t>& arr, float mean)
{
    float sum = 0;
    for (auto a:arr)
        sum = sum + (a - mean) *
                    (a - mean);

    return sqrt(sum / (float)(arr.size() - 1));
}

// Function to find t-test of
// two set of statistical data.
float tTest(const vector<uint32_t>& arr1,
            const vector<uint32_t>& arr2)
{
    float mean1 = Mean(arr1);
    float mean2 = Mean(arr2);
    float sd1 = standardDeviation(arr1, mean1);
    float sd2 = standardDeviation(arr2, mean2);

    // Formula to find t-test
    // of two set of data.
    float t_test = (mean1 - mean2) / sqrt((sd1 * sd1)
                                          / (float)arr1.size() + (sd2 * sd2) / (float)arr2.size());
    return t_test;
}


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
    const string countColName="count";

// //load kmers,count from kmc DB for sample and control
     for(const auto& filename:allSamples) {
         kDataFrame* currentFrame=kDataFrame::load(filename);
         any totalCountAny=kProcessor::aggregate(currentFrame,(uint64_t)0,
                                                 [countColName](kDataFrameIterator& it, any v) -> any {
                 uint32_t count0=it.getCount();
                 return (any)(any_cast<uint64_t>(v) + (uint64_t)count0);
         });
         uint64_t totalCount=  any_cast<uint64_t>(totalCountAny);

         kProcessor::transformInPlace(currentFrame,  [=](kDataFrameIterator& it) -> void {
             uint32_t count0=it.getCount();
             double normalized = (double)count0*(100000000.0) / totalCount;
             it.setCount((uint32_t)normalized);
         });
         cout<< "Filename loaded "<<filename<<endl;

         kFrames.push_back(currentFrame);
    }



     uint64_t kSize=kFrames.back()->getkSize();


     int chunkSize = 1000;
     kDataFrame * genesFrame = kDataFrameFactory::createBtree(kSize);
     kmerDecoder * KMERS = new Kmers(genes_file, chunkSize, kSize,genesFrame->KD->hash_mode);
     kProcessor::index(KMERS, genes_file+".names", genesFrame);

     //   kProcessor::createColorColumn(genesFrame);
     kFrames.push_back(genesFrame);
     requiredIndices={kFrames.size()-1};
     string colorColumn="color"+to_string(kFrames.size()-1);
     cout<<"Load "<<genes_file<< " kmers: "<<kFrames.back()->size()<<endl;

     kDataFrame* res= kProcessor::innerJoin(kFrames, requiredIndices);

     cout<<"Joined "<<res->size()<<" kmers"<<endl;



     res=kProcessor::filter(res,[&](kDataFrameIterator& r) -> bool {
         for(unsigned i=0; i < allDatasets ;i++ ){
             uint32_t count;
             r.getColumnValue<vectorColumn<uint32_t> >("count"+ to_string(i),count);
             if(count>0)
                 return true;
         }
         return false;
     });


     cout<<"Filtered "<<res->size()<<endl;

     vector<uint32_t> samplesCounts;
     vector<uint32_t> controlCounts;
     samplesCounts.reserve(samplesInput.size());
     controlCounts.reserve(controlInput.size());
     // calculate pvalues for kmers

    const string ttestColName="t-test";
    res->addColumn(ttestColName, new vectorColumn<float>(res->size()));
     kProcessor::transformInPlace(res,  [&](kDataFrameIterator& it) -> void {
         samplesCounts.clear();
         controlCounts.clear();
         unsigned i=0;
         for(;i < nSamples ; i++)
         {
             uint32_t count;
             string colName = "count"+to_string(i);
             it.getColumnValue<vectorColumn<uint32_t> >(colName,count);
             samplesCounts.push_back(count);
         }
         for(;i < allDatasets ; i++)
         {
             uint32_t count;
             string colName = "count"+to_string(i);
             it.getColumnValue<vectorColumn<uint32_t> >(colName,count);
             controlCounts.push_back(count);
         }
         float ttest= tTest(samplesCounts,controlCounts);
         it.setColumnValue<vectorColumn<float> >(ttestColName, ttest);
     });

     // if we calculate average we dont need to save all the kmers

     auto ttestByGene=new unordered_map<string ,vector<float> >();
     any genesGatherAny=kProcessor::aggregate(res, ttestByGene, [=](kDataFrameIterator& it, any v) -> any {
         auto dict=any_cast<unordered_map<string ,vector<float>>*>(v);
         float ttest;
         it.getColumnValue<vectorColumn<float> >(ttestColName, ttest);
         vector<string> color;
         it.getColumnValue<deduplicatedColumn<StringColorColumn>>(colorColumn,color);
         for(auto c: color)
         {
             (*dict)[c].push_back(ttest);
         }
         return (any)(dict);
     });


     ofstream output(outputFilename.c_str());

     for(auto k:*ttestByGene)
     {
         if(!k.second.empty()) {
             sort(k.second.begin(), k.second.end());
             float median = k.second[k.second.size() / 2];
             output<<k.first<<"\t"<<median<<"\n";
         }

     }
     output.close();


    return;
}

int main(int argc, char *argv[]){
    CLI::App app;
    string genes_file;
    string samplesFileLst;
    string controlFileLst;
    vector<string> samples;
    vector<string> control;
    string outputFilename;


    app.add_option("-g,--genes_file", genes_file,
                   "Genes File fasta")->required();
    app.add_option("-s,--sample_file", samplesFileLst,
                   "Sample KMC db ")->required();
    app.add_option("-c,--control", controlFileLst,
                   "Control KMC db")->required();

    app.add_option("-o,--output", outputFilename,
                   "Output file of median change per genes")->required();


    CLI11_PARSE(app, argc, argv);

    string fileName;

    ifstream samplesInput(samplesFileLst);
    while(samplesInput >> fileName)
        samples.push_back(fileName);

    samplesInput.close();

    ifstream controlInput(controlFileLst);
    while(controlInput >> fileName)
        control.push_back(fileName);

    controlInput.close();

    differntialExpression(genes_file,samples,control,outputFilename);
    return 0;
}
