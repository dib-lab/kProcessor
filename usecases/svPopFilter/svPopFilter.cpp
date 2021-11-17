#include <iostream>
#include <string>
#include "CLI11.hpp"
#include "kDataFrame.hpp"
#include <vector>
#include <algorithm>
#include "algorithms.hpp"
#include <any>
#include "Utils/utils.hpp"
using namespace std;



int main(int argc, char *argv[]){
    CLI::App app;
    string kframePath;
    string bcalm_input;
    string outputFilename;


    app.add_option("-k,--kframe", kframePath,
                   "Kdataframe path")->required();
    app.add_option("-s,--sv", bcalm_input,
                   "structural variants in vcf format ")->required();
    app.add_option("-s,--sv", bcalm_input,
                   "structural variants in vcf format ")->required();
    app.add_option("-o,--output", outputFilename,
                   "Output file of median change per genes")->required();



    CLI11_PARSE(app, argc, argv);


    kDataFrame* kframe= kDataFrame::load(kframePath);

    uint32_t kSize=kframe->ksize();



/// calculating count Histogram

    unordered_map<uint32_t,uint32_t> hist;
    uint32_t maximum=0;
    for(auto k:*kframe){
        uint32_t count=k.getCount();
        hist[count]++;
        maximum=max(maximum,count);
    }

    uint32_t histSize=maximum;
    cout<<"Maximum is "<<histSize<<endl;
    uint64_t* histArr= new uint64_t [histSize+1];
    memset(histArr,0,histSize *sizeof(uint64_t));
    for(auto h:hist)
    {
        histArr[h.first-1]=h.second;
    }

    double fp,fn;
    int kThreshold=kProcessor::utils::cleaning_pick_kmer_threshold(histArr, histSize,
                                 NULL, NULL,
                                 &fp, &fn);
    cout<<"Kmer threshold = "<<kThreshold<<endl;


    uint32_t minKeepTip=2*kSize;
    uint32_t numRemovedLowCov=0, numRemovedTips=0, numRemoved=0,numRemovedBoth=0,  numRemovedKmers=0;

    kmerDecoder *KD = kmerDecoder::getInstance(bcalm_input, 1000, kframe->KD->slicing_mode, kframe->KD->hash_mode, {{"kSize", kframe->ksize()}});

    vector<uint32_t> untigsCount(1000000);
    while (!KD->end()) {
        KD->next_chunk();
        for (const auto &seq : *KD->getKmers()) {
            size_t pos= seq.first.find("L:") +1;
            bool isTip= seq.first.find("L:",pos) == std::string::npos;
            uint32_t unitigLen = kSize+ seq.second.size() -1;

            bool removeTip= isTip && unitigLen < minKeepTip;
            size_t top=0;
            for (const auto &kmer : seq.second) {
                untigsCount[top++]=kframe->getCount(kmer.hash);
            }
            sort(untigsCount.begin(),untigsCount.begin()+top);
            uint32_t median=untigsCount[top/2];
            bool removeLowCov= median < kThreshold;


            numRemovedTips += removeTip && ! removeLowCov;
            numRemovedLowCov += !removeTip &&  removeLowCov;
            numRemovedBoth += removeLowCov && removeTip;

            if(removeLowCov || removeTip)
            {
                numRemoved++;
                numRemovedKmers+=seq.second.size();
                for (const auto &kmer : seq.second) {
                    kframe->erase(kmer.hash);
                }
            }
        }
    }

    cout<<"Number of Removed Low Coverage Unitigs = "<<numRemovedLowCov<<endl<<
    "Number of Removed Tips = "<<numRemovedTips<<endl<<
    "Number of Removed Both = "<<numRemovedBoth<<endl<<
    "Number of Removed = "<<numRemoved<<endl<<
    "Number of Removed kmers = "<<numRemovedKmers<<endl;

    cout<<"Final kdataframe size = "<<kframe->size()<<endl;

    kframe->save(outputFilename);



    return 0;
}
