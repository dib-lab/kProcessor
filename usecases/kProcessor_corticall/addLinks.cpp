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
    vector<string> readsPath;
    string outputFilename;


    app.add_option("-k,--kframe", kframePath,
                   "Kdataframe path")->required();
    app.add_option("-s,--single-reads", readsPath,
                   "Single end reads in fastq format ")->required();

    app.add_option("-o,--output", outputFilename,
                   "Output prefix for result kDataframe")->required();



    CLI11_PARSE(app, argc, argv);


    kDataFrame* kframe= kDataFrame::load(kframePath);

    uint32_t kSize=kframe->ksize();

    uint32_t chunk_size=10000;

    vector<uint32_t> alignedSegmentID(10000);


    for( auto filename:readsPath){
        cout<<"Processing "<<filename<<endl;
        kmerDecoder *KD = kmerDecoder::getInstance(filename, chunk_size, kframe->KD->slicing_mode, kframe->KD->hash_mode, {{"kSize", kframe->ksize()}});

        while (!KD->end()) {
            KD->next_chunk();
            for (const auto &seq : *KD->getKmers()) {
                if(alignedSegmentID.size() < seq.second.size())
                    alignedSegmentID.resize(seq.second.size() +100);
                uint32_t currSegment=1;
                bool aligned=false;
                deque<uint32_t> gapStart;
                deque<uint32_t> gapEnd;
                for (unsigned i=0; i < seq.second.size() ; i++) {
                    if(kframe->kmerExist(seq.second[i].hash))
                    {
                        alignedSegmentID[i]=currSegment;
                        aligned=true;
                        if(i!=0 && alignedSegmentID[i-1]==0 && currSegment!=1)
                            gapEnd.push_back(i-1);
                    }
                    else{
                        if(i!=0 && alignedSegmentID[i-1]!=0){
                            currSegment++;
                            gapStart.push_back(i);
                        }
                        alignedSegmentID[i]=0;
                    }
                }
                if(alignedSegmentID[seq.second.size()-1]==0)
                    gapStart.pop_back();
                if(aligned && (gapStart.size()>0||gapEnd.size()>0)){
                    for (unsigned i=0; i < seq.second.size() ; i++)
                    {
                        cout<<alignedSegmentID[i]<<" ";
                    }
                    cout<<endl;
                    for(unsigned i=0;i<gapStart.size();i++)
                    {
                        cout<<gapStart[i]<<" ";
                    }
                    cout<<endl;

                    for(unsigned i=0;i<gapEnd.size();i++)
                    {
                        cout<<gapEnd[i]<<" ";
                    }
                    cout<<endl;
                }
                if(aligned)
                {
                    if(currSegment>1)
                    {

                    }

                }
            }
            break;
        }
        delete KD;


    }

    kframe->save(outputFilename);



    return 0;
}
