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
    string outPath;
    string fq1;
    string fq2;
    CLI::App app;

    app.add_option("-i,--input", inputPath,
                   "Prefix of source kdataframe")->required();

    app.add_option("-f,--fastq_pair1", fq1,
                   "path for first pair fastq")->required();

    app.add_option("-r,--fastq-pair2", fq2,
                   "path for second pair fastq")->required();




    CLI11_PARSE(app, argc, argv);



    kDataFrame* kframe=kDataFrame::load(inputPath);

    uint32_t chunk_size=10000;
    kmerDecoder *KD1 = kmerDecoder::getInstance(fq1, chunk_size, kframe->KD->slicing_mode, kframe->KD->hash_mode, {{"kSize", kframe->ksize()}});
    kmerDecoder *KD2 = kmerDecoder::getInstance(fq2, chunk_size, kframe->KD->slicing_mode, kframe->KD->hash_mode, {{"kSize", kframe->ksize()}});

    unordered_map<uint32_t,uint32_t> partitionsMap;

    partitionsMap[0]=0;
    uint32_t lastPartition=1;
    string partionColName="partition";
    auto pColumn=new vectorColumn<uint32_t>(kframe->size()+1);
    kframe->addColumn(partionColName,pColumn);
    uint32_t numberOfMatchingReads=0;
    while (!KD1->end()) {
        KD1->next_chunk();
        KD2->next_chunk();
        auto fqChunk1=KD1->getKmers()->begin();
        auto fqChunk2=KD2->getKmers()->begin();

        while(fqChunk1!=KD1->getKmers()->end() && fqChunk2!=KD2->getKmers()->end())
        {
            set<uint32_t> paritions;
            bool readFound=false;
            for(auto kmer: fqChunk1->second)
            {
                if(kframe->kmerExist(kmer.hash)){
                    readFound=true;
                    uint32_t p=kframe->getKmerColumnValue<uint32_t,vectorColumn<uint32_t>>(partionColName,kmer.str);
                    if(p!=0)
                        paritions.insert(p);
                }
            }
            for(auto kmer: fqChunk2->second)
            {
                if(kframe->kmerExist(kmer.hash)){
                    readFound=true;
                    uint32_t p=kframe->getKmerColumnValue<uint32_t,vectorColumn<uint32_t>>(partionColName,kmer.hash);
                    if(p!=0)
                        paritions.insert(p);
                }
            }
            if(readFound){
                numberOfMatchingReads++;
                uint32_t currPartiton;
                if(paritions.size() ==0)
                {
                    currPartiton=lastPartition++;
                    partitionsMap[currPartiton]=currPartiton;
                }
                else if(paritions.size() ==1)
                {
                    currPartiton=*(paritions.begin());
                }
                else
                {
                   // cout<<"Joining Partitions"<<endl;
                    currPartiton=*(paritions.begin());
                    for(auto p:paritions)
                        partitionsMap[p]=currPartiton;
                    kProcessor::transformInPlace(kframe,[&](kDataFrameIterator& it)
                    {
                        uint32_t p;
                        it.getColumnValue<uint32_t, vectorColumn<uint32_t> >(partionColName,p);
                        it.setColumnValue<uint32_t, vectorColumn<uint32_t> >(partionColName,partitionsMap[p]);

                    });
                }
              //  cout<<currPartiton<<endl;

                for(auto kmer: fqChunk1->second)
                {
                    if(kframe->kmerExist(kmer.hash)){
                        kframe->setKmerColumnValue<uint32_t,vectorColumn<uint32_t>>(partionColName,kmer.hash,currPartiton);
                    }
                }
                for(auto kmer: fqChunk2->second)
                {
                    if(kframe->kmerExist(kmer.hash)){
                        kframe->setKmerColumnValue<uint32_t,vectorColumn<uint32_t>>(partionColName,kmer.hash,currPartiton);
                    }
                }


            }
            fqChunk1++;
            fqChunk2++;


        }

    }
    uint32_t numPartitons=0;
    for(auto p:partitionsMap)
    {
        if(p.first==p.second)
            numPartitons++;
    }



    numPartitons--;//zero- zero

    unordered_map<uint32_t,uint32_t> hist;
    for(auto v:pColumn->dataV)
        hist[v]++;




    cout<<"Number of reads "<<numberOfMatchingReads<<endl;
    cout<<"Number of paritions "<<numPartitons<<endl;

    cout<<"Histogram"<<endl;
    for(auto v:hist)
        cout<<v.first<<"\t"<<v.second<<endl;
    return 0;

}
