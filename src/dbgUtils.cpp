//
// Created by Mostafa on 10/13/2021.
//

#include "dbgUtils.hpp"


namespace kProcessor {
    Assembler::Assembler()
    {
        kframe=nullptr;
    }
    Assembler::Assembler(kDataFrame* kframe)
    {
        if(kframe->KD->hash_mode!= TwoBits_hasher)
            throw std::logic_error("KFrame should be using two bit hasher");
        this->kframe=kframe;
    }
    bool Assembler::edgeExists(string kmer){
        return kframe->kmerExist(kmer);
    }
    vector<string> Assembler::getLinks(kDataFrameIterator& it){
        return vector<string>();
    }
    string Assembler::getContig(string seed){
        string nucls="ACGT";
        string contig=seed;
        //move to left
        bool moreWork=true;
        string currKmer=seed;
        const uint32_t kmerSize=seed.size();
        while(moreWork)
        {
            deque<char> edges;
            for(auto n: nucls){
                string newKmer=currKmer.substr(0,kmerSize-1)+n;
                if(edgeExists(newKmer))
                    edges.push_back(n);
            }
            if(edges.size()==1)
            {
                moreWork=true;
                contig+=edges[0];
                currKmer=currKmer.substr(0,kmerSize-1)+edges[0];
            }
            else{
                moreWork=false;
            }
        }
        currKmer=seed;
        string leftPart="";
        moreWork=true;
        while(moreWork)
        {
            deque<char> edges;
            for(auto n: nucls){
                string newKmer=n+currKmer.substr(1,kmerSize);
                if(edgeExists(newKmer))
                    edges.push_back(n);
            }
            if(edges.size()==1)
            {
                moreWork=true;
                leftPart+=edges[0];
                currKmer=edges[0]+currKmer.substr(1,kmerSize);
            }
            else{
                moreWork=false;
            }
        }
        reverse(leftPart.begin(),leftPart.end());
        return leftPart + contig;
    }

}