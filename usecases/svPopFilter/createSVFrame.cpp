#include <iostream>
#include <string>
#include "CLI11.hpp"
#include "kDataFrame.hpp"
#include <vector>
#include <algorithm>
#include "algorithms.hpp"
#include <any>
#include "Utils/utils.hpp"
#include <kseq/kseq.h>
using namespace std;

std::vector<std::string> tokenize(const std::string& s, char c) {
    auto end = s.cend();
    auto start = end;

    std::vector<std::string> v;
    for( auto it = s.cbegin(); it != end; ++it ) {
        if( *it != c ) {
            if( start == end )
                start = it;
            continue;
        }
        if( start != end ) {
            v.emplace_back(start, it);
            start = end;
        }
    }
    if( start != end )
        v.emplace_back(start, end);
    return v;
}

inline void printFasta(ofstream& out,string& header,string& seq)
{
    out<<">"<<header;
    for(unsigned i=0;i<seq.size();i++)
    {
        if(i%81==0)
            out<<"\n";
        out<<seq[i];
    }
    out<<endl;
}
int main(int argc, char *argv[]){
    CLI::App app;
    uint32_t kSize;
    string vcf_input;
    string refName;
    string outputFilename;
    string tempDir="";


    app.add_option("-k,--kSize", kSize,
                   "Kmer size")->required();
    app.add_option("-s,--sv", vcf_input,
                   "structural variants in vcf format ")->required();
    app.add_option("-g,--ref", refName,
                   "Reference Name in vcf format ")->required();


    app.add_option("-o,--output", outputFilename,
                   "Output prefix")->required();
    app.add_option("-t,--temp-dir", tempDir,
                   "Temporary Directory ");


    CLI11_PARSE(app, argc, argv);

    kseq_t *kseqObj{};
    auto fp = gzopen(refName.c_str(), "r");
    kseqObj = kseq_init(fp);
    KS_FULL_COMMENT = true;
    int size=0;

    unordered_map<string,string> refs;

    while(kseq_read(kseqObj)>=0){
        refs[kseqObj->name.s]=kseqObj->seq.s;
    }

    ifstream vcfFile(vcf_input);
    string line;
    string contigsFileName=tempDir+"contigs.fa";
    ofstream contigsFile(contigsFileName);
    ofstream namesFile(contigsFileName+".names");
    uint32_t numVariants=0;
    unordered_map<string,uint32_t> variantToId;
    while(getline(vcfFile,line))
    {
        if(line[0] == '#')
            continue;
        vector<string> cols=tokenize(line,'\t');
        string ref=cols[0];
        uint32_t pos=atoi(cols[1].c_str());
        string id=cols[2];
        string refSeq=cols[3];
        string sampleSeq=cols[4];
        string p=cols[7];
        vector<string> params= tokenize(p,';');
        string svType="UNKWON";
        for(auto s:params)
        {
            if(s.size()>8 && s.substr(0,7)== "SVTYPE=")
                svType=s.substr(7);
        }
        string eref=refs[ref].substr(pos-1,refSeq.size());
        string newSeq;
        uint32_t start=0;
        if(pos>=50)
            start=pos-50;
        if(svType== "INS")
        {

            uint32_t end=min((uint32_t)refs[ref].size(),pos+50);
            newSeq=refs[ref].substr(start,pos-1-start)+ sampleSeq +refs[ref].substr(pos,end-pos);

        }
        else if (svType == "DEL")
        {
            uint32_t afterDel=pos+refSeq.size();
            uint32_t end=min((uint32_t)refs[ref].size(), afterDel +50 );
            newSeq=refs[ref].substr(start,pos-1-start)+ sampleSeq +refs[ref].substr(afterDel-1,end-afterDel);
        }
        printFasta(contigsFile,id,newSeq);
        namesFile<<id<<"\t"<<id<<endl;
        variantToId[id]=numVariants;
        numVariants++;

        //cout<<ref<<"\t"<<pos<<"\t"<<refSeq<<"\t"<<sampleSeq<<"\t"<<svType<<endl;

    }
    cout<<"variants parsed"<<endl;
    contigsFile.close();
    namesFile.close();

    int chunkSize = 1000;
    kDataFrame* svkmers = new kDataFramePHMAP(kSize, integer_hasher);
    kmerDecoder *KD_KMERS = kProcessor::initialize_kmerDecoder(contigsFileName, chunkSize, "kmers", {{"k_size", kSize}});
    KD_KMERS->setHashingMode(TwoBits_hasher,kSize);
    kProcessor::index(KD_KMERS, contigsFileName+".names", svkmers);
    svkmers->save(outputFilename);

    cout<<"Indexing finished"<<endl;
    //delete svkmers;






    return 0;
}
