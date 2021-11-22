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
    string kframePath;
    string vcf_input;
    string refName;
    string outputFilename;
    string tempDir="";
    string fq1;
    string fq2;


    app.add_option("-k,--kframe", kframePath,
                   "Kdataframe path")->required();
    app.add_option("-s,--sv", vcf_input,
                   "structural variants in vcf format ")->required();
    app.add_option("-r,--ref", refName,
                   "Reference Name in vcf format ")->required();

    app.add_option("-1,--pair1", fq1,
                   "Reads fq format pair 1 ")->required();

    app.add_option("-2,--pair2", fq2,
                   "Reads fq format pair 2")->required();



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
    cout<<"variants parsed parsed"<<endl;
    contigsFile.close();
    namesFile.close();

    kDataFrame* kframe= kDataFrame::load(kframePath);
    uint32_t kSize=kframe->ksize();
    delete kframe;
    int chunkSize = 1000;
    kDataFrame* svkmers = new kDataFramePHMAP(kSize, integer_hasher);
    kmerDecoder *KD_KMERS = kProcessor::initialize_kmerDecoder(contigsFileName, chunkSize, "kmers", {{"k_size", kSize}});
    kProcessor::index(KD_KMERS, contigsFileName+".names", svkmers);
    svkmers->save(tempDir + "SV");
    //delete svkmers;

//    vector<string> joiningInput;
//    joiningInput.push_back(kframePath);
//    joiningInput.push_back(tempDir+"SV");
//
//    kDataFrame* joined=kProcessor::parallelJoin(joiningInput,{0,1});
//    joined->save(outputFilename);


    // filter svs here

    // extract reads

    kseq_t *kseqObjFq1{};
    kseq_t *kseqObjFq2{};
    auto fp1 = gzopen(fq1.c_str(), "r");
    auto fp2 = gzopen(fq2.c_str(), "r");
    kseqObjFq1 = kseq_init(fp1);
    kseqObjFq2 = kseq_init(fp2);

    auto kmerDecoder=svkmers->getkmerDecoder();
    deque<tuple<string,string,string> > readsBuffer1;
    deque<tuple<string,string,string> > readsBuffer2;
    vector<deque<uint32_t> > readsBufferIDS(numVariants);

    unsigned readsProcessed=0;
    while(kseq_read(kseqObjFq1)>=0 && kseq_read(kseqObjFq2)>=0){
        string read1=kseqObjFq1->seq.s;
        string read2=kseqObjFq2->seq.s;
        vector<kmer_row> kmers;
        kmerDecoder->seq_to_kmers(read1,kmers);
        kmerDecoder->seq_to_kmers(read2,kmers);

        string header1=kseqObjFq1->name.s;
        string header2=kseqObjFq2->name.s;

        string qual1= kseqObjFq1->qual.s;
        string qual2= kseqObjFq2->qual.s;

        vector<uint32_t> numMatches(numVariants,0);
        for(auto k:kmers)
        {
            if(svkmers->kmerExist(k.hash) && kframe->kmerExist(k.hash) )
            {
                vector<string> colors=svkmers->getKmerColumnValue<vector<string> ,deduplicatedColumn<StringColorColumn>>("color",k.hash);
                for(auto c:colors)
                {
                    uint32_t i=variantToId[c];
                    numMatches[i]++;
                }
            }
        }
        bool readAdded=false;
        for(unsigned  i=0;i<numMatches.size();i++)
        {
            double perc=(double)numMatches[i]/(double)kmers.size();
            if(perc>=0.5)
            {
                if(!readAdded){
                    readAdded=true;
                    readsBuffer1.push_back(make_tuple(header1,read1,qual1));
                    readsBuffer2.push_back(make_tuple(header2,read2,qual2));
                }
                readsBufferIDS[i].push_back(readsBuffer1.size()-1);
            }
        }
        if(readsProcessed++ %100000 == 0)
        {
            cout<<"Processed "<<readsProcessed<< "reads and saved "<<readsBufferIDS.size()<<"reads"<<endl;
        }

    }
    for(unsigned i=0; i<numVariants; i++)
    {
        ofstream readsOut1(outputFilename+"."+ to_string(i)+".1.fq");
        for(auto r: readsBufferIDS[i])
        {
            readsOut1<< get<0>(readsBuffer1[r])<<"\n"
            <<get<1>(readsBuffer1[r])<<"\n+\n"
            <<get<2>(readsBuffer1[r])<<"\n";
        }
        readsOut1.close();

        ofstream readsOut2(outputFilename+"."+ to_string(i)+".2.fq");
        for(auto r: readsBufferIDS[i])
        {
            readsOut2<< get<0>(readsBuffer2[r])<<"\n"
            <<get<1>(readsBuffer2[r])<<"\n+\n"
            <<get<2>(readsBuffer2[r])<<"\n";
        }
        readsOut2.close();

    }






    return 0;
}
