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
    uint32_t ksize;
    string vcfFrame;
    string outputFilename;

    double matchPerc=0.7;
    string tempDir="";
    vector<string> fq;


    app.add_option("-k,--ksize", ksize,
                   "ksize")->required();
    app.add_option("-s,--sv", vcfFrame,
                   "structural variants kdataframe ")->required();


    app.add_option("-r,--reads", fq,
                   "Reads pair1:pair2. can be repeated ")->required();


    app.add_option("-p,--percent", matchPerc,
                   "Matching percent for the read default 0.7");

    app.add_option("-o,--output", outputFilename,
                   "Output prefix")->required();
    app.add_option("-t,--temp-dir", tempDir,
                   "Temporary Directory ");


    CLI11_PARSE(app, argc, argv);




    kDataFrame* svkmers = kDataFrame::load(vcfFrame);
    kmerDecoder *kmerDecoder = svkmers->getkmerDecoder();
    //delete svkmers;

//    vector<string> joiningInput;
//    joiningInput.push_back(kframePath);
//    joiningInput.push_back(tempDir+"SV");
//
//    kDataFrame* joined=kProcessor::parallelJoin(joiningInput,{0,1});
//    joined->save(outputFilename);


    // filter svs here

    // extract reads

    auto namesMap=((deduplicatedColumn<StringColorColumn>*)svkmers->columns["color"])->values->namesMap;
    uint32_t numVariants=namesMap.size();

    unordered_map<string,uint32_t> variantToId;
    for(auto v:namesMap)
    {
        variantToId[v.second]=v.first;
    }
    cout<<"NUm of variants "<<numVariants<<endl;
    unsigned readsProcessed=0;

    deque<tuple<string,string,string> > readsBuffer1;
    deque<tuple<string,string,string> > readsBuffer2;
    vector<deque<uint32_t> > readsBufferIDS(numVariants);
    
    for(auto file:fq){

        vector<string> files= tokenize(file,':');
        cout<<"Processing "<<files[0]<<" and "<<files[1]<<endl;

        kseq_t *kseqObjFq1{};
        kseq_t *kseqObjFq2{};

        auto fp1 = gzopen(files[0].c_str(), "r");
        auto fp2 = gzopen(files[1].c_str(), "r");
        kseqObjFq1 = kseq_init(fp1);
        kseqObjFq2 = kseq_init(fp2);


        while(kseq_read(kseqObjFq1)>=0 && kseq_read(kseqObjFq2)>=0){
            string read1=kseqObjFq1->seq.s;
            string read2=kseqObjFq2->seq.s;
	    // cout<<read1<<endl<<read2<<endl;
            vector<kmer_row> kmers;
	    vector<kmer_row> kmers2;
	    kmerDecoder->seq_to_kmers(read1,kmers);
            kmerDecoder->seq_to_kmers(read2,kmers2);
	    kmers.insert(std::end(kmers), std::begin(kmers2), std::end(kmers2));

            vector<uint32_t> numMatches(numVariants,0);
            for(auto k:kmers)
            {
	      //  cout<<k.str<<"\t"<<svkmers->kmerExist(k.str)<<"\t"<< kframe->kmerExist(k.str)<<endl;
                if(svkmers->kmerExist(k.hash))
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
		
                if(perc>=0.6)
                {
		  //  cout<<read1<<endl;
		  //	  cout<<i<<" "<<perc<<endl;
                    if(!readAdded){
		      	string header1=kseqObjFq1->name.s;
			string header2=kseqObjFq2->name.s;

			string qual1= kseqObjFq1->qual.s;
			string qual2= kseqObjFq2->qual.s;

                        readAdded=true;
                        readsBuffer1.push_back(make_tuple(header1,read1,qual1));
                        readsBuffer2.push_back(make_tuple(header2,read2,qual2));
                    }
                    readsBufferIDS[i].push_back(readsBuffer1.size()-1);
                }
            }
            if(++readsProcessed %100000 == 0)
            {
                cout<<"Processed "<<readsProcessed<< "reads and saved "<<readsBuffer1.size()*2<<"reads"<<endl;
            }

        }
    }
    for(unsigned i=0; i<numVariants; i++)
    {
        ofstream readsOut1(outputFilename+"."+ to_string(i)+".1.fq");
        for(auto r: readsBufferIDS[i])
        {
            readsOut1<<"@"<< get<0>(readsBuffer1[r])<<"\n"
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
