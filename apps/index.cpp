//
// Created by mostafa on 4/6/20.
//

#include <iostream>
#include <string>
#include "kDataFrame.hpp"
#include <vector>
#include "algorithms.hpp"
#include "CLI11.hpp"
#include <typeinfo>

using namespace std;

int main(int argc, char *argv[])
{
    CLI::App app;
    string inputPath;
    uint64_t  q;
    string tmpFolder="";
    string outPath;
    uint32_t num_vectors=20;
    uint32_t vector_size=1000000;

    app.add_option("-i,--input", inputPath,
                   "File containig a list of kDataframe paths")->required();
    app.add_option("-q", q,
                   "Q size of the result MQF frame")->required();
    app.add_option("-t,--tempFolder", tmpFolder,
                   "Path for Temporary Folder");
    app.add_option("-o,--output", outPath,
                   "Output Path")->required();
    app.add_option("-n,--numVectors", num_vectors,
                   "Number of vectors in the output mixvectors column");
    app.add_option("-s,--vectorSize", vector_size,
                   "size of vectors in the output mixvectors column");


    CLI11_PARSE(app, argc, argv);




    vector<string> filenames;
    vector<kDataFrame*> frames;
    string sample;
    ifstream input(inputPath);


    while(input>>sample)
    {
        filenames.push_back(sample);
        frames.push_back(kDataFrame::load(sample));
        if(dynamic_cast<kDataFrameBMQF*>(frames.back()))
        {
            ((kDataFrameBMQF*)frames.back())->deleteMemoryBuffer();
        }
	cerr<<"sample "<<sample<<" loaded"<<endl; 
    }
    uint64_t  kSize=frames[0]->getkSize();
    for(auto f:frames)
    {
        if(f->getkSize()!=kSize)
        {
            cerr<<"All Kdataframes should have the same kSize "<<endl;
            return -1;
        }
    }

    kDataFrame* output= new kDataFrameMQF(kSize,q,integer_hasher);
    kProcessor::indexPriorityQueue(frames,tmpFolder,output,num_vectors,vector_size);
    cout<<"Indexing Finished"<<endl;

    // output->save(outPath);
    //    cout<<"Saving Finished"<<endl;
    // delete output;
    // output=kDataFrame::load(outPath);
    uint64_t testedKmers=0;
    uint64_t failedKmers=0;
    uint64_t notFoundKmers =0;
    unordered_map<uint32_t ,uint32_t > sizeHist;
     for(int i=0;i<filenames.size();i++)
     {
         ifstream inp(filenames[i]+".testkmers");
         string kmer;
         uint64_t count;
         while(inp>>kmer>>count)
         {
     	  if(notFoundKmers+failedKmers==1000)
     	    break;
             testedKmers++;
             vector<uint32_t> colors=output->getKmerDefaultColumnValue<vector<uint32_t >, mixVectors >(kmer);
             sizeHist[colors.size()]++;
             if(colors.size()==0)
             {
                 cout<<filenames[i]<<" KMER "<<kmer<<endl;
                 notFoundKmers++;
                 continue;
             }
             auto colorIt=find(colors.begin(),colors.end(),i);
             if(colorIt==colors.end())
             {
                 cerr<<"Error detected in sample #"<<i <<" "<<
                     filenames[i]<<" at kmer "<<kmer<<" Combination got "<<endl;
                 for(auto c: colors)
                     cerr<<c <<" ";
                 cerr<<endl;

                 failedKmers++;
             }

         }
         inp.close();
     }
     cout<<"Numbers of tested Kmers = "<<testedKmers<<endl;
     cout<<"Numbers of non found kmers = "<<notFoundKmers<<endl;
     cout<<"Numbers of wrong combination = "<<failedKmers<<endl;

     //    cout<<"sizeHist"<<endl;
     //    for(auto s:sizeHist)
     //        cout<<s.first<<" : " <<s.second<<endl;

    output->save(outPath);
    cout<<"Saving Finished"<<endl;
    delete output;
    for(auto f:frames)
        delete f;

    return 0;
}
