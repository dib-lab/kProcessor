//
// Created by mostafa on 4/6/20.
//

#include <iostream>
#include <string>
#include "kDataFrame.hpp"
#include <vector>
#include "algorithms.hpp"
#include "CLI11.hpp"
#include <omp.h>

using namespace std;

int main(int argc, char *argv[])
{
    CLI::App app;
    string inputPath;
    uint64_t  q;
    string tmpFolder="";
    string outPath;
    bool checkIndex=false;
    int nThreads=1;
    app.add_option("-i,--input", inputPath,
                   "File containig a list of kDataframe indexes paths")->required();
    app.add_option("-t,--tempFolder", tmpFolder,
                   "Path for Temporary Folder");
    app.add_option("-n,--numThreads", nThreads,
                   "Number of threads");
    app.add_option("-o,--output", outPath,
                   "Output Path")->required();
    app.add_flag("-c,--checkResults", checkIndex,
                 "Check if the index created successfully");


    CLI11_PARSE(app, argc, argv);

    omp_set_num_threads(nThreads);

    vector<string> filenames;

    vector<uint32_t> kmersToKeep;
    ifstream input(inputPath);
    string sample;
    uint32_t ii=0;

    while(input>>sample)
    {
        filenames.push_back(sample);
        kmersToKeep.push_back(ii++);
    }
    for(auto sample:filenames)
    {

    }
    kDataFrame* output=kProcessor::parallelJoin(filenames,kmersToKeep,nThreads);
    cout<<"Merging finished "<<endl;
    cout<<"Final number of kmers "<<output->size()<<endl;

    if(checkIndex)
    {
        ii=0;
        int Errors=100;

#pragma omp parallel for
        for(unsigned i=0;i<filenames.size();i++)
        {
            string sample=filenames[i];
            cout<<"Testing "<<sample<<endl;
            kDataFrame* kf=kDataFrame::load(sample);
            string colorColumnName="color"+to_string(i);
            string sampleColor="color";
            for(auto k:*kf)
            {
                vector<uint32_t> colorsCorrect;

                k.getColumnValue<vector<uint32_t >, deduplicatedColumn<prefixTrie,phmap::flat_hash_map<uint32_t,uint32_t>> >(sampleColor,colorsCorrect);

                vector<uint32_t> colorsQuered=
                        output->getKmerColumnValue<vector<uint32_t >, deduplicatedColumn<prefixTrie,phmap::flat_hash_map<uint32_t,uint32_t>> >(colorColumnName,k.getHashedKmer());
                if(colorsQuered!=colorsCorrect)
                {
                    vector<uint32_t> colorsCorrect;

                    k.getColumnValue<vector<uint32_t >, deduplicatedColumn< prefixTrie,phmap::flat_hash_map<uint32_t,uint32_t>> >(sampleColor,colorsCorrect);

                    vector<uint32_t> colorsQuered=
                            output->getKmerColumnValue<vector<uint32_t >, deduplicatedColumn< prefixTrie,phmap::flat_hash_map<uint32_t,uint32_t>> >(colorColumnName,k.getHashedKmer());
                    cout<<"Error Found at sample "<<sample<< " at kmer "<<k.getKmer()<<endl;
                    cout<<"Expected color is ";
                    for(auto c:colorsCorrect)
                    {
                        cout<<c<<" ";
                    }
                    cout<<endl<<"Found Color ";
                    for(auto c:colorsQuered)
                    {
                        cout<<c<<" ";
                    }
                    cout<<endl;
                    Errors--;
                    if(Errors<=0)
                        break ;
                }
            }
            cout<<"Finished "<<sample<<endl;
            delete kf;
        }
    }
    output->save(outPath);
    delete output;
    return 0;


    // output->save(outPath);
    //    cout<<"Saving Finished"<<endl;
    // delete output;
    // output=kDataFrame::load(outPath);
//    uint64_t testedKmers=0;
//    uint64_t failedKmers=0;
//    uint64_t notFoundKmers =0;
//    unordered_map<uint32_t ,uint32_t > sizeHist;
//     for(int i=0;i<filenames.size();i++)
//     {
//         ifstream inp(filenames[i]+".testkmers");
//         string kmer;
//         uint64_t count;
//         while(inp>>kmer>>count)
//         {
//     	  if(notFoundKmers+failedKmers==1000)
//     	    break;
//             testedKmers++;
//             vector<uint32_t> colors=output->getKmerDefaultColumnValue<vector<uint32_t >, mixVectors >(kmer);
//             sizeHist[colors.size()]++;
//             if(colors.size()==0)
//             {
//                 cout<<filenames[i]<<" KMER "<<kmer<<endl;
//                 notFoundKmers++;
//                 continue;
//             }
//             auto colorIt=find(colors.begin(),colors.end(),i);
//             if(colorIt==colors.end())
//             {
//                 cerr<<"Error detected in sample #"<<i <<" "<<
//                     filenames[i]<<" at kmer "<<kmer<<" Combination got "<<endl;
//                 for(auto c: colors)
//                     cerr<<c <<" ";
//                 cerr<<endl;
//
//                 failedKmers++;
//             }
//
//         }
//         inp.close();
//     }
//     cout<<"Numbers of tested Kmers = "<<testedKmers<<endl;
//     cout<<"Numbers of non found kmers = "<<notFoundKmers<<endl;
//     cout<<"Numbers of wrong combination = "<<failedKmers<<endl;

     //    cout<<"sizeHist"<<endl;
     //    for(auto s:sizeHist)
     //        cout<<s.first<<" : " <<s.second<<endl;


    return 0;
}
