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

    app.add_option("-i,--input", inputPath,
                   "File containig a list of kDataframe indexes paths")->required();
//    app.add_option("-q", q,
//                   "Q size of the result MQF frame")->required();
    app.add_option("-t,--tempFolder", tmpFolder,
                   "Path for Temporary Folder");
    app.add_option("-o,--output", outPath,
                   "Output Path")->required();


    CLI11_PARSE(app, argc, argv);




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
    kProcessor::parallelJoin(filenames,kmersToKeep,1);
    cout<<"Merging finished"<<endl;

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
