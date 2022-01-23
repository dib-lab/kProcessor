//
// Created by mostafa on 4/6/20.
//

#include <iostream>
#include <string>
#include "kDataFrame.hpp"
#include <vector>
#include "algorithms.hpp"
#include "omp.h"

using namespace std;

int main(int argc, char *argv[]) {
    string inputIndex=argv[1];
    string outputIndex=argv[2];
    int nThreads=atoi(argv[3]);
    omp_set_num_threads(nThreads);
    kDataFrame* KF=kDataFrame::load(inputIndex);
    auto prevColor=(deduplicatedColumn< mixVectors>*)KF->columns["color"];
    auto newColor=new deduplicatedColumn<prefixTrie>();
    newColor->index=prevColor->index;
    newColor->values=new prefixTrie(prevColor->values);
    KF->columns["color"]=newColor;
    KF->save(outputIndex);
    return 0;
//    string inputList = argv[1];
//    string framePath = argv[2];
//    int numThreads=atoi(argv[3]);
//    string outputColumn = argv[4];
//    omp_set_num_threads(numThreads);
//    vector<string> filenames;
//    vector<kDataFrame *> frames;
//
//    string sample;
//    ifstream input(inputList);
//    while (input >> sample) {
//        filenames.push_back(sample);
//    }
//
//    kDataFrame *indexFrame = kDataFrame::load(framePath);
//    ((mixVectors *) indexFrame->getDefaultColumn())->explainSize();
//    //((mixVectors*)indexFrame->getDefaultColumn())->sortColors();
//
//    uint64_t testedKmers = 0;
//    uint64_t failedKmers = 0;
//    uint64_t notFoundKmers = 0;
//
//    mixVectors *qColumn = ((mixVectors *) indexFrame->getDefaultColumn());
//    prefixTrie *pColumn = new prefixTrie(qColumn);
//    pColumn->serialize(outputColumn);
//
//
//    pColumn->explainSize();
//    indexFrame->changeDefaultColumnType(pColumn);
//   // indexFrame->save(framePath + ".prefixTrie");
//    //delete indexFrame;
//    //indexFrame = kDataFrame::load(framePath + ".prefixTrie");
//
//    for (int i = 0; i < filenames.size(); i++) {
//        //  cerr<<"Testing "<<filenames[i]+".testkmers"<<endl;
//        ifstream inp(filenames[i] + ".testkmers");
//        string kmer;
//        uint64_t count;
//        while (inp >> kmer >> count && failedKmers<1000) {
//            testedKmers++;
//            vector<uint32_t> colors = indexFrame->getKmerDefaultColumnValue<vector<uint32_t>, prefixTrie>(
//                    kmer);
//
//            if (colors.size() == 0) {
//                notFoundKmers++;
//                continue;
//            }
//            auto colorIt = find(colors.begin(), colors.end(), i);
//            if (colorIt == colors.end()) {
//                cerr << "Error detected in sample #" << i << " " <<
//                     filenames[i] << " at kmer " << kmer << " Combination got " << endl;
//                for (auto c: colors)
//                    cerr << c << " ";
//                cerr << endl;
//
//                indexFrame->changeDefaultColumnType(qColumn);
//                vector<uint32_t> colors = indexFrame->getKmerDefaultColumnValue<vector<uint32_t>, mixVectors>(
//                        kmer);
//                cerr << "It should be ";
//                for (auto c: colors)
//                    cerr << c << " ";
//                cerr << endl;
//                failedKmers++;
//                indexFrame->changeDefaultColumnType(pColumn);
//            }
//
//        }
//        inp.close();
//    }
//    cout << "Numbers of tested Kmers = " << testedKmers << endl;
//    cout << "Numbers of non found kmers = " << notFoundKmers << endl;
//    cout << "Numbers of wrong combination = " << failedKmers << endl;
//    return 0;
}
