#include <iostream>
#include <string>
#include "CLI11.hpp"
#include "kDataFrame.hpp"
#include <cstdlib>
#include "set"
#include <fstream>
#include <algorithm>
#include "algorithms.hpp"
using namespace std;


int main(int argc, char *argv[]){
//    CLI::App app;
//    string input_file;
//
//
//    app.add_option("-i,--input", input_file,
//                   "KdataFrame")->required();
//
//
//
//    CLI11_PARSE(app, argc, argv);
//    //kDataFrame* frame=kDataFrame::load(input_file);
//    kDataFrame* frame =new kDataFrameMQF();
//    kProcessor::loadFromKMC(frame, "input.unitigs");
//
//    cout<<"size "<<frame->size()<<endl;
//    frame->addColumn("visited",new vectorColumn<bool>(frame->size()));
//
//    for(auto kmer:*frame)
//    {
//        cout<<kmer.getKmer()<<endl;
//        if(!frame->getKmerColumnValue<bool,vectorColumn<bool> >("visited",kmer.kmer)) {
//            dbgIterator it = frame->getDBGIterator(kmer.kmer);
//            bool moreKmers=true;
//            while(moreKmers)
//            {
//                cout<<" "<<it.currentKmer<<endl;
//                frame->setKmerColumnValue<bool,vectorColumn<bool> >("visited",it.currentKmer,true);
//                moreKmers=false;
//                for(unsigned int i=0;i<it.nextFwdKmers.size();i++)
//                {
//                    if(!frame->getKmerColumnValue<bool,vectorColumn<bool> >("visited",it.nextFwdKmers[i]))
//                    {
//                        moreKmers=true;
//                        it.nextFWD(i);
//                        break;
//                    }
//                }
//            }
//        }
//
//    }
//    cout<<"Navigate dbg is completed"<<endl;






    return 0;
}
