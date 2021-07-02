#include <iostream>
#include <string>
#include "CLI11.hpp"
#include <vector>
#include <stdint.h>
#include <gqf.h>
#include "algorithms.hpp"
#include "Utils/kmer.h"
#include "Utils/utils.hpp"
#include <cmath>
#include "kDataFrame.hpp"

using namespace std;


int main(int argc, char *argv[]){
  if(argc!=2){
    cout<<"Help: ./dumpMain <kdataframe prefix>"<<endl;
    return 0;
    }
  string input=argv[1];
  kDataFrame* KF=kDataFrame::load(input);
  auto it=KF->begin();
  while(it!=KF->end())
  {
    cout<<it.getKmer()<<"\t"<<it.getkmerOrder()<<"\n";
    it++;
  }
  return 0;
}
