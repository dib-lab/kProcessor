#include <iostream>
#include <string>
#include "ThirdParty/CLI11.hpp"
#include <vector>
#include <stdint.h>
#include <gqf.h>
#include "KmerCounter/KmerCounter.hpp"
#include "KmerCounter/kmer.h"
#include "Utils/utils.hpp"
#include <cmath>

using namespace std;


int estimateMemory_main(int argc, char *argv[]){
  CLI::App app;
  vector<string> input_files;
  double fpr=0;
  uint64_t tagSize=0;
  int k;

  app.add_option("-i,--input", input_files,
   "Ntcard Files containing the approximate histogram of the kmers.")->required()
  ->check(CLI::ExistingFile);

  app.add_option("-r,--fpr",fpr,
  "False Positive Rate of MQF. use 0 for exact counting and less than 1 for probalistic counting. Default 0");
  app.add_option("-k,--kmer-length",k,"kmer length")->required();
  app.add_option("-t,--tag-size",tagSize,"Tag number of bits. Default 0");
  CLI11_PARSE(app, argc, argv);

  uint64_t res_noSlots, res_fixedSizeCounter, res_slotSize, res_memory;
  uint64_t slotSize;

  if(fpr==0){
    slotSize=2*k;
  }
  else if(fpr<1){
    slotSize=-(uint64_t)(log2(fpr))+1;
  }

  estimateMemRequirement(input_files[0],
      slotSize, tagSize,
     &res_noSlots, &res_fixedSizeCounter, &res_slotSize, &res_memory);

  cout<<"Number Slots = "<<res_noSlots<<endl
      <<"Fixed size counters= "<<res_fixedSizeCounter<<endl
      <<"Slot Size= "<<res_slotSize<<endl
      <<"Memory= "<<res_memory<<"KB"<<endl;
  return 0;
}
