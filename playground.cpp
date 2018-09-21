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



int playground_main(int argc, char *argv[]){
  CLI::App app;
  string input_file="";

  app.add_option("-i,--input", input_file,
   "MQF file")->required()->check(CLI::ExistingFile);


  CLI11_PARSE(app, argc, argv);

  QF qf;
  qf_deserialize(&qf,input_file.c_str());
  double res=(double)qf.metadata->noccupied_slots/(double)qf.metadata->xnslots;
  cout<<qf_space(&qf)<<" "<<qf.metadata->noccupied_slots<<" "<<qf.metadata->maximum_occupied_slots<<endl;
  QF memqf;
  qf_init(&memqf, qf.metadata->nslots, qf.metadata->key_bits  , 0,qf.metadata->fixed_counter_size, true, "", 2038074761);
    qf_migrate(&qf,&memqf);
  
  cout<<"Loaded"<<endl;
  cout<<"Space= "<<qf_space(&qf)<<" occpuied slots= "<<qf.metadata->noccupied_slots<<" max solts= "<<qf.metadata->maximum_occupied_slots<<  " slots used by function = "<<slotsUsedInCounting(&qf)<<endl;
  cout<<"Space= "<<qf_space(&memqf)<<" occpuied slots= "<<memqf.metadata->noccupied_slots<<" max solts= "<<memqf.metadata->maximum_occupied_slots<<  " slots used by function = "<<slotsUsedInCounting(&memqf)<<endl;
  cout<<qf_space(&memqf)<<endl;
  cout<<memqf.metadata->noccupied_slots<<endl;
  cout<<"finish"<<endl;
  return 0;
}
