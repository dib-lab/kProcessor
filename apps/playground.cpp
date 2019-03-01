#include <iostream>
#include <string>
#include "CLI11.hpp"
#include <vector>
#include <stdint.h>
#include <gqf.hpp>
#include "KmerCounter/KmerCounter.hpp"
#include "KmerCounter/kmer.h"
#include "Utils/utils.hpp"
#include "kDataFrame.hpp"
#include <cmath>

using namespace std;


int playground_main(int argc, char *argv[]){
  CLI::App app;
  string input_file="";

  app.add_option("-i,--input", input_file,
   "MQF file")->required();


  CLI11_PARSE(app, argc, argv);
  kDataFrameMQF* indexFrame=(kDataFrameMQF*)kDataFrame::load(input_file, "MQF");

  cout<<indexFrame->getCounter("TCCCACATGGATATAGACAAA")<<endl;
  indexFrame->setCounter("TCCCACATGGATATAGACAAA",86505);
  cout<<indexFrame->getCounter("TCCCACATGGATATAGACAAA")<<endl;
  return 0;
}
