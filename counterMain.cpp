#include <iostream>
#include <string>
#include "CLI11.hpp"
#include <vector>
#include <stdint.h>
#include <gqf.h>
#include "KmerCounter/KmerCounter.hpp"
#include "Utils/utils.hpp"

using namespace std;


int KmerCounter_main(int argc, char *argv[]){
  CLI::App app;
  vector<string> input_files;
  string outputMQF;
  string outputKmers="";
  uint64_t nslots=32768;
  uint64_t fixed_size_counter=1;
  double accuracy=1;
  int noThreads=1;
  uint64_t maxMemory=0;
  int k;

  app.add_option("-i,--input", input_files,
   "Sequence Files to count. can be fasta,fastq or BAM")->required()
  ->check(CLI::ExistingFile)->group("I/O");
  app.add_option("-o,--output", outputMQF,
   "Output MQF filename")->required()->group("I/O");
  app.add_option("-u,--output-kmers", outputKmers,
    "Output in the format of Kmer\tCount. Available only in Exact Counting")->group("I/O");

  app.add_option("-k,--kmer-length",k,"kmer length")->required()->group("MQF Options");
  app.add_option("-s,--no-slots",nslots,"Number of slots in MQF. Should be of power of two")->group("MQF Options");
  app.add_option("-f,--fixed-size-counter",fixed_size_counter,
  "Number of bits in Fixed-size counter size in MQF. Default 1")->group("MQF Options");
  app.add_option("-a,--accuracy",accuracy,
  "Accuracy of MQF. use 1 for exact counting and less than 1 for probalistic counting. Default 1")
  ->group("MQF Options");

  app.add_option("-t,--threads", noThreads,
   "Number of threads used in kmer counting. Default 1")->group("Misc");
  app.add_option("-m,--max-memory", maxMemory,
   "Max Memory allocated by KmerCounter. Default Unlimited")->group("Misc");



  CLI11_PARSE(app, argc, argv);

  QF qf;
  qf_init(&qf, nslots, 2*k+15, 0,fixed_size_counter, true, "", 2038074761);

  for(auto file: input_files)
    loadIntoMQF(file,k,noThreads,&qf);

  if(outputKmers!=""){
    dumpMQF(&qf,k,outputKmers);
  }
  qf_serialize(&qf,outputMQF.c_str());



  return 0;
}
