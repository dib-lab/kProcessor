#ifndef KmerCounter_HPP
#define KmerCounter_HPP
#include <stdint.h>
#include <string>
#include <gqf.h>


void loadIntoMQF(std::string sequenceFilename, int k,int noThreads,QF * memoryMQF);

void dumpMQF(QF * memoryMQF,int ksize,std::string outputFilename);

#endif
