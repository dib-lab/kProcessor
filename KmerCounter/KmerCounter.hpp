#ifndef KmerCounter_HPP
#define KmerCounter_HPP
#include <stdint.h>
#include <string>
#include <gqf.h>
#include "../HashUtils/hashutil.h"


void loadIntoMQF(std::string sequenceFilename, int k,int noThreads,Hasher *hasher,QF * memoryMQF);

void dumpMQF(QF * memoryMQF,int ksize,std::string outputFilename);

#endif
