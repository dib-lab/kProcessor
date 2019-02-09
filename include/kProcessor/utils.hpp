#ifndef UTILS_HPP
#define UTILS_HPP
#include <stdint.h>
#include <string>
#include <gqf.hpp>
#include <seqan/seq_io.h>

using namespace seqan;
// void removeReadsWithN(std::string inputFilename,std::string outputFilename)
// {
//   SeqFileIn seqFileIn(inputFilename.c_str());
//   SeqFileOut seqFileOut(outputFilename.c_str());
//   CharString id;
//   CharString quals;
//   std::string readT;
//
//   while(!atEnd(seqFileIn))
//   {
//     bool hasN=false;
//     readRecord(id, readT,quals, seqFileIn);
//     for(uint64_t i=0;i<readT.size();i++)
//     {
//       if(readT[i]=='N')
//       {
//         hasN=true;
//         break;
//       }
//     }
//     if(!hasN)
//     {
//       writeRecord(seqFileOut,id,readT,quals);
//     }
//   }
//
//
// }



#endif
