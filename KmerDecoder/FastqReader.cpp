#include "FastqReader.hpp"
#include <seqan/seq_io.h>
#include <vector>
#include <utility>

FastqReader::FastqReader(string path){
    seqIn=(new seqan::SeqFileIn(path.c_str()));
    chunkSize=10000;
}
pair<string,string> FastqReader::readSeq(){
  // string read,id;
  // seqan::readRecord(id, read, *seqIn);
  // return make_pair(read,"");
  return make_pair(" "," ");
}

void FastqReader::readNSeq(vector<pair<string,string> >* res, uint64_t N){
  if(N==0)
    N=chunkSize;
  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::CharString> reads;
  seqan::readRecords(ids, reads, *seqIn,N);
  //cout<<string((char*)toCString(reads[0]))<<endl;
  uint64_t i=0;
  for(auto read:reads){
    (*res)[i++]=make_pair(string((char*)toCString(read)),"");
  }
}

bool FastqReader::isEOF(){
  return atEnd(*seqIn);
}
