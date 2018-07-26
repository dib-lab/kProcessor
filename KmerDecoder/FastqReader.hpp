#include "KmerDecoder.hpp"
#include <vector>
#include <seqan/seq_io.h>
#include <utility>
using namespace std;


class FastqReader: KmerDecoder{
private:
    seqan::SeqFileIn *seqIn;
public:
  FastqReader(string path);
  pair<string,string> readSeq();
  void readNSeq(vector<pair<string,string> >* res,uint64_t N=0);
  bool isEOF();
};
