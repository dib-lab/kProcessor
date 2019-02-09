#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <deque>
#include <utility>
#include <string>
using namespace std;

class KmerDecoder{
protected:
  uint64_t chunkSize;
public:
  KmerDecoder(){
  };
  pair<string,string> readSeq();
  void readNSeq(deque<pair<string,string> > *res,uint64_t N=0);
  bool isEOF();
};
