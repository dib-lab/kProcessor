#include "KmerDecoder.hpp"
#include <vector>
#include <deque>
#include <seqan/seq_io.h>
#include <utility>
using namespace std;


class FastqReader: KmerDecoder{
private:
    seqan::SeqFileIn *seqIn;
public:
  FastqReader(string path);
  pair<string,string> readSeq();
  void readNSeq(deque<pair<string,string> >* res,uint64_t N=0);
  bool isEOF();
};

struct file_pointer {
	char* part{nullptr};
	char* part_buffer{nullptr};
	uint64_t size{0};
	uint64_t part_filled{0};
};


class FastqReaderSqueker: KmerDecoder{
private:
    FILE *seqIn;
    file_pointer* fp;
    char* part;
    void parseReads(deque<pair<string,string> >* res);
public:
  FastqReaderSqueker(string path);
  pair<string,string> readSeq();
  void readNSeq(deque<pair<string,string> >* res,uint64_t N=0);
  bool isEOF();
};
