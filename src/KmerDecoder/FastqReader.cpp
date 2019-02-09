#include "KmerDecoder/FastqReader.hpp"
#include <seqan/seq_io.h>
#include <vector>
#include <deque>
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

void FastqReader::readNSeq(deque<pair<string,string> >* res, uint64_t N){
  if(N==0)
    N=chunkSize;
  seqan::StringSet<seqan::CharString> ids;
  seqan::StringSet<seqan::CharString> reads;
  seqan::readRecords(ids, reads, *seqIn,N);
  res->clear();
  //cout<<string((char*)toCString(reads[0]))<<endl;
  uint64_t i=0;
  for(auto read:reads){
    res->push_back(make_pair(string((char*)toCString(read)),""));
  }
}

bool FastqReader::isEOF(){
  return atEnd(*seqIn);
}


FastqReaderSqueker::FastqReaderSqueker(string path){
  seqIn = fopen(path.c_str(), "rb");
  chunkSize=10000;
  uint32_t OVERHEAD_SIZE = 65535;
  fp=new file_pointer;
  fp->part_buffer = new char[OVERHEAD_SIZE];
  uint64_t part_size = 1ULL << 23;
  part = (char *)malloc((part_size + OVERHEAD_SIZE)*sizeof(char));
}

pair<string,string> FastqReaderSqueker::readSeq(){
  // string read,id;
  // seqan::readRecord(id, read, *seqIn);
  // return make_pair(read,"");
  return make_pair(" "," ");
}

/* move the pointer to the end of the next newline. */
bool skip_next_eol(char *part, int64_t &pos, int64_t max_pos)
{
	int64_t i;
	for(i = pos; i < max_pos-2; ++i)
		if((part[i] == '\n' || part[i] == '\r') && !(part[i+1] == '\n' ||
																								 part[i+1] == '\r'))
			break;

	if(i >= max_pos-2)
		return false;
	pos = i+1;

	return true;
}
void FastqReaderSqueker::parseReads(deque<pair<string,string> >* res)
{
  auto fs = fp->part;
  auto fe = fp->part;
  auto end = fs + fp->size;

  while (fs && fs!=end) {
    fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore the first line
    fs++; // increment the pointer

    fe = static_cast<char*>(memchr(fs, '\n', end-fs)); // read the read
    string read(fs, fe-fs);
    //(*res)[readID++]=make_pair(read,"");
    res->push_back(make_pair(read,""));
    fs = ++fe;		// increment the pointer
    fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one line
    fs++; // increment the pointer
    fs = static_cast<char*>(memchr(fs, '\n', end-fs)); // ignore one more line
    fs++; // increment the pointer
  }
}

void FastqReaderSqueker::readNSeq(deque<pair<string,string> >* res, uint64_t N){


  if(N==0)
  N=chunkSize;
  uint64_t readID=0;
  char *& _part = (fp->part);
  uint64_t& _size = fp->size;
  char*& part_buffer = (fp->part_buffer);
  uint64_t& part_filled = fp->part_filled;

  uint32_t OVERHEAD_SIZE = 65535;
  uint64_t part_size = 1ULL << 23;
  res->clear();
  memcpy(part, part_buffer, part_filled);
  if(isEOF())
    return;
  uint64_t readed = 0;
  readed = fread(part+part_filled, 1, part_size, seqIn);
  int64_t total_filled = part_filled + readed;
  int64_t i;
  if(part_filled >= OVERHEAD_SIZE)
  {
    cout << "Error: Wrong input file!\n";
    exit(EXIT_FAILURE);
  }
  if(isEOF())
  {
    _part = part;
    _size = total_filled;
    part = NULL;
    parseReads(res);
    return;
  }
  // Looking for a FASTQ record at the end of the area
  {
    int64_t line_start[9];
    int32_t j;
    i = total_filled - OVERHEAD_SIZE / 2;
    for(j = 0; j < 9; ++j)
    {
      if(!skip_next_eol(part, i, total_filled))
      break;
      line_start[j] = i;
    }
    _part = part;
    if(j < 9)
    _size = 0;
    else
    {
      int k;
      for(k = 0; k < 4; ++k)
      {
        if(part[line_start[k]+0] == '@' && part[line_start[k+2]+0] == '+')
        {
          if(part[line_start[k+2]+1] == '\n' || part[line_start[k+2]+1] == '\r')
          break;
          if(line_start[k+1]-line_start[k] == line_start[k+3]-line_start[k+2] &&
            memcmp(part+line_start[k]+1, part+line_start[k+2]+1,
              line_start[k+3]-line_start[k+2]-1) == 0)
              break;
            }
          }
          if(k == 4)
          _size = 0;
          else
          _size = line_start[k];
        }
      }

      copy(_part+_size, _part+total_filled, part_buffer);
      part_filled = total_filled - _size;
      parseReads(res);
}


bool FastqReaderSqueker::isEOF(){
  return feof(seqIn) != 0;
}
