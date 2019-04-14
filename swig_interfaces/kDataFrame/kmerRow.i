class kmerRow{
public:
  string kmer;
  uint64_t hashedKmer;
  uint64_t count;
  kmerRow(){
    kmer="";
    hashedKmer=0;
    count=0;
  }
  kmerRow(string kmer,uint64_t hashedKmer,uint64_t count)
  {
    this->kmer=kmer;
    this->hashedKmer=hashedKmer;
    this->count=count;
  }
  kmerRow(const kmerRow& other)
  {
    kmer=other.kmer;
    hashedKmer=other.hashedKmer;
    count=other.count;
  }

  kmerRow copy(const kmerRow& other)
  {
    return * (new kmerRow(other));
  }

 bool operator == (kmerRow& other)
  {
    return hashedKmer==other.hashedKmer;
  }

  bool operator < (kmerRow& other)
  {
    return hashedKmer<other.hashedKmer;
  }
  bool operator > ( kmerRow& other)
  {
    return hashedKmer>other.hashedKmer;
  }
  bool operator < (const kmerRow& other)
  {
    return hashedKmer<other.hashedKmer;
  }
  bool operator > (const kmerRow& other)
  {
    return hashedKmer>other.hashedKmer;
  }

};