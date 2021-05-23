class kDataFrameBMQF: public kDataFrame{
public:

  kDataFrameBMQF(uint64_t kSize);
  kDataFrameBMQF(readingModes RM, hashingModes HM, map<string, int> params);
  kDataFrameBMQF(uint64_t ksize,uint8_t q,uint8_t fixedCounterSize,uint8_t tagSize,double falsePositiveRate);
  kDataFrameBMQF(bufferedMQF* bufferedmqf,uint64_t ksize,double falsePositiveRate);

  ~kDataFrameBMQF(){
    bufferedMQF_destroy(bufferedmqf);
    delete bufferedmqf;
  }

  void reserve (uint64_t n);

  kDataFrame* getTwin();

  bool insert(string kmer,uint64_t count);
  bool insert(string kmer);
  bool insert(uint64_t kmer,uint64_t count);
  bool insert(uint64_t kmer);
  bool setCount(string kmer,uint64_t count);
  bool setCount(uint64_t kmer, uint64_t count);

  uint64_t getCount(string kmer);
  uint64_t getCount(uint64_t kmer);


  bool erase(string kmer);
  bool erase(uint64_t kmer);

  uint64_t size();

  uint64_t max_size();
  float load_factor();
  float max_load_factor();


  bufferedMQF* getBMQF(){
    return bufferedmqf;
  }

  void save(string filePath);
  static kDataFrame* load(string filePath);

  kDataFrameIterator begin();
  kDataFrameIterator end();
  kDataFrameIterator find(string kmer);
};