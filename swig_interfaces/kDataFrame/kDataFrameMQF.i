class kDataFrameMQF: public kDataFrame{

public:
  kDataFrameMQF();
  kDataFrameMQF(uint64_t kSize);

  ~kDataFrameMQF(){
    qf_destroy(mqf);
    delete mqf;
  }

  void reserve (uint64_t n);

  kDataFrame* getTwin();

  bool setCount(string kmer,uint64_t count);
  bool insert(string kmer,uint64_t count);
  bool insert(string kmer);
  uint64_t count(string kmer);
  bool erase(string kmer);

  uint64_t size();
  uint64_t max_size();
  float load_factor();
  float max_load_factor();

  QF* getMQF(){
    return mqf;
  }

  void save(string filePath);
  static kDataFrame* load(string filePath);

  kDataFrameIterator begin();
  kDataFrameIterator end();
};