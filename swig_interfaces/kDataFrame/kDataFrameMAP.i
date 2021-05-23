class kDataFrameMAP : public kDataFrame
{
public:
  kDataFrameMAP();
  kDataFrameMAP(uint64_t ksize);
  kDataFrameMAP(readingModes RM, hashingModes HM, map<string, int> params);
  kDataFrame* getTwin();
  void reserve (uint64_t n);

  bool setCount(string kmer, uint64_t count);
  bool insert(string kmer);
  bool insert(string kmer, uint64_t count);
  uint64_t getCount(string kmer);
  bool erase(string kmer);

  uint64_t size();
  uint64_t max_size();
  float load_factor();
  float max_load_factor();
  kDataFrameIterator begin();
  kDataFrameIterator end();

  uint64_t bucket(string kmer);
  void save(string filePath);
  static kDataFrame *load(string filePath);

    ~kDataFrameMAP() {
        this->MAP.clear();
    }
};