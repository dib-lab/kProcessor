class kDataFrame{
public:
  kDataFrame();
  kDataFrame(uint8_t kSize);
  virtual ~kDataFrame(){
  }
  virtual kDataFrame* getTwin()=0;
  virtual void reserve (uint64_t n )=0;
  virtual bool insert(string kmer)=0;
  virtual bool insert(string kmer,uint64_t N)=0;  
  bool insert(kmerRow k);
  virtual bool setCount(string kmer,uint64_t N)=0;
  virtual uint64_t count(string kmer)=0;
  virtual bool erase(string kmer)=0;
  virtual uint64_t size()=0;
  virtual uint64_t max_size()=0;
  bool empty();
  virtual float load_factor()=0;
  virtual float max_load_factor()=0;
  virtual kDataFrameIterator begin()=0;
  virtual kDataFrameIterator end()=0;
  virtual void save(string filePath)=0;
  Hasher* getHasher(){
    return hasher;
  };
  static kDataFrame* load(string filePath);
  uint64_t getkSize(){return kSize;}
  void setkSize(uint64_t k){kSize=k;}
};
