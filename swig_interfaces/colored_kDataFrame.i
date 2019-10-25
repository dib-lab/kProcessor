class colored_kDataFrame{
public:
  colored_kDataFrame();
  void addNewColor(uint32_t color, vector<uint32_t> & samplesIds);
  void setKmerColor(string kmer,uint32_t color);
  uint32_t getKmerColor(string kmer);
  vector<uint32_t> getSamplesIDForKmer(string kmer);
  vector<uint32_t> getSamplesIDForColor(uint32_t color);

  void colorKmer(string kmer,vector<uint32_t> & samplesIds);

  void setColorTable(colorTable* table);
  void setkDataFrame(kDataFrame* f);
  void save(string prefix);
  static colored_kDataFrame* load(string prefix);
  uint64_t getkSize();

  // Converting phmap to unordered_map, mainly for the python interface
  unordered_map<int, string> names_map();
  unordered_map<string, int> inverse_names_map();

  // Get the kDataFrame of the colored_kDataFrame
  kDataFrame * getkDataFrame();
};
