namespace kProcessor{

kDataFrame* kFrameUnion(const vector<kDataFrame*>& input);
kDataFrame* kFrameIntersect(const vector<kDataFrame*>& input);
kDataFrame* kFrameDiff(const vector<kDataFrame*>& input);
void countKmersFromString(kDataFrame * kframe, std::map<std::string, int> params, string sequence);
void countKmersFromFile(kDataFrame * kframe, std::map<std::string, int> params, string filename, int chunk_size = 1000);
colored_kDataFrame * index(kDataFrame *kframe, std::map<std::string, int> params, string filename, int chunks, string names_fileName);

//void parseSequences(kmerDecoder * KD, kDataFrame* output);
//void countKmersFromString(kmerDecoder *KD, string sequence,kDataFrame* output);
//colored_kDataFrame *index(kmerDecoder *KD, string names_fileName, kDataFrame *frame);
//kmerDecoder* initialize_kmerDecoder(std::string filename, int chunkSize, std::string mode, std::map<std::string, int> params);
//kmerDecoder* initialize_kmerDecoder(std::string mode, std::map<std::string, int> params);
//kmerDecoder* initialize_kmerDecoder(int kSize, int hash_mode);
//void kmerDecoder_setHashing(kmerDecoder * KD, int hash_mode, bool canonical);

}