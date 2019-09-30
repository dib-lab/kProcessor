namespace kProcessor{
void parseSequences(kmerDecoder * KD, kDataFrame* output);
void parseSequencesFromString(kmerDecoder *KD, string sequence,kDataFrame* output);
kDataFrame* kFrameUnion(const vector<kDataFrame*>& input);
kDataFrame* kFrameIntersect(const vector<kDataFrame*>& input);
kDataFrame* kFrameDiff(const vector<kDataFrame*>& input);
kmerDecoder* initialize_kmerDecoder(std::string filename, int chunkSize, std::string mode, std::map<std::string, int> params);
kmerDecoder* initialize_kmerDecoder(std::string mode, std::map<std::string, int> params);
kmerDecoder* initialize_kmerDecoder(int kSize, int hash_mode);
void kmerDecoder_setHashing(kmerDecoder * KD, int hash_mode, bool canonical);
colored_kDataFrame *index(kmerDecoder *KD, string names_fileName, kDataFrame *frame);
}