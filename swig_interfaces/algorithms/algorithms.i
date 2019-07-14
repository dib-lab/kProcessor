namespace kProcessor{
void parseSequences(string seqFileName,int nThreads,kDataFrame* output);
kDataFrame* kFrameUnion(const vector<kDataFrame*>& input);
kDataFrame* kFrameIntersect(const vector<kDataFrame*>& input);
kDataFrame* kFrameDiff(const vector<kDataFrame*>& input);
kmerDecoder* initialize_kmerDecoder(std::string filename, int chunkSize, std::string mode, std::map<std::string, int> params);
colored_kDataFrame *index(kmerDecoder *KD, string names_fileName, uint64_t Q = 27);
}