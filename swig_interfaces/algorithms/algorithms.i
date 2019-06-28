namespace kProcessor{
void parseSequences(string seqFileName,int nThreads,kDataFrame* output);
kDataFrame* kFrameUnion(const vector<kDataFrame*>& input);
kDataFrame* kFrameIntersect(const vector<kDataFrame*>& input);
kDataFrame* kFrameDiff(const vector<kDataFrame*>& input);
colored_kDataFrame* index(string seqFileName,string namesFiles,uint64_t kSize);
colored_kDataFrame* index(string seqFileName, string namesFiles, uint64_t kSize, uint64_t Q);
}