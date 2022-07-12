namespace kProcessor{

kDataFrame* kFrameUnion(const vector<kDataFrame*>& input);
kDataFrame* kFrameIntersect(const vector<kDataFrame*>& input);
kDataFrame* kFrameDiff(const vector<kDataFrame*>& input);
void countKmersFromFile(kDataFrame * kframe, string filename, int chunk_size = 1000);
void countKmersFromString(string sequence,kDataFrame* output);

colored_kDataFrame * index(kDataFrame * frame, string filename, int chunk_size, string names_fileName);
}