class dbgIterator{
public:
    kDataFrame* frame;
    string currentKmer;
    vector<string> nextFwdKmers;
    vector<string> nextRevKmers;

    dbgIterator();
    dbgIterator(kDataFrame*,string kmer);

    dbgIterator(const dbgIterator& other);
    dbgIterator& operator= (const dbgIterator& other);

    void nextFWD(uint32_t index);
    void nextREV(uint32_t index);
};