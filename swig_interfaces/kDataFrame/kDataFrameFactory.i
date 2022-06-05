class kDataFrameFactory {
public:
    static kDataFrame* loadMQF(string filePath);
    static kDataFrame* createMQF(uint32_t kSize, uint32_t numKmers = 10000);
    static kDataFrame* createMQF(kDataFrame* kframe);

    static kDataFrame* loadPHMAP(string filePath);
    static kDataFrame* createPHMAP(uint32_t kSize, uint32_t numKmers = 10000);
    static kDataFrame* createPHMAP(uint64_t ksize, hashingModes hash_mode);


    static kDataFrame* loadMAP(string filePath);
    static kDataFrame* createMAP(uint32_t kSize, uint32_t numKmers = 10000);

    static kDataFrame* loadBtree(string filePath);
    static kDataFrame* createBtree(uint32_t kSize, uint32_t numKmers = 10000);

    static kDataFrame* loadBMQF(string filePath);
    static kDataFrame* createBMQF(uint32_t kSize, string filePath, uint32_t numKmers = 10000);
    static kDataFrame* createBMQF(kDataFrame* kframe, string filePath);


    static kDataFrame* loadBlight(string filePath);
    static kDataFrame* createBlight(uint32_t kSize, string filePath);

    static kDataFrame* loadSSHASH(string filePath);
    static kDataFrame* createSSHASH(uint32_t kSize, string filePath);

};


class kDataFrameUtility {
public:
    static void deleteMemoryBufferBMQF(kDataFrame* frame);
};