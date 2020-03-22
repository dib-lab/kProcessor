#include <parallel_hashmap/phmap.h>
#include "colored_kDataFrame.hpp"

using namespace std;

// ----------------- kDataFrame Batch Query -----------------

class kf_batchQuery {
private:
    kmerDecoder *KD;
    kDataFrame *kFrame;
    unordered_map<string, vector<uint32_t>> results;

    void flush_results();

public:
    kf_batchQuery(colored_kDataFrame *ckFrame, string filename, map<string, int> parse_params,
                  int chunk_size = 500);

    kf_batchQuery(kDataFrame *kFrame, string filename, map<string, int> parse_params,
                  int chunk_size = 500);

    void next();

    bool end();

    unordered_map<string, vector<uint32_t>> get_transcripts();

};

// ----------------- colored kDataFrame Batch Query -----------------

class ckf_batchQuery {
private:
    kmerDecoder *KD;
    colored_kDataFrame *ckFrame;
    unordered_map<string, vector<vector<uint32_t>>> results;

    void flush_results();

public:
    ckf_batchQuery(colored_kDataFrame *ckFrame, string filename, map<string, int> parse_params,
                   int chunk_size = 500);

    void next();

    bool end();

    unordered_map<string, vector<vector<uint32_t>>> get_transcripts();

};