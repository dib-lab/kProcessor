#include "batchQuery.hpp"
#include "algorithms.hpp"

// ----------------- colored kDataFrame Batch Query -----------------


ckf_batchQuery::ckf_batchQuery(colored_kDataFrame *ckFrame, std::string filename,
                               std::map<std::string, int> parse_params, int chunk_size) {

    // parse_params["mode"] = 1 > Default: Kmers
    // parse_params["mode"] = 2 > Skipmers
    // parse_params["mode"] = 3 > Minimizers

    // Initialize kmerDecoder
    std::string mode = "kmers";
    bool check_mode = (parse_params.find("mode") != parse_params.end());
    if (check_mode) {
        if (parse_params["mode"] == 2) mode = "skipmers";
        else if (parse_params["mode"] == 3) mode = "minimizers";
    }

    this->KD = kProcessor::initialize_kmerDecoder(filename, chunk_size, mode, parse_params);
    this->ckFrame = ckFrame;
}

void ckf_batchQuery::next() {
    this->KD->next_chunk();
}

bool ckf_batchQuery::end() {
    return this->KD->end();
}

void ckf_batchQuery::flush_results() {
    this->results.clear();
}

std::unordered_map<std::string, vector<vector<uint32_t>>> ckf_batchQuery::get_transcripts() {
    this->flush_results();
    for (const auto &seq : *KD->getKmers()) {

        std::vector<vector<uint32_t>> temp_tr_ids;

        for (const auto &kmer : seq.second) {
            std::vector<uint32_t> temp_kmer_sources;
            this->ckFrame->getKmerSource(kmer.hash, temp_kmer_sources);
            temp_tr_ids.emplace_back(temp_kmer_sources);
        }
        this->results[seq.first] = temp_tr_ids;

    }
    return this->results;
}



// ----------------- kDataFrame Batch Query -----------------

kf_batchQuery::kf_batchQuery(colored_kDataFrame *ckFrame, std::string filename,
                             std::map<std::string, int> parse_params, int chunk_size) {

    // parse_params["mode"] = 1 > Default: Kmers
    // parse_params["mode"] = 2 > Skipmers
    // parse_params["mode"] = 3 > Minimizers

    // Initialize kmerDecoder

    std::string mode = "kmers";
    bool check_mode = (parse_params.find("mode") != parse_params.end());
    if (check_mode) {
        if (parse_params["mode"] == 2) mode = "skipmers";
        else if (parse_params["mode"] == 3) mode = "minimizers";
    }

    this->KD = kProcessor::initialize_kmerDecoder(filename, chunk_size, mode, parse_params);

    this->kFrame = ckFrame->getkDataFrame();

}

kf_batchQuery::kf_batchQuery(kDataFrame *kFrame, std::string filename,
                             std::map<std::string, int> parse_params, int chunk_size) {

    std::string mode = "kmers";
    bool check_mode = (parse_params.find("mode") != parse_params.end());

    if (check_mode) {
        if (parse_params["mode"] == 2) mode = "skipmers";
        else if (parse_params["mode"] == 3) mode = "minimizers";
    }

    this->KD = kProcessor::initialize_kmerDecoder(filename, chunk_size, mode, parse_params);
    this->kFrame = kFrame;
}

void kf_batchQuery::next() {
    this->KD->next_chunk();
}

bool kf_batchQuery::end() {
    return this->KD->end();
}

void kf_batchQuery::flush_results() {
    this->results.clear();
}

std::unordered_map<std::string, vector<uint32_t>> kf_batchQuery::get_transcripts() {
    this->flush_results();

    for (const auto &seq : *KD->getKmers()) {
        std::vector<uint32_t> temp_counts;

        for (const auto &kmer : seq.second) {
            temp_counts.emplace_back(this->kFrame->getCount(kmer.hash));
        }
        this->results[seq.first] = temp_counts;
    }

    return this->results;
}