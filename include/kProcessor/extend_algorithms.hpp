#include "algorithms.hpp"

namespace kProcessor{

uint64_t aggregate_count(kDataFrame *kf, const string &column_name);
kDataFrame *filter_zeroCounts(kDataFrame *res, uint32_t nAllSamples, const string & counts_col_prefix);
kDataFrame *transform_normalize(kDataFrame *kf, const string &column_name, uint64_t totalCount);
kDataFrame *transform_foldchange(kDataFrame *res, const string &output_col_name, const string &counts_col_prefix, uint32_t nTestSamples, uint32_t nAllSamples);
unordered_map<uint32_t, vector<double>> aggregate_foldChangeByGene(kDataFrame *res, const string &colorColumn);

}