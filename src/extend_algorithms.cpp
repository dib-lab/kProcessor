#include "extend_algorithms.hpp"

namespace kProcessor {

    uint64_t aggregate_count(kDataFrame *kf, const string &column_name) {

        any totalCountAny = kProcessor::aggregate(kf, (uint64_t) 0, [column_name](kmerRow it, any v) -> any {
            uint32_t count0;
            it.getColumnValue<uint32_t, vectorColumn<uint32_t> >(column_name, count0);
            return (any) (any_cast<uint64_t>(v) + (uint64_t) count0);
        });

        return any_cast<uint64_t>(totalCountAny);
    }

    kDataFrame *filter_zeroCounts(kDataFrame *res, uint32_t nAllSamples, const string &counts_col_prefix) {
        res = kProcessor::filter(res, [=](kmerRow r) -> bool {
            for (unsigned i = 0; i < nAllSamples; i++) {
                uint32_t count;
                r.getColumnValue<uint32_t, vectorColumn<uint32_t> >(counts_col_prefix + ".0", count);
                if (count > 0)
                    return true;
            }
            return false;
        });

        return res;
    }

    kDataFrame *transform_normalize(kDataFrame *kf, const string &column_name, uint64_t totalCount) {

        kf = kProcessor::transform(kf, [=, & column_name, & totalCount](kDataFrameIterator& it) -> void {
            uint32_t count0;
            it.getColumnValue<uint32_t, vectorColumn<uint32_t> >(column_name, count0);
            double normalized = (double) count0 * (100000000.0) / totalCount;
            it.setColumnValue<uint32_t, vectorColumn<uint32_t> >(column_name, (uint32_t) normalized);
        });

        return kf;
    }

    kDataFrame *transform_foldchange(kDataFrame *res, const string &output_col_name, const string &counts_col_prefix,
                                     uint32_t nTestSamples, uint32_t nAllSamples) {

        uint32_t nControlSamples = nAllSamples - nTestSamples;
        res->addColumn(output_col_name, new vectorColumn<double>(res->size()));

        res = kProcessor::transform(res, [=](kDataFrameIterator& it) -> kmerRow {
            unsigned i = 0;
            uint32_t sampleSum = 0;
            for (; i < nTestSamples; i++) {
                uint32_t count;
                string colName = "count." + to_string(i);
                it.getColumnValue<uint32_t, vectorColumn<uint32_t> >(colName, count);
                sampleSum += count;
            }
            uint32_t controlSum = 0;
            for (; i < nAllSamples; i++) {
                uint32_t count;
                string colName = counts_col_prefix + "." + to_string(i);
                it.getColumnValue<uint32_t, vectorColumn<uint32_t> >(colName, count);
                controlSum += count;
            }
            double sampleAVG = (double) sampleSum / (double) nAllSamples;
            double controlAVG = (double) controlSum / (double) nControlSamples;
            double foldChange = sampleAVG / controlAVG;
            it.setColumnValue<double, vectorColumn<double> >(output_col_name, foldChange);

        });

        return res;
    }


    unordered_map<uint32_t, vector<double>> aggregate_foldChangeByGene(kDataFrame *res, const string &colorColumn) {

        auto foldChangeByGene = new unordered_map<string, vector<double> >();
        any genesGatherAny = kProcessor::aggregate(res, foldChangeByGene, [=](kmerRow it, any v) -> any {
            auto dict = any_cast<unordered_map<string, vector<double>> *>(v);
            double foldChange;
            it.getColumnValue<double, vectorColumn<double> >("foldChange", foldChange);
            vector<string> color;
            it.getColumnValue<vector<string>, deduplicatedColumn<StringColorColumn> >(colorColumn,
                                                                                                          color);
            for (auto c: color) {
                (*dict)[c].push_back(foldChange);
            }
            return (any) (dict);
        });

        return any_cast<unordered_map<uint32_t, vector<double>>>(foldChangeByGene);
    }
}