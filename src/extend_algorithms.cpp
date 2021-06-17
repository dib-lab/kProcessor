#include "extend_algorithms.hpp"

namespace kProcessor {

    uint64_t aggregate_count(kDataFrame *kf, const string &column_name) {

        any totalCountAny = kProcessor::aggregate(kf, (uint64_t) 0, [column_name](kmerRow it, any v) -> any {
            any count0;
            it.getColumnValue(column_name, count0);
            return (any) (any_cast<uint64_t>(v) + any_cast<uint64_t>(count0));
        });

        return any_cast<uint64_t>(totalCountAny);
    }

    kDataFrame *filter_zeroCounts(kDataFrame *res, uint32_t nAllSamples, const string &counts_col_prefix) {
        res = kProcessor::filter(res, [=](kmerRow r) -> bool {
            for (unsigned i = 0; i < nAllSamples; i++) {
                any count;
                r.getColumnValue(counts_col_prefix + ".0", count);
                if (any_cast<uint32_t>(count) > 0)
                    return true;
            }
            return false;
        });

        return res;
    }

    kDataFrame *transform_normalize(kDataFrame *kf, const string &column_name, uint64_t totalCount) {

        kf = kProcessor::transform(kf, [=, & column_name, & totalCount](kmerRow it) -> kmerRow {
            any count0;
            it.getColumnValue(column_name, count0);
            double normalized = (double) any_cast<uint32_t>(count0)  * (100000000.0) / totalCount;
            it.setColumnValue(column_name, (uint32_t) normalized);
            return it;
        });

        return kf;
    }

    kDataFrame *transform_foldchange(kDataFrame *res, const string &output_col_name, const string &counts_col_prefix,
                                     uint32_t nTestSamples, uint32_t nAllSamples) {

        uint32_t nControlSamples = nAllSamples - nTestSamples;
        res->addColumn(output_col_name, new vectorColumn<double>(res->size()));

        res = kProcessor::transform(res, [=](kmerRow it) -> kmerRow {
            unsigned i = 0;
            uint32_t sampleSum = 0;
            for (; i < nTestSamples; i++) {
                any count;
                string colName = "count." + to_string(i);
                it.getColumnValue(colName, count);
                sampleSum += any_cast<uint32_t>(count);
            }
            uint32_t controlSum = 0;
            for (; i < nAllSamples; i++) {
                any count;
                string colName = counts_col_prefix + "." + to_string(i);
                it.getColumnValue(colName, count);
                controlSum += any_cast<uint32_t>(count);
            }
            double sampleAVG = (double) sampleSum / (double) nAllSamples;
            double controlAVG = (double) controlSum / (double) nControlSamples;
            double foldChange = sampleAVG / controlAVG;
            it.setColumnValue(output_col_name, foldChange);
            return it;
        });

        return res;
    }


    unordered_map<uint32_t, vector<double>> aggregate_foldChangeByGene(kDataFrame *res, const string &colorColumn) {

        auto foldChangeByGene = new unordered_map<uint32_t, vector<double> >();
        any genesGatherAny = kProcessor::aggregate(res, foldChangeByGene, [=](kmerRow it, any v) -> any {
            auto dict = any_cast<unordered_map<uint32_t, vector<double>> *>(v);
            any foldChange_any;
            it.getColumnValue("foldChange", foldChange_any);
            any color_any;
            it.getColumnValue(colorColumn,color_any);
            double foldChange=any_cast<double>(foldChange_any);
            vector<uint32_t> color=any_cast<vector<uint32_t> >(color_any);
            for (auto c: color) {
                (*dict)[c].push_back(foldChange);
            }
            return (any) (dict);
        });

        return any_cast<unordered_map<uint32_t, vector<double>>>(foldChangeByGene);
    }
}