#include "algorithms.hpp"

namespace kProcessor {

    uint64_t aggregate_count(kDataFrame* kf, const string& countColName) {
        any totalCountAny = kProcessor::aggregate(kf, (uint64_t)0, [countColName](kDataFrameIterator& it, any v) -> any {
            uint32_t count0;

            it.getColumnValue<uint32_t, vectorColumn<uint32_t> >(countColName, count0);
            return (any)(any_cast<uint64_t>(v) + (uint64_t)count0);
            });

        return any_cast<uint64_t>(totalCountAny);
    }

    void transform_normalize(kDataFrame* kf, const string& countColName, uint64_t totalCount) {

        kProcessor::transformInPlace(kf, [=](kDataFrameIterator& it) -> void {
            uint32_t count0;
            it.getColumnValue<uint32_t, vectorColumn<uint32_t> >(countColName, count0);
            double normalized = (double)count0 * (100000000.0) / totalCount;
            it.setColumnValue<uint32_t, vectorColumn<uint32_t> >(countColName, (uint32_t)normalized);
            });

    }

    kDataFrame* filter_zeroCounts(kDataFrame* res, uint32_t allDatasets) {
        res = kProcessor::filter(res, [=](kDataFrameIterator& r) -> bool {
            for (unsigned i = 0; i < allDatasets;i++) {
                uint32_t count;
                r.getColumnValue<uint32_t, vectorColumn<uint32_t> >("count" + to_string(i), count);
                if (count > 0)
                    return true;
            }
            return false;
            });
        return res;
    }

    void transform_foldchange(kDataFrame* res, uint32_t nSamples, uint32_t nControl, uint32_t allDatasets, const string& foldChangeColName) {
        kProcessor::transformInPlace(res, [=](kDataFrameIterator& it) -> void {
            unsigned i = 0;
            uint32_t sampleSum = 0;
            for (;i < nSamples; i++)
            {
                uint32_t count;
                string colName = "count" + to_string(i);
                it.getColumnValue<uint32_t, vectorColumn<uint32_t> >(colName, count);
                sampleSum += count;
            }
            uint32_t controlSum = 0;
            for (;i < allDatasets; i++)
            {
                uint32_t count;
                string colName = "count" + to_string(i);
                it.getColumnValue<uint32_t, vectorColumn<uint32_t> >(colName, count);
                controlSum += count;
            }
            double sampleAVG = (double)sampleSum / (double)nSamples;
            double controlAVG = (double)controlSum / (double)nControl;
            double foldChange = sampleAVG / controlAVG;
            it.setColumnValue<double, vectorColumn<double> >(foldChangeColName, foldChange);

            });
    }

    inline void _dedup(kDataFrameIterator &it, string colorColumnName, vector<string> &color){
        it.getColumnValue<vector<string>, deduplicatedColumn<vector<string>, StringColorColumn>>(colorColumnName, color);
    }

    unordered_map<string, vector<double>> aggregate_foldChangeByGene(kDataFrame* res, const string foldChangeColName, const string colorColumnName) {

        auto foldChangeByGene = unordered_map<string, vector<double>>();

        any genesGatherAny = kProcessor::aggregate(res, &foldChangeByGene, [=](kDataFrameIterator& it, any v) -> any {
            auto dict = any_cast<unordered_map<string, vector<double>>*>(v);
            double foldChange;
            it.getColumnValue<double, vectorColumn<double> >(foldChangeColName, foldChange);
            vector<string> color;
            _dedup(it, colorColumnName, color);
            it.getColumnValue<vector<string>, deduplicatedColumn<vector<string>, StringColorColumn>>(colorColumnName, color);
            for (auto c : color)
            {
                (*dict)[c].push_back(foldChange);
            }
            return (any)(dict);
            });

        return foldChangeByGene;
    }

}