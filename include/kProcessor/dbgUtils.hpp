//
// Created by Mostafa on 10/13/2021.
//
#ifndef dbgUtils_HPP
#define dbgUtils_HPP
#include "kDataFrame.hpp"

namespace kProcessor {
    class Assembler{
    private:
        kDataFrame* kframe;
        vector<string> links;
    protected:
        virtual bool edgeExists(string kmer);
        virtual vector<string> getLinks(kDataFrameIterator& it);
    public:
        Assembler();
        Assembler(kDataFrame* kframe);
        string getContig(string seed);

    };
}

#endif
