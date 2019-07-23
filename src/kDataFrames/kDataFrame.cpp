#include "kDataFrame.hpp"
#include "Utils/kmer.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <limits>
#include <sstream>

using namespace std;

inline bool fileExists(const std::string &name) {
    ifstream f(name.c_str());
    return f.good();
}

kDataFrame::kDataFrame() {
    kSize = 31;
}

kDataFrame::kDataFrame(uint8_t k_size) {
    kSize = k_size;

}

bool kDataFrame::empty() {
    return this->size() == 0;
}

bool kDataFrame::insert(kmerRow k) {
    return this->insert(k.kmer, k.count);
}


kDataFrame *kDataFrame::load(string filePath) {
    if (fileExists(filePath + ".mqf"))
        return kDataFrameMQF::load(filePath);
    else if (fileExists(filePath + ".phmap"))
        return kDataFrameMAP::load(filePath);
    else
        throw std::runtime_error("Could not open kDataFrame file");
}
