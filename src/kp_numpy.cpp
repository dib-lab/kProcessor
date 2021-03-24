#include "kp_numpy.hpp"
#include "kDataFrame.hpp"


void kDataFramePHMAP::to_numpy(unsigned long long *rangevec, int rows, int cols) {

    int i = 0;
    
    for (auto itr = this->MAP.begin(); itr != this->MAP.end(); ++itr) {
        rangevec[i++] = itr->first;
        rangevec[i++] = itr->second;
    }

}