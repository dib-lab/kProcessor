#ifndef _defaultColumn_H_
#define _defaultColumn_H_

#include <vector>
#include <stdint.h>
#include <iostream>

using namespace std;


template<typename  T>
class vectorColumn{
public:
    vector<T> data;
    vectorColumn(){

    }
    ~vectorColumn(){

    }
    uint32_t  insertAndGetIndex(T item);
    T getWithIndex(uint32_t index);


};



#endif