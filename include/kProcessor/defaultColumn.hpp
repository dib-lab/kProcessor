#ifndef _defaultColumn_H_
#define _defaultColumn_H_

#include <vector>
#include <stdint.h>
#include <iostream>

using namespace std;

//to add new column create a column class and add the name in getContainer Bynamefn
class Column{
public:
    Column(){}
    virtual ~Column(){}

    static Column* getContainerByName(size_t name);

    virtual void serialize(string filename)=0;
    virtual void deserialize(string filename)=0;


};

template<typename  T>
class vectorColumn: public Column{
public:
    vector<T> dataV;
    vectorColumn(){

    }
    vectorColumn(uint32_t size){
        dataV=vector<T>(size);
    }
    ~vectorColumn(){

    }
    uint32_t  insertAndGetIndex(T item);
    T getWithIndex(uint32_t index);

    void insert(T item,uint32_t index);
    T get(uint32_t index);

    void serialize(string filename);
    void deserialize(string filename);


};



#endif