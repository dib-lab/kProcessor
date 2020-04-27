#ifndef _defaultColumn_H_
#define _defaultColumn_H_

#include <vector>
#include <stdint.h>
#include <iostream>
#include <parallel_hashmap/phmap.h>
#include <stack>

using phmap::flat_hash_map;
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

class colorNode{
public:
    flat_hash_map<uint32_t, colorNode*> edges;
    uint32_t currColor;
    colorNode() {
        currColor = 0;
    }
    ~colorNode()
    {
//        for(auto e:edges)
//            delete e.second;
    }
};
class colorIndex{
public:
    colorNode* root;
    uint32_t lastColor;
    colorIndex()
    {
        root=new colorNode();
        lastColor=0;
    }
    ~colorIndex();

    bool hasColorID(vector<uint32_t>& v);
    uint32_t getColorID(vector<uint32_t>& v);
    void serialize(string fileName);
    void deserialize(string filename);

    void populateColors(vector<vector<uint32_t > >& colors);

};

class stringColorIndex{
public:
    flat_hash_map<string,uint32_t> colors;
    uint32_t lastColor;
    stringColorIndex()
    {

        lastColor=0;
    }
    ~stringColorIndex(){

    }
    inline string toString(vector<uint32_t>& v)
    {
        string res="";
        for(auto c:v)
        {
            res+=to_string(c);
            res+=";";
        }
        return res;

    }
    bool hasColorID(vector<uint32_t>& v);
    uint32_t getColorID(vector<uint32_t>& v);
    void serialize(string fileName);
    void deserialize(string filename);

    void populateColors(vector<vector<uint32_t > >& colors);

};



class colorColumn: public Column{
public:
    vector<vector<uint32_t > > colors;
    colorIndex colorInv;
    uint64_t  noSamples;
    colorColumn(){
        colors.push_back(vector<uint32_t >());
    }
    colorColumn(uint64_t noSamples){
        this->noSamples=noSamples;
        colors.push_back(vector<uint32_t >());
    }

    ~colorColumn(){

    }
    uint32_t  insertAndGetIndex(vector<uint32_t > item);
    vector<uint32_t > getWithIndex(uint32_t index);

    void insert(vector<uint32_t > item,uint32_t index);
    vector<uint32_t > get(uint32_t index);

    void serialize(string filename);
    void deserialize(string filename);

    void populateColors();



};

class StringColorColumn: public Column{
public:
    vector<vector<uint32_t > > colors;
    flat_hash_map<uint32_t,string> namesMap;

    StringColorColumn(){
        colors.push_back(vector<uint32_t > ());
    }

    ~StringColorColumn(){

    }

    vector<string > getWithIndex(uint32_t index);
    //uint32_t  insertAndGetIndex(vector<uint32_t > item);
//    void insert(vector<uint32_t > item,uint32_t index);
//    vector<uint32_t > get(uint32_t index);

    void serialize(string filename);
    void deserialize(string filename);

    void populateColors();



};


#endif
