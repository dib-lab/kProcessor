#ifndef _defaultColumn_H_
#define _defaultColumn_H_

#include <vector>
#include <stdint.h>
#include <iostream>
#include <parallel_hashmap/phmap.h>

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

struct colorNode{
    flat_hash_map<uint32_t, colorNode*> edges;
    uint32_t currColor;
    colorNode() {
        currColor = 0;
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
    bool hasColorID(vector<uint32_t>& v)
    {
        colorNode* currNode=root;
        int i=0;
        while(i<v.size())
        {
            auto it=currNode->edges.find(v[i]);
            if(it==currNode->edges.end())
                break;
            currNode=it->second;
            i++;
        }
        if(i!=v.size())
            return false;
        return currNode->currColor!=0;
    }
    uint32_t getColorID(vector<uint32_t>& v)
    {
        colorNode* currNode=root;
        int i=0;
        while(i<v.size())
        {
            auto it=currNode->edges.find(v[i]);
            if(it==currNode->edges.end())
                break;
            currNode=it->second;
            i++;
        }
        for(;i<v.size();i++)
        {
            currNode->edges[v[i]]=new colorNode();
            currNode=currNode->edges[v[i]];
        }
        if(currNode->currColor==0)
            currNode->currColor=++lastColor;
        return currNode->currColor;
    }
};



class colorColumn: public Column{
public:
    vector<vector<uint32_t > > colors;
    colorIndex colorInv;
    colorColumn(){

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

#endif