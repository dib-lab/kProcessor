#ifndef _defaultColumn_H_
#define _defaultColumn_H_

#include <vector>
#include <stdint.h>
#include <iostream>
#include <parallel_hashmap/phmap.h>
#include "sdsl/vectors.hpp"
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
    map<uint32_t, colorNode*> edges;
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
    uint32_t noSamples;
    colorIndex()
    {
        root=new colorNode();
        lastColor=0;
        noSamples=0;
    }
    ~colorIndex();

    bool hasColorID(vector<uint32_t>& v);
    uint32_t getColorID(vector<uint32_t>& v);
    void serialize(string fileName);
    void deserialize(string filename);

    void populateColors(vector<vector<uint32_t > >& colors);

    void optimize();
    void stats();

};

class stringColorIndex{
public:
    flat_hash_map<string,uint32_t> colors;
    uint32_t lastColor;
    uint32_t noSamples;
    stringColorIndex()
    {

        lastColor=0;
        noSamples=0;
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
  void optimize(){}

};



class inExactColorIndex{
public:
    flat_hash_map<uint64_t,uint32_t> colors;
    uint32_t lastColor;
    uint32_t noSamples;
    inExactColorIndex()
    {

        lastColor=0;
        noSamples=0;
    }
    ~inExactColorIndex(){

    }
    bool hasColorID(vector<uint32_t>& v);
    uint32_t getColorID(vector<uint32_t>& v);
    void serialize(string fileName){}
    void deserialize(string filename){}

    void populateColors(vector<vector<uint32_t > >& colors){}
    void optimize(){}

};

#define NUM_VECTORS 20
#define VECTOR_SIZE 100000000

class insertColorColumn: public Column{
public:

    vector<sdsl::int_vector<> > colors;
    vector<uint32_t> colorsTop;
    vector<uint32_t> vecCount;
    inExactColorIndex colorInv;
    uint64_t  noSamples;
    uint64_t noColors;
    string tmpFolder;
    insertColorColumn(){
        colorsTop=vector<uint32_t>(NUM_VECTORS);
        vecCount=vector<uint32_t>(NUM_VECTORS);
        for(int i=0;i<NUM_VECTORS;i++) {
            colors.push_back(sdsl::int_vector<>(VECTOR_SIZE));
            colorsTop[i]=0;
            vecCount[i]=0;
        }
        tmpFolder="";
        noColors=0;
    }
    insertColorColumn(uint64_t noSamples,string tmp){
        colorsTop=vector<uint32_t>(NUM_VECTORS);
        vecCount=vector<uint32_t>(NUM_VECTORS);
        this->noSamples=noSamples;
        for(int i=0;i<NUM_VECTORS;i++) {
            colors.push_back(sdsl::int_vector<>(VECTOR_SIZE));
            colorsTop[i]=0;
            vecCount[i]=0;
        }
        tmpFolder=tmp;
        noColors=0;

    }

    ~insertColorColumn(){

    }
    uint32_t  insertAndGetIndex(vector<uint32_t >& item);
    vector<uint32_t > getWithIndex(uint32_t index);

    void insert(vector<uint32_t >& item,uint32_t index);
    vector<uint32_t >& get(uint32_t index);

    void serialize(string filename);
    void deserialize(string filename);

    void populateColors();

    uint32_t getNumColors(){
        return colors.size();
    }



};
class vectorBase{
public:
    uint32_t beginID;
    vectorBase(){
        beginID=0;
    }
    vectorBase(uint32_t b){
        beginID=b;
    }
    virtual ~vectorBase(){}
    virtual uint32_t size()=0;
    void save(string filename);
    void load(string filename);

    virtual void serialize(ofstream& of)=0;
    virtual void deserialize(ifstream& f)=0;

    virtual uint64_t sizeInBytes()=0;

    virtual vector<uint32_t> get(uint32_t index)=0;
    virtual void set(uint32_t index,vector<uint32_t>& v)=0;

    friend inline bool operator< (const vectorBase& lhs, const vectorBase& rhs){
        return lhs.beginID < rhs.beginID ;
    }

};

class vectorOfVectors: public vectorBase{
public:
  typedef  sdsl::enc_vector<> vectype;
    vectype  vecs;
    vectype starts;

    vectorOfVectors()
    {

    }
    vectorOfVectors(uint32_t beginId)
            :vectorBase(beginId)
    {
      //  starts.resize(noColors);
    //    vecs.resize(noColors);
    }
    vectorOfVectors(uint32_t beginId,uint32_t noColors)
            :vectorBase(beginId)
    {
            starts=vectype(sdsl::int_vector<>(noColors));
      //        starts=sdsl::enc_vector(sdsl::int_vector(noColors));
        //  starts.resize(noColors);
        //    vecs.resize(noColors);
    }
    ~vectorOfVectors(){

    }
    vector<uint32_t> get(uint32_t index){
        uint32_t start=starts[index];
        uint32_t end=vecs.size();
        if(index< starts.size()-1)
        {
            end=starts[index+1];
        }
        vector<uint32_t> res(end-start);
        for(int i=start;i<end;i++)
            res[i-start]=vecs[i];
        return res;
    };
    void set(uint32_t index,vector<uint32_t>& v)
    {
        throw logic_error("set is not supported");
//        vecs[index]=sdsl::int_vector<>(v.size());
//        for(unsigned int i=0;i<v.size();i++)
//            vecs[index][i]=v[i];
       // sdsl::util::bit_compress(vecs[index]);
    };

    uint32_t size()override {
        return starts.size();
    }
    void loadFromInsertOnly(string path,sdsl::int_vector<>& idsMap);
    void serialize(ofstream& of);
    void deserialize(ifstream& iif);
    uint64_t sizeInBytes(){
        uint64_t res=0;
        res+=sdsl::size_in_bytes(vecs);
        res+=sdsl::size_in_bytes(starts);

        return res;
    }
};


class constantVector: public vectorBase{
public:
    uint32_t  noColors;
    constantVector()
    {
        noColors=0;
    }
    constantVector(uint32_t noColors)
    {
        this->noColors=noColors;
    }

    vector<uint32_t> get(uint32_t index){
        //there is no color id =0 but there is a sample equal zero
        vector<uint32_t> res={index-1};
        return res;
    }
    void set(uint32_t index,vector<uint32_t>& v)
    {
    }

    uint32_t size()override {
        return noColors;
    }
    void serialize(ofstream& of){}
    void deserialize(ifstream& iif){}
    uint64_t sizeInBytes(){

        return 4;
    }
};

class fixedSizeVector: public vectorBase{
public:
    typedef  sdsl::vlc_vector<> vectype;
    vectype vec;
    uint32_t colorsize;
    fixedSizeVector()
    {

    }
    fixedSizeVector(uint32_t beginId,uint32_t colorsize)
            :vectorBase(beginId)
    {
        this->colorsize=colorsize;
        // vec.resize(noColors*size);
    }
    void loadFromInsertOnly(string path,sdsl::int_vector<>& idsMap);
    vector<uint32_t> get(uint32_t index){
        vector<uint32_t> res(colorsize);
        for(unsigned int i=0;i<colorsize;i++)
        {
            res[i]=vec[(index*colorsize)+i];
        }
        return res;
    }
    void set(uint32_t index,vector<uint32_t>& v)
    {
        throw logic_error("set is not supported");
//        if(v.size()!=colorsize)
//            cout<<"error "<<index<<" "<<colorsize<<" "<<v.size()<<endl;
//        for(unsigned int i=0;i<colorsize;i++)
//        {
//            vec[(index*colorsize)+i]=v[i];
//        }
    }

    uint32_t size()override {
        return vec.size()/colorsize;
    }
    void serialize(ofstream& of);
    void deserialize(ifstream& iif);
    uint64_t sizeInBytes(){

        return sdsl::size_in_bytes(vec)+4;
    }
};
class RLEfixedSizeVector: public vectorBase{
public:
    typedef  sdsl::enc_vector<> vectype;
    vectype vec;
    vectype starts;
    uint32_t colorsize;
    uint32_t numColors;
    RLEfixedSizeVector()
    {

    }
    RLEfixedSizeVector(fixedSizeVector* fv,sdsl::int_vector<>& idsMap);
    void loadFromInsertOnly(string path,sdsl::int_vector<>& idsMap);
    vector<uint32_t> get(uint32_t index);
    void set(uint32_t index,vector<uint32_t>& v)
    {
        throw logic_error("set is not supported");
    }

    uint32_t size()override {
        return numColors;
    }
    void serialize(ofstream& of);
    void deserialize(ifstream& iif);
    uint64_t sizeInBytes(){

        return sdsl::size_in_bytes(vec)+sdsl::size_in_bytes(starts)+8;
    }
};
class queryColorColumn: public Column{
public:
    deque<vectorBase*> colors;
    uint64_t  noSamples;
    sdsl::int_vector<> idsMap;
    uint64_t numColors;
    queryColorColumn(){
        colors.push_back(new vectorOfVectors());
    }
    queryColorColumn(uint64_t noSamples){
        this->noSamples=noSamples;
        colors.push_back(new vectorOfVectors());
    }
    queryColorColumn(uint64_t noSamples,uint64_t noColors,string tmpFolder);
    queryColorColumn(insertColorColumn* col);
    ~queryColorColumn(){
        for(auto v:colors)
            delete v;
    }
    uint32_t  insertAndGetIndex(vector<uint32_t >& item);
    vector<uint32_t > getWithIndex(uint32_t index);

    void insert(vector<uint32_t >& item,uint32_t index);
    vector<uint32_t > get(uint32_t index);

    void serialize(string filename);
    void deserialize(string filename);
    void optimizeRLE();
    void optimize(insertColorColumn* col);
    void optimize2();
    void optimize3(insertColorColumn* col);

    uint32_t getNumColors(){
        uint32_t res=0;
        for(int i=1;i<colors.size();i++)
            res+=colors[i]->size();
        return res;
    }
    uint64_t sizeInBytes();



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
