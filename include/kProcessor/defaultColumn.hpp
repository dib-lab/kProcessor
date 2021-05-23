#ifndef _defaultColumn_H_
#define _defaultColumn_H_

#include <vector>
#include <stdint.h>
#include <iostream>
#include <parallel_hashmap/phmap.h>
#include "sdsl/vectors.hpp"
#include "sdsl/bp_support.hpp"
#include <stack>
#include <unordered_set>

using phmap::flat_hash_map;
using namespace std;

//to add new column create a column class and add the name in getContainer Bynamefn
class Column{
public:
    Column(){}
    virtual ~Column(){}

    virtual Column* getTwin()=0;
    virtual void setSize(uint32_t size)=0;
    virtual void resize(uint32_t size)=0;
    static Column* getContainerByName(size_t name);

    virtual void serialize(string filename)=0;
    virtual void deserialize(string filename)=0;

    virtual void setValueFromColumn(Column* Container, uint32_t inputOrder,uint32_t outputOrder){

    }


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

    void resize(uint32_t size)
    {
        dataV.resize(size);
    }
    void insert(T item,uint32_t index);
    T get(uint32_t index);

    void serialize(string filename);
    void deserialize(string filename);
    Column* getTwin();
    void setSize(uint32_t size);
    void setValueFromColumn(Column* Container, uint32_t inputOrder,uint32_t outputOrder);



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

    Column* getTwin();
    void setSize(uint32_t size);
    void resize(uint32_t size);



};


class _vectorBaseIterator{
    public:
    _vectorBaseIterator(){}

    virtual _vectorBaseIterator& operator ++ (int)=0;
    virtual _vectorBaseIterator* clone()=0;
    virtual bool operator ==(const _vectorBaseIterator& other)=0;
    virtual bool operator !=(const _vectorBaseIterator& other)=0;
    virtual vector<uint32_t> operator*()=0;
    virtual uint32_t getID()=0;
    virtual ~_vectorBaseIterator(){};

};

class vectorBaseIterator{
public:
    _vectorBaseIterator* iterator;
    vectorBaseIterator(){
        iterator=NULL;
    }

    vectorBaseIterator(const vectorBaseIterator& other){
        if(other.iterator!=NULL){
            iterator=other.iterator->clone();
        }
        else{
            iterator=NULL;
        }
    }
    vectorBaseIterator& operator= (const vectorBaseIterator& other){
        if(other.iterator!=NULL){
            iterator=other.iterator->clone();
        }
        else{
            iterator=NULL;
        }
        return *this;
    }
    vectorBaseIterator(_vectorBaseIterator* it){
        iterator=it;
    }
/// Increment the iterator to the next kmer
    vectorBaseIterator& operator ++ (int){
        (*iterator)++;
        return *this;
    }

    /// Increment the iterator to the next kmer (Implemented mainly for python interface)
    vectorBaseIterator& next(){
        (*iterator)++;
        return *this;
    }
    /// Compare the position of each iterator in the underlying datastructure.
    /*! returns True when current and other points to the same kmer */
    bool operator ==(const vectorBaseIterator& other)
    {
        return *iterator == *other.iterator;
    }
    /// Compare the position of each iterator in the underlying datastructure.
    /*! returns True when current and other points to different kmers */
    bool operator !=(const vectorBaseIterator& other)
    {
        return *iterator != *other.iterator;
    }
    vector<uint32_t> operator*(){
        return *(*iterator);
    }

    uint32_t getID(){
        return iterator->getID();
    }
    ~vectorBaseIterator(){
        delete iterator;
    }
};



class vectorBase{
    protected:
        vectorBaseIterator* endIterator;
    public:
    uint32_t beginID;
    vectorBase(){
        beginID=0;
    }
    vectorBase(uint32_t b){
        beginID=b;
    }
    virtual ~vectorBase(){delete endIterator;}
    virtual uint32_t size()=0;
    void save(string filename);
    void load(string filename);

    virtual void sort(sdsl::int_vector<>& idsMap)=0;

    virtual void serialize(ofstream& of)=0;
    virtual void deserialize(ifstream& f)=0;

    virtual uint64_t sizeInBytes()=0;
    virtual double sizeInMB()=0;
    virtual void explainSize()=0;
    virtual uint64_t numIntegers()=0;

    virtual vector<uint32_t> get(uint32_t index)=0;
    virtual void set(uint32_t index,vector<uint32_t>& v)=0;

    friend inline bool operator< (const vectorBase& lhs, const vectorBase& rhs){
        return lhs.beginID < rhs.beginID ;
    }

    ///Returns an iterator at the beginning of the vectorBaseIterator.
    virtual vectorBaseIterator begin()=0;
    ///Returns an iterator at the end of the vectorBaseIterator.
    virtual vectorBaseIterator end(){
        return *endIterator;
    }

};

class vectorOfVectors: public vectorBase{
public:
  typedef  sdsl::enc_vector<> vectype;
    vectype  vecs;
    vectype starts;

    vectorOfVectors();
    vectorOfVectors(uint32_t beginId);
    vectorOfVectors(uint32_t beginId,uint32_t noColors);
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
        for(unsigned int i=start;i<end;i++)
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
    void explainSize() override {
        cout<<"Vector of Vectors "<<size()<<" colors in"<< sizeInBytes()/(1024.0*1024.0)<<"MB"<<endl;
    }
    uint64_t numIntegers()override {
        return vecs.size()+starts.size();
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
    double sizeInMB(){
        return sdsl::size_in_mega_bytes(vecs)+sdsl::size_in_mega_bytes(starts);
    }
    void sort(sdsl::int_vector<>& idsMap);
    vectorBaseIterator begin() override ;
    vectorBaseIterator end();



};

class vectorOfVectorsIterator: public _vectorBaseIterator{
public:
    vectorOfVectors::vectype::iterator vecsIt;
    vectorOfVectors::vectype::iterator startsIt;
    vectorOfVectors* origin;
    vectorOfVectorsIterator(vectorOfVectors* o){
        origin=o;
        vecsIt=o->vecs.begin();
        startsIt=o->starts.begin();
    }

    _vectorBaseIterator& operator ++ (int){
        uint32_t end=origin->vecs.size();
        if((startsIt+1)!=origin->starts.end())
            end=*(startsIt+1);
        int diff=end-*startsIt;
        startsIt++;
        vecsIt+=diff;
        return *this;
    }
    _vectorBaseIterator* clone(){
        vectorOfVectorsIterator* newIt=new vectorOfVectorsIterator(origin);
        newIt->vecsIt=vecsIt;
        newIt->startsIt=startsIt;
        newIt->origin=origin;
        return newIt;
    }
    bool operator ==(const _vectorBaseIterator& other){
        return startsIt == ((vectorOfVectorsIterator*)&other)->startsIt;
    };
    bool operator !=(const _vectorBaseIterator& other){
        return startsIt != ((vectorOfVectorsIterator*)&other)->startsIt;
    };
    vector<uint32_t> operator*(){
        auto tmpIt=vecsIt;
        uint32_t end=origin->vecs.size();
        if((startsIt+1)!=origin->starts.end())
            end=*(startsIt+1);
        int diff=end-*startsIt;
        vector<uint32_t> res(diff);
        for(int i=0;i<diff;i++) {
            res[i] = *tmpIt;
            tmpIt++;
        }
        return res;
    };
    uint32_t getID()
    {
        return (startsIt-origin->starts.begin())+origin->beginID;
    }
    ~vectorOfVectorsIterator(){};

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
    ~constantVector(){}

    vector<uint32_t> get(uint32_t index){
        //there is no color id =0 but there is a sample equal zero
        vector<uint32_t> res={index-1};
        return res;
    }
    void sort(sdsl::int_vector<>& idsMap)
    {

    }
    void set(uint32_t index,vector<uint32_t>& v)
    {
    }
    void explainSize(){
        cout<<"constant "<<size()<<" colors in"<< sizeInBytes()/(1024.0*1024.0)<<"MB"<<endl;
    }
    uint64_t numIntegers(){
        return 0;
    }
    uint32_t size()override {
        return noColors;
    }
    void serialize(ofstream& of){}
    void deserialize(ifstream& iif){}
    uint64_t sizeInBytes(){

        return 4;
    }
    vectorBaseIterator begin(){
        return vectorBaseIterator();
    }
    double sizeInMB(){
        return 0.0;
    }


};



class fixedSizeVectorIterator;


class fixedSizeVector: public vectorBase{
public:
    typedef  sdsl::vlc_vector<> vectype;
    vectype vec;
    uint32_t colorsize;
    fixedSizeVector();
    ~fixedSizeVector(){}

    fixedSizeVector(uint32_t beginId,uint32_t colorsize);
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
    void sort(sdsl::int_vector<>& idsMap);
    void explainSize(){
        cout<<"fixed size("<<colorsize<<") of "<<size()<<" colors in "<< sizeInBytes()/(1024.0*1024.0)<<"MB"<<endl;
    }
    uint64_t numIntegers(){
        return vec.size();
    }
    uint64_t sizeInBytes(){

        return sdsl::size_in_bytes(vec)+4;
    }
    double sizeInMB(){
        return sdsl::size_in_mega_bytes(vec);
    }
    vectorBaseIterator begin();
    vectorBaseIterator end();

};

class fixedSizeVectorIterator: public _vectorBaseIterator{
public:
    fixedSizeVector::vectype::iterator it;
    fixedSizeVector* origin;
    fixedSizeVectorIterator(fixedSizeVector* o){
        origin=o;
        it=o->vec.begin();
    }

    _vectorBaseIterator& operator ++ (int){
        it+=origin->colorsize;
        return *this;
    }
    _vectorBaseIterator* clone(){
        fixedSizeVectorIterator* newIt=new fixedSizeVectorIterator(origin);
        newIt->it=it;
        newIt->origin=origin;
        return newIt;
    }
    bool operator ==(const _vectorBaseIterator& other){
        return it == ((fixedSizeVectorIterator*)&other)->it;
    };
    bool operator !=(const _vectorBaseIterator& other){
        return it != ((fixedSizeVectorIterator*)&other)->it;
    };
    vector<uint32_t> operator*(){
        vector<uint32_t> res(origin->colorsize);
        auto tmpIt=it;
        for(unsigned int i=0;i<origin->colorsize;i++) {
            res[i] = *tmpIt;
            tmpIt++;
        }
        return res;
    }
    uint32_t getID()
    {
        return (it-origin->vec.begin())/origin->colorsize+origin->beginID;
    }
    ~fixedSizeVectorIterator(){};

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
    void sort(sdsl::int_vector<>& idsMap)
    {

    }
    void explainSize(){
        cout<<"RLE fixed size("<<colorsize<<") of "<<size()<<" colors in "<< sizeInBytes()/(1024.0*1024.0)<<"MB"<<endl;
    }
    uint64_t numIntegers(){
        return numColors*colorsize;
    }
    vectorBaseIterator begin(){
        return vectorBaseIterator();
    }
    vectorBaseIterator end(){
        return vectorBaseIterator();
    }
    double sizeInMB(){
        return sdsl::size_in_mega_bytes(vec)+sdsl::size_in_mega_bytes(starts);
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
    void sortColors();
    void optimizeRLE();
    void optimize(insertColorColumn* col);
    void optimize2();
    void optimize3(insertColorColumn* col);

    uint32_t size(){
        uint32_t res=0;
        for(unsigned int i=1;i<colors.size();i++)
            res+=colors[i]->size();
        return res;
    }
    uint32_t numIntegers(){
        uint32_t res=0;
        for(unsigned int i=1;i<colors.size();i++)
            res+=colors[i]->numIntegers();
        return res;
    }
    uint64_t sizeInBytes();
    void explainSize();

    Column* getTwin();
    void setSize(uint32_t size);
    void resize(uint32_t size);


};


class prefixTrieQueryColorColumn: public Column{
private:
    unordered_map<uint64_t ,vector<uint32_t > > nodesCache;
    deque<sdsl::int_vector<>*>  unCompressedEdges;
public:
    typedef  sdsl::vlc_vector<> vectype;
    deque<vectype*>  edges;
    deque<sdsl::bit_vector*> tree;
    deque<sdsl::bp_support_sada<>*> bp_tree;
    sdsl::int_vector<64> starts;
    uint64_t  noSamples;
    sdsl::int_vector<64> idsMap;
    uint64_t numColors;
    prefixTrieQueryColorColumn(){

    }
    prefixTrieQueryColorColumn(queryColorColumn* col);
    ~prefixTrieQueryColorColumn(){
        for(auto t:tree)
            delete t;
        for(auto b:bp_tree)
            delete b;
        for(auto e:edges)
            delete e;

    }
    uint32_t  insertAndGetIndex(vector<uint32_t >& item);
    vector<uint32_t > getWithIndex(uint32_t index);

    void insert(vector<uint32_t >& item,uint32_t index);
    //vector<uint32_t > get(uint32_t index);

    void serialize(string filename);
    void deserialize(string filename);
    void shorten(deque<uint32_t> & input,deque<uint32_t> & output);
   // void _shorten(vector<uint32_t> & input,vector<uint32_t> & output, vector<uint32_t> & remaining);
    uint32_t getNumColors();
    uint64_t sizeInBytes();
    void explainSize();
    void exportTree(string prefix,int tree);
    Column* getTwin();
    void setSize(uint32_t size);
    void resize(uint32_t size);


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
    vector<uint32_t > get(uint32_t index);

    void serialize(string filename);
    void deserialize(string filename);

    void populateColors();


    Column* getTwin();
    void setSize(uint32_t size);
    void resize(uint32_t size);


};

template<typename  T, typename ColumnType>
class deduplicatedColumn: public Column{
public:
    vector<uint32_t> index;
    ColumnType* values;
    deduplicatedColumn(){

    }
    deduplicatedColumn(uint32_t size){
        index=vector<uint32_t>(size);
    }
    ~deduplicatedColumn(){

    }


    T get(uint32_t index);

    void serialize(string filename);
    void deserialize(string filename);
    Column* getTwin();
    void setSize(uint32_t size);
    void resize(uint32_t size);
    void setValueFromColumn(Column* Container, uint32_t inputOrder,uint32_t outputOrder);



};



#endif
