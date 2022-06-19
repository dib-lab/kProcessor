#include "kDataFrame.hpp"
#include "Utils/kmer.h"
#include <iostream>
#include <math.h>
#include "defaultColumn.hpp"

using namespace std;

inline bool fileExists(const std::string &name) {
    ifstream f(name.c_str());
    return f.good();
}

kDataFrame::kDataFrame() {
    kSize = 31;
    lastKmerOrder=1;
}

kDataFrame::kDataFrame(uint8_t k_size) {
    kSize = k_size;
    lastKmerOrder=1;
}

kDataFrame::~kDataFrame(){
    delete KD;
    for(auto c:columns)
        delete c.second;
    delete endIterator;
}

bool kDataFrame::empty() {
    return this->size() == 0;
}


void kDataFrame::save(string filePath)
{
    ofstream out(filePath+".multiColumn");
    out<<"isStatic=\t"<<true<<endl;
    out<<"numColumns=\t"<<columns.size()<<endl;
    for(auto c:columns)
    {
	    string suffix=".multiColumn."+c.first;
        string filename=filePath+ suffix;
        size_t columnType=typeid(*(c.second)).hash_code();
        out<<c.first<<"\t"<<columnType<<"\t"<<suffix<<endl;
        c.second->serialize(filename);
    }
    out<<"default\t"<<0<<"\tNULL"<<endl;
    out<<"lastKmerOrder\t"<<lastKmerOrder<<"\tNULL"<<endl;



    out.close();
    this->serialize(filePath);
}

kDataFrame * kDataFrame::load(string filePath) {
    kDataFrame* res;
    if (fileExists(filePath + ".mqf"))
        res=kDataFrameFactory::loadMQF(filePath);
    else if (fileExists(filePath + ".map"))
        res=kDataFrameFactory::loadMAP(filePath);
    else if (fileExists(filePath + ".phmap"))
        res=kDataFrameFactory::loadPHMAP(filePath);
    else if (fileExists(filePath + ".btree"))
        res=kDataFrameFactory::loadBtree(filePath);
    else if (fileExists(filePath+ ".bmqf"))
        res=kDataFrameFactory::loadBMQF(filePath);
    else if (fileExists(filePath+ ".blight.gz"))
        res=kDataFrameFactory::loadBlight(filePath);
    else if (fileExists(filePath+ ".sshash"))
        res=kDataFrameFactory::loadSSHASH(filePath);
    else
        throw std::runtime_error("Could not open kDataFrame file");

    if (fileExists(filePath +".multiColumn") ){
        ifstream inp(filePath + ".multiColumn");
        string tmp;
        string name, path;
        uint64_t type;
        bool tmpB;
        inp >> tmp >> tmpB;
        uint32_t numColumns;
        inp >> tmp >> numColumns;
        for (uint32_t i = 0; i < numColumns; i++) {
            inp >> name >> type >> path;
            Column *c = Column::getContainerByName(type);
            c->deserialize(filePath+path);
            res->columns[name] = c;
	    if(name== "count")
	      res->countColumn=(vectorColumn<unsigned int>*)c;
        if(name== "color") {
            if(type == typeid(StringColorColumn).hash_code() )
            {
                res->colorColumn = ((StringColorColumn *) c)->colors;
            }
            else
                res->colorColumn = (deduplicatedColumn<queryColorColumn> *) c;
        }
        }
        while(inp >> name >> type >> path)
        {
            if(name=="lastKmerOrder")
            {
                res->lastKmerOrder=type;
            }
        }


       
    }

    return res;

}


void kDataFrameIterator::setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName)
{
    std::uint64_t inputKmerOrder=input->getkmerOrder(getHashedKmer());

    origin->columns[outputColName]->setValueFromColumn(input->columns[inputColName],inputKmerOrder,getOrder());
}
void kDataFrame::setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName, std::uint64_t kmer)
{
    std::uint64_t inputKmerOrder=input->getkmerOrder(kmer);
    std::uint64_t outputKmerOrder=getkmerOrder(kmer);
    columns[outputColName]->setValueFromColumn(input->columns[inputColName],inputKmerOrder,outputKmerOrder);
}
void kDataFrame::setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName, string kmer)
{
    std::uint64_t inputKmerOrder=input->getkmerOrder(kmer);
    std::uint64_t outputKmerOrder=getkmerOrder(kmer);
    columns[outputColName]->setValueFromColumn(input->columns[inputColName],inputKmerOrder,outputKmerOrder);
}


//template int kDataFrame::getKmerColumnValue<int, vectorColumn<int> >(string columnName,string kmer);
//template double kDataFrame::getKmerColumnValue<double, vectorColumn<double> >(string columnName,string kmer);
//template bool kDataFrame::getKmerColumnValue<bool, vectorColumn<bool> >(string columnName,string kmer);
//
//template void kDataFrame::setKmerColumnValue<int, vectorColumn<int>  >(string columnName,string kmer, int value);
//template void kDataFrame::setKmerColumnValue<double, vectorColumn<double>  >(string columnName,string kmer, double value);
//template void kDataFrame::setKmerColumnValue<bool, vectorColumn<bool>  >(string columnName,string kmer, bool value);
//
//

void kDataFrame::reserve (std::uint64_t n )
{
    this->_reserve(n);
    for(auto c : columns)
        c.second->resize(n);
}

void kDataFrame::addColumn(string columnName,Column* ptr)
{
  columns[columnName]=ptr;
  if(columnName== "count")
    this->countColumn=(vectorColumn<unsigned int>*)ptr;
}
void kDataFrame::removeColumn(string columnName)
{
  auto it=columns.find(columnName);
  delete it->second;
  columns.erase(it);
}

void kDataFrame::addCountColumn()
{
  if(columns.find("count")==columns.end())  
    addColumn("count",new vectorColumn<uint32_t>(size()+1));
  countColumn=(vectorColumn<uint32_t>*)columns["count"];    
}
bool kDataFrame::setCount(const string &kmer, std::uint64_t N)
{
    uint32_t order=this->getkmerOrder(kmer);
    if(order==0)
    {
        this->insert(kmer);
        order=this->getkmerOrder(kmer);
    }
    countColumn->insert(N,order);

    return true;
}
bool kDataFrame::setCount(std::uint64_t kmer,std::uint64_t N)
{
    uint32_t order=this->getkmerOrder(kmer);
    if(order==0)
    {
        this->insert(kmer);
        order=this->getkmerOrder(kmer);
    }
    countColumn->insert(N,order);
    return true;
}
std::uint64_t kDataFrame::getCount(const string &kmer)
{
    uint32_t o=this->getkmerOrder(kmer);
    if(o==0)
        return 0;
    return countColumn->get(o);
}
std::uint64_t kDataFrame::getCount(std::uint64_t kmer)
{
    uint32_t o=this->getkmerOrder(kmer);
    if(o==0)
        return 0;
    return countColumn->get(o);
}
void kDataFrame::incrementCount(std::uint64_t kmer)
{
    uint32_t order=this->getkmerOrder(kmer);
    uint32_t count=0;
    if(order==0)
    {
        this->insert(kmer);
        order=this->getkmerOrder(kmer);
    }
    else{
        count=countColumn->get(order);
    }
    countColumn->insert(count+1,order);
}
void kDataFrame::incrementCount(const string kmer)
{
    uint32_t order=this->getkmerOrder(kmer);
    uint32_t count=0;
    if(order==0)
    {
        this->insert(kmer);
        order=this->getkmerOrder(kmer);
    }
    else{
        count=countColumn->get(order);
    }
    countColumn->insert(count+1,order);
}

void kDataFrame::addColorColumn(deduplicatedColumn<queryColorColumn>*  col){
    columns["color"]=col;
    colorColumn=col;
}
std::uint32_t kDataFrame::getColorID(const string &kmer){
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    if(kmerOrder==0)
        throw std::logic_error("kmer not found!");
    return colorColumn->index[kmerOrder];
}
std::uint32_t kDataFrame::getColorID(std::uint64_t kmer){
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    if(kmerOrder==0)
        throw std::logic_error("kmer not found!");
    return colorColumn->index[kmerOrder];
}

vector<uint32_t> kDataFrame::getColor(const string &kmer){
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    if(kmerOrder==0)
        throw std::logic_error("kmer not found!");
    return colorColumn->get(kmerOrder);
}
vector<uint32_t> kDataFrame::getColor(std::uint64_t kmer){
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    if(kmerOrder==0)
        throw std::logic_error("kmer not found!");
    return colorColumn->get(kmerOrder);
}



kDataFrameIterator kDataFrame::end(){
//    kDataFrameBMQFIterator* it=new kDataFrameBMQFIterator(bufferedmqf,kSize,KD);
//    it->endIterator();
//    return (kDataFrameIterator(it,(kDataFrame*)this));
    return *endIterator;
}


dbgIterator::dbgIterator(){
    frame = nullptr;
}
dbgIterator::dbgIterator(kDataFrame* f,string kmer){
    frame=f;
    currentKmer = kmer;
    generateNextKmers();
}

dbgIterator::dbgIterator(const dbgIterator& other){
    frame=other.frame;
    currentKmer=other.currentKmer;
    generateNextKmers();
}
dbgIterator& dbgIterator::operator= (const dbgIterator& other){
    frame=other.frame;
    currentKmer=other.currentKmer;
    generateNextKmers();
    return *this;
}


void dbgIterator::generateNextKmers(){
    nextFwdKmers.clear();
    nextRevKmers.clear();
    char possibleNuc[]= {'A','C','G','T'};
    string suffix=currentKmer.substr(1,currentKmer.size()-1);
    for(auto c:possibleNuc)
    {
        string candidate=suffix+c;
        if(candidate==currentKmer)
            continue;
        if(frame->kmerExist(candidate)  )
            nextFwdKmers.push_back(candidate);
    }
    uint32_t k=frame->ksize();
    uint64_t current=kmer::str_to_int(currentKmer);
    uint64_t reverse=kmer::reverse_complement(current,k);
    string reverseKmer=kmer::int_to_str(reverse,k);
    suffix=reverseKmer.substr(1,reverseKmer.size()-1);
    for(auto c:possibleNuc)
    {
        string candidate=suffix+c;
        if(frame->kmerExist(candidate) )
            nextRevKmers.push_back(candidate);
    }
}
void dbgIterator::nextFWD(uint32_t index){
    currentKmer=nextFwdKmers[index];
    generateNextKmers();
}
void dbgIterator::nextREV(uint32_t index){
    currentKmer=nextRevKmers[index];
    generateNextKmers();
}


dbgIterator kDataFrame::getDBGIterator(const string &kmer)
{
    if(!this->kmerExist(kmer))
        throw std::logic_error("Kmer not found in the frame");
    return dbgIterator(this,kmer);
}

vector<string> kDataFrame::getColumnNames()
{
    vector<string> res(columns.size());
    unsigned i=0;
    for(auto c :columns)
        res[i++]=c.first;
    return res;
}


bool kDataFrameIterator::setCount(std::uint64_t count){
    return origin->setCount(iterator->getHashedKmer(),count);
  }
std::uint64_t kDataFrameIterator::getCount(){

    uint32_t o = iterator->getOrder();
    return origin->countColumn->get(o);
}

vector<std::uint32_t> kDataFrameIterator::getColor(){

    uint32_t o = iterator->getOrder();
    return origin->colorColumn->get(o);
}
std::uint32_t kDataFrameIterator::getColorID(){

    uint32_t o = iterator->getOrder();
    return origin->colorColumn->index[o];

}