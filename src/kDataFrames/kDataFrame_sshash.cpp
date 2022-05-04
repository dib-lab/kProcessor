#include "kDataframes/kDataFrameSSHash.hpp"
#include <iostream>
#include "builder/build.cpp"
#include "lookup.cpp"


/*
 *****************************
 *** kDataFrameSSHAshIterator ***
 *****************************
 */


kDataFrame_sshashIterator::kDataFrame_sshashIterator(uint32_t kmerID, kDataFrame_sshash *origin, uint64_t kSize)
        : _kDataFrameIterator(kSize) {

    this->origin = origin;
    this->kmerID= kmerID;

    iterator=new sshash::dictionary::iterator(&origin->dict);
    if(kmerID==origin->dict.size())
    {

        *iterator = origin->dict.at(kmerID-1);
        iterator->next();
    }
    else{

        *iterator = origin->dict.at(kmerID);
        currKmer=iterator->next();
    }
}

kDataFrame_sshashIterator::kDataFrame_sshashIterator(const kDataFrame_sshashIterator &other) :
        _kDataFrameIterator(other.kSize) {


    this->origin = other.origin;
    kmerID= other.kmerID;
    iterator=new sshash::dictionary::iterator(&origin->dict);
    if(kmerID==origin->dict.size())
    {

        *iterator = origin->dict.at(kmerID-1);
        iterator->next();
    }
    else{

        *iterator = origin->dict.at(kmerID);
        currKmer=iterator->next();
    }
}

_kDataFrameIterator *kDataFrame_sshashIterator::clone() {
    return new kDataFrame_sshashIterator(*this);
}

kDataFrame_sshashIterator &kDataFrame_sshashIterator::operator++(int) {
    if(iterator->has_next())
        currKmer=iterator->next();
    kmerID++;
    return *this;
}

uint64_t kDataFrame_sshashIterator::getHashedKmer() {
    //return origin->getHasher()->hash(iterator->first);
    return currKmer.first;

}

string kDataFrame_sshashIterator::getKmer() {
    return currKmer.second;
    // return iterator->first;
}

uint64_t kDataFrame_sshashIterator::getCount() {
    return origin->getkmerOrder(currKmer.second);
}

uint64_t kDataFrame_sshashIterator::getOrder() {
    return origin->getkmerOrder(currKmer.second);
}

bool kDataFrame_sshashIterator::setOrder(uint64_t count) {
    throw logic_error("kmer_sshash is static");
}

void kDataFrame_sshashIterator::endIterator() {
    *iterator=origin->dict.at(origin->dict.size()-1);
    iterator->next();
    kmerID=origin->dict.size();
}


bool kDataFrame_sshashIterator::operator==(const _kDataFrameIterator &other) {

    return kmerID == ((kDataFrame_sshashIterator *) &other)->kmerID;
}

bool kDataFrame_sshashIterator::operator!=(const _kDataFrameIterator &other) {
    return kmerID != ((kDataFrame_sshashIterator *) &other)->kmerID;
}

kDataFrame_sshashIterator::~kDataFrame_sshashIterator() {
    delete iterator;
}
string gen_random_sshash(const int len) {
    srand (time(NULL));
    static const char alphanum[] =
            "0123456789"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            "abcdefghijklmnopqrstuvwxyz";
    string s="_sshash.";
    for (int i = 0; i < len; ++i) {
        s+= alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    return s;
}
/*
 **********************
 *** kDataFrame_sshash ***
 **********************
 */


// kDataFrame_sshash _____________________________


template<typename T,typename Container>
T kDataFrame_sshash::getKmerColumnValue(const string& columnName,string kmer)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    return ((Container*)columns[columnName])->get(kmerOrder);
}
template<typename T,typename Container>
void kDataFrame_sshash::setKmerColumnValue(const string& columnName,string kmer,T value)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    ((Container*)columns[columnName])->insert(value,kmerOrder);
}



kDataFrame_sshash::kDataFrame_sshash()
{
    KD= nullptr;
    endIterator= nullptr;
}

kDataFrame_sshash::kDataFrame_sshash(uint64_t ksize,string input_fasta_file) {
    this->class_name = "sshash"; // Temporary until resolving #17
    this->kSize = ksize;


    sshash::build_configuration build_config;
    build_config.k = ksize;
    build_config.m = 10;
    build_config.canonical_parsing=true;
    dict.build(input_fasta_file, build_config);



    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrame_sshashIterator(dict.size() ,this, kSize),
            (kDataFrame *) this);

    KD= nullptr;
    // this->hasher = (new IntegerHasher(ksize));
}


bool kDataFrame_sshash::kmerExist(string kmerS) {
    return dict.is_member(kmerS.c_str());
}
bool kDataFrame_sshash::kmerExist(uint64_t kmerS) {
    throw logic_error("kDataFrame_sshash doesnt implement kmerExit for uint64_T input. use string version instead");
}



uint32_t kDataFrame_sshash::insert(const string &kmerS) {
    throw logic_error("kDataFrame_sshash is static. Insertion is not allowed");
    return true;
}



uint32_t kDataFrame_sshash::insert(uint64_t kmer) {
    throw logic_error("kDataFrame_sshash is static. Insertion is not allowed");
    return true;
}

bool  kDataFrame_sshash::setOrder(const string &kmer, std::uint64_t count){
    throw logic_error("kDataFrame_sshash is static. Insertion is not allowed");
    return true;
}
bool  kDataFrame_sshash::setOrder(std::uint64_t kmer, std::uint64_t count){
    throw logic_error("kDataFrame_sshash is static. Insertion is not allowed");
    return true;
}


uint64_t kDataFrame_sshash::getkmerOrder(const string &kmerS) {
    return dict.lookup(kmerS.c_str())+1;
}

uint64_t kDataFrame_sshash::getkmerOrder(uint64_t kmerS) {
   return 0;
    // return dict.lookup_uint64(kmerS);
}




bool kDataFrame_sshash::erase(const string &kmerS) {
    throw logic_error("kDataFrame_sshash is static. Deletion is not allowed");
    return true;
}

bool kDataFrame_sshash::erase(uint64_t kmer) {
    throw logic_error("kDataFrame_sshash is static. Deletion is not allowed");
    return true;
}

uint64_t kDataFrame_sshash::size() {
    return (uint64_t) dict.size();
}

uint64_t kDataFrame_sshash::max_size() {
    return (uint64_t) dict.size();
}

float kDataFrame_sshash::load_factor() {
    return 1;
}

float kDataFrame_sshash::max_load_factor() {
    return 1;
}


void kDataFrame_sshash::serialize(string filePath) {

    // Write the kmerSize
    ofstream file(filePath + ".extra");
    file << kSize << endl;
   // file << this->KD->hash_mode << endl;
    file.close();
    filePath += ".sshash";
    {
        essentials::save(dict, filePath.c_str());
    }

}

kDataFrame *kDataFrame_sshash::load(string filePath) {

    // Load kSize
    ifstream file(filePath + ".extra");
    uint64_t kSize;
    file >> kSize;
    file.close();
    filePath += ".sshash";

    kDataFrame_sshash *KMAP = new kDataFrame_sshash();
    KMAP->kSize = kSize;
    {
        essentials::load(KMAP->dict, filePath.c_str());
    }
    if(KMAP->endIterator != nullptr)
        delete KMAP->endIterator;




    KMAP->endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrame_sshashIterator(KMAP->dict.size() ,KMAP, kSize),
            (kDataFrame *) KMAP);

    return KMAP;
}

kDataFrame *kDataFrame_sshash::getTwin() {
    return ((kDataFrame *) new kDataFrame_sshash(kSize));
}

void kDataFrame_sshash::_reserve(uint64_t n) {
    throw logic_error("kDataFrame_sshash is static. _reserve is not allowed");

}
void kDataFrame_sshash::_reserve(vector<uint64_t> countHistogram) {
    throw logic_error("kDataFrame_sshash is static. _reserve is not allowed");


}
kDataFrameIterator kDataFrame_sshash::begin() {

    return (kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrame_sshashIterator(0, this, kSize),
            (kDataFrame *) this));
}

//kDataFrameIterator kDataFrameBlight::end() {
//    return *(new kDataFrameIterator(
//            (_kDataFrameIterator *) new kDataFrameBlightIterator(MAP.end(), this, kSize),
//            (kDataFrame *) this));
//}
kDataFrameIterator kDataFrame_sshash::find(const string &kmer) {
    throw logic_error("not implemented yet");
}
kDataFrameIterator kDataFrame_sshash::find(uint64_t kmer) {
    throw logic_error("not implemented yet");
}

kDataFrame *kDataFrame_sshash::clone() {
    throw logic_error("not implemented yet");
}


kDataFrame *kDataFrameFactory::loadSSHASH(string filePath) {
    return kDataFrame_sshash::load(filePath);
}

kDataFrame *kDataFrameFactory::createSSHASH(uint32_t kSize, string filePath) {
    return new kDataFrame_sshash(kSize,filePath);
}


