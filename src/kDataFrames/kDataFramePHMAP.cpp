#include "kDataFrame.hpp"
#include "parallel_hashmap/phmap_dump.h"
#include <iostream>
#include <fstream>
#include "Utils/kmer.h"

/*
 *****************************
 *** kDataFramePHMAPIterator ***
 *****************************
 */

kDataFramePHMAPIterator::kDataFramePHMAPIterator(itType it, kDataFramePHMAP *origin,
                                             uint64_t kSize)
        : _kDataFrameIterator(kSize) {
    iterator = it;
    this->origin = origin;
    this->KD = this->origin->KD;
}

kDataFramePHMAPIterator::kDataFramePHMAPIterator(const kDataFramePHMAPIterator &other) :
        _kDataFrameIterator(other.kSize) {
    iterator = other.iterator;
    this->origin = other.origin;
    this->KD = other.KD;
}

_kDataFrameIterator *kDataFramePHMAPIterator::clone() {
    return new kDataFramePHMAPIterator(*this);
}

kDataFramePHMAPIterator &kDataFramePHMAPIterator::operator++(int) {
    iterator++;
    return *this;
}

uint64_t kDataFramePHMAPIterator::getHashedKmer() {
    //return origin->getHasher()->hash(iterator->first);
    return iterator->first;

}

string kDataFramePHMAPIterator::getKmer() {
    return origin->getkmerDecoder()->ihash_kmer(iterator->first);
    // return iterator->first;
}

uint64_t kDataFramePHMAPIterator::getCount() {
    return iterator->second;
}
uint64_t kDataFramePHMAPIterator::getOrder() {
    return iterator->second;
}


bool kDataFramePHMAPIterator::setOrder(uint64_t count) {
    iterator->second = count;
}


bool kDataFramePHMAPIterator::operator==(const _kDataFrameIterator &other) {
    return iterator == ((kDataFramePHMAPIterator *) &other)->iterator;
}

bool kDataFramePHMAPIterator::operator!=(const _kDataFrameIterator &other) {
    return iterator != ((kDataFramePHMAPIterator *) &other)->iterator;
}

kDataFramePHMAPIterator::~kDataFramePHMAPIterator() {
}

/*
 **********************
 *** kDataFramePHMAP ***
 **********************
 */

kDataFramePHMAP::kDataFramePHMAP(uint64_t ksize) {
    this->class_name = "PHMAP"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(kSize, TwoBits_hasher);
//    hasher = new wrapperHasher<flat_hash_map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
    this->MAP = MapType(1000000);
    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);
    // this->hasher = (new IntegerHasher(ksize));
}
kDataFramePHMAP::kDataFramePHMAP(uint64_t ksize,uint64_t nKmers) {
    this->class_name = "PHMAP"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(kSize, TwoBits_hasher);
//    hasher = new wrapperHasher<flat_hash_map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
    this->MAP = MapType(nKmers);
    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);
    // this->hasher = (new IntegerHasher(ksize));
}

kDataFramePHMAP::kDataFramePHMAP(uint64_t ksize, hashingModes hash_mode) {
    this->class_name = "PHMAP"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(kSize, hash_mode);
//    hasher = new wrapperHasher<flat_hash_map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
    this->MAP = MapType(1000);
    // this->hasher = (new IntegerHasher(ksize));
    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);
}

kDataFramePHMAP::kDataFramePHMAP(readingModes RM, hashingModes hash_mode, map<string, int> params) {
    this->class_name = "PHMAP"; // Temporary until resolving #17
    KD = kmerDecoder::getInstance(RM, hash_mode, params);
    this->kSize = KD->get_kSize();
    this->MAP = MapType(1000);
}

kDataFramePHMAP::kDataFramePHMAP(uint64_t ksize,vector<uint64_t> kmersHistogram) {
    this->class_name = "PHMAP"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(kSize, TwoBits_hasher);
//    hasher = new wrapperHasher<flat_hash_map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);

    uint64_t countSum=0;
    for(auto h:kmersHistogram)
      countSum+=h;
    reserve(countSum);
//    hasher = new wrapperHasher<std::map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
//    this->MAP = std::map<uint64_t, uint64_t>(1000);
    // this->hasher = (new IntegerHasher(ksize));
    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);
}

kDataFramePHMAP::kDataFramePHMAP() {
    this->class_name = "PHMAP"; // Temporary until resolving #17
    this->kSize = 23;
    this->MAP = MapType(1000);
    KD = new Kmers(kSize, TwoBits_hasher);
    // hasher=new wrapperHasher<flat_hash_map<uint64_t,uint64_t>::hasher >(MAP.hash_function(),kSize);
    // this->hasher = (new IntegerHasher(23));
    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);
}

bool kDataFramePHMAP::kmerExist(string kmerS) {
    return (this->MAP.find(KD->hash_kmer(kmerS)) == this->MAP.end()) ? 0 : 1;
}
bool kDataFramePHMAP::kmerExist(uint64_t kmer) {
    return (this->MAP.find(kmer) == this->MAP.end()) ? 0 : 1;

}




uint32_t kDataFramePHMAP::insert(const string &kmerS) {
    auto it=this->MAP.find(KD->hash_kmer(kmerS));
    if(it!=this->MAP.end())
        return it->second;
    this->MAP[KD->hash_kmer(kmerS)] = lastKmerOrder++;
    return lastKmerOrder-1;
}



uint32_t kDataFramePHMAP::insert(uint64_t kmer) {
    auto it=this->MAP.find(kmer);
    if(it!=this->MAP.end())
        return it->second;
    this->MAP.insert(make_pair(kmer,lastKmerOrder++));
    //this->MAP[kmer] = lastKmerOrder++;
    return lastKmerOrder-1;
}


bool kDataFramePHMAP::setOrder(const string &kmerS, uint64_t tag) {
    this->MAP[KD->hash_kmer(kmerS)] = tag;
    return true;
}

bool kDataFramePHMAP::setOrder(uint64_t kmerS, uint64_t tag) {
    this->MAP[kmerS] = tag;
    return true;
}

uint64_t kDataFramePHMAP::getkmerOrder(const string &kmerS) {
    MapType::const_iterator got = this->MAP.find(KD->hash_kmer(kmerS));
    if ( got == this->MAP.end() )
        return 0;
    else
        return got->second;
}

uint64_t kDataFramePHMAP::getkmerOrder(uint64_t kmerS) {
    MapType::iterator got = this->MAP.find(kmerS);
    if ( got == this->MAP.end() )
        return 0;
    else
        return got->second;
}

uint64_t kDataFramePHMAP::bucket(string kmerS) {
    return 1;
    // return this->MAP.bucket(kmer::str_to_canonical_int(kmerS));
}


bool kDataFramePHMAP::erase(const string &kmerS) {
    return this->MAP.erase(KD->hash_kmer(kmerS));
}

bool kDataFramePHMAP::erase(uint64_t kmer) {
    return this->MAP.erase(kmer);
}

uint64_t kDataFramePHMAP::size() {
    return (uint64_t) this->MAP.size();
}

uint64_t kDataFramePHMAP::max_size() {
    return (uint64_t) this->MAP.max_size();
}

float kDataFramePHMAP::load_factor() {
    return this->MAP.load_factor();
}

float kDataFramePHMAP::max_load_factor() {
    return this->MAP.max_load_factor();
}


void kDataFramePHMAP::serialize(string filePath) {

    // Write the kmerSize
    ofstream file(filePath + ".extra");
    file << kSize << endl;
    file << this->KD->hash_mode << endl;
    file.close();
    filePath += ".phmap";
    {   
        phmap::BinaryOutputArchive ar_out(filePath.c_str());
        this->MAP.dump(ar_out);
    }
    
}

kDataFrame *kDataFramePHMAP::load(string filePath) {

    // Load kSize
    ifstream file(filePath + ".extra");
    uint64_t kSize;
    int hashing_mode;
    file >> kSize;
    file >> hashing_mode;
    file.close();
    filePath += ".phmap";
    hashingModes hash_mode = static_cast<hashingModes>(hashing_mode);
    kDataFramePHMAP *KMAP = new kDataFramePHMAP(kSize, hash_mode);
    {
        phmap::BinaryInputArchive ar_in(filePath.c_str());
        KMAP->MAP.load(ar_in);
    }
    if(KMAP->endIterator != nullptr)
        delete KMAP->endIterator;
    KMAP->endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(KMAP->MAP.end(), KMAP, kSize),
            (kDataFrame *) KMAP);

    return KMAP;
}

kDataFrame *kDataFramePHMAP::getTwin() {
    return ((kDataFrame *) new kDataFramePHMAP(kSize, this->KD->hash_mode));
}

void kDataFramePHMAP::_reserve(uint64_t n) {
    this->MAP.reserve(n);
    if(endIterator != nullptr)
        delete endIterator;
    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);

}
void kDataFramePHMAP::_reserve(vector<uint64_t> countHistogram) {
    uint64_t countSum=0;
    for(auto h:countHistogram)
      countSum+=h;
    reserve(countSum);


}
kDataFrameIterator kDataFramePHMAP::begin() {
    return (kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.begin(), this, kSize),
            (kDataFrame *) this));
}

kDataFrameIterator kDataFramePHMAP::end() {
    ((kDataFramePHMAPIterator*)endIterator->iterator)->iterator=MAP.end();
    return *endIterator;
}
kDataFrameIterator kDataFramePHMAP::find(const string &kmer) {
  return (kDataFrameIterator(
          (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.find(KD->hash_kmer(kmer)), this, kSize),
          (kDataFrame *) this));
}

kDataFrameIterator kDataFramePHMAP::find(uint64_t kmer) {
    return (kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.find((kmer)), this, kSize),
            (kDataFrame *) this));
}

kDataFramePHMAP::MapType *kDataFramePHMAP::getMap()  {
    return &MAP;
}
