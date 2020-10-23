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

kDataFramePHMAPIterator::kDataFramePHMAPIterator(flat_hash_map<uint64_t, uint64_t>::iterator it, kDataFramePHMAP *origin,
                                             uint64_t kSize)
        : _kDataFrameIterator(kSize) {
    iterator = it;
    this->origin = origin;
    this->KD = (new Kmers(kSize));
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
    return KD->ihash_kmer(iterator->first);
    // return iterator->first;
}

uint64_t kDataFramePHMAPIterator::getCount() {
    return iterator->second;
}

bool kDataFramePHMAPIterator::setCount(uint64_t count) {
    iterator->second = count;
}


bool kDataFramePHMAPIterator::operator==(const _kDataFrameIterator &other) {
    return iterator == ((kDataFramePHMAPIterator *) &other)->iterator;
}

bool kDataFramePHMAPIterator::operator!=(const _kDataFrameIterator &other) {
    return iterator != ((kDataFramePHMAPIterator *) &other)->iterator;
}

kDataFramePHMAPIterator::~kDataFramePHMAPIterator() {
delete KD;
}

/*
 **********************
 *** kDataFramePHMAP ***
 **********************
 */

kDataFramePHMAP::kDataFramePHMAP(uint64_t ksize) {
    this->class_name = "PHMAP"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(kSize, 2);
//    hasher = new wrapperHasher<flat_hash_map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
    this->MAP = flat_hash_map<uint64_t, uint64_t>(1000);
    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);
    // this->hasher = (new IntegerHasher(ksize));
}

kDataFramePHMAP::kDataFramePHMAP(uint64_t ksize, int mode) {
    this->class_name = "PHMAP"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(kSize, mode);
//    hasher = new wrapperHasher<flat_hash_map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
    this->MAP = flat_hash_map<uint64_t, uint64_t>(1000);
    // this->hasher = (new IntegerHasher(ksize));
    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);
}

kDataFramePHMAP::kDataFramePHMAP(uint64_t ksize,vector<uint64_t> kmersHistogram) {
    this->class_name = "PHMAP"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(kSize, 2);
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
    this->MAP = flat_hash_map<uint64_t, uint64_t>(1000);
    KD = new Kmers(kSize, 2);
    // hasher=new wrapperHasher<flat_hash_map<uint64_t,uint64_t>::hasher >(MAP.hash_function(),kSize);
    // this->hasher = (new IntegerHasher(23));
    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);
}

inline bool kDataFramePHMAP::kmerExist(string kmerS) {
    return (this->MAP.find(KD->hash_kmer(kmerS)) == this->MAP.end()) ? 0 : 1;
}


bool kDataFramePHMAP::insert(string kmerS, uint64_t count) {
    this->MAP[KD->hash_kmer(kmerS)] += count;
    return true;
}

bool kDataFramePHMAP::insert(string kmerS) {
    this->MAP[KD->hash_kmer(kmerS)] += 1;
    return true;
}


bool kDataFramePHMAP::insert(uint64_t kmer, uint64_t count) {
    this->MAP[kmer] += count;
    return true;
}

bool kDataFramePHMAP::insert(uint64_t kmer) {
    this->MAP[kmer] += 1;
    return true;
}


bool kDataFramePHMAP::setCount(string kmerS, uint64_t tag) {
    this->MAP[KD->hash_kmer(kmerS)] = tag;
    return true;
}

bool kDataFramePHMAP::setCount(uint64_t kmerS, uint64_t tag) {
    this->MAP[kmerS] = tag;
    return true;
}

uint64_t kDataFramePHMAP::getCount(string kmerS) {
    phmap::flat_hash_map<uint64_t,uint64_t>::const_iterator got = this->MAP.find(KD->hash_kmer(kmerS));
    if ( got == this->MAP.end() )
        return 0;
    else
        return got->second;
}

uint64_t kDataFramePHMAP::getCount(uint64_t kmerS) {
    phmap::flat_hash_map<uint64_t,uint64_t>::const_iterator got = this->MAP.find(kmerS);
    if ( got == this->MAP.end() )
        return 0;
    else
        return got->second;
}

uint64_t kDataFramePHMAP::bucket(string kmerS) {
    return 1;
    // return this->MAP.bucket(kmer::str_to_canonical_int(kmerS));
}


bool kDataFramePHMAP::erase(string kmerS) {
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
    
    kDataFramePHMAP *KMAP = new kDataFramePHMAP(kSize, hashing_mode);
    {
        phmap::BinaryInputArchive ar_in(filePath.c_str());
        KMAP->MAP.load(ar_in);
    }
    if(endIterator != nullptr)
        delete endIterator;
    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);

    return KMAP;
}

kDataFrame *kDataFramePHMAP::getTwin() {
    return ((kDataFrame *) new kDataFramePHMAP(kSize));
}

void kDataFramePHMAP::reserve(uint64_t n) {
    this->MAP.reserve(n);
    if(endIterator != nullptr)
        delete endIterator;
    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);

}
void kDataFramePHMAP::reserve(vector<uint64_t> countHistogram) {
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

//kDataFrameIterator kDataFramePHMAP::end() {
//    return *(new kDataFrameIterator(
//            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
//            (kDataFrame *) this));
//}
kDataFrameIterator kDataFramePHMAP::find(string kmer) {
  return (kDataFrameIterator(
          (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.find(kmer::str_to_canonical_int(kmer)), this, kSize),
          (kDataFrame *) this));
}
