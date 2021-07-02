#include "kDataFrame.hpp"
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/archives/binary.hpp>
#include <iostream>
#include <fstream>
#include "Utils/kmer.h"

/*
 *****************************
 *** kDataFrameMAPIterator ***
 *****************************
 */

kDataFrameMAPIterator::kDataFrameMAPIterator(std::map<uint64_t, uint64_t>::iterator it, kDataFrameMAP *origin,
                                             uint64_t kSize)
        : _kDataFrameIterator(kSize) {
    iterator = it;
    this->origin = origin;
}

kDataFrameMAPIterator::kDataFrameMAPIterator(const kDataFrameMAPIterator &other) :
        _kDataFrameIterator(other.kSize) {
    iterator = other.iterator;
    this->origin = other.origin;
}

_kDataFrameIterator *kDataFrameMAPIterator::clone() {
    return new kDataFrameMAPIterator(*this);
}

kDataFrameMAPIterator &kDataFrameMAPIterator::operator++(int) {
    iterator++;
    return *this;
}

uint64_t kDataFrameMAPIterator::getHashedKmer() {
    //return origin->getHasher()->hash(iterator->first);
    return iterator->first;

}

string kDataFrameMAPIterator::getKmer() {
    return origin->getkmerDecoder()->ihash_kmer(iterator->first);
    // return iterator->first;
}

uint64_t kDataFrameMAPIterator::getCount() {
    return iterator->second;
}
uint64_t kDataFrameMAPIterator::getOrder() {
    return iterator->second;
}

bool kDataFrameMAPIterator::setOrder(uint64_t count) {
    iterator->second = count;
}


bool kDataFrameMAPIterator::operator==(const _kDataFrameIterator &other) {
    return iterator == ((kDataFrameMAPIterator *) &other)->iterator;
}

bool kDataFrameMAPIterator::operator!=(const _kDataFrameIterator &other) {
    return iterator != ((kDataFrameMAPIterator *) &other)->iterator;
}

kDataFrameMAPIterator::~kDataFrameMAPIterator() {

}

/*
 **********************
 *** kDataFrameMAP ***
 **********************
 */

kDataFrameMAP::kDataFrameMAP(uint64_t ksize) {
    this->class_name = "MAP"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(ksize, TwoBits_hasher);
//    hasher = new wrapperHasher<std::map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
//    this->MAP = std::map<uint64_t, uint64_t>(1000);
    // this->hasher = (new IntegerHasher(ksize));
    endIterator=new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);
}

kDataFrameMAP::kDataFrameMAP(uint64_t ksize,uint64_t nKmers) {
    this->class_name = "MAP"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(ksize, TwoBits_hasher);
//    hasher = new wrapperHasher<std::map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
//    this->MAP = std::map<uint64_t, uint64_t>(1000);
    // this->hasher = (new IntegerHasher(ksize));
    endIterator=new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);
}
kDataFrameMAP::kDataFrameMAP(uint64_t ksize,vector<uint64_t> kmersHistogram) {
    this->class_name = "MAP"; // Temporary until resolving #17
    KD = new Kmers(ksize, TwoBits_hasher);
    this->kSize = ksize;
    uint64_t countSum=0;
    for(auto h:kmersHistogram)
      countSum+=h;
    reserve(countSum);
    endIterator=new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);
//    hasher = new wrapperHasher<std::map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
//    this->MAP = std::map<uint64_t, uint64_t>(1000);
    // this->hasher = (new IntegerHasher(ksize));
}
kDataFrameMAP::kDataFrameMAP() {
    this->class_name = "MAP"; // Temporary until resolving #17
    this->kSize = 23;
    KD = new Kmers(this->kSize, TwoBits_hasher);
    endIterator=new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);
//    this->MAP = std::map<uint64_t, uint64_t>(1000);
    // hasher=new wrapperHasher<std::map<uint64_t,uint64_t>::hasher >(MAP.hash_function(),kSize);
    // this->hasher = (new IntegerHasher(23));
}

bool kDataFrameMAP::kmerExist(string kmerS) {
    return (this->MAP.find(kmer::str_to_canonical_int(kmerS)) == this->MAP.end()) ? 0 : 1;
}
bool kDataFrameMAP::kmerExist(uint64_t kmer) {
    return (this->MAP.find(kmer) == this->MAP.end()) ? 0 : 1;
}



bool kDataFrameMAP::insert(const string &kmerS) {
    auto it=this->MAP.find(KD->hash_kmer(kmerS));
    if(it==this->MAP.end())
        return false;
    this->MAP[KD->hash_kmer(kmerS)] = lastKmerOrder++;
    return true;
}



bool kDataFrameMAP::insert(uint64_t kmer) {
    auto it=this->MAP.find(kmer);
    if(it==this->MAP.end())
        return false;
    this->MAP[kmer] = lastKmerOrder++;
    return true;
}


bool kDataFrameMAP::setCount(const string &kmerS, uint64_t tag) {
    this->MAP[KD->hash_kmer(kmerS)] = tag;
    return true;
}

bool kDataFrameMAP::setOrder(uint64_t kmerS, uint64_t tag) {
    this->MAP[kmerS] = tag;
    return true;
}

uint64_t kDataFrameMAP::getkmerOrder(const string &kmerS) {
    auto it= this->MAP.find(KD->hash_kmer(kmerS));
    if(it==this->MAP.end())
        return 0;

    return it->second;
}

uint64_t kDataFrameMAP::getkmerOrder(uint64_t kmerS) {
    auto it= this->MAP.find(kmerS);
    if(it==this->MAP.end())
        return 0;

    return it->second;
}

uint64_t kDataFrameMAP::bucket(string kmerS) {
    return 1;
    // return this->MAP.bucket(kmer::str_to_canonical_int(kmerS));
}


bool kDataFrameMAP::erase(const string &kmerS) {
    return this->MAP.erase(KD->hash_kmer(kmerS));
}

bool kDataFrameMAP::erase(uint64_t kmer) {
    return this->MAP.erase(kmer);
}

uint64_t kDataFrameMAP::size() {
    return (uint64_t) this->MAP.size();
}

uint64_t kDataFrameMAP::max_size() {
    return (uint64_t) this->MAP.max_size();
}

float kDataFrameMAP::load_factor() {
    return 0.5;
 //   return (float)this->MAP.size()/(float)this->MAP.capacity();
}

float kDataFrameMAP::max_load_factor() {
    return 1;
  //  return this->MAP.max_load_factor();
}


void kDataFrameMAP::serialize(string filePath) {

    // Write the kmerSize
    ofstream file(filePath + ".extra");
    file << kSize << endl;
    file << 2 << endl;
    file.close();
    std::ofstream os(filePath + ".map", std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(this->MAP);

}

kDataFrame *kDataFrameMAP::load(string filePath) {

    // Load kSize
    ifstream file(filePath + ".extra");
    uint64_t kSize;
    int hashing_mode;
    file >> kSize;
    file >> hashing_mode;

    if(hashing_mode != 2){
        std::cerr << "Error: In the kDataFrameMAP, hashing must be 2:TwoBitsRepresentation mode" << endl;
        exit(1);
    }
    file.close();
    // Initialize kDataFrameMAP
    kDataFrameMAP *KMAP = new kDataFrameMAP(kSize);

    // Load the hashMap into the kDataFrameMAP
    std::ifstream os(filePath + ".map", std::ios::binary);
    cereal::BinaryInputArchive iarchive(os);
    iarchive(KMAP->MAP);
    delete KMAP->endIterator;
    KMAP->endIterator=new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameMAPIterator(KMAP->MAP.end(), KMAP, kSize),
            (kDataFrame *) KMAP);

    return KMAP;
}

kDataFrame *kDataFrameMAP::getTwin() {
    return ((kDataFrame *) new kDataFrameMAP(kSize));
}

void kDataFrameMAP::reserve(uint64_t n) {
//    this->MAP.reserve(n);
}
void kDataFrameMAP::reserve(vector<uint64_t> countHistogram) {
//    this->MAP.reserve(n);
}

kDataFrameIterator kDataFrameMAP::begin() {
    return (kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameMAPIterator(MAP.begin(), this, kSize),
            (kDataFrame *) this));
}

kDataFrameIterator kDataFrameMAP::end() {
    ((kDataFrameMAPIterator*)endIterator->iterator)->iterator=MAP.end();
    return *endIterator;
}
kDataFrameIterator kDataFrameMAP::find(const string &kmer) {
  return ( kDataFrameIterator(
          (_kDataFrameIterator *) new kDataFrameMAPIterator(MAP.find(KD->hash_kmer(kmer)), this, kSize),
          (kDataFrame *) this));
}

kDataFrameIterator kDataFrameMAP::find(uint64_t kmer) {
    return ( kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameMAPIterator(MAP.find(kmer), this, kSize),
            (kDataFrame *) this));
}
