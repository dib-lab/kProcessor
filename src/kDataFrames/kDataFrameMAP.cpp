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
    return kmer::int_to_str(iterator->first, this->kSize);
    // return iterator->first;
}

uint64_t kDataFrameMAPIterator::getKmerCount() {
    return iterator->second;
}

bool kDataFrameMAPIterator::setKmerCount(uint64_t count) {
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
    KD = new Kmers(ksize, 2);
//    hasher = new wrapperHasher<std::map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
//    this->MAP = std::map<uint64_t, uint64_t>(1000);
    // this->hasher = (new IntegerHasher(ksize));
}
kDataFrameMAP::kDataFrameMAP(uint64_t ksize,vector<uint64_t> kmersHistogram) {
    this->class_name = "MAP"; // Temporary until resolving #17
    KD = new Kmers(ksize, 2);
    this->kSize = ksize;
    uint64_t countSum=0;
    for(auto h:kmersHistogram)
      countSum+=h;
    reserve(countSum);
//    hasher = new wrapperHasher<std::map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
//    this->MAP = std::map<uint64_t, uint64_t>(1000);
    // this->hasher = (new IntegerHasher(ksize));
}
kDataFrameMAP::kDataFrameMAP() {
    this->class_name = "MAP"; // Temporary until resolving #17
    this->kSize = 23;
    KD = new Kmers(this->kSize, 2);
//    this->MAP = std::map<uint64_t, uint64_t>(1000);
    // hasher=new wrapperHasher<std::map<uint64_t,uint64_t>::hasher >(MAP.hash_function(),kSize);
    // this->hasher = (new IntegerHasher(23));
}

inline bool kDataFrameMAP::kmerExist(string kmerS) {
    return (this->MAP.find(kmer::str_to_canonical_int(kmerS)) == this->MAP.end()) ? 0 : 1;
}


bool kDataFrameMAP::insert(string kmerS, uint64_t count) {
    this->MAP[KD->hash_kmer(kmerS)] += count;
    return true;
}

bool kDataFrameMAP::insert(string kmerS) {
    this->MAP[KD->hash_kmer(kmerS)] = 1;
    return true;
}

bool kDataFrameMAP::insert(uint64_t kmer, uint64_t count) {
    this->MAP[kmer] += count;
    return true;
}

bool kDataFrameMAP::insert(uint64_t kmer) {
    this->MAP[kmer] = 1;
    return true;
}


bool kDataFrameMAP::setCount(string kmerS, uint64_t tag) {
    this->MAP[KD->hash_kmer(kmerS)] = tag;
    return true;
}

bool kDataFrameMAP::setCount(uint64_t kmerS, uint64_t tag) {
    this->MAP[kmerS] = tag;
    return true;
}

uint64_t kDataFrameMAP::count(string kmerS) {
    return this->MAP[KD->hash_kmer(kmerS)];
}

uint64_t kDataFrameMAP::count(uint64_t kmerS) {
    return this->MAP[kmerS];
}

uint64_t kDataFrameMAP::bucket(string kmerS) {
    return 1;
    // return this->MAP.bucket(kmer::str_to_canonical_int(kmerS));
}


bool kDataFrameMAP::erase(string kmerS) {
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
    return 1;
//    return this->MAP.load_factor();
}

float kDataFrameMAP::max_load_factor() {
    return 1;
    //    return this->MAP.max_load_factor();
}


void kDataFrameMAP::save(string filePath) {

    // Write the kmerSize
    ofstream file(filePath + ".extra");
    file << kSize << endl;

    std::ofstream os(filePath + ".map", std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(this->MAP);

}

kDataFrame *kDataFrameMAP::load(string filePath) {

    // Load kSize
    ifstream file(filePath + ".extra");
    uint64_t kSize;
    file >> kSize;

    // Initialize kDataFrameMAP
    kDataFrameMAP *KMAP = new kDataFrameMAP(kSize);

    // Load the hashMap into the kDataFrameMAP
    std::ifstream os(filePath + ".map", std::ios::binary);
    cereal::BinaryInputArchive iarchive(os);
    iarchive(KMAP->MAP);

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
    return *(new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameMAPIterator(MAP.begin(), this, kSize),
            (kDataFrame *) this));
}

kDataFrameIterator kDataFrameMAP::end() {
    return *(new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this));
}
