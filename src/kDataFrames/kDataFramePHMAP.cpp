#include "kDataFrame.hpp"
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/archives/binary.hpp>
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
}

kDataFramePHMAPIterator::kDataFramePHMAPIterator(const kDataFramePHMAPIterator &other) :
        _kDataFrameIterator(other.kSize) {
    iterator = other.iterator;
    this->origin = other.origin;
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
    return kmer::int_to_str(iterator->first, this->kSize);
    // return iterator->first;
}

uint64_t kDataFramePHMAPIterator::getKmerCount() {
    return iterator->second;
}

bool kDataFramePHMAPIterator::setKmerCount(uint64_t count) {
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
    this->class_name = "MAP"; // Temporary until resolving #17
    this->kSize = ksize;
    hasher = new wrapperHasher<flat_hash_map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
    this->MAP = flat_hash_map<uint64_t, uint64_t>(1000);
    // this->hasher = (new IntegerHasher(ksize));
}

kDataFramePHMAP::kDataFramePHMAP() {
    this->class_name = "MAP"; // Temporary until resolving #17
    this->kSize = 23;
    this->MAP = flat_hash_map<uint64_t, uint64_t>(1000);
    // hasher=new wrapperHasher<flat_hash_map<uint64_t,uint64_t>::hasher >(MAP.hash_function(),kSize);
    // this->hasher = (new IntegerHasher(23));
}

inline bool kDataFramePHMAP::kmerExist(string kmerS) {
    return (this->MAP.find(kmer::str_to_canonical_int(kmerS)) == this->MAP.end()) ? 0 : 1;
}


bool kDataFramePHMAP::insert(string kmerS, uint64_t count) {
    this->MAP[kmer::str_to_canonical_int(kmerS)] += count;
    return true;
}

bool kDataFramePHMAP::insert(string kmerS) {
    this->MAP[kmer::str_to_canonical_int(kmerS)] = 1;
    return true;
}


bool kDataFramePHMAP::setCount(string kmerS, uint64_t tag) {
    this->MAP[kmer::str_to_canonical_int(kmerS)] = tag;
    return true;
}

uint64_t kDataFramePHMAP::count(string kmerS) {
    return this->MAP[kmer::str_to_canonical_int(kmerS)];
}

uint64_t kDataFramePHMAP::bucket(string kmerS) {
    return 1;
    // return this->MAP.bucket(kmer::str_to_canonical_int(kmerS));
}


bool kDataFramePHMAP::erase(string kmerS) {
    return this->MAP.erase(kmer::str_to_canonical_int(kmerS));
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


void kDataFramePHMAP::save(string filePath) {

    // Write the kmerSize
    ofstream file(filePath + ".extra");
    file << kSize << endl;

    std::ofstream os(filePath + ".phmap", std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(this->MAP);

}

kDataFrame *kDataFramePHMAP::load(string filePath) {

    // Load kSize
    ifstream file(filePath + ".extra");
    uint64_t kSize;
    file >> kSize;

    // Initialize kDataFramePHMAP
    kDataFramePHMAP *KMAP = new kDataFramePHMAP(kSize);

    // Load the hashMap into the kDataFramePHMAP
    std::ifstream os(filePath + ".phmap", std::ios::binary);
    cereal::BinaryInputArchive iarchive(os);
    iarchive(KMAP->MAP);

    return KMAP;
}

kDataFrame *kDataFramePHMAP::getTwin() {
    return ((kDataFrame *) new kDataFramePHMAP(kSize));
}

void kDataFramePHMAP::reserve(uint64_t n) {
    this->MAP.reserve(n);
}

kDataFrameIterator kDataFramePHMAP::begin() {
    return *(new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.begin(), this, kSize),
            (kDataFrame *) this));
}

kDataFrameIterator kDataFramePHMAP::end() {
    return *(new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this));
}
