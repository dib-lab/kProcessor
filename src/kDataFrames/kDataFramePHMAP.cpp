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

kDataFramePHMAPIterator::kDataFramePHMAPIterator(KDATAFRAME_PHMAP::iterator it, kDataFramePHMAP *origin,
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
    this->MAP = KDATAFRAME_PHMAP(1000);
    // this->hasher = (new IntegerHasher(ksize));
}

kDataFramePHMAP::kDataFramePHMAP(uint64_t ksize, hashingModes hash_mode) {
    this->class_name = "PHMAP"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(kSize, hash_mode);
//    hasher = new wrapperHasher<flat_hash_map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
    this->MAP = KDATAFRAME_PHMAP(1000);
    // this->hasher = (new IntegerHasher(ksize));
}

kDataFramePHMAP::kDataFramePHMAP(readingModes RM, hashingModes hash_mode, map<string, int> params) {
    this->class_name = "PHMAP"; // Temporary until resolving #17
    KD = kmerDecoder::getInstance(RM, hash_mode, params);
    this->kSize = KD->get_kSize();
    this->MAP = KDATAFRAME_PHMAP(1000);
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
}

kDataFramePHMAP::kDataFramePHMAP() {
    this->class_name = "PHMAP"; // Temporary until resolving #17
    this->kSize = 23;
    this->MAP = KDATAFRAME_PHMAP(1000);
    KD = new Kmers(kSize, TwoBits_hasher);
    // hasher=new wrapperHasher<flat_hash_map<uint64_t,uint64_t>::hasher >(MAP.hash_function(),kSize);
    // this->hasher = (new IntegerHasher(23));
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
    KDATAFRAME_PHMAP::const_iterator got = this->MAP.find(KD->hash_kmer(kmerS));
    if ( got == this->MAP.end() )
        return 0;
    else
        return got->second;
}

uint64_t kDataFramePHMAP::getCount(uint64_t kmerS) {
    KDATAFRAME_PHMAP::const_iterator got = this->MAP.find(kmerS);
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


void kDataFramePHMAP::save(string filePath) {

    // Write the kmerSize
    ofstream file(filePath + ".extra");
    file << kSize << endl;
    file << this->KD->hash_mode << endl;
    file << this->KD->slicing_mode << endl;
    file << this->KD->params_to_string() << endl;

    filePath += ".phmap";
    {   
        phmap::BinaryOutputArchive ar_out(filePath.c_str());
        this->MAP.phmap_dump(ar_out);
    }
    
}

kDataFrame *kDataFramePHMAP::load(string filePath) {

    // Load kSize
    ifstream file(filePath + ".extra");
    uint64_t filekSize;

    int hashing_mode, reading_mode;

    string KD_params_string;

    file >> filekSize;
    file >> hashing_mode;
    file >> reading_mode;
    file >> KD_params_string;

    hashingModes hash_mode = static_cast<hashingModes>(hashing_mode);
    readingModes slicing_mode = static_cast<readingModes>(reading_mode);
    map<string, int> kmerDecoder_params = kmerDecoder::string_to_params(KD_params_string);

    filePath += ".phmap";
    
    kDataFramePHMAP *KMAP = new kDataFramePHMAP(slicing_mode, hash_mode, kmerDecoder_params);
    {
        phmap::BinaryInputArchive ar_in(filePath.c_str());
        KMAP->MAP.phmap_load(ar_in);
    }

    return KMAP;
}

kDataFrame *kDataFramePHMAP::getTwin() {
    return ((kDataFrame *) new kDataFramePHMAP(this->KD->slicing_mode, this->KD->hash_mode, this->KD->string_to_params(this->KD->params_to_string())));
}

void kDataFramePHMAP::reserve(uint64_t n) {
    this->MAP.reserve(n);
}

void kDataFramePHMAP::reserve(vector<uint64_t> countHistogram) {
    uint64_t countSum=0;
    for(auto h:countHistogram)
      countSum+=h;
    reserve(countSum);

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
kDataFrameIterator kDataFramePHMAP::find(string kmer) {
  return *(new kDataFrameIterator(
          (_kDataFrameIterator *) new kDataFramePHMAPIterator(MAP.find(kmer::str_to_canonical_int(kmer)), this, kSize),
          (kDataFrame *) this));
}
