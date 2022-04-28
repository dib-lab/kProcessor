#include "kDataframes/kDataFrameSTL.hpp"
#include "parallel_hashmap/phmap_dump.h"
#include <iostream>
#include <fstream>
#include "Utils/kmer.h"
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/archives/binary.hpp>

template
class kDataFramePHMAPIterator;

template
class kDataFramePHMAP;

template
class kDataFrameMAP;
template
class kDataFrameMAPIterator;

template
class kDataFrameBtree;
template
class kDataFrameBtreeIterator;




/*
 *****************************
 *** kDataFrameSTLIterator ***
 *****************************
 */

template<class MapType>
kDataFrameSTLIterator<MapType>::kDataFrameSTLIterator(itType it, kDataFrameSTL<MapType> *origin,
                                                      uint64_t kSize)
        : _kDataFrameIterator(kSize) {
    iterator = it;
    this->origin = origin;
    this->KD = this->origin->KD;
}

template<class MapType>
kDataFrameSTLIterator<MapType>::kDataFrameSTLIterator(const kDataFrameSTLIterator<MapType> &other) :
        _kDataFrameIterator(other.kSize) {
    iterator = other.iterator;
    this->origin = other.origin;
    this->KD = other.KD;
}

template<class MapType>
_kDataFrameIterator *kDataFrameSTLIterator<MapType>::clone() {
    return new kDataFrameSTLIterator(*this);
}

template<class MapType>
kDataFrameSTLIterator<MapType> &kDataFrameSTLIterator<MapType>::operator++(int) {
    iterator++;
    return *this;
}

template<class MapType>
uint64_t kDataFrameSTLIterator<MapType>::getHashedKmer() {
    //return origin->getHasher()->hash(iterator->first);
    return iterator->first;

}

template<class MapType>
string kDataFrameSTLIterator<MapType>::getKmer() {
    return origin->getkmerDecoder()->ihash_kmer(iterator->first);
    // return iterator->first;
}

template<class MapType>
uint64_t kDataFrameSTLIterator<MapType>::getCount() {
    return iterator->second;
}

template<class MapType>
uint64_t kDataFrameSTLIterator<MapType>::getOrder() {
    return iterator->second;
}

template<class MapType>
bool kDataFrameSTLIterator<MapType>::setOrder(uint64_t count) {
    iterator->second = count;
}

template<class MapType>
bool kDataFrameSTLIterator<MapType>::operator==(const _kDataFrameIterator &other) {
    return iterator == ((kDataFrameSTLIterator *) &other)->iterator;
}

template<class MapType>
bool kDataFrameSTLIterator<MapType>::operator!=(const _kDataFrameIterator &other) {
    return iterator != ((kDataFrameSTLIterator *) &other)->iterator;
}


template<class MapType>
kDataFrameSTLIterator<MapType>::~kDataFrameSTLIterator() {
}

/*
 **********************
 *** kDataFrameSTL ***
 **********************
 */

template<class MapType>
kDataFrameSTL<MapType>::kDataFrameSTL(uint64_t ksize) {
    this->class_name = "STL"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(kSize, TwoBits_hasher);
    endIterator= nullptr;
    reserve(1000000);
//    hasher = new wrapperHasher<flat_hash_map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
//    this->MAP = MapType(1000000);
//    endIterator = new kDataFrameIterator(
//            (_kDataFrameIterator *) new kDataFrameSTLIterator(MAP.end(), this, kSize),
//            (kDataFrame *) this);
    // this->hasher = (new IntegerHasher(ksize));
}


template<class MapType>
kDataFrameSTL<MapType>::kDataFrameSTL(uint64_t ksize, uint64_t nKmers) {
    this->class_name = "STL"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(kSize, TwoBits_hasher);
    endIterator= nullptr;
    reserve(nKmers);
//    hasher = new wrapperHasher<flat_hash_map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
//    this->MAP = MapType(nKmers);
//    endIterator = new kDataFrameIterator(
//            (_kDataFrameIterator *) new kDataFrameSTLIterator(MAP.end(), this, kSize),
//            (kDataFrame *) this);
    // this->hasher = (new IntegerHasher(ksize));
}

template<class MapType>
kDataFrameSTL<MapType>::kDataFrameSTL(uint64_t ksize, hashingModes hash_mode) {
    this->class_name = "STL"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(kSize, hash_mode);
    endIterator= nullptr;
    reserve(1000);
//    hasher = new wrapperHasher<flat_hash_map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
//    this->MAP = MapType(1000);
//    // this->hasher = (new IntegerHasher(ksize));
//    endIterator = new kDataFrameIterator(
//            (_kDataFrameIterator *) new kDataFrameSTLIterator(MAP.end(), this, kSize),
//            (kDataFrame *) this);
}

template<class MapType>
kDataFrameSTL<MapType>::kDataFrameSTL(readingModes RM, hashingModes hash_mode, map<string, int> params) {
    this->class_name = "STL"; // Temporary until resolving #17
    KD = kmerDecoder::getInstance(RM, hash_mode, params);
    this->kSize = KD->get_kSize();
    //this->MAP = MapType(1000);
    endIterator= nullptr;
    reserve(1000);
}

template<class MapType>
kDataFrameSTL<MapType>::kDataFrameSTL(uint64_t ksize, vector<uint64_t> kmersHistogram) {
    this->class_name = "STL"; // Temporary until resolving #17
    this->kSize = ksize;
    KD = new Kmers(kSize, TwoBits_hasher);
    endIterator= nullptr;
//    hasher = new wrapperHasher<flat_hash_map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);

    uint64_t countSum = 0;
    for (auto h:kmersHistogram)
        countSum += h;
    reserve(countSum);
//    hasher = new wrapperHasher<std::map<uint64_t, uint64_t>::hasher>(MAP.hash_function(), ksize);
//    this->MAP = std::map<uint64_t, uint64_t>(1000);
    // this->hasher = (new IntegerHasher(ksize));
//    endIterator = new kDataFrameIterator(
//            (_kDataFrameIterator *) new kDataFrameSTLIterator<MapType>(MAP.end(), this, kSize),
//            (kDataFrame *) this);
}

template<class MapType>
kDataFrameSTL<MapType>::kDataFrameSTL() {
    this->class_name = "STL"; // Temporary until resolving #17
    this->kSize = 23;
    this->MAP = MapType();
    KD = new Kmers(kSize, TwoBits_hasher);
    endIterator= nullptr;
    // hasher=new wrapperHasher<flat_hash_map<uint64_t,uint64_t>::hasher >(MAP.hash_function(),kSize);
    // this->hasher = (new IntegerHasher(23));
//    endIterator = new kDataFrameIterator(
//            (_kDataFrameIterator *) new kDataFrameSTLIterator<MapType>(MAP.end(), this, kSize),
//            (kDataFrame *) this);
}

template<class MapType>
bool kDataFrameSTL<MapType>::kmerExist(string kmerS) {
    return (this->MAP.find(KD->hash_kmer(kmerS)) == this->MAP.end()) ? 0 : 1;
}

template<class MapType>
bool kDataFrameSTL<MapType>::kmerExist(uint64_t kmer) {
    return (this->MAP.find(kmer) == this->MAP.end()) ? 0 : 1;

}


template<class MapType>
uint32_t kDataFrameSTL<MapType>::insert(const string &kmerS) {
    auto it = this->MAP.find(KD->hash_kmer(kmerS));
    if (it != this->MAP.end())
        return it->second;
    this->MAP[KD->hash_kmer(kmerS)] = lastKmerOrder++;
    return lastKmerOrder - 1;
}


template<class MapType>
uint32_t kDataFrameSTL<MapType>::insert(uint64_t kmer) {
    auto it = this->MAP.find(kmer);
    if (it != this->MAP.end())
        return it->second;
    this->MAP.insert(make_pair(kmer, lastKmerOrder++));
    //this->MAP[kmer] = lastKmerOrder++;
    return lastKmerOrder - 1;
}

template<class MapType>
bool kDataFrameSTL<MapType>::setOrder(const string &kmerS, uint64_t tag) {
    this->MAP[KD->hash_kmer(kmerS)] = tag;
    return true;
}

template<class MapType>
bool kDataFrameSTL<MapType>::setOrder(uint64_t kmerS, uint64_t tag) {
    this->MAP[kmerS] = tag;
    return true;
}

template<class MapType>
uint64_t kDataFrameSTL<MapType>::getkmerOrder(const string &kmerS) {
    typename MapType::const_iterator got = this->MAP.find(KD->hash_kmer(kmerS));
    if (got == this->MAP.end())
        return 0;
    else
        return got->second;
}

template<class MapType>
uint64_t kDataFrameSTL<MapType>::getkmerOrder(uint64_t kmerS) {
    typename MapType::iterator got = this->MAP.find(kmerS);
    if (got == this->MAP.end())
        return 0;
    else
        return got->second;
}

template<class MapType>
uint64_t kDataFrameSTL<MapType>::bucket(string kmerS) {
    return 1;
    // return this->MAP.bucket(kmer::str_to_canonical_int(kmerS));
}


template<class MapType>
bool kDataFrameSTL<MapType>::erase(const string &kmerS) {
    return this->MAP.erase(KD->hash_kmer(kmerS));
}

template<class MapType>
bool kDataFrameSTL<MapType>::erase(uint64_t kmer) {
    return this->MAP.erase(kmer);
}

template<class MapType>
uint64_t kDataFrameSTL<MapType>::size() {
    return (uint64_t) this->MAP.size();
}

template<class MapType>
uint64_t kDataFrameSTL<MapType>::max_size() {
    return (uint64_t) this->MAP.max_size();
}

template<class MapType>
float kDataFrameSTL<MapType>::load_factor() {
    return this->MAP.load_factor();
}

template<class MapType>
float kDataFrameSTL<MapType>::max_load_factor() {
    return this->MAP.max_load_factor();
}


template<>
float kDataFrameMAP::load_factor() {
    return 0.5;
}

template<>
float kDataFrameMAP::max_load_factor() {
    return 1.0;
}

template<>
float kDataFrameBtree::load_factor() {
    return 0.5;
}

template<>
float kDataFrameBtree::max_load_factor() {
    return 1.0;
}


template<class MapType>
void kDataFrameSTL<MapType>::serialize(string filePath) {

    ofstream file(filePath + ".extra");
    file << kSize << endl;
    file << this->KD->hash_mode << endl;
    file.close();
    std::ofstream os(filePath + ".map", std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(this->MAP);

}

template<class MapType>
kDataFrame *kDataFrameSTL<MapType>::load(string filePath) {

    // Load kSize
    ifstream file(filePath + ".extra");
    uint64_t kSize;
    int hashing_mode;
    file >> kSize;
    file >> hashing_mode;


    file.close();
    // Initialize kDataFrameMAP
    kDataFrameSTL<MapType> *KMAP = new kDataFrameSTL<MapType>(kSize);

    // Load the hashMap into the kDataFrameMAP
    std::ifstream os(filePath + ".map", std::ios::binary);
    cereal::BinaryInputArchive iarchive(os);
    iarchive(KMAP->MAP);
    delete KMAP->endIterator;
    KMAP->endIterator = new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameSTLIterator<MapType>(KMAP->MAP.end(), KMAP, kSize),
            (kDataFrame *) KMAP);

    return KMAP;
}

template<>
void kDataFramePHMAP::serialize(string filePath) {

    ofstream file(filePath + ".extra");
    file << kSize << endl;
    file << this->KD->hash_mode << endl;
    file.close();
    filePath += ".phmap";
    {
        phmap::BinaryOutputArchive ar_out(filePath.c_str());
        this->MAP.phmap_dump(ar_out);
    }


}

template<>
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
        KMAP->MAP.phmap_load(ar_in);
    }
    if (KMAP->endIterator != nullptr)
        delete KMAP->endIterator;
    KMAP->endIterator = new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFramePHMAPIterator(KMAP->MAP.end(), KMAP, kSize),
            (kDataFrame *) KMAP);

    return KMAP;
}

template<>
void kDataFrameMAP::serialize(string filePath) {

    ofstream file(filePath + ".extra");
    file << kSize << endl;
    file << this->KD->hash_mode << endl;
    file.close();
    std::ofstream os(filePath + ".map", std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(this->MAP);


}

template<>
kDataFrame *kDataFrameMAP::load(string filePath) {
    // Load kSize
    ifstream file(filePath + ".extra");
    uint64_t kSize;
    int hashing_mode;
    file >> kSize;
    file >> hashing_mode;


    file.close();
    // Initialize kDataFrameMAP
    kDataFrameMAP *KMAP = new kDataFrameMAP(kSize);

    // Load the hashMap into the kDataFrameMAP
    std::ifstream os(filePath + ".map", std::ios::binary);
    cereal::BinaryInputArchive iarchive(os);
    iarchive(KMAP->MAP);
    delete KMAP->endIterator;
    KMAP->endIterator = new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameMAPIterator(KMAP->MAP.end(), KMAP, kSize),
            (kDataFrame *) KMAP);

    return KMAP;
}

template<>
void kDataFrameBtree::serialize(string filePath) {

    ofstream file(filePath + ".extra");
    file << kSize << endl;
    file << this->KD->hash_mode << endl;
    file.close();
    std::ofstream os(filePath + ".btree", std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(this->MAP);


}

template<>
kDataFrame *kDataFrameBtree::load(string filePath) {
    // Load kSize
    ifstream file(filePath + ".extra");
    uint64_t kSize;
    int hashing_mode;
    file >> kSize;
    file >> hashing_mode;


    file.close();
    // Initialize kDataFrameMAP
    kDataFrameBtree *KMAP = new kDataFrameBtree(kSize);

    // Load the hashMap into the kDataFrameMAP
    std::ifstream os(filePath + ".btree", std::ios::binary);
    cereal::BinaryInputArchive iarchive(os);
    iarchive(KMAP->MAP);
    delete KMAP->endIterator;
    KMAP->endIterator = new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameBtreeIterator(KMAP->MAP.end(), KMAP, kSize),
            (kDataFrame *) KMAP);

    return KMAP;
}


template<class MapType>
kDataFrame *kDataFrameSTL<MapType>::getTwin() {
    return ((kDataFrame *) new kDataFrameSTL(kSize, this->KD->hash_mode));
}

template<class MapType>
void kDataFrameSTL<MapType>::_reserve(uint64_t n) {
    this->MAP.reserve(n);
    if (endIterator != nullptr)
        delete endIterator;
    endIterator = new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameSTLIterator<MapType>(MAP.end(), this, kSize),
            (kDataFrame *) this);

}

template<class MapType>
void kDataFrameSTL<MapType>::_reserve(vector<uint64_t> countHistogram) {
    uint64_t countSum = 0;
    for (auto h:countHistogram)
        countSum += h;
    reserve(countSum);


}

template<>
void kDataFrameMAP::_reserve(uint64_t n) {

    if (endIterator != nullptr)
        delete endIterator;
    endIterator = new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameMAPIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);

}

template<>
void kDataFrameBtree::_reserve(uint64_t n) {

    if (endIterator != nullptr)
        delete endIterator;
    endIterator = new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameBtreeIterator(MAP.end(), this, kSize),
            (kDataFrame *) this);

}

template<class MapType>
kDataFrameIterator kDataFrameSTL<MapType>::begin() {
    return (kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameSTLIterator<MapType>(MAP.begin(), this, kSize),
            (kDataFrame *) this));
}

template<class MapType>
kDataFrameIterator kDataFrameSTL<MapType>::end() {
    ((kDataFrameSTLIterator<MapType> *) endIterator->iterator)->iterator = MAP.end();
    return *endIterator;
}

template<class MapType>
kDataFrameIterator kDataFrameSTL<MapType>::find(const string &kmer) {
    return (kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameSTLIterator<MapType>(MAP.find(KD->hash_kmer(kmer)), this, kSize),
            (kDataFrame *) this));
}

template<class MapType>
kDataFrameIterator kDataFrameSTL<MapType>::find(uint64_t kmer) {
    return (kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameSTLIterator<MapType>(MAP.find((kmer)), this, kSize),
            (kDataFrame *) this));
}

template<class MapType>
MapType *kDataFrameSTL<MapType>::getMap() {
    return &MAP;
}

template<class MapType>
kDataFrame *kDataFrameSTL<MapType>::clone() {
    kDataFrameSTL *newOne = new kDataFrameSTL(kSize);
    newOne->MAP = MAP;
    for (auto c:columns) {
        newOne->columns[c.first] = c.second->clone();
        if (c.first == "count")
            newOne->countColumn = (vectorColumn<uint32_t> *) newOne->columns[c.first];
    }
    delete newOne->endIterator;
    newOne->endIterator = new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameSTLIterator<MapType>(newOne->MAP.end(), newOne, kSize),
            (kDataFrame *) newOne);
    return newOne;
}


kDataFrame *kDataFrameFactory::loadPHMAP(string filePath) {
    return kDataFramePHMAP::load(filePath);
}

kDataFrame *kDataFrameFactory::createPHMAP(uint32_t kSize, uint32_t numKmers) {
    return new kDataFramePHMAP(kSize,numKmers);
}

kDataFrame *kDataFrameFactory::createPHMAP(uint64_t ksize, hashingModes hash_mode) {
    return new kDataFramePHMAP(ksize,hash_mode);
}

kDataFrame *kDataFrameFactory::loadMAP(string filePath) {
    return kDataFrameMAP::load(filePath);
}

kDataFrame *kDataFrameFactory::createMAP(uint32_t kSize, uint32_t numKmers) {
    return new kDataFrameMAP(kSize,numKmers);
}

kDataFrame *kDataFrameFactory::loadBtree(string filePath) {
    return kDataFrameBtree::load(filePath);
}

kDataFrame *kDataFrameFactory::createBtree(uint32_t kSize, uint32_t numKmers) {
    return new kDataFrameBtree(kSize,numKmers);
}
