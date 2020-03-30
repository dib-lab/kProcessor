#include "kDataFrame.hpp"
#include "Utils/kmer.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <limits>
#include <sstream>


/*
 *****************************
 *** kDataFrameMQFIterator ***
 *****************************
 */


kDataFrameMQFIterator::kDataFrameMQFIterator(QF *mqf, uint64_t kSize, kmerDecoder *KD)
        : _kDataFrameIterator(kSize) {
    qfi = new QFi();
    qf_iterator(mqf, qfi, 0);
    this->KD = KD;
}
kDataFrameMQFIterator::kDataFrameMQFIterator(QFi *mqfIt, uint64_t kSize, kmerDecoder *KD)
        : _kDataFrameIterator(kSize) {
    qfi = mqfIt;
    this->KD = KD;
}
kDataFrameMQFIterator::kDataFrameMQFIterator(const kDataFrameMQFIterator &other) :
        _kDataFrameIterator(other.kSize) {
    qfi = new QFi();
    qfi->qf = other.qfi->qf;
    qfi->run = other.qfi->run;
    qfi->current = other.qfi->current;
    qfi->cur_start_index = other.qfi->cur_start_index;
    qfi->cur_length = other.qfi->cur_length;
    qfi->num_clusters = other.qfi->num_clusters;
    qfi->c_info = other.qfi->c_info;
    KD = other.KD;
}

_kDataFrameIterator *kDataFrameMQFIterator::clone() {
    return new kDataFrameMQFIterator(*this);
}

kDataFrameMQFIterator &kDataFrameMQFIterator::operator++(int) {
    qfi_next(qfi);
    return *this;
}

uint64_t kDataFrameMQFIterator::getHashedKmer() {
    uint64_t key, value, count;
    qfi_get(qfi, &key, &value, &count);
    return key;

}

string kDataFrameMQFIterator::getKmer() {
    return KD->ihash_kmer(getHashedKmer());
}

uint64_t kDataFrameMQFIterator::getCount() {
    uint64_t key, value, count;
    qfi_get(qfi, &key, &value, &count);
    return count;
}

bool kDataFrameMQFIterator::setCount(uint64_t count) {
    uint64_t key, value, currentCount;
    qfi_get(qfi, &key, &value, &currentCount);
    if (currentCount > count) {
        qf_remove(qfi->qf, key, currentCount - count, false, false);
    } else {
        qf_insert(qfi->qf, key, count - currentCount, false, false);
    }
    return true;
}

void kDataFrameMQFIterator::endIterator() {
    qfi->current = qfi->qf->metadata->xnslots;
}

bool kDataFrameMQFIterator::operator==(const _kDataFrameIterator &other) {
    if (qfi->current >= qfi->qf->metadata->xnslots &&
        ((kDataFrameMQFIterator *) &other)->qfi->current >=
        ((kDataFrameMQFIterator *) &other)->qfi->qf->metadata->xnslots)
        return true;

    return qfi->current == ((kDataFrameMQFIterator *) &other)->qfi->current;
}

bool kDataFrameMQFIterator::operator!=(const _kDataFrameIterator &other) {
    if (qfi->current >= qfi->qf->metadata->xnslots &&
        ((kDataFrameMQFIterator *) &other)->qfi->current >=
        ((kDataFrameMQFIterator *) &other)->qfi->qf->metadata->xnslots)
        return false;
    if (qfi->current >= qfi->qf->metadata->xnslots &&
        ((kDataFrameMQFIterator *) &other)->qfi->current <
        ((kDataFrameMQFIterator *) &other)->qfi->qf->metadata->xnslots)
        return true;
    if (qfi->current < qfi->qf->metadata->xnslots &&
        ((kDataFrameMQFIterator *) &other)->qfi->current >=
        ((kDataFrameMQFIterator *) &other)->qfi->qf->metadata->xnslots)
        return true;
    return qfi->current != ((kDataFrameMQFIterator *) &other)->qfi->current;
}

kDataFrameMQFIterator::~kDataFrameMQFIterator() {
    delete qfi;
}

kDataFrameIterator kDataFrameMQF::begin() {
    return (kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameMQFIterator(mqf, kSize, KD),
            (kDataFrame *) this));
}

kDataFrameIterator kDataFrameMQF::end() {
    kDataFrameMQFIterator *it = new kDataFrameMQFIterator(mqf, kSize, KD);
    it->endIterator();
    return (kDataFrameIterator(it, (kDataFrame *) this));
}
kDataFrameIterator kDataFrameMQF::find(string kmer) {
    QFi* mqfIt = new QFi();
    uint64_t hash = KD->hash_kmer(kmer) % mqf->metadata->range;
    qfi_find(mqf,mqfIt,hash);
    kDataFrameMQFIterator *it = new kDataFrameMQFIterator(mqfIt, kSize, KD);
    return (kDataFrameIterator(it, (kDataFrame *) this));
}

/*
 *********************
 *** kDataFrameMQF ***
 *********************
 */


kDataFrameMQF::kDataFrameMQF() : kDataFrame() {
    this->class_name = "MQF"; // Temporary until resolving #17
    mqf = new QF();
    qf_init(mqf, (1ULL << 16), 2 * kSize, 0, 2, 32, true, "", 2038074761);
    KD = (new Kmers(kSize));
    falsePositiveRate = 0;
    hashbits = 2 * kSize;
    range = (1ULL << hashbits);
}

kDataFrameMQF::kDataFrameMQF(uint64_t ksize, uint8_t q, int mode) : kDataFrame(ksize) {

    this->class_name = "MQF"; // Temporary until resolving #17

    // Mode 0: Murmar Hashing | Irreversible
    // Mode 1: Integer Hashing | Reversible | Full Hashing
    // Mode 2: TwoBitsHashing | Not considered hashing, just store the two bits representation

    KD = new Kmers(ksize, mode);


//    switch (mode) {
//        case 0:
//            hasher = (new MumurHasher(2038074761));
//            break;
//        case 1:
//            hasher = (new IntegerHasher(kSize));
//            break;
//        case 2:
//            hasher = (new QHasher(kSize, q));
//            break;
//        case 3:
//            hasher = (new TwoBitsHasher(kSize));
//            break;
//        default:
//            hasher = (new IntegerHasher(kSize));
//            break;
//    }

    mqf = new QF();
    qf_init(mqf, (1ULL << q), 2 * ksize, 0, 2, 32, true, "", 2038074761);
    this->falsePositiveRate = falsePositiveRate;
    hashbits = 2 * kSize;
    range = (1ULL << hashbits);

}

kDataFrameMQF::kDataFrameMQF(uint64_t ksize, uint8_t q, uint8_t fixedCounterSize, uint8_t tagSize,
                             double falsePositiveRate) :
        kDataFrame(ksize) {
    this->class_name = "MQF"; // Temporary until resolving #17
    mqf = new QF();
    qf_init(mqf, (1ULL << q), 2 * ksize, tagSize, fixedCounterSize, 32, true, "", 2038074761);
    this->falsePositiveRate = falsePositiveRate;
    if (falsePositiveRate == 0) {
        KD = (new Kmers(kSize,1));
    } else if (falsePositiveRate < 1) {
        KD = (new Kmers(kSize, 0));
    }
    hashbits = 2 * kSize;
    range = (1ULL << hashbits);

}

kDataFrameMQF::kDataFrameMQF(uint64_t ksize, int mode) :
    kDataFrame(ksize) {
    this->class_name = "MQF"; // Temporary until resolving #17
    this->falsePositiveRate = 0.0;
    KD = new Kmers(ksize, mode);
    hashbits = 2 * kSize;
    range = (1ULL << hashbits);
    mqf = NULL;
    reserve(1000000);
}

kDataFrameMQF::kDataFrameMQF(uint64_t ksize) :
        kDataFrame(ksize) {
    this->class_name = "MQF"; // Temporary until resolving #17
    this->falsePositiveRate = 0.0;
    KD = (new Kmers(kSize));
    hashbits = 2 * kSize;
    range = (1ULL << hashbits);
    mqf = NULL;
    reserve(1000000);
}

kDataFrameMQF::kDataFrameMQF(QF *mqf, uint64_t ksize, double falsePositiveRate) :
        kDataFrame(ksize) {
    this->class_name = "MQF"; // Temporary until resolving #17
    this->mqf = mqf;
    this->falsePositiveRate = falsePositiveRate;
    if (falsePositiveRate == 0) {
        KD = (new Kmers(kSize));
    } else if (falsePositiveRate < 1) {
        KD = (new Kmers(kSize, 0));
    }
    hashbits = this->mqf->metadata->key_bits;
    hashbits = 2 * kSize;
    range = (1ULL << hashbits);
}

kDataFrame *kDataFrameMQF::getTwin() {
    uint64_t q = log2(mqf->metadata->nslots);
    return ((kDataFrame *) new kDataFrameMQF(kSize, q, mqf->metadata->fixed_counter_size,
                                             mqf->metadata->label_bits, falsePositiveRate));
}

void kDataFrameMQF::reserve(uint64_t n) {
    QF *old = mqf;
    mqf = new QF();
    uint64_t q = (uint64_t) ceil(log2((double) n * 1.4));
//    std::cerr << "[DEBUG] Q: " << q << std::endl;
    qf_init(mqf, (1ULL << q), hashbits, 0, 2, 32, true, "", 2038074761);
    if (old != NULL) {
        qf_migrate(old, mqf);
        qf_destroy(old);
        delete old;
    }
}
void kDataFrameMQF::reserve(vector<uint64_t> countHistogram) {
    QF *old = mqf;
    mqf = new QF();
    uint64_t nSlots;
    uint64_t fixedCounterSize;
    uint64_t memory;
    kDataFrameMQF::estimateParameters(countHistogram, 2 * getkSize(), 0,
                                      &nSlots, &fixedCounterSize, &memory);
//    std::cerr << "[DEBUG] Q: " << q << std::endl;
    qf_init(mqf, nSlots, 2 * getkSize(), 0, fixedCounterSize, 0, true, "", 2038074761);
    if (old != NULL) {
        qf_migrate(old, mqf);
        qf_destroy(old);
        delete old;
    }
}
kDataFrameMQF::kDataFrameMQF(uint64_t ksize, vector<uint64_t> countHistogram, uint8_t tagSize, double falsePositiveRate)
        :
        kDataFrame(ksize) {
    this->class_name = "MQF"; // Temporary until resolving #17
    uint64_t nSlots;
    uint64_t fixedCounterSize;
    uint64_t memory;
    kDataFrameMQF::estimateParameters(countHistogram, 2 * ksize, tagSize,
                                      &nSlots, &fixedCounterSize, &memory);
    qf_init(mqf, nSlots, 2 * ksize, tagSize, fixedCounterSize, 32, true, "", 2038074761);
}
kDataFrameMQF::kDataFrameMQF(uint64_t ksize, vector<uint64_t> countHistogram)
        :
        kDataFrameMQF(ksize,countHistogram,0,0.0) {
}
uint64_t kDataFrameMQF::estimateMemory(uint64_t nslots, uint64_t slotSize, uint64_t fcounter, uint64_t tagSize) {
    uint64_t SLOTS_PER_BLOCK_2 = 64;
    uint64_t xnslots = nslots + 10 * sqrt((double) nslots);
    uint64_t nblocks = (xnslots + SLOTS_PER_BLOCK_2 - 1) / SLOTS_PER_BLOCK_2;
    uint64_t blocksize = 17;

    return ((nblocks) * (blocksize + 8 * (slotSize + fcounter + tagSize))) / 1024;

}

void kDataFrameMQF::estimateParameters(vector<uint64_t> countHistogram,
                                       uint64_t numHashBits, uint64_t tagSize,
                                       uint64_t *res_noSlots, uint64_t *res_fixedSizeCounter, uint64_t *res_memory) {

    uint64_t noDistinctKmers = countHistogram[0];
    *res_memory = numeric_limits<uint64_t>::max();
    for (int i = 8; i < 64; i++) {
        uint64_t noSlots = (1ULL) << i;
        if (noSlots < noDistinctKmers)
            continue;
        bool moreWork = false;
        uint64_t slotSize = numHashBits - log2((double) noSlots);
        for (uint64_t fixedSizeCounter = 1; fixedSizeCounter < slotSize; fixedSizeCounter++) {
            if (isEnough(countHistogram, noSlots, fixedSizeCounter, slotSize)) {
                uint64_t tmpMem = estimateMemory(noSlots, slotSize, fixedSizeCounter, tagSize);
                if (*res_memory > tmpMem) {
                    *res_memory = tmpMem;
                    *res_fixedSizeCounter = fixedSizeCounter;
                    *res_noSlots = noSlots;
                    moreWork = true;
                } else {
                    break;
                }
            }

        }
        if (!moreWork && *res_memory != numeric_limits<uint64_t>::max())
            break;
    }
    if (*res_memory == numeric_limits<uint64_t>::max()) {
        throw std::overflow_error("Data limits exceeds MQF capabilities(> uint64). Check if ntcard file is corrupted");
    }

}

bool
kDataFrameMQF::isEnough(vector<uint64_t> histogram, uint64_t noSlots, uint64_t fixedSizeCounter, uint64_t slotSize) {
    // cout<<"noSlots= "<<noSlots<<endl
    //     <<"fcounter= "<<fixedSizeCounter<<endl
    //     <<"slot size= "<<numHashBits<<endl;

    noSlots = (uint64_t) ((double) noSlots * 0.90);
    for (uint64_t i = 1; i < 1000; i++) {
        uint64_t usedSlots = 1;

        if (i > ((1ULL) << fixedSizeCounter) - 1) {
            uint64_t nSlots2 = 0;
            __uint128_t capacity;
            do {
                nSlots2++;
                capacity = ((__uint128_t) (1ULL) << (nSlots2 * slotSize + fixedSizeCounter)) - 1;
                //  cout<<"slots num "<<nSlots2<<" "<<capacity<<endl;
            } while ((__uint128_t) i > capacity);
            usedSlots += nSlots2;
        }
        //cout<<"i= "<<i<<"->"<<usedSlots<<" * "<<histogram[i]<<endl;
        if (noSlots >= (usedSlots * histogram[i])) {
            noSlots -= (usedSlots * histogram[i]);
        } else {
            //  cout<<"failed"<<endl<<endl;
            return false;
        }

    }
    //cout<<"success"<<endl<<endl;
    return true;
}


bool kDataFrameMQF::setCount(string kmer, uint64_t count) {
    uint64_t hash = KD->hash_kmer(kmer) % mqf->metadata->range;
    uint64_t currentCount = qf_count_key(mqf, hash);
    if (currentCount > count) {
        qf_remove(mqf, hash, currentCount - count, false, false);
    } else {
        try {
            qf_insert(mqf, hash, count - currentCount, false, false);
        }
        catch (overflow_error &e) {
            reserve(mqf->metadata->nslots);
            return setCount(kmer, count);
        }
    }
    return true;
}

bool kDataFrameMQF::setCount(uint64_t kmer, uint64_t count) {
    uint64_t hash = kmer % mqf->metadata->range;
    uint64_t currentCount = qf_count_key(mqf, hash);
    if (currentCount > count) {
        qf_remove(mqf, hash, currentCount - count, false, false);
    } else {
        try {
            qf_insert(mqf, hash, count - currentCount, false, false);
        }
        catch (overflow_error &e) {
            reserve(mqf->metadata->nslots);
            return setCount(kmer, count);
        }
    }
    return true;
}

bool kDataFrameMQF::insert(string kmer, uint64_t count) {
    uint64_t hash = KD->hash_kmer(kmer) % mqf->metadata->range;
    try {
        qf_insert(mqf, hash, count, true, true);
    }
    catch (overflow_error &e) {
        reserve(mqf->metadata->nslots);
        return insert(kmer, count);
    }
    return true;
}

bool kDataFrameMQF::insert(string kmer) {
    if (load_factor() > 0.9)
        reserve(mqf->metadata->nslots);
    uint64_t hash = KD->hash_kmer(kmer) % mqf->metadata->range;
    // cout << "Inserting kmer: " << kmer << ", Hash: " << hash << endl;
    try {
        qf_insert(mqf, hash, 1, false, false);
    }
    catch (overflow_error &e) {
        reserve(mqf->metadata->nslots);
        return insert(kmer);
    }
    return true;
}


bool kDataFrameMQF::insert(uint64_t kmer, uint64_t count) {
    uint64_t hash = kmer % mqf->metadata->range;
    try {
        qf_insert(mqf, hash, count, true, true);
    }
    catch (overflow_error &e) {
        reserve(mqf->metadata->nslots);
        return insert(kmer, count);
    }
    return true;
}

bool kDataFrameMQF::insert(uint64_t kmer) {
    if (load_factor() > 0.9)
        reserve(mqf->metadata->nslots);
    uint64_t hash = kmer % mqf->metadata->range;
    // cout << "Inserting kmer: " << kmer << ", Hash: " << hash << endl;
    try {
        qf_insert(mqf, hash, 1, false, false);
    }
    catch (overflow_error &e) {
        reserve(mqf->metadata->nslots);
        return insert(kmer);
    }
    return true;
}

uint64_t kDataFrameMQF::getCount(string kmer) {
    uint64_t hash = KD->hash_kmer(kmer) % mqf->metadata->range;
    return qf_count_key(mqf, hash);
}

uint64_t kDataFrameMQF::getCount(uint64_t kmer) {
    uint64_t hash = kmer % mqf->metadata->range;
    return qf_count_key(mqf, hash);
}


bool kDataFrameMQF::erase(string kmer) {
    uint64_t hash = KD->hash_kmer(kmer) % mqf->metadata->range;
    uint64_t currentCount = qf_count_key(mqf, hash);

    qf_remove(mqf, hash, currentCount, false, false);
    return true;
}

bool kDataFrameMQF::erase(uint64_t kmer) {
    uint64_t hash = kmer % mqf->metadata->range;
    uint64_t currentCount = qf_count_key(mqf, hash);

    qf_remove(mqf, hash, currentCount, false, false);
    return true;
}

uint64_t kDataFrameMQF::size() {
    return mqf->metadata->ndistinct_elts;
}

uint64_t kDataFrameMQF::max_size() {
    return mqf->metadata->xnslots;
}

float kDataFrameMQF::load_factor() {
    return (float) qf_space(mqf) / 100.0;
}

float kDataFrameMQF::max_load_factor() {
    return 0.9;
}


void kDataFrameMQF::serialize(string filePath) {
    //filePath += ".mqf";
    ofstream file(filePath + ".extra");
    file << kSize << endl;
    file << this->KD->hash_mode << endl;
    file.close();
    // uint64_t legendSize=tagsLegend.size();
    // file<<legendSize<<endl;
    // auto it = tagsLegend.begin();
    // while(it==tagsLegend.end())
    // {
    //   file<<it->first<<" "<<it->second<<endl;
    //   it++;
    // }
    // file.close();
    qf_serialize(mqf, (filePath + ".mqf").c_str());
}

kDataFrame *kDataFrameMQF::load(string filePath) {
    ifstream file(filePath + ".extra");
    uint64_t filekSize, hashing_mode;
    file >> filekSize;
    file >> hashing_mode;
    double flasePositiveRate;
    flasePositiveRate = (hashing_mode == 1) ? 0 : 0.5;
    file.close();
    QF *mqf = new QF();
    qf_deserialize(mqf, (filePath + ".mqf").c_str());
    return new kDataFrameMQF(mqf, filekSize, flasePositiveRate);
}

void kDataFrameMQF::preprocessKmerOrder()
{
  qf_ComputeItemsOrder(mqf);
}
uint64_t kDataFrameMQF::getkmerOrder(string kmer)
{
  uint64_t hash = KD->hash_kmer(kmer) % mqf->metadata->range;
  return itemOrder(mqf,hash);
}
