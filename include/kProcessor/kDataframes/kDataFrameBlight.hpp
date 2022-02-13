#ifndef _kDataFRAMEBlight_H_
#define _kDataFRAMEBlight_H_

#include "kDataFrame.hpp"
#include "blight.h"

class kDataFrameBlight;


class kDataFrameBlightIterator : public _kDataFrameIterator {
private:
    kmer_Set_Light_iterator iterator;
    kDataFrameBlight *origin;
public:
    kDataFrameBlightIterator(kmer_Set_Light_iterator it, kDataFrameBlight *origin, std::uint64_t kSize);

    kDataFrameBlightIterator(const kDataFrameBlightIterator &);

    kDataFrameBlightIterator &operator++(int);

    _kDataFrameIterator *clone();

    std::uint64_t getHashedKmer();

    string getKmer();

    uint64_t getOrder();
    uint64_t getCount();

    bool setOrder(uint64_t count);

    void endIterator();

    bool operator==(const _kDataFrameIterator &other);

    bool operator!=(const _kDataFrameIterator &other);

    ~kDataFrameBlightIterator();
};


// kDataFrameBlight _____________________________

class kDataFrameBlight : public kDataFrame {
private:
    kmer_Set_Light* blight_index;

public:
    kDataFrameBlight();

    kDataFrameBlight(uint64_t ksize)
    {
        blight_index=new kmer_Set_Light(ksize);
        kSize=ksize;
    }

    kDataFrameBlight(std::uint64_t ksize,string input_fasta_file);
    kDataFrame *getTwin();

    void _reserve(std::uint64_t n);
    void _reserve (vector<std::uint64_t> countHistogram);

    bool kmerExist(string kmer);
    bool kmerExist(uint64_t kmer);

    bool setOrder(const string &kmer, std::uint64_t count);
    bool setOrder(std::uint64_t kmer, std::uint64_t count);

    uint32_t insert(const string &kmer) override;
    uint32_t insert(std::uint64_t kmer) override;

    std::uint64_t getkmerOrder(const string &kmer);
    std::uint64_t getkmerOrder(std::uint64_t kmer);

    bool erase(const string &kmer);
    bool erase(std::uint64_t kmer);

    std::uint64_t size();

    std::uint64_t max_size();

    float load_factor();

    float max_load_factor();

    kDataFrameIterator begin();

    // kDataFrameIterator end();
    kDataFrameIterator find(const string &kmer);
    kDataFrameIterator find(uint64_t kmer);


    void serialize(string filePath);

    static kDataFrame *load(string filePath);

    ~kDataFrameBlight() = default;

    kDataFrame *clone() override;

    template<typename T,typename Container>
    T getKmerColumnValue(const string& columnName,string kmer);

    template<typename T,typename Container>
    void setKmerColumnValue(const string& columnName,string kmer, T value);


};


#endif
