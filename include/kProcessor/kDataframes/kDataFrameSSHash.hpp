#ifndef _kDataFRAMESSSHASH_H_
#define _kDataFRAMESSSHASH_H_
#include "kDataFrame.hpp"
#include "dictionary.hpp"


class kDataFrame_sshash;

class kDataFrame_sshashIterator : public _kDataFrameIterator {
private:
    sshash::dictionary::iterator* iterator;
    kDataFrame_sshash *origin;
    std::pair<uint64_t, std::string> currKmer;
    uint32_t kmerID;
public:
    kDataFrame_sshashIterator();
    kDataFrame_sshashIterator(uint32_t kmerID, kDataFrame_sshash *origin, std::uint64_t kSize);

    kDataFrame_sshashIterator(const kDataFrame_sshashIterator &);

    kDataFrame_sshashIterator &operator++(int);

    _kDataFrameIterator *clone();

    std::uint64_t getHashedKmer();

    string getKmer();

    uint64_t getOrder();
    uint64_t getCount();

    bool setOrder(uint64_t count);

    void endIterator();

    bool operator==(const _kDataFrameIterator &other);

    bool operator!=(const _kDataFrameIterator &other);

    ~kDataFrame_sshashIterator();
};

class kDataFrame_sshash : public kDataFrame {


public:
    sshash::dictionary dict;
    kDataFrame_sshash();

    kDataFrame_sshash(uint64_t ksize)
    {
        kSize=ksize;
        KD= nullptr;
    }

    kDataFrame_sshash(std::uint64_t ksize,string input_fasta_file,std::uint64_t minimizer=9);
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

    ~kDataFrame_sshash() = default;

    kDataFrame *clone() override;

    template<typename T,typename Container>
    T getKmerColumnValue(const string& columnName,string kmer);

    template<typename T,typename Container>
    void setKmerColumnValue(const string& columnName,string kmer, T value);






};

#endif
