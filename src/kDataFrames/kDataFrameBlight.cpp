#include "kDataFrame.hpp"
#include <iostream>
#include <fstream>
#include "Utils/kmer.h"
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


template<typename T,typename Container>
T kDataFrameBlight::getKmerColumnValue(const string& columnName,string kmer)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    return ((Container*)columns[columnName])->get(kmerOrder);
}
template<typename T,typename Container>
void kDataFrameBlight::setKmerColumnValue(const string& columnName,string kmer,T value)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    ((Container*)columns[columnName])->insert(value,kmerOrder);
}



/*
 *****************************
 *** kDataFrameBlightIterator ***
 *****************************
 */

kDataFrameBlightIterator::kDataFrameBlightIterator(kmer_Set_Light_iterator it, kDataFrameBlight *origin, uint64_t kSize)
        : _kDataFrameIterator(kSize) {
    iterator = it;
    this->origin = origin;
}

kDataFrameBlightIterator::kDataFrameBlightIterator(const kDataFrameBlightIterator &other) :
        _kDataFrameIterator(other.kSize) {
    iterator = other.iterator;
    this->origin = other.origin;
}

_kDataFrameIterator *kDataFrameBlightIterator::clone() {
    return new kDataFrameBlightIterator(*this);
}

kDataFrameBlightIterator &kDataFrameBlightIterator::operator++(int) {
    iterator.next();
    if(iterator.kmer_id > iterator.index_ptr->number_kmer+1)
        iterator.kmer_id = iterator.index_ptr->number_kmer+1;
    return *this;
}

uint64_t kDataFrameBlightIterator::getHashedKmer() {
    //return origin->getHasher()->hash(iterator->first);
    return iterator.get_kmer();

}

string kDataFrameBlightIterator::getKmer() {
    return iterator.get_kmer_str();
    // return iterator->first;
}

uint64_t kDataFrameBlightIterator::getCount() {
    return origin->getkmerOrder(iterator.get_kmer_str());
}

uint64_t kDataFrameBlightIterator::getOrder() {
    return origin->getkmerOrder(iterator.get_kmer_str());
}

bool kDataFrameBlightIterator::setOrder(uint64_t count) {
    throw logic_error("kmerBlight is static");
}


bool kDataFrameBlightIterator::operator==(const _kDataFrameIterator &other) {

    return iterator.kmer_id == ((kDataFrameBlightIterator *) &other)->iterator.kmer_id;
}

bool kDataFrameBlightIterator::operator!=(const _kDataFrameIterator &other) {
    return iterator.kmer_id != ((kDataFrameBlightIterator *) &other)->iterator.kmer_id;
}

kDataFrameBlightIterator::~kDataFrameBlightIterator() {

}
string gen_randomBlight(const int len) {
    srand (time(NULL));
    static const char alphanum[] =
            "0123456789"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            "abcdefghijklmnopqrstuvwxyz";
    string s="blight.";
    for (int i = 0; i < len; ++i) {
        s+= alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    return s;
}
/*
 **********************
 *** kDataFrameBlight ***
 **********************
 */
kDataFrameBlight::kDataFrameBlight()
{

}

kDataFrameBlight::kDataFrameBlight(uint64_t ksize,string input_fasta_file) {
    this->class_name = "Blight"; // Temporary until resolving #17
    this->kSize = ksize;

    int core_number(31);
    int minimizer_size(10);
    int file_number_exponent(4);
    int subsampling_bits(0);

    blight_index=  new kmer_Set_Light(ksize, core_number, minimizer_size, file_number_exponent, subsampling_bits);
 //   string workingDirectory="./workdir"+ gen_randomBlight(10);
    string workingDirectory="./workdir";
 //   string mkdirCommand="mkdir -p "+workingDirectory;
//    system(mkdirCommand.c_str());
    blight_index->construct_index(input_fasta_file, workingDirectory);

    kmer_Set_Light_iterator it(blight_index);
    it.kmer_id=it.index_ptr->number_kmer+1;

    endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameBlightIterator(it, this, kSize),
            (kDataFrame *) this);
    // this->hasher = (new IntegerHasher(ksize));
}


bool kDataFrameBlight::kmerExist(string kmerS) {
    return blight_index->get_presence_query(kmerS)[0];
}
bool kDataFrameBlight::kmerExist(uint64_t kmerS) {
    throw logic_error("kDataFrameBlight doesnt implement kmerExit for uint64_T input. use string version instead");
}



uint32_t kDataFrameBlight::insert(const string &kmerS) {
    throw logic_error("kDataFrameBlight is static. Insertion is not allowed");
    return true;
}



uint32_t kDataFrameBlight::insert(uint64_t kmer) {
    throw logic_error("kDataFrameBlight is static. Insertion is not allowed");
    return true;
}

bool  kDataFrameBlight::setOrder(const string &kmer, std::uint64_t count){
    throw logic_error("kDataFrameBlight is static. Insertion is not allowed");
    return true;
}
bool  kDataFrameBlight::setOrder(std::uint64_t kmer, std::uint64_t count){
    throw logic_error("kDataFrameBlight is static. Insertion is not allowed");
    return true;
}


uint64_t kDataFrameBlight::getkmerOrder(const string &kmerS) {
    return blight_index->get_hashes_query(kmerS)[0];
}

uint64_t kDataFrameBlight::getkmerOrder(uint64_t kmerS) {
    return kmerS;
}




bool kDataFrameBlight::erase(const string &kmerS) {
    throw logic_error("kDataFrameBlight is static. Deletion is not allowed");
    return true;
}

bool kDataFrameBlight::erase(uint64_t kmer) {
    throw logic_error("kDataFrameBlight is static. Deletion is not allowed");
    return true;
}

uint64_t kDataFrameBlight::size() {
    return (uint64_t) this->blight_index->number_kmer;
}

uint64_t kDataFrameBlight::max_size() {
    return (uint64_t) this->blight_index->number_kmer;
}

float kDataFrameBlight::load_factor() {
    return 1;
}

float kDataFrameBlight::max_load_factor() {
    return 1;
}


void kDataFrameBlight::serialize(string filePath) {

    // Write the kmerSize
    ofstream file(filePath + ".extra");
    file << kSize << endl;
   // file << this->KD->hash_mode << endl;
    file.close();
    filePath += ".blight.gz";
    {
        blight_index->dump_disk(filePath);
    }

}

kDataFrame *kDataFrameBlight::load(string filePath) {

    // Load kSize
    ifstream file(filePath + ".extra");
    uint64_t kSize;
    file >> kSize;
    file.close();
    filePath += ".blight.gz";

    kDataFrameBlight *KMAP = new kDataFrameBlight();
    KMAP->kSize = kSize;
    {
        cout<<filePath<<endl;
        KMAP->blight_index= new kmer_Set_Light(filePath);

    }
    if(KMAP->endIterator != nullptr)
        delete KMAP->endIterator;



    kmer_Set_Light_iterator it(KMAP->blight_index);
    it.kmer_id=it.index_ptr->number_kmer+1;

    KMAP->endIterator= new kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameBlightIterator(it, KMAP, kSize),
            (kDataFrame *) KMAP);

    return KMAP;
}

kDataFrame *kDataFrameBlight::getTwin() {
    return ((kDataFrame *) new kDataFrameBlight(kSize));
}

void kDataFrameBlight::_reserve(uint64_t n) {
    throw logic_error("kDataFrameBlight is static. _reserve is not allowed");

}
void kDataFrameBlight::_reserve(vector<uint64_t> countHistogram) {
    throw logic_error("kDataFrameBlight is static. _reserve is not allowed");


}
kDataFrameIterator kDataFrameBlight::begin() {

    return (kDataFrameIterator(
            (_kDataFrameIterator *) new kDataFrameBlightIterator(kmer_Set_Light_iterator(blight_index), this, kSize),
            (kDataFrame *) this));
}

//kDataFrameIterator kDataFrameBlight::end() {
//    return *(new kDataFrameIterator(
//            (_kDataFrameIterator *) new kDataFrameBlightIterator(MAP.end(), this, kSize),
//            (kDataFrame *) this));
//}
kDataFrameIterator kDataFrameBlight::find(const string &kmer) {
    throw logic_error("not implemented yet");
}
kDataFrameIterator kDataFrameBlight::find(uint64_t kmer) {
    throw logic_error("not implemented yet");
}

kDataFrame *kDataFrameBlight::clone() {
    throw logic_error("not implemented yet");
}

kDataFrame *kDataFrameFactory::loadBlight(string filePath) {
    return kDataFrameBlight::load(filePath);
}

kDataFrame *kDataFrameFactory::createBlight(uint32_t kSize, string filePath) {
    return new kDataFrameBlight(kSize,filePath);
}

