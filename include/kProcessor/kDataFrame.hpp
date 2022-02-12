#ifndef _kDataFRAME_H_
#define _kDataFRAME_H_

#include <vector>

#include <iostream>
#include <parallel_hashmap/phmap.h>
#include <parallel_hashmap/btree.h>
#include "kmerDecoder.hpp"
#include <any>
//#include "dictionary.hpp"


#include "defaultColumn.hpp"
#include <cstdint>

using phmap::flat_hash_map;
using namespace std;

class kDataFrame;
class kDataFrameMQF;
class kDataFrameBMQF;

#define kDataFramePHMAP kDataFrameSTL<phmap::parallel_flat_hash_map<std::uint64_t,std::uint32_t,std::hash<uint64_t>,std::equal_to<uint64_t>,std::allocator<std::pair<const uint64_t, uint64_t>>,4,std::mutex>>
#define kDataFramePHMAPIterator kDataFrameSTLIterator<phmap::parallel_flat_hash_map<std::uint64_t,std::uint32_t,std::hash<uint64_t>,std::equal_to<uint64_t>,std::allocator<std::pair<const uint64_t, uint64_t>>,4,std::mutex>>
#define kDataFrameMAP kDataFrameSTL<std::map<std::uint64_t, std::uint32_t>>
#define kDataFrameMAPIterator kDataFrameSTLIterator<std::map<std::uint64_t, std::uint32_t>>
#define kDataFrameBtree kDataFrameSTL<phmap::btree_map<uint64_t,uint32_t>>
#define kDataFrameBtreeIterator kDataFrameSTLIterator<phmap::btree_map<uint64_t,uint32_t>>




class kmerRow{
public:
  string kmer;
  uint64_t hashedKmer;
  uint64_t count;
  uint64_t order;
  kDataFrame * origin;
  kmerRow(){
    kmer="";
    hashedKmer=0;
    count=0;
    order=0;
    origin=nullptr;
  }
  kmerRow(string kmer,std::uint64_t hashedKmer,std::uint64_t count,std::uint64_t order,kDataFrame* o)
  {
    this->kmer=kmer;
    this->hashedKmer=hashedKmer;
    this->count=count;
    this->order=order;
    this->origin = o;
  }
  kmerRow(const kmerRow& other)
  {
    kmer=other.kmer;
    hashedKmer=other.hashedKmer;
    count=other.count;
    origin = other.origin ;
    order=0;
  }

  static kmerRow copy(const kmerRow& other)
  {
    return * (new kmerRow(other));
  }

  template<typename T,typename Container>
  void getColumnValue(const string& colName, T& res);

  template<typename T,typename Container>
  void setColumnValue(const string& colName, T value);

  bool operator==(kmerRow &other) const
  {
    return hashedKmer == other.hashedKmer;
  }

  bool operator < (kmerRow& other) const
  {
    return hashedKmer<other.hashedKmer;
  }
  bool operator > ( kmerRow& other) const
  {
    return hashedKmer>other.hashedKmer;
  }
  bool operator < (const kmerRow& other) const
  {
    return hashedKmer<other.hashedKmer;
  }
  bool operator > (const kmerRow& other) const
  {
    return hashedKmer>other.hashedKmer;
  }

};

class _kDataFrameIterator{
protected:
    std::uint64_t kSize;
public:
  _kDataFrameIterator()= default;
  explicit _kDataFrameIterator(std::uint64_t k):kSize(k){}
  virtual uint64_t getOrder()=0;
  virtual _kDataFrameIterator& operator ++ (int)=0;
  virtual _kDataFrameIterator* clone()=0;
  virtual std::uint64_t getHashedKmer()=0;
  virtual string getKmer()=0;
  virtual std::uint64_t getCount()=0;
  virtual bool setOrder(std::uint64_t count)=0;
  virtual bool operator ==(const _kDataFrameIterator& other)=0;
  virtual bool operator !=(const _kDataFrameIterator& other)=0;
  virtual ~_kDataFrameIterator()= default;

};

class kDataFrameIterator{
private:
  kDataFrame* origin;
public:
    _kDataFrameIterator* iterator;
    using iterator_category = std::forward_iterator_tag;
  kDataFrameIterator(){
    iterator=nullptr;
  }
  kDataFrameIterator(const kDataFrameIterator& other){
    if(other.iterator!=nullptr){
      iterator=other.iterator->clone();
      origin=other.origin;
    }
    else{
      iterator=nullptr;
      origin=nullptr;
    }
  }
  kDataFrameIterator& operator= (const kDataFrameIterator& other){
    if(other.iterator!=nullptr){
      iterator=other.iterator->clone();
      origin=other.origin;
    }
    else{
      iterator=nullptr;
      origin=nullptr;
    }
    return *this;
  }
  kDataFrameIterator(_kDataFrameIterator* it,kDataFrame* o){
    this->origin=o;
    iterator=it;
  }
/// Increment the iterator to the next kmer
    kDataFrameIterator& operator * (){
        return *this;
    }

/// Increment the iterator to the next kmer
  kDataFrameIterator& operator ++ (int){
    (*iterator)++;
    return *this;
  }
    kDataFrameIterator& operator ++ (){
        (*iterator)++;
        return *this;
    }

  /// Increment the iterator to the next kmer (Implemented mainly for python interface)
  kDataFrameIterator& next(){
    (*iterator)++;
    return *this;
  }

// /// Increment the iterator to the next kmer
//   kDataFrameIterator operator ++ (int){
//     kDataFrameIterator temp=*this;
//     (*iterator)++;
//     return temp;
//   }
// /// Compare the position of each iterator in the underlying datastructure.
// /*! returns True when other points to kmer points to a further position than the current */
//   bool operator <(const kDataFrameIterator& other){
//     return *iterator < *other.iterator;
//   }
//   /// Compare the position of each iterator in the underlying datastructure.
//   /*! returns True when other points to kmer points to a nerarer position than the current */
//   bool operator >(const kDataFrameIterator& other){
//     return *iterator > *other.iterator;
//   }
//   /// Compare the position of each iterator in the underlying datastructure.
//   /*! returns True when other points to kmer points to a further or equal position than the current */
//   bool operator <=(const kDataFrameIterator& other){
//     return *iterator <= *other.iterator;
//   }
//   /// Compare the position of each iterator in the underlying datastructure.
//   /*! returns True when other points to kmer points to a nearer or equal position than the current */
//   bool operator >=(const kDataFrameIterator& other){
//     return *iterator >= *other.iterator;
//   }
  /// Compare the position of each iterator in the underlying datastructure.
  /*! returns True when current and other points to the same kmer */
  bool operator ==(const kDataFrameIterator& other)
  {
    return *iterator == *other.iterator;
  }
  /// Compare the position of each iterator in the underlying datastructure.
  /*! returns True when current and other points to different kmers */
  bool operator !=(const kDataFrameIterator& other)
  {
    return *iterator != *other.iterator;
  }
  /// Returns the hash value of the current kmer
  std::uint64_t getHashedKmer(){
    return iterator->getHashedKmer();
  };
  /// Returns the current kmer order
    uint64_t getOrder(){
        return iterator->getOrder();
    }
  /// Returns the current kmer
  string getKmer(){
    return iterator->getKmer();
  }
  /// Returns the count of the current kmer
  std::uint64_t getCount();
  /// sets the count of the current kmer
  bool setCount(std::uint64_t count);
  
  bool setOrder(std::uint64_t count){
    return iterator->setOrder(count);
  }
  kmerRow getKmerRow(){
    return kmerRow(iterator->getKmer(),
                   iterator->getHashedKmer(),
                   getCount(),
                   iterator->getOrder(),
                   origin
                  );
  }
  template<typename T,typename Container>
  void getColumnValue(const string& colName, T& res);
  template<typename T,typename Container>
  void setColumnValue(const string& colName, T value);


  void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName);

  ~kDataFrameIterator(){
    delete iterator;
  }
};



class dbgIterator{
public:
    kDataFrame* frame;
    string currentKmer;
    vector<string> nextFwdKmers;
    vector<string> nextRevKmers;

    dbgIterator();
    dbgIterator(kDataFrame*,string kmer);

    dbgIterator(const dbgIterator& other);
    dbgIterator& operator= (const dbgIterator& other);

    void nextFWD(uint32_t index);
    void nextREV(uint32_t index);
private:
    void generateNextKmers();

};

class kDataFrame{
protected:
  std::uint64_t kSize;
  string class_name; // Default = MQF, change if MAP. Temporary until resolving #17



  unordered_map<uint64_t,uint32_t> orderCheckpoints;
  uint32_t lastCheckpoint;
  vectorColumn<uint32_t>* countColumn;
  kDataFrameIterator* endIterator;
  friend class kDataFrameIterator;
public:
    uint64_t lastKmerOrder;
    unordered_map<string, Column*> columns;
    typedef kDataFrameIterator iterator;
    typedef kmerRow value_type;
    kmerDecoder * KD;
  virtual string get_class_name(){ return class_name;}  // Temporary until resolving #17
  kDataFrame();
  explicit kDataFrame(uint8_t kSize);


  virtual kDataFrame* clone()=0;
  virtual ~kDataFrame();
/// creates a new kDataframe using the same parameters as the current kDataFrame.
/*! It is like clone but without copying the data */
  virtual kDataFrame* getTwin()=0;

  /// request a capacity change so that the kDataFrame can approximately hold at least n kmers
  void reserve (std::uint64_t n );
/// request a capacity change so that the kDataFrame can approximately hold at least n kmers
  virtual void _reserve (std::uint64_t n )=0;
/// request a capacity change so that the kDataFrame can approximately hold kmers with countHistogram distribution
  virtual void _reserve (vector<std::uint64_t> countHistogram)=0;
/// insert the kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
  virtual uint32_t insert(const string &kmer)=0;
/// insert the hashed kmer one time in the kDataFrame, or increment the kmer count if it is already exists.
/*! Returns bool value indicating whether the kmer is inserted or not*/
  virtual uint32_t insert(std::uint64_t kmer)=0;
  /// insert the kmer in the kmer row time in the kDataFrame, or increment the kmer count with the count in the row if it is already exists.
  /*! Returns bool value indicating whether the kmer is inserted or not*/
  uint32_t insert(kmerRow k);
  kDataFrame::iterator insert(kDataFrame::iterator& it,kmerRow k);
/// set the kmer's count to N time in the kDataFrame
/*! Returns bool value indicating whether the kmer is inserted or not.
The difference between setCount and insert is that setCount set the count to N no matter the previous kmer count was*/
  void addCountColumn();
  virtual bool setCount(const string &kmer, std::uint64_t N);
  virtual bool setCount(std::uint64_t kmer,std::uint64_t N);
  std::uint64_t getCount(const string &kmer);
  std::uint64_t getCount(std::uint64_t kmer);
  void incrementCount(std::uint64_t kmer);
  void incrementCount(const string kmer);
/// returns the count of the kmer in the kDataFrame, i.e. the number of times the kmer is inserted in the kdataFrame.
  virtual std::uint64_t getkmerOrder(const string &kmer)=0;
  virtual std::uint64_t getkmerOrder(std::uint64_t kmer)=0;
  virtual bool setOrder(const string &kmer, std::uint64_t N)=0;
  virtual  bool setOrder(std::uint64_t kmer,std::uint64_t N)=0;
// Removes  a kmer from the kDataFrame
/*! Returns bool value indicating whether the kmer is erased or not*/
  virtual bool erase(const string &kmer)=0;
  virtual bool erase(std::uint64_t kmer)=0;

/// Returns the number of kmers in the kDataframe.
  virtual std::uint64_t size()=0;
/// Returns the maximum number of kmers that the kDataframe can hold.
  virtual std::uint64_t max_size()=0;
/// Test whether the kDataFrame is empty.
/*! Returns a bool value indicating whether the kDataFrame is empty, i.e. whether its size is 0.*/
  bool empty();
/// Returns the current load factor in the kDataFrame.
  virtual float load_factor()=0;
/// Returns the current maximum load factor for the kDataFrame.
  virtual float max_load_factor()=0;

///Returns an iterator at the beginning of the kDataFrame.
  virtual kDataFrameIterator begin()=0;
///Returns an iterator at the end of the kDataFrame.
  virtual kDataFrameIterator end();
///Returns an iterator at the specific kmer.
  virtual kDataFrameIterator find(const string &kmer)=0;
  virtual kDataFrameIterator find(uint64_t kmer)=0;


  dbgIterator getDBGIterator(const string &kmer);
  virtual void serialize(string filePath)=0;
/// Returns the kmerDecoder used by kDataframe.
  kmerDecoder* getkmerDecoder() const{
    return KD;
  };



  static kDataFrame* load(string filePath);
  void save(string filePath);

  vector<string> getColumnNames();
  std::uint64_t getkSize() const{return kSize;}

  // duplicate for easier name in python, getkSize won't be wrapped
  std::uint64_t ksize() const{return kSize;}

  void setkSize(std::uint64_t k){
    kSize=k;
    if(KD!= nullptr)
    {
        kmerDecoder* newKD=new Kmers(kSize,KD->hash_mode);
        delete KD;
        KD=newKD;
    }else{
        KD= new Kmers(kSize);
    }
  }


  void addColumn(string columnName, Column*);
  void removeColumn(string columnName);
  template<typename T,typename Container>
  T getKmerColumnValue(const string& columnName,string kmer);

  template<typename T,typename Container>
  void setKmerColumnValue(const string& columnName,string kmer, T value);


  template<typename T,typename Container>
  T getKmerColumnValue(const string& columnName,uint64_t kmer);

  template<typename T,typename Container>
  void setKmerColumnValue(const string& columnName,uint64_t kmer, T value);


  template<typename T,typename Container>
  T getKmerColumnValueByOrder(const string& columnName,uint64_t kmerOrder);

  template<typename T,typename Container>
  void setKmerColumnValueByOrder(const string& columnName,uint64_t kmerOrder, T value);




  virtual bool kmerExist(string kmer)=0;
  virtual bool kmerExist(uint64_t kmer)=0;

  void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName, std::uint64_t kmer);
  void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName, string kmer);

};

template<typename T,typename Container>
T kDataFrame::getKmerColumnValue(const string& columnName,string kmer)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    if(kmerOrder==0)
        throw std::logic_error("kmer not found!");
    return ((Container*)columns[columnName])->get(kmerOrder);
}
template<typename T,typename Container>
void kDataFrame::setKmerColumnValue(const string& columnName,string kmer,T value)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    if(kmerOrder==0)
        throw std::logic_error("kmer not found!");
    ((Container*)columns[columnName])->insert(value,kmerOrder);
}


template<typename T,typename Container>
T kDataFrame::getKmerColumnValue(const string& columnName,uint64_t kmer)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    if(kmerOrder==0)
        throw std::logic_error("kmer not found!");
    return ((Container*)columns[columnName])->get(kmerOrder);
}
template<typename T,typename Container>
void kDataFrame::setKmerColumnValue(const string& columnName,uint64_t kmer,T value)
{
    std::uint64_t kmerOrder=getkmerOrder(kmer);
    if(kmerOrder==0)
        throw std::logic_error("kmer not found!");
    ((Container*)columns[columnName])->insert(value,kmerOrder);
}

template<typename T,typename Container>
T kDataFrame::getKmerColumnValueByOrder(const string& columnName,uint64_t kmerOrder)
{
    return ((Container*)columns[columnName])->get(kmerOrder);
}
template<typename T,typename Container>
void kDataFrame::setKmerColumnValueByOrder(const string& columnName,uint64_t kmerOrder,T value)
{
    ((Container*)columns[columnName])->insert(value,kmerOrder);
}



template<typename T,typename Container>
void kmerRow::getColumnValue(const string& colName, T& res)
{
    res= origin->getKmerColumnValueByOrder<T,Container>(colName, order);
}

template<typename T,typename Container>
void kmerRow::setColumnValue(const string& colName, T value)
{
    origin->setKmerColumnValueByOrder<T,Container>(colName, order,value);
}


template<typename T,typename Container>
void kDataFrameIterator::getColumnValue(const string& colName, T& res)
{
    res= origin->getKmerColumnValueByOrder<T,Container>(colName, iterator->getOrder());
}

template<typename T,typename Container>
void kDataFrameIterator::setColumnValue(const string& colName, T value)
{
    origin->setKmerColumnValueByOrder<T,Container>(colName, iterator->getOrder(),value);
}







template <class MapType>
        class kDataFrameSTL;

template <class MapType>
class kDataFrameSTLIterator : public _kDataFrameIterator {
private:
    kDataFrameSTL<MapType> *origin;
    kmerDecoder * KD;
public:
    typedef typename MapType::iterator itType;
    itType iterator;
    kDataFrameSTLIterator(itType, kDataFrameSTL<MapType> *origin, std::uint64_t kSize);

    kDataFrameSTLIterator(const kDataFrameSTLIterator<MapType> &);

    kDataFrameSTLIterator &operator++(int);

    _kDataFrameIterator *clone();

    std::uint64_t getHashedKmer();

    string getKmer();
    uint64_t getOrder();
    std::uint64_t getCount();

    bool setOrder(std::uint64_t count);

    void endIterator();

    bool operator==(const _kDataFrameIterator &other);

    bool operator!=(const _kDataFrameIterator &other);

    ~kDataFrameSTLIterator();
};




// kDataFrameSTL _____________________________


template <class MapType>
class kDataFrameSTL : public kDataFrame {
public:
    typedef  MapType  MAPType;
private:
    MapType MAP;
public:
    kDataFrameSTL();

    explicit kDataFrameSTL(uint64_t ksize);
    kDataFrameSTL(std::uint64_t kSize,uint64_t nKmers);
    kDataFrameSTL(readingModes RM, hashingModes hash_mode, map<string, int> params);
    kDataFrameSTL(uint64_t ksize, hashingModes hash_mode);
    kDataFrameSTL(uint64_t kSize,vector<uint64_t> kmersHistogram);


    kDataFrame *clone() override;


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

    kDataFrameIterator end();
    kDataFrameIterator find(const string &kmer);
    kDataFrameIterator find(uint64_t kmer);

    std::uint64_t bucket(string kmer);

    void serialize(string filePath);

    static kDataFrame *load(string filePath);

    ~kDataFrameSTL() {
        this->MAP.clear();
    }
    MapType *getMap();


};

template <>
void kDataFramePHMAP::serialize(string filePath);
template <>
kDataFrame * kDataFramePHMAP::load(string filePath);

template <>
void kDataFrameMAP::_reserve(std::uint64_t n);
template<>
float kDataFrameMAP::load_factor();
template<>
float kDataFrameMAP::max_load_factor();
template <>
void kDataFrameMAP::serialize(string filePath);
template <>
kDataFrame * kDataFrameMAP::load(string filePath);


template <>
void kDataFrameBtree::_reserve(std::uint64_t n);
template<>
float kDataFrameBtree::load_factor();
template<>
float kDataFrameBtree::max_load_factor();
template <>
void kDataFrameBtree::serialize(string filePath);
template <>
kDataFrame * kDataFrameBtree::load(string filePath);





class kDataFrame_sshash;


class kDataFrameFactory{
public:
    static kDataFrame* loadMQF(string filePath);
    static kDataFrame* createMQF(uint32_t kSize,uint32_t numKmers=10000);
    static kDataFrame* createMQF(kDataFrame* kframe);

    static kDataFrame* loadPHMAP(string filePath);
    static kDataFrame* createPHMAP(uint32_t kSize,uint32_t numKmers=10000);
    static kDataFrame* createPHMAP(uint64_t ksize, hashingModes hash_mode);


    static kDataFrame* loadMAP(string filePath);
    static kDataFrame* createMAP(uint32_t kSize,uint32_t numKmers=10000);

    static kDataFrame* loadBtree(string filePath);
    static kDataFrame* createBtree(uint32_t kSize,uint32_t numKmers=10000);

    static kDataFrame* loadBMQF(string filePath);
    static kDataFrame* createBMQF(uint32_t kSize,string filePath,uint32_t numKmers=10000);
    static kDataFrame* createBMQF(kDataFrame* kframe,string filePath);


    static kDataFrame* loadBlight(string filePath);
    static kDataFrame* createBlight(uint32_t kSize,string filePath);

    static kDataFrame* loadSSHASH(string filePath);
    static kDataFrame* createSSHASH(uint32_t kSize,string filePath);

};





#endif
