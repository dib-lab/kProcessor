#ifndef COLORTABLE_H__
#define COLORTABLE_H__
#include <vector>
#include "sdsl/vectors.hpp"
#include <memory>

using namespace std;
class colorTable{
public:
  colorTable();
  virtual ~colorTable();
  static colorTable* load(string folderName);
  virtual bool getSamples(uint64_t colorID,vector<uint32_t>& res)=0;
  virtual bool setColor(uint64_t colorID,vector<uint32_t>& v)=0;
  virtual void save(string folderName)=0;
  uint64_t numSamples;
  uint32_t numColors;

};

class colorTableInv{
public:
  colorTableInv();
  virtual ~colorTableInv();
  virtual uint64_t getColorId(vector<uint32_t>& )=0;
  virtual void setColorId(uint64_t colorID,vector<uint32_t>& v)=0;
};

class stringColorTableInv: public colorTableInv{
private:
  unordered_map<string,uint64_t> table;
  string getKey(vector<uint32_t>&);
public:
  stringColorTableInv();
  ~stringColorTableInv();
  uint64_t getColorId(vector<uint32_t>& );
  void setColorId(uint64_t colorID,vector<uint32_t>& v);
};

// this class is copied and edited from mantis project
//https://github.com/splatlab/mantis
class BitVectorsTable: public colorTable{
public:
  static const uint64_t NUM_BV_BUFFER=20000000;
  BitVectorsTable(){}
  BitVectorsTable(uint64_t numSamples);
  BitVectorsTable(string folderName);
  //BitVectorsTable(vector<string> fileNames,uint64_t numSamples);
  virtual ~BitVectorsTable();
  bool getSamples(uint64_t colorID,vector<uint32_t>& res)override;
  bool setColor(uint64_t colorID,vector<uint32_t>& v)override;
  void save(string folderName)override;
private:
  std::vector<sdsl::rrr_vector<63> > eqclasses;
  sdsl::bit_vector current;
  uint64_t nCurrentColors;
};

// class samplesCombination{
// public:
//   samplesCombination(){}
//   ~samplesCombination(){}
//   bool serialize(ofstream& out);
//   //virtual bool deserialize(ifstream& input)=0;
//   bool getSamples(uint32_t numSamples,vector<uint32_t>& res);
//   uint64_t sizeInBytes();
//   static unique_ptr<samplesCombination> load(ifstream& input);
// };
// class sdslEncVector: samplesCombination{
// public:
//   sdslEncVector(vector<uint32_t> samples){};
//   ~sdslEncVector(){};
//   bool serialize(ofstream& out);
//   bool deserialize(ifstream& out);
//   bool getSamples(uint32_t numSamples,vector<uint32_t>& res);
//   uint64_t sizeInBytes();
// private:
//   sdsl:enc_vector<> data;
// }
// class sdslRRRBitVector: samplesCombination{
// public:
//   sdslRRRBitVector(vector<uint32_t> samples,uint32_t numSamples);
//   ~sdslRRRBitVector();
//   bool serialize(ofstream& out);
//   bool deserialize(ifstream& out);
//   bool getSamples(uint32_t numSamples,vector<uint32_t>& res);
//   uint64_t sizeInBytes();
// private:
//   sdsl::rrr_vector<63> rrr_bits;
// }
//
// class pointer: public samplesCombination{
// public:
//   pointer(uint32_t value): value(value){};
//   pointer(ifstream& out);
//   ~pointer(){};
//   bool serialize(ofstream& out);
//   bool getSamples(uint32_t numSamples,vector<uint32_t>& res);
//   uint64_t sizeInBytes();
//   uint32_t value;
//
// };

// template<int size=8>
// class samplesList: samplesCombination{
// public:
//   samplesList(vector<uint32_t> samples);
//   ~samplesList();
//   bool serialize(ofstream& out);
//   bool deserialize(ifstream& out);
//   bool getSamples(uint32_t numSamples,vector<uint32_t>& res);
//   uint64_t sizeInBytes();
// private:
//   uint8_t data[size];
// }
//
// class HybridTable: public colorTable{
// public:
//   HybridTable(uint64_t numSamples,uint64_t numColors);
//   HybridTable(string folderName);
//   ~HybridTable(){
//
//   }
//   bool getSamples(uint64_t colorID,uint8_t colorType,vector<uint32_t>& res);
//   bool setColor(uint64_t colorID,uint8_t colorType,vector<uint32_t>& v);
//   void save(string folderName)override;
// private:
//   vector<unique_ptr<samplesCombination> > combinations;
//   BitVectorsTable bitvectors;
//   uint64_t lastBVpointer;
// };
//
//

#endif
