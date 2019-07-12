#ifndef _coloredkDataFRAME_H_
#define _coloredkDataFRAME_H_

#include "kDataFrame.hpp"
#include "colorTable.hpp"
#include <parallel_hashmap/phmap.h>
using namespace std;
using phmap::flat_hash_map;

class colored_kDataFrame{
private:
  kDataFrame* frame;
  colorTable* colors;
  colorTableInv* colorsInv;
  uint64_t nextAvailableColor;
public:
  flat_hash_map<uint32_t,string> namesMap;
  flat_hash_map<string,uint32_t> namesMapInv;
  colored_kDataFrame();
  void addNewColor(uint32_t color, vector<uint32_t> & samplesIds);
  void setKmerColor(string kmer,uint32_t color);
  uint32_t getKmerColor(string kmer);
  void getSamplesIDForKmer(string kmer,vector<uint32_t>& result);
  void getSamplesIDForColor(uint32_t color,vector<uint32_t>& result);

  void colorKmer(string kmer,vector<uint32_t> & samplesIds);

  void setColorTable(colorTable* table);
  void setkDataFrame(kDataFrame* f);
  void save(string prefix);
  static colored_kDataFrame* load(string prefix);
  uint64_t getkSize();


};

#endif
