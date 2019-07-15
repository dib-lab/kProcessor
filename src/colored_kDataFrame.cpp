#include "colored_kDataFrame.hpp"

using namespace std;


colored_kDataFrame::colored_kDataFrame()
{
  frame=new kDataFrameMQF();
  colors=new BitVectorsTable();
  nextAvailableColor=1;
}

void colored_kDataFrame::addNewColor(uint32_t color, vector<uint32_t> & samplesIds)
{
  colors->setColor(color,samplesIds);
}
void colored_kDataFrame::setKmerColor(string kmer,uint32_t color)
{
  frame->setCount(kmer,color);
}
uint32_t colored_kDataFrame::getKmerColor(string kmer)
{
  return frame->count(kmer);
}
void colored_kDataFrame::getSamplesIDForKmer(string kmer,vector<uint32_t>& result)
{
  uint32_t color=getKmerColor(kmer);
  result.clear();
  if(color==0)
  {
    return;
  }
  colors->getSamples(color,result);
}
void colored_kDataFrame::getSamplesIDForColor(uint32_t color,vector<uint32_t>& result)
{
  colors->getSamples(color,result);
}
void colored_kDataFrame::colorKmer(string kmer,vector<uint32_t> & samplesIds){
  uint64_t color=colorsInv->getColorId(samplesIds);
  if(color==0)
  {
    color=nextAvailableColor++;
    colorsInv->setColorId(color,samplesIds);
    colors->setColor(color,samplesIds);
  }
  setKmerColor(kmer,color);
}
void colored_kDataFrame::setColorTable(colorTable* table){
  colors=table;
}
void colored_kDataFrame::setkDataFrame(kDataFrame* f)
{
  frame=f;
}
uint64_t colored_kDataFrame::getkSize(){
  return frame->getkSize();
}
void colored_kDataFrame::save(string prefix)
{
  frame->save(prefix);
  colors->save(prefix);
  ofstream namesMapOut(prefix+".namesMap");
  namesMapOut<<namesMap.size()<<endl;
  for(auto it:namesMap)
  {
    namesMapOut<<it.first<<" "<<it.second<<endl;
  }
  namesMapOut.close();
}
colored_kDataFrame* colored_kDataFrame::load(string prefix)
{
  colored_kDataFrame* res=new colored_kDataFrame();
  res->frame=kDataFrame::load(prefix);
  res->colors=colorTable::load(prefix);
  ifstream namesMapIn(prefix+".namesMap");
  uint64_t size;
  namesMapIn>>size;
  for(int i=0;i<size;i++)
  {
    uint32_t color;
    string name;
    namesMapIn>>color>>name;
    res->namesMap[color]=name;
    res->namesMapInv[name]=color;
  }
  return res;
}

unordered_map<int, string> colored_kDataFrame::names_map(){
    unordered_map<int, string> names_map;
    for(const auto &item : this->namesMap){
        names_map[item.first] = item.second;
    }
    return names_map;
}

unordered_map<string, int> colored_kDataFrame::inverse_names_map(){
    unordered_map<string, int> inv_names_map;
    for(const auto &item : this->namesMapInv){
        inv_names_map[item.first] = item.second;
    }
    return inv_names_map;
}