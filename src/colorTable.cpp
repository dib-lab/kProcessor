#include "colorTable.hpp"
#include <fstream>
#include <iostream>
#include "json.hpp"
#include <memory>
#include "Utils/utils.hpp"
using nlohmann::json;
using namespace std;

colorTable* colorTable::load(string prefix){
  string configFilename=prefix+"colorTableConfig.json";
  ifstream configStream(configFilename);
  json config;
  configStream>>config;
  configStream.close();
  string type=config["Type"];
  if(type=="BitVectorsTable")
    return new BitVectorsTable(prefix);
  throw std::logic_error("Unknown color table file type ="+ type);
}
BitVectorsTable::BitVectorsTable(uint64_t numSamples)
{
  current=sdsl::bit_vector(numSamples*NUM_BV_BUFFER);
  nCurrentColors=0;
  numColors=0;
  this->numSamples=numSamples;
}
BitVectorsTable::BitVectorsTable(string prefix){
  string configFilename=prefix+"colorTableConfig.json";
  ifstream configStream(configFilename);
  json config;
  configStream>>config;
  configStream.close();
  numSamples=config["numSamples"];
  int counter=0;
  string fileName=prefix +"."+ std::to_string(counter) + ".cls.eqclass_rrr";
  std::map<int, std::string> sorted_files;
  while(kProcessor::utils::FileExists(fileName))
  {

    sorted_files[counter] = fileName;
    counter++;
    fileName=prefix +"."+ std::to_string(counter) + ".cls.eqclass_rrr";
  }
  eqclasses.reserve(sorted_files.size());
  uint32_t num_serializations=0;
  for (auto file : sorted_files) {
    sdsl::rrr_vector<63> rrr_bits;
    sdsl::load_from_file(rrr_bits, file.second.c_str());
    eqclasses.push_back(rrr_bits);
    num_serializations++;
  }
  uint64_t total = 0;
  for (uint32_t i = 0; i < num_serializations; i++)
    total += eqclasses[i].size();


  numColors=total/numSamples;
  if(numColors!=config["numColors"])
  {
    cerr<<"Bitvector color table is not loaded correctly. Loaded Num Colors is "
    <<numColors<<" and it is expected to be "<<config["numColors"]<<endl;
  }
}

// BitVectorsTable::BitVectorsTable(vector<string> eqclass_files,uint64_t num_Samples): colorTable(){
//   numSamples=num_Samples;
//
//   std::map<int, std::string> sorted_files;
//   for (std::string file : eqclass_files) {
//     int id = std::stoi(first_part(last_part(file, '/'), '_'));
//     sorted_files[id] = file;
//   }
//   eqclasses.reserve(sorted_files.size());
//   uint32_t num_serializations=0;
//   for (auto file : sorted_files) {
//     eqclasses.push_back(BitVectorRRR(file.second));
//     num_serializations++;
//   }
//
//   uint64_t total = 0;
//   for (uint32_t i = 0; i < num_serializations; i++)
//     total += eqclasses[i].bit_size();
//   numColors=total/numSamples;
// }
colorTable::colorTable(){
  numSamples=0;
  numColors=0;
}
colorTable::~colorTable(){

}
BitVectorsTable::~BitVectorsTable(){
  eqclasses.clear();
}
bool BitVectorsTable::getSamples(uint64_t colorID,vector<uint32_t>& res)
{
  if(nCurrentColors!=0)
  {
    if(nCurrentColors<NUM_BV_BUFFER)
    {
      current.resize(nCurrentColors*numSamples);
    }
    eqclasses.push_back(sdsl::rrr_vector<63>(current));
    current=sdsl::bit_vector(current.size());
    nCurrentColors=0;
  }
  res.clear();
  uint32_t sum=0;
  uint64_t start_idx = (colorID - 1);
  uint64_t bucket_idx = start_idx / NUM_BV_BUFFER;
  uint64_t bucket_offset = (start_idx % NUM_BV_BUFFER) * numSamples;
  uint32_t top=0;
  for (uint32_t w = 0; w <= numSamples / 64; w++) {
    // if(sum>classesThreshold)
    //   break;
    uint64_t len = std::min((uint64_t)64, numSamples - w * 64);
    uint64_t wrd = eqclasses.at(bucket_idx).get_int(bucket_offset, len);
    for (uint32_t i = 0, sCntr = w * 64; i < len; i++, sCntr++)
        if ((wrd >> i) & 0x01){
            sum += 1;
            res.push_back(sCntr);
          }
    bucket_offset += len;
  }

  return true;

}
bool BitVectorsTable::setColor(uint64_t colorID,vector<uint32_t>& samples){
  if(colorID!=numColors+1)
  {
    cerr<<"Bitvectors supports only appendig technique when adding colors"<<endl;
    return false;
  }
  uint i=nCurrentColors*numSamples;
  for(auto s:samples)
  {
    current[s+i]=1;
  }
  nCurrentColors++;
  if(nCurrentColors==NUM_BV_BUFFER){
    eqclasses.push_back(sdsl::rrr_vector<63>(current));
    current=sdsl::bit_vector(current.size());
  }
  numColors++;
  return true;
}
void BitVectorsTable::save(string prefix)
{
  if(nCurrentColors!=0)
  {
    if(nCurrentColors<NUM_BV_BUFFER)
    {
      current.resize(nCurrentColors*numSamples);
    }
    eqclasses.push_back(sdsl::rrr_vector<63>(current));
    current=sdsl::bit_vector(current.size());
    nCurrentColors=0;
  }
  json config={
    {"Type","BitVectorsTable"},
    {"numSamples" ,numSamples},
    {"numColors", numColors}
  };
  string configFilename=prefix+"colorTableConfig.json";
  ofstream configStream(configFilename);
  configStream<<config;
  configStream.close();
  int currentBuffer=0;
  for(auto color: eqclasses)
  {
    std::string bv_file(prefix +"."+ std::to_string(currentBuffer++)
    + ".cls.eqclass_rrr");
    sdsl::store_to_file(color, bv_file.c_str());
  }
}

colorTableInv::colorTableInv(){

}
colorTableInv::~colorTableInv(){

}

stringColorTableInv::stringColorTableInv(){

}
stringColorTableInv::~stringColorTableInv(){

}
uint64_t stringColorTableInv::getColorId(vector<uint32_t>& combination)
{

  return table[getKey(combination)];

}
string stringColorTableInv::getKey(vector<uint32_t>& combination)
{
  string key="";
  for(auto i:combination)
  {
    key+=std::to_string(i)+";";
  }
  return key;
}


void stringColorTableInv::setColorId(uint64_t colorID,vector<uint32_t>& combination)
{  
  table[getKey(combination)]=colorID;
}

//
// static unique_ptr<samplesCombination> samplesCombination::load(ifstream& input){
//   uint8_t code;
//   input.read((char*)(&code),sizeof(code));
//   if(code==0)
//   {
//     return make_unique<pointer>(input);
//   }else{
//     throw runtime_error("Unknown combination type");
//   }
// }
//
// bool pointer::serialize(ofstream& out){
//   uint8_t code=0;
//   out.write((char*)(&code),sizeof(code));
//   out.write((char*)(&value),sizeof(value));
// }
// pointer::pointer(ifstream& input){
//   input.read((char*)(&value),sizeof(value));
// }
// bool pointer::getSamples(uint32_t numSamples,vector<uint32_t>& res){
//   throw runtime_error("get samples shouldnt be called on a pointer object");
//   return false;
// }
// uint64_t pointer::sizeInBytes(){
//   return sizeof(value);
// }
//
//
// HybridTable::HybridTable(uint64_t numSamples,uint64_t numColors){
//   bitvectors=BitVectorsTable(numSamples);
//   combinations=vector<unique_ptr<samplesCombination> >(numColors);
//   lastBVpointer=0;
// }
// HybridTable::HybridTable(string prefix){
//   string configFilename=prefix+"colorTableConfig.json";
//   ifstream configStream(configFilename);
//   json config;
//   configStream>>config;
//   configStream.close();
//   numSamples=config["numSamples"];
//   numColors=config["numColors"];
//   string bvFileName=config["BitVectorTable"];
//   bitvectors=BitVectorsTable(bvFileName);
//   combinations.reserve(numColors);
//   string hybrid=config["hybridTable"];
//   ifstream input(hybrid,ios::in | ios::binary);
//   for(int i=0;i<numColors;i++)
//   {
//     combinations.push_back(samplesCombination::load(input));
//   }
// }
// bool HybridTable::getSamples(uint64_t colorID,uint8_t colorType,vector<uint32_t>& res){
//   if(colorType==0)
//   {
//     return bitvectors.getSamples(((pointer*)&combinations[colorID])->value,res);
//   }
//   else{
//     return combinations[colorID]->getSamples(numSamples,res);
//   }
// }
// bool HybridTable::setColor(uint64_t colorID,uint8_t colorType,vector<uint32_t>& samples){
//   if(true || colorType==0)
//   {
//     combinations[colorID]=unique_ptr<samplesCombination>(new pointer(lastBVpointer));
//     bitvectors.setColor(lastBVpointer,samples);
//     lastBVpointer++;
//   }
//   else{
//   //  combinations[colorID]=sdslEncVector(samples);
//   }
// }
// void HybridTable::save(string prefix){
//   string bvFileName=prefix+".hybrid.bv.";
//   string hybridTableFileName=prefix+"eqclass.hybrid";
//   bitvectors.save(bvFileName);
//   json config={
//     {"Type","HybridTable"},
//     {"numSamples" ,numSamples},
//     {"numColors", numColors},
//     {"BitVectorTable",bvFileName},
//     {"hybridTable",hybridTableFileName}
//   };
//   string configFilename=prefix+"colorTableConfig.json";
//   ofstream configStream(configFilename);
//   configStream<<config;
//   configStream.close();
//   ofstream out(hybridTableFileName, ios::out | ios::binary);
//   for(int i=0;i<combinations.size();i++){
//     combinations[i]->serialize(out);
//   }
//   out.close();
// }
