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


class intVectorsTable: public colorTable{
public:
  intVectorsTable(){}  
  intVectorsTable(string folderName);
  //BitVectorsTable(vector<string> fileNames,uint64_t numSamples);
  virtual ~intVectorsTable();
  bool getSamples(uint64_t colorID,vector<uint32_t>& res)override;
  bool setColor(uint64_t colorID,vector<uint32_t>& v)override;
  void save(string folderName)override;
};
