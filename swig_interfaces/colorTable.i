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


class intVectorsTable: public colorTable{
    public:
        intVectorsTable(){}
        intVectorsTable(string folderName);
        //BitVectorsTable(vector<string> fileNames,uint64_t numSamples);
        virtual ~intVectorsTable();
        bool getSamples(uint64_t colorID,vector<uint32_t>& res)override;
        bool setColor(uint64_t colorID,vector<uint32_t>& v)override;
        void save(string folderName)override;
        private:
        flat_hash_map<uint64_t, std::vector<uint32_t> > colors;
};