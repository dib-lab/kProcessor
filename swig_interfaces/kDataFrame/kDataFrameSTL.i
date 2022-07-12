#define kDataFramePHMAP kDataFrameSTL<phmap::parallel_flat_hash_map<std::uint64_t,std::uint32_t,std::hash<uint64_t>,std::equal_to<uint64_t>,std::allocator<std::pair<const uint64_t, uint32_t>>,4,std::mutex>>
#define kDataFramePHMAPIterator kDataFrameSTLIterator<phmap::parallel_flat_hash_map<std::uint64_t,std::uint32_t,std::hash<uint64_t>,std::equal_to<uint64_t>,std::allocator<std::pair<const uint64_t, uint32_t>>,4,std::mutex>>
#define kDataFrameMAP kDataFrameSTL<std::map<std::uint64_t, std::uint32_t>>
#define kDataFrameMAPIterator kDataFrameSTLIterator<std::map<std::uint64_t, std::uint32_t>>
#define kDataFrameBtree kDataFrameSTL<phmap::btree_map<uint64_t,uint32_t>>
#define kDataFrameBtreeIterator kDataFrameSTLIterator<phmap::btree_map<uint64_t,uint32_t>>


template <class MapType>
class kDataFrameSTL;

template <class MapType>



// kDataFrameSTL _____________________________


template <class MapType>
class kDataFrameSTL : public kDataFrame {
public:
        typedef  MapType  MAPType;
        kDataFrameSTL();

        explicit kDataFrameSTL(uint64_t ksize);
        kDataFrameSTL(std::uint64_t kSize, uint64_t nKmers);
        kDataFrameSTL(readingModes RM, hashingModes hash_mode, map<string, int> params);
        kDataFrameSTL(uint64_t ksize, hashingModes hash_mode);
        kDataFrameSTL(uint64_t kSize, vector<uint64_t> kmersHistogram);


        kDataFrame* clone() override;


        kDataFrame* getTwin();

        void _reserve(std::uint64_t n);
        void _reserve(vector<std::uint64_t> countHistogram);

        bool kmerExist(string kmer);
        bool kmerExist(uint64_t kmer);

        bool setOrder(const string& kmer, std::uint64_t count);
        bool setOrder(std::uint64_t kmer, std::uint64_t count);

        uint32_t insert(const string& kmer) override;
        uint32_t insert(std::uint64_t kmer) override;

        std::uint64_t getkmerOrder(const string& kmer);
        std::uint64_t getkmerOrder(std::uint64_t kmer);

        bool erase(const string& kmer);
        bool erase(std::uint64_t kmer);

        std::uint64_t size();

        std::uint64_t max_size();

        float load_factor();

        float max_load_factor();

        kDataFrameIterator begin();

        kDataFrameIterator end();
        kDataFrameIterator find(const string& kmer);
        kDataFrameIterator find(uint64_t kmer);

        std::uint64_t bucket(string kmer);

        void serialize(string filePath);

        static kDataFrame* load(string filePath);

        ~kDataFrameSTL() {
                this->MAP.clear();
        }
        MapType* getMap();


};

template <>
void kDataFramePHMAP::serialize(string filePath);
template <>
kDataFrame* kDataFramePHMAP::load(string filePath);

template <>
void kDataFrameMAP::_reserve(std::uint64_t n);
template<>
float kDataFrameMAP::load_factor();
template<>
float kDataFrameMAP::max_load_factor();
template <>
void kDataFrameMAP::serialize(string filePath);
template <>
kDataFrame* kDataFrameMAP::load(string filePath);


template <>
void kDataFrameBtree::_reserve(std::uint64_t n);
template<>
float kDataFrameBtree::load_factor();
template<>
float kDataFrameBtree::max_load_factor();
template <>
void kDataFrameBtree::serialize(string filePath);
template <>
kDataFrame* kDataFrameBtree::load(string filePath);