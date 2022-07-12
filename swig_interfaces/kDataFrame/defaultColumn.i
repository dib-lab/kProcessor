class Column {
public:
    Column() {}
    virtual ~Column() {}

    virtual Column* getTwin() = 0;
    virtual void resize(uint32_t size) = 0;
    virtual uint32_t size() = 0;
    static Column* getContainerByName(size_t name);

    virtual void serialize(string filename) = 0;
    virtual void deserialize(string filename) = 0;
    virtual Column* clone() = 0;

    virtual void setValueFromColumn(Column* Container, uint32_t inputOrder, uint32_t outputOrder) {

    }


};



template<typename ColumnType, typename indexType = vector<uint32_t> >
class deduplicatedColumn : public Column {
public:
    typedef typename ColumnType::dataType dataType;
    indexType index;
    ColumnType* values;
    deduplicatedColumn() {
    }
    deduplicatedColumn(uint32_t size);
    ~deduplicatedColumn() {
        delete values;
    }

    uint32_t  size()
    {
        return index.size();
    }
    void insert(dataType item, uint32_t index);
    dataType get(uint32_t index);

    void serialize(string filename);
    void deserialize(string filename);
    Column* getTwin();
    void resize(uint32_t size);
    void setValueFromColumn(Column* Container, uint32_t inputOrder, uint32_t outputOrder);

    Column* clone() override;


};

class inExactColorIndex {
public:
    flat_hash_map<uint64_t, uint32_t> colors;
    uint32_t lastColor;
    uint32_t noSamples;
    inExactColorIndex()
    {

        lastColor = 0;
        noSamples = 0;
    }
    ~inExactColorIndex() {

    }
    bool hasColorID(vector<uint32_t>& v);
    uint32_t getColorID(vector<uint32_t>& v);
    void serialize(string fileName) {}
    void deserialize(string filename) {}

    void populateColors(vector<vector<uint32_t > >& colors) {}
    void optimize() {}

};


template<typename  T>
class vectorColumn : public Column {
public:
    typedef T dataType;
    vector<T> dataV;
    vectorColumn() {
        dataV = vector<T>(1000);
    }
    vectorColumn(uint32_t size) {
        dataV.resize(size + 1);
    }
    ~vectorColumn() {

    }

    uint32_t  insertAndGetIndex(T item);
    T getWithIndex(uint32_t index);
    uint32_t size() {
        return dataV.size();
    }
    void resize(uint32_t size)
    {
        dataV.resize(size + 1);
    }
    void insert(T item, uint32_t index);
    T get(uint32_t index);

    void serialize(string filename);
    void deserialize(string filename);
    Column* getTwin();
    void setValueFromColumn(Column* Container, uint32_t inputOrder, uint32_t outputOrder);

    Column* clone() override;

};




class inExactColorIndex {
public:
    flat_hash_map<uint64_t, uint32_t> colors;
    uint32_t lastColor;
    uint32_t noSamples;
    inExactColorIndex()
    {

        lastColor = 0;
        noSamples = 0;
    }
    ~inExactColorIndex() {

    }
    bool hasColorID(vector<uint32_t>& v);
    uint32_t getColorID(vector<uint32_t>& v);
    void serialize(string fileName) {}
    void deserialize(string filename) {}

    void populateColors(vector<vector<uint32_t > >& colors) {}
    void optimize() {}

};



class queryColorColumn : public Column {
public:
    typedef vector<uint32_t> dataType;
    uint64_t  noSamples;

    uint64_t numColors;
    queryColorColumn() {}
    ~queryColorColumn() override = default;
    uint32_t  insertAndGetIndex(vector<uint32_t >& item);
    virtual vector<uint32_t > getWithIndex(uint32_t index) = 0;
    uint32_t getNumColors()
    {
        return numColors;
    }
    virtual void insert(vector<uint32_t >& item, uint32_t index) = 0;
    virtual vector<uint32_t > get(uint32_t index) = 0;

    virtual void serialize(string filename) = 0;
    virtual void deserialize(string filename) = 0;

    virtual uint32_t size() = 0;
    virtual uint64_t sizeInBytes() = 0;
    virtual void explainSize() = 0;

    virtual Column* getTwin() = 0;

    virtual void resize(uint32_t size) = 0;

};


class StringColorColumn;
class mixVectors : public queryColorColumn {
public:
    typedef vector<uint32_t> dataType;
    deque<vectorBase*> colors;
    sdsl::int_vector<> idsMap;
    mixVectors() {
        colors.push_back(new vectorOfVectors());
    }
    mixVectors(uint64_t noSamples) {
        this->noSamples = noSamples;
        colors.push_back(new vectorOfVectors());
    }
    mixVectors(insertColorColumn* col);

    mixVectors(flat_hash_map<uint64_t, std::vector<uint32_t>>& colors, uint32_t noSamples, uint32_t num_vectors = 20, uint32_t vector_size = 1000000);

    ~mixVectors() {
        for (auto v : colors)
            delete v;
    }

    uint32_t  insertAndGetIndex(vector<uint32_t >& item);
    vector<uint32_t > getWithIndex(uint32_t index);

    void insert(vector<uint32_t >& item, uint32_t index);
    vector<uint32_t > get(uint32_t index);

    void serialize(string filename);
    void deserialize(string filename);
    void sortColors();

    //    void optimize2();
    //    void optimize3(insertColorColumn* col);

    uint32_t size() {
        return numColors;
        uint32_t res = 0;
        for (unsigned int i = 1;i < colors.size();i++)
            res += colors[i]->size();
        return res;
    }
    uint32_t numIntegers() {
        uint32_t res = 0;
        for (unsigned int i = 1;i < colors.size();i++)
            res += colors[i]->numIntegers();
        return res;
    }
    uint64_t sizeInBytes();
    void explainSize();

    Column* getTwin();

    void resize(uint32_t size);

    Column* clone() override;


};

