class kDataFrameIterator {

public:
    _kDataFrameIterator* iterator;
    using iterator_category = std::forward_iterator_tag;
    kDataFrameIterator() {
        iterator = nullptr;
    }
    kDataFrameIterator(const kDataFrameIterator& other) {
        if (other.iterator != nullptr) {
            iterator = other.iterator->clone();
            origin = other.origin;
        }
        else {
            iterator = nullptr;
            origin = nullptr;
        }
    }
    
    /// Returns the hash value of the current kmer
    std::uint64_t getHashedKmer() {
        return iterator->getHashedKmer();
    };
    /// Returns the current kmer order
    uint64_t getOrder() {
        return iterator->getOrder();
    }
    /// Returns the current kmer
    string getKmer() {
        return iterator->getKmer();
    }
    /// Returns the count of the current kmer
    std::uint64_t getCount();
    /// sets the count of the current kmer
    bool setCount(std::uint64_t count);

    bool setOrder(std::uint64_t count) {
        return iterator->setOrder(count);
    }
    kmerRow getKmerRow() {
        return kmerRow(iterator->getKmer(),
            iterator->getHashedKmer(),
            getCount(),
            iterator->getOrder(),
            origin
        );
    }
    // template<typename T, typename Container>
    // void getColumnValue(const string& colName, T& res);
    template<typename T, typename Container>
    void setColumnValue(const string& colName, T value);


    void setKmerColumnValueFromOtherColumn(kDataFrame* input, string inputColName, string outputColName);

    ~kDataFrameIterator() {
        delete iterator;
    }
};


// %extend kDataFrameIterator {
//     %template(getColumnValue_deduplicate) getColumnValue<vector<string>, deduplicatedColumn<vector<string>, StringColorColumn>>;
// }