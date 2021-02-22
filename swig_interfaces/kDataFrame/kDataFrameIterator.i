
class kDataFrameIterator{
        private:
        kDataFrame* origin;
        _kDataFrameIterator* iterator;
        public:
        using iterator_category = std::forward_iterator_tag;
        kDataFrameIterator(){
            iterator=NULL;
        }
        kDataFrameIterator(const kDataFrameIterator& other){
            if(other.iterator!=NULL){
                iterator=other.iterator->clone();
                origin=other.origin;
            }
            else{
                iterator=NULL;
                origin=NULL;
            }
        }
        kDataFrameIterator& operator= (const kDataFrameIterator& other){
            if(other.iterator!=NULL){
                iterator=other.iterator->clone();
                origin=other.origin;
            }
            else{
                iterator=NULL;
                origin=NULL;
            }
            return *this;
        }
        kDataFrameIterator(_kDataFrameIterator* it,kDataFrame* o){
            this->origin=o;
            iterator=it;
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

        /// Returns the current kmer
        string getKmer(){
            return iterator->getKmer();
        }
        /// Returns the count of the current kmer
        std::uint64_t getCount(){
            return iterator->getCount();
        }
        /// sets the count of the current kmer
        bool setCount(std::uint64_t count){
            return iterator->setCount(count);
        }
        kmerRow operator*(){
            return kmerRow(iterator->getKmer(),
                           iterator->getHashedKmer(),
                           iterator->getCount(),
                           origin
            );
        }
        ~kDataFrameIterator(){
            delete iterator;
        }
};