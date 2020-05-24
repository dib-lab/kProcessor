using namespace std;
class Hasher{
public:
	virtual uint64_t hash(string key){return 0;};
	virtual uint64_t hash(uint64_t key){return 0;};
	virtual string Ihash(uint64_t key){
		throw logic_error("Reverese Hash function is not/ cannot be implemented for this hash function.");
	}
	virtual Hasher* clone(){return this;};
};

// TwoBitsHasher
class TwoBitsHasher: public Hasher{
        private:
        uint64_t kSize;
        public:
        TwoBitsHasher(uint64_t kSize);
        Hasher* clone() { return new TwoBitsHasher(kSize);}
        uint64_t hash(string key);
        uint64_t hash(uint64_t key);
        string Ihash(uint64_t key);
};
