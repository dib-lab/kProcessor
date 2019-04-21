class kmer {
public:
	static uint64_t str_to_int(string str);
	static uint64_t str_to_canonical_int(string str);
	static string canonicalKmer(string k);
	static string int_to_str(uint64_t kmer, uint32_t K);
	static uint64_t reverse_complement(uint64_t kmer, uint32_t K);
private:
	kmer();
};
