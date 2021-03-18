#include "algorithms.hpp"
#include <iostream>
#include "Utils/kmer.h"
#include <fstream>
#include <limits>
#include <omp.h>
#include <stdexcept>
#include <math.h>
#include <deque>
#include <gqf.h>
#include <string>
#include <queue>
#include <functional>
#include <limits>
#include <parallel_hashmap/phmap.h>
#include "ntcard.hpp"
#include <cstdio>
#include <stack>
#include <chrono>
#include "defaultColumn.hpp"
#include "kmc_file.h"


using namespace std::chrono;


using std::string;
using std::vector;
using std::cerr;
using std::cout;

using phmap::flat_hash_map;

#define QBITS_LOCAL_QF 16

namespace kProcessor {
    // TO BE REMOVED TODO V2
    /*
    static inline void insertToLevels(uint64_t item, QF *local, QF *main, QF *diskMQF = NULL) {
        if (!qf_insert(main, item, 1,
                       true, false)) {
            qf_insert(local, item, 1,
                      false, false);
            // check of the load factor of the local QF is more than 50%
            if (qf_space(local) > 50) {

                {
                    if (main->metadata->noccupied_slots + local->metadata->noccupied_slots
                        < main->metadata->maximum_occupied_slots) {
                        qf_migrate(local, main);
                    } else if (diskMQF != NULL) {
                        SEQAN_OMP_PRAGMA(critical){
                            qf_general_lock(main, true);
                            qf_migrate(main, diskMQF);
                            qf_reset(main);
                            qf_general_unlock(main);
                            qf_migrate(local, main);
                        }
                    } else {
                        throw overflow_error("memory MQF doesn't have enough space to migrate");
                    }

                }
                qf_reset(local);
            }
        } else {
            if (qf_space(main) > 90) {
                SEQAN_OMP_PRAGMA(critical)
                {
                    if (qf_space(main) > 90) {
                        if (diskMQF != NULL) {
                            qf_general_lock(main, true);
                            qf_migrate(main, diskMQF);
                            qf_reset(main);
                            qf_general_unlock(main);
                        } else {
                            throw overflow_error("memory MQF doesn't have enough space");
                        }
                    }
                }
            }
        }
    }
    */

// TO BE REMOVED TODO V2
/*
    void loadIntoMQF(string sequenceFilename, unsigned int ksize, int noThreads, Hasher *hasher, QF *memoryMQF,
                     QF *diskMQF) {
        FastqReaderSqueker reader(sequenceFilename);
        omp_set_num_threads(noThreads);
        QF *localMQF;
        bool moreWork = true;
        uint64_t numReads = 0;
        deque<pair<string, string> > reads;
        string read, tag;
#pragma omp parallel private(reads, localMQF, read, tag) shared(reader, moreWork, numReads)  firstprivate(ksize, noThreads, memoryMQF, diskMQF)
        {
            auto localHasher = hasher->clone();
            localMQF = new QF();
            reads = deque<pair<string, string> >(15000);
            qf_init(localMQF, (1ULL << QBITS_LOCAL_QF), memoryMQF->metadata->key_bits,
                    0, memoryMQF->metadata->fixed_counter_size, 0, true, "", 2038074761);
            while (moreWork) {
                SEQAN_OMP_PRAGMA(critical)
                {
                    reader.readNSeq(&reads, 15000);
                    numReads += 15000;
                    bool tmp = !reader.isEOF();
                    moreWork = tmp;
                }

                for (unsigned int j = 0; j < reads.size(); j++) {
                    read = reads[j].first;
                    start_read:
                    if (read.size() < ksize) {
                        continue;
                    }

                    uint64_t first = 0;
                    uint64_t first_rev = 0;
                    uint64_t item = 0;

                    for (unsigned int i = 0; i < ksize; i++) {
                        //First kmer
                        uint8_t curr = kmer::map_base(read[i]);
                        if (curr > DNA_MAP::G) {
                            // 'N' is encountered

                            read = read.substr(i + 1, read.length());

                            //continue;
                            goto start_read;
                        }
                        first = first | curr;
                        first = first << 2;
                    }
                    first = first >> 2;
                    first_rev = kmer::reverse_complement(first, ksize);


                    if (kmer::compare_kmers(first, first_rev))
                        item = first;
                    else
                        item = first_rev;

                    item = localHasher->hash(item) % memoryMQF->metadata->range;
                    insertToLevels(item, localMQF, memoryMQF, diskMQF);

                    uint64_t next = (first << 2) & BITMASK(2 * ksize);
                    uint64_t next_rev = first_rev >> 2;

                    for (uint32_t i = ksize; i < length(read); i++) {
                        //next kmers
                        //cout << "K: " << read.substr(i-K+1,K) << endl;
                        uint8_t curr = kmer::map_base(read[i]);
                        if (curr > DNA_MAP::G) {
                            // 'N' is encountered
                            //continue;
                            //read = read.substr(i+1, length(read));
                            read = read.substr(i + 1, read.length());
                            //erase(read,0,i+1);

                            goto start_read;
                        }
                        next |= curr;
                        uint64_t tmp = kmer::reverse_complement_base(curr);
                        tmp <<= (ksize * 2 - 2);
                        next_rev = next_rev | tmp;
                        if (kmer::compare_kmers(next, next_rev))
                            item = next;
                        else
                            item = next_rev;


                        item = localHasher->hash(item) % memoryMQF->metadata->range;
                        insertToLevels(item, localMQF, memoryMQF, diskMQF);
                        next = (next << 2) & BITMASK(2 * ksize);
                        next_rev = next_rev >> 2;
                    }
                }

            }
            //    #pragma omp critical
            {
                qf_migrate(localMQF, memoryMQF);
            }
            qf_destroy(localMQF);
        }
        if (diskMQF != NULL) {
            qf_migrate(memoryMQF, diskMQF);
        }

    }
*/

    void dumpMQF(QF *MQF, int ksize, std::string outputFilename) {
        IntegerHasher Ihasher(BITMASK(2 * ksize));
        ofstream output(outputFilename.c_str());
        QFi qfi;
        qf_iterator(MQF, &qfi, 0);
        do {
            uint64_t key, value, count;
            qfi_get(&qfi, &key, &value, &count);
            string kmer = Ihasher.Ihash(key);
            output << kmer << " " << count << "\n";
        } while (!qfi_next(&qfi));
    }

    bool isEnough(vector<uint64_t> histogram, uint64_t noSlots, uint64_t fixedSizeCounter, uint64_t slotSize) {
        // cout<<"noSlots= "<<noSlots<<endl
        //     <<"fcounter= "<<fixedSizeCounter<<endl
        //     <<"slot size= "<<numHashBits<<endl;

        noSlots = (uint64_t) ((double) noSlots * 0.90);
        for (uint64_t i = 1; i < 1000; i++) {
            uint64_t usedSlots = 1;

            if (i > ((1ULL) << fixedSizeCounter) - 1) {
                uint64_t nSlots2 = 0;
                __uint128_t capacity;
                do {
                    nSlots2++;
                    capacity = ((__uint128_t) (1ULL) << (nSlots2 * slotSize + fixedSizeCounter)) - 1;
                    //  cout<<"slots num "<<nSlots2<<" "<<capacity<<endl;
                } while ((__uint128_t) i > capacity);
                usedSlots += nSlots2;
            }
            //cout<<"i= "<<i<<"->"<<usedSlots<<" * "<<histogram[i]<<endl;
            if (noSlots >= (usedSlots * histogram[i])) {
                noSlots -= (usedSlots * histogram[i]);
            } else {
                //  cout<<"failed"<<endl<<endl;
                return false;
            }

        }
        //cout<<"success"<<endl<<endl;
        return true;
    }


    void estimateMemRequirement(std::string ntcardFilename,
                                uint64_t numHashBits, uint64_t tagSize,
                                uint64_t *res_noSlots, uint64_t *res_fixedSizeCounter, uint64_t *res_memory) {
        uint64_t noDistinctKmers = 0, totalNumKmers=0;
        vector<uint64_t> histogram(1000, 0);
        ifstream ntcardFile(ntcardFilename);
        string f;
        uint64_t count;
        while (ntcardFile >> f >> count) {
            if (count == numeric_limits<uint64_t>::max())
                continue;
            if (f == "F0")
                noDistinctKmers = count;
            else if (f == "F1")
                totalNumKmers = count;
            else {
                f = f.substr(1, f.size());
                int n = atoi(f.c_str());
                histogram[n] = count;
            }
        }
        *res_memory = numeric_limits<uint64_t>::max();
        for (int i = 8; i < 64; i++) {
            uint64_t noSlots = (1ULL) << i;
            if (noSlots < noDistinctKmers)
                continue;
            bool moreWork = false;
            uint64_t slotSize = numHashBits - log2((double) noSlots);
            for (uint64_t fixedSizeCounter = 1; fixedSizeCounter < slotSize; fixedSizeCounter++) {
                if (isEnough(histogram, noSlots, fixedSizeCounter, slotSize)) {
                    uint64_t tmpMem = estimateMemory(noSlots, slotSize, fixedSizeCounter, tagSize);
                    if (*res_memory > tmpMem) {
                        *res_memory = tmpMem;
                        *res_fixedSizeCounter = fixedSizeCounter;
                        *res_noSlots = noSlots;
                        moreWork = true;
                    } else {
                        break;
                    }
                }

            }
            if (!moreWork && *res_memory != numeric_limits<uint64_t>::max())
                break;
        }
        if (*res_memory == numeric_limits<uint64_t>::max()) {
            throw std::overflow_error(
                    "Data limits exceeds MQF capabilities(> uint64). Check if ntcard file is corrupted");
        }


    }

    kDataFrame *transform(kDataFrame *input, function<kmerRow (kmerRow i)> fn) {
        kDataFrame *res = input->getTwin();
        kDataFrameIterator it = input->begin();
        while (it != input->end()) {
            kmerRow newkmer = fn(*it);
            res->insert(newkmer);
            it++;
        }
        for(auto col: input->columns)
        {
            string newColName= col.first;
            Column* column=col.second->getTwin();
            column->setSize(res->size());
            res->addColumn(newColName, column);
        }
        for(auto kmer:*res)
        {
            for (auto col: input->columns) {
                string newColName = col.first;
                res->setKmerColumnValueFromOtherColumn(input,col.first, newColName,kmer.kmer);
            }
        }
        delete input;
        return res;

    }
    kDataFrame* innerJoin(vector<kDataFrame *> input, vector<uint32_t> kmersToKeep) {
        kDataFrame *res = input[0]->getTwin();
        uint64_t numKmers = 0;
        for (auto kframe:input) {
            numKmers += kframe->size();
        }
        res->reserve((uint64_t) ((double) numKmers * 0.75));
        merge(input, res, [&](vector<kmerRow> &input) -> kmerRow {
            kmerRow res;
            bool exists=true;
            for (auto i : kmersToKeep ) {
                if (input[i].count == 0) {
                    exists = false;
                    break;
                }
            }
            if(exists) {
                uint32_t i=kmersToKeep[0];
                return kmerRow(input[i].kmer,input[i].hashedKmer,input[i].count,nullptr);
            }
            return kmerRow();
        });
        for(unsigned int i=0; i <  input.size() ; i++)
        {
            for(auto col: input[i]->columns)
            {
                string newColName= col.first+"."+to_string(i);
                Column* column=col.second->getTwin();
                column->setSize(res->size());
                res->addColumn(newColName, column);
            }
        }
        for(auto kmer:*res)
        {
            for(unsigned int i=0; i <  input.size() ; i++)
            {
                if(input[i]->kmerExist(kmer.kmer)) {
                    for (auto col: input[i]->columns) {
                        string newColName = col.first + "." + to_string(i);
                        res->setKmerColumnValueFromOtherColumn(input[i],col.first, newColName,kmer.kmer);
                    }
                }
            }
        }
        return res;
    }

    kDataFrame *filter(kDataFrame *input, function<bool (kmerRow i)> fn) {
        kDataFrame *res = input->getTwin();
        kDataFrameIterator it = input->begin();
        while (it != input->end()) {
            if (fn(*it))
                res->insert(*it);
            it++;
        }
        for(auto col: input->columns)
        {
            string newColName= col.first;
            Column* column=col.second->getTwin();
            column->setSize(res->size());
            res->addColumn(newColName, column);
        }
        for(auto kmer:*res)
        {
            for (auto col: input->columns) {
                string newColName = col.first;
                res->setKmerColumnValueFromOtherColumn(input,col.first, newColName,kmer.kmer);
            }
        }
        delete input;
        return res;

    }

    any aggregate(kDataFrame *input, any initial, function<any (kmerRow i, any v)> fn) {
        kDataFrameIterator it = input->begin();
        while (it != input->end()) {
            initial = fn(*it, initial);
            it++;
        }
        return initial;
    }

    void parseSequences(kmerDecoder *KD, kDataFrame *output) {
        if (KD->get_kSize() != (int) output->getkSize()) {
            std::cerr << "kmerDecoder kSize must be equal to kDataFrame kSize" << std::endl;
            exit(1);
        }

        while (!KD->end()) {
            KD->next_chunk();
            for (const auto &seq : *KD->getKmers()) {
                for (const auto &kmer : seq.second) {
                    output->insert(kmer.hash);
                }
            }
        }
    }


    void countKmersFromFile(kDataFrame *kframe, std::map<std::string, int> parse_params, string filename, int chunk_size) {
        // parse_params["mode"] = 1 > Default: Kmers
        // parse_params["mode"] = 2 > Skipmers
        // parse_params["mode"] = 3 > Minimizers

        // Initialize kmerDecoder
        //  vector<uint64_t> countHistogram= estimateKmersHistogram(filename, kframe->getkSize() ,1);
        kframe->reserve(100000);
        std::string mode = "kmers";
        bool check_mode = (parse_params.find("mode") != parse_params.end());
        if (check_mode) {
            if (parse_params["mode"] == 2) mode = "skipmers";
            else if (parse_params["mode"] == 3) mode = "minimizers";
        }

        parse_params["k_size"] = kframe->ksize();
        parse_params["k"] = kframe->ksize();
        kmerDecoder *KD = initialize_kmerDecoder(filename, chunk_size, mode, parse_params);

        // Clone the hashing

        kmerDecoder_setHashing(KD, kframe->KD->hash_mode, kframe->KD->canonical);

        // Processing

        while (!KD->end()) {
            KD->next_chunk();
            for (const auto &seq : *KD->getKmers()) {
                for (const auto &kmer : seq.second) {
                    kframe->insert(kmer.hash);
                }
            }
        }
        delete KD;


    }

    void countKmersFromString(kmerDecoder *KD, string sequence, kDataFrame *output) {
        if (KD->get_kSize() != (int) output->getkSize()) {
            std::cerr << "kmerDecoder kSize must be equal to kDataFrame kSize" << std::endl;
            exit(1);
        }

        std::vector<kmer_row> kmers;
        KD->seq_to_kmers(sequence, kmers);

        for (const auto &kmer : kmers) {
            output->insert(kmer.hash);
        }

    }

    void countKmersFromString(kDataFrame *frame, std::map<std::string, int> parse_params, string sequence) {

        // parse_params["mode"] = 1 > Default: Kmers
        // parse_params["mode"] = 2 > Skipmers
        // parse_params["mode"] = 3 > Minimizers

        // Initialize kmerDecoder
        std::string mode = "kmers";
        bool check_mode = (parse_params.find("mode") != parse_params.end());
        if (check_mode) {
            if (parse_params["mode"] == 2) mode = "skipmers";
            else if (parse_params["mode"] == 3) mode = "minimizers";
        }

        parse_params["k_size"] = frame->ksize();
        parse_params["k"] = frame->ksize();
        kmerDecoder *KD = initialize_kmerDecoder(mode, parse_params);

        // Clone the hashing

        kmerDecoder_setHashing(KD, frame->KD->hash_mode, frame->KD->canonical);

        if (KD->get_kSize() != (int) frame->getkSize()) {
            std::cerr << "kmerDecoder kSize must be equal to kDataFrame kSize" << std::endl;
            exit(1);
        }

        std::vector<kmer_row> kmers;
        KD->seq_to_kmers(sequence, kmers);

        for (const auto &kmer : kmers) {
            frame->insert(kmer.hash);
        }

    }


    struct CustomKmerRow {
        bool operator()(const pair<kmerRow, int> &lhs, const pair<kmerRow, int> &rhs) {
            return lhs.first.hashedKmer > rhs.first.hashedKmer;
        }
    };

    inline void terminate_if_kDataFrameMAP(const vector<kDataFrame *> &input) {
        for (const auto &kFrame : input) {
            if (kFrame->get_class_name() == "MAP") {
                throw logic_error("can't apply set functions kDataFramMAP on kDataFrameMAP.");
            }
        }
    }

    void merge(const vector<kDataFrame *> &input, kDataFrame *res,function<kmerRow (vector<kmerRow> &i)> fn) {
        terminate_if_kDataFrameMAP(input);
        priority_queue<pair<kmerRow, int>, vector<pair<kmerRow, int> >, CustomKmerRow> Q;
        vector<kDataFrameIterator> iterators(input.size());
        for (unsigned int i = 0; i < input.size(); i++) {
            iterators[i] = input[i]->begin();
            if (iterators[i] != input[i]->end()) {
                //  cout<<i<<" "<<(*iterators[i]).hashedKmer<<endl;
                Q.push(make_pair(*iterators[i], i));
            }
        }
        vector<kmerRow> current(input.size());
        while (Q.size() > 0) {
            for (unsigned int i = 0; i < current.size(); i++)
                current[i] = kmerRow();
            pair<kmerRow, int> top = Q.top();
            Q.pop();
            iterators[top.second]++;
            current[top.second] = top.first;
            //  cout<<top.first.hashedKmer<<" "<<top.second<<endl;
            if (iterators[top.second] != input[top.second]->end())
                Q.push(make_pair(*iterators[top.second], top.second));
            while (Q.size() > 0 && top.first.hashedKmer == Q.top().first.hashedKmer) {
                top = Q.top();
                Q.pop();
                current[top.second] = top.first;
                iterators[top.second]++;
                if (iterators[top.second] != input[top.second]->end())
                    Q.push(make_pair(*iterators[top.second], top.second));
            }

            kmerRow newRow = fn(current);
            if(newRow.count > 0)
                res->insert(newRow);
        }

    }

    kDataFrame *kFrameUnion(const vector<kDataFrame *> &input) {
        kDataFrame *res = input[0]->getTwin();
        uint64_t numKmers = 0;
        for (auto kframe:input) {
            numKmers += kframe->size();
        }
        res->reserve((uint64_t) ((double) numKmers * 0.75));
        merge(input, res, [](vector<kmerRow> &input) -> kmerRow {
            kmerRow res;
            for (unsigned int i = 0; i < input.size(); i++) {
                if (input[i].count != 0) {
                    res.kmer = input[i].kmer;
                    res.hashedKmer = input[i].hashedKmer;
                }
                res.count += input[i].count;
            }
            return res;
        });
        return res;
    }

    kDataFrame *kFrameIntersect(const vector<kDataFrame *> &input) {
        kDataFrame *res = input[0]->getTwin();
        uint64_t numKmers = numeric_limits<uint64_t>::max();
        for (auto kframe:input) {
            numKmers = min(numKmers, (uint64_t) kframe->size());
        }
        res->reserve((uint64_t) ((double) numKmers * 1.2));
        merge(input, res, [](vector<kmerRow> &input) -> kmerRow {
            kmerRow res = input[0];
            for (unsigned int i = 1; i < input.size(); i++) {
                //cout<<input[i].kmer<<endl;
                res.count = min(res.count, input[i].count);
            }
            if (res.count == 0)
                return kmerRow();
            //cout<<res.count<<endl;
            return res;
        });
        //cout<<"Size "<<res->size()<<endl;
        return res;
    }

    kDataFrame *kFrameDiff(const vector<kDataFrame *> &input) {
        kDataFrame *res = input[0]->getTwin();
        merge(input, res, [](vector<kmerRow> &input) -> kmerRow {
            kmerRow res = input[0];
            bool found = false;
            for (unsigned int i = 1; i < input.size(); i++) {
                if (input[i].count > 0) {
                    found = true;
                    break;
                }
            }
            if (found)
                return kmerRow();
            return res;
        });
        return res;
    }


    void kmerDecoder_setHashing(kmerDecoder *KD, int hash_mode, bool canonical) {
        KD->setHashingMode(hash_mode, canonical);
    }

    void kmerDecoder_setHashing(kDataFrame *KF, int hash_mode, bool canonical) {
        KF->KD->setHashingMode(hash_mode, canonical);
    }


    kmerDecoder *initialize_kmerDecoder(std::string filename, int chunkSize, std::string mode,
                                        std::map<std::string, int> parse_params) {

        std::string func_name = "wrong parameters in initialize_kmerDecoder() : \n";

        // for avoiding case sensitivity issues.
        transform(mode.begin(), mode.end(), mode.begin(), ::tolower);

        if (mode == "kmers") {
            if (parse_params.find("k_size") != parse_params.end()) {
                return new Kmers(filename, chunkSize, parse_params["k_size"]);
            } else {
                std::cerr << func_name << "kmerDecoder Kmers parameters {k_size} validation failed" << std::endl;
                exit(1);
            }
        } else if (mode == "skipmers") {
            bool check_k = (parse_params.find("k_size") != parse_params.end());
            bool check_m = (parse_params.find("m") != parse_params.end());
            bool check_n = (parse_params.find("n") != parse_params.end());
            bool check_orf = (parse_params.find("orf") != parse_params.end());

            if (check_k && check_m && check_n) {
                if (check_orf)
                    return new Skipmers(filename, chunkSize, parse_params["m"], parse_params["n"],
                                        parse_params["k_size"], parse_params["orf"]);
                return new Skipmers(filename, chunkSize, parse_params["m"], parse_params["n"], parse_params["k_size"]);
            } else {
                std::cerr << func_name << "kmerDecoder Skipmers parameters {k_size, m, n} validation failed"
                          << std::endl;
                exit(1);
            }
        } else if (mode == "minimizers") {
            bool check_k = (parse_params.find("k_size") != parse_params.end());
            bool check_w = (parse_params.find("w") != parse_params.end());

            if (check_k && check_w) {
                return new Minimizers(filename, chunkSize, parse_params["k_size"], parse_params["w"]);
            } else {
                std::cerr << func_name << "kmerDecoder Skipmers parameters {k_size, w} validation failed" << std::endl;
                exit(1);
            }

        } else {
            std::cerr << func_name << "supported kmerDecoder modes: {kmers, skipmers, minimizers}" << std::endl;
            exit(1);
        }
    }

    kmerDecoder *initialize_kmerDecoder(std::string mode, std::map<std::string, int> parse_params) {

        std::string func_name = "wrong parameters in initialize_kmerDecoder() : \n";

        // for avoiding case sensitivity issues.
        transform(mode.begin(), mode.end(), mode.begin(), ::tolower);

        if (mode == "kmers") {
            if (parse_params.find("k_size") != parse_params.end()) {
                return new Kmers(parse_params["k_size"]);
            } else {
                std::cerr << func_name << "kmerDecoder Kmers parameters {k_size} validation failed" << std::endl;
                exit(1);
            }
        } else if (mode == "skipmers") {
            bool check_k = (parse_params.find("k_size") != parse_params.end());
            bool check_m = (parse_params.find("m") != parse_params.end());
            bool check_n = (parse_params.find("n") != parse_params.end());
            bool check_orf = (parse_params.find("orf") != parse_params.end());

            if (check_k && check_m && check_n) {
                if (check_orf)
                    return new Skipmers(parse_params["m"], parse_params["n"], parse_params["k_size"],
                                        parse_params["orf"]);
                return new Skipmers(parse_params["m"], parse_params["n"], parse_params["k_size"]);
            } else {
                std::cerr << func_name << "kmerDecoder Skipmers parameters {k_size, m, n} validation failed"
                          << std::endl;
                exit(1);
            }
        } else if (mode == "minimizers") {
            bool check_k = (parse_params.find("k_size") != parse_params.end());
            bool check_w = (parse_params.find("w") != parse_params.end());

            if (check_k && check_w) {
                return new Minimizers(parse_params["k_size"], parse_params["w"]);
            } else {
                std::cerr << func_name << "kmerDecoder Skipmers parameters {k_size, w} validation failed" << std::endl;
                exit(1);
            }

        } else {
            std::cerr << func_name << "supported kmerDecoder modes: {kmers, skipmers, minimizers}" << std::endl;
            exit(1);
        }
    }

    kmerDecoder *initialize_kmerDecoder(int kSize, int hash_mode) {
        // Mode 0: Murmar Hashing | Irreversible
        // Mode 1: Integer Hashing | Reversible | Full Hashing
        // Mode 2: TwoBitsHashing | Not considered hashing, just store the two bits representation
        return new Kmers(kSize, hash_mode);
    }

    void index(kmerDecoder *KD, string names_fileName, kDataFrame *frame) {

        if (KD->get_kSize() != (int) frame->ksize()) {
            std::cerr << "kmerDecoder kSize must be equal to kDataFrame kSize" << std::endl;
            exit(1);
        }


        flat_hash_map<string, string> namesMap;
        flat_hash_map<string, uint64_t> tagsMap;
        flat_hash_map<string, uint64_t> groupNameMap;
        flat_hash_map<uint64_t, std::vector<uint32_t>> *legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
        flat_hash_map<uint64_t, uint64_t> colorsCount;
        uint64_t readID = 0, groupID = 1;
        ifstream namesFile(names_fileName.c_str());
        string seqName, groupName;
        string line;
        priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
        flat_hash_map<string, uint64_t> groupCounter;
//        while (namesFile >> seqName >> groupName) {
        while (std::getline(namesFile, line)) {
            std::vector<string> tokens;
            std::istringstream iss(line);
            std::string token;
            while (std::getline(iss, token, '\t'))   // but we can specify a different one
                tokens.push_back(token);
            seqName = tokens[0];
            groupName = tokens[1];
            namesMap.insert(make_pair(seqName, groupName));
            auto it = groupNameMap.find(groupName);
            groupCounter[groupName]++;
            if (it == groupNameMap.end()) {
                groupNameMap.insert(make_pair(groupName, groupID));
                tagsMap.insert(make_pair(to_string(groupID), groupID));
                vector<uint32_t> tmp;
                tmp.clear();
                tmp.push_back(groupID);
                legend->insert(make_pair(groupID, tmp));
                colorsCount.insert(make_pair(groupID, 0));
                groupID++;
            }
        }


        flat_hash_map<uint64_t, string> inv_groupNameMap;
            for (auto &_ : groupNameMap)
                inv_groupNameMap[_.second] = _.first;

        vector<kDataFrameMQF *> frames;
//        int currIndex = 0;
        string kmer;
//        uint64_t tagBits = 0;
//        uint64_t maxTagValue = (1ULL << tagBits) - 1;
        //  kDataFrame *frame;
//        int kSize = KD->get_kSize();


//        uint64_t lastTag = 0;
        readID = 0;

        while (!KD->end()) {
            KD->next_chunk();

            flat_hash_map<uint64_t, uint64_t> convertMap;

            for (const auto &seq : *KD->getKmers()) {
                string readName = seq.first;

                auto it = namesMap.find(readName);
                if (it == namesMap.end()) {
                    cerr << "WARNING: " << "read " << readName << "doesn't have a group. Please, check the names file." << endl;
                    continue;
                }
                string groupName = it->second;

                uint64_t readTag = groupNameMap.find(groupName)->second;


                convertMap.clear();
                convertMap.insert(make_pair(0, readTag));
                convertMap.insert(make_pair(readTag, readTag));
                //    cout<<readName<<"   "<<seq.size()<<endl;
                for (const auto &kmer : seq.second) {
                    uint64_t currentTag = frame->getCount(kmer.hash);
                    auto itc = convertMap.find(currentTag);
                    if (itc == convertMap.end()) {
                        vector<uint32_t> colors = legend->find(currentTag)->second;
                        auto tmpiT = find(colors.begin(), colors.end(), readTag);
                        if (tmpiT == colors.end()) {
                            colors.push_back(readTag);
                            sort(colors.begin(), colors.end());
                        }

                        string colorsString = to_string(colors[0]);
                        for (unsigned int k = 1; k < colors.size(); k++) {
                            colorsString += ";" + to_string(colors[k]);
                        }

                        auto itTag = tagsMap.find(colorsString);
                        if (itTag == tagsMap.end()) {
                            uint64_t newColor;
                            if (freeColors.size() == 0) {
                                newColor = groupID++;
                            } else {
                                newColor = freeColors.top();
                                freeColors.pop();
                            }

                            tagsMap.insert(make_pair(colorsString, newColor));
                            legend->insert(make_pair(newColor, colors));
                            itTag = tagsMap.find(colorsString);
                            colorsCount[newColor] = 0;
                            // if(groupID>=maxTagValue){
                            //   cerr<<"Tag size is not enough. ids reached "<<groupID<<endl;
                            //   return -1;
                            // }
                        }
                        uint64_t newColor = itTag->second;

                        convertMap.insert(make_pair(currentTag, newColor));
                        itc = convertMap.find(currentTag);
                    }

                    if (itc->second != currentTag) {

                        colorsCount[currentTag]--;
                        if (colorsCount[currentTag] == 0 && currentTag != 0) {
                            auto _invGroupNameIT = inv_groupNameMap.find(currentTag);
                            if (_invGroupNameIT == inv_groupNameMap.end()){
                                freeColors.push(currentTag);

                                vector<uint32_t> colors = legend->find(currentTag)->second;
                                string colorsString = to_string(colors[0]);
                                for (unsigned int k = 1; k < colors.size(); k++) {
                                    colorsString += ";" + to_string(colors[k]);
                                }
                                tagsMap.erase(colorsString);
                                
                                legend->erase(currentTag);
                                if (convertMap.find(currentTag) != convertMap.end())
                                    convertMap.erase(currentTag);
                            }

                        }
                        colorsCount[itc->second]++;
                    }

                    frame->setCount(kmer.hash, itc->second);
                    if (frame->getCount(kmer.hash) != itc->second) {
                        //frame->setC(kmer,itc->second);
                        cout << "Error Founded " << kmer.str << " from sequence " << readName << " expected "
                             << itc->second << " found " << frame->getCount(kmer.hash) << endl;
                        return;
                    }
                }
                readID += 1;
                groupCounter[groupName]--;
                if (colorsCount[readTag] == 0) {                   
                    if (groupCounter[groupName] == 0) {
                        freeColors.push(readTag);
                        legend->erase(readTag);
                    }
                }
            }
        }
        StringColorColumn *col = new StringColorColumn();
        frame->changeDefaultColumnType((Column *) col);
        col->colors = vector<vector<uint32_t> >(legend->size());
        col->colors.push_back(vector<uint32_t>());
        for (auto it : *legend) {
            col->colors[it.first] = it.second;
        }


        for (auto iit = namesMap.begin(); iit != namesMap.end(); iit++) {
            uint32_t sampleID = groupNameMap[iit->second];
            col->namesMap[sampleID] = iit->second;
        }

    }

    void index(kDataFrame *frame, std::map<std::string, int> parse_params, string filename, int chunk_size,
               string names_fileName) {

        // parse_params["mode"] = 1 > Default: Kmers
        // parse_params["mode"] = 2 > Skipmers
        // parse_params["mode"] = 3 > Minimizers

        // Initialize kmerDecoder
        std::string mode = "kmers";
        bool check_mode = (parse_params.find("mode") != parse_params.end());
        if (check_mode) {
            if (parse_params["mode"] == 2) mode = "skipmers";
            else if (parse_params["mode"] == 3) mode = "minimizers";
        }

        parse_params["k_size"] = frame->ksize();
        parse_params["k"] = frame->ksize();
        kmerDecoder *KD = initialize_kmerDecoder(filename, chunk_size, mode, parse_params);

        // Clone the hashing

        kmerDecoder_setHashing(KD, frame->KD->hash_mode, frame->KD->canonical);

        // Processing

        if (KD->get_kSize() != (int) frame->ksize()) {
            std::cerr << "kmerDecoder kSize must be equal to kDataFrame kSize" << std::endl;
            exit(1);
        }
        index(KD, names_fileName, frame);


    }

    void indexPriorityQueue(vector<kDataFrame *> &input, string tmpFolder, kDataFrame *output) {
        insertColorColumn *colors = new insertColorColumn(input.size(), tmpFolder);
        output->changeDefaultColumnType(colors);
        for (unsigned int i = 0; i < input.size(); i++) {
            vector<uint32_t> tmp = {i};
            colors->insertAndGetIndex(tmp);
        }


        auto compare = [](tuple<uint64_t, uint64_t, kDataFrameIterator *, kDataFrameIterator *> lhs,
                          tuple<uint64_t, uint64_t, kDataFrameIterator *, kDataFrameIterator *> rhs) {
            if (get<0>(lhs) == get<0>(rhs))
                return get<1>(lhs) > get<1>(rhs);
            return get<0>(lhs) > get<0>(rhs);
        };

        priority_queue<tuple<uint64_t, uint64_t, kDataFrameIterator *, kDataFrameIterator *>, vector<tuple<uint64_t, uint64_t, kDataFrameIterator *, kDataFrameIterator *> >, decltype(compare)> nextKmer(
                compare);

        for (unsigned int i = 0; i < input.size(); i++) {
            kDataFrameIterator *it = new kDataFrameIterator(input[i]->begin());
            kDataFrameIterator *itend = new kDataFrameIterator(input[i]->end());
            nextKmer.push(make_tuple(it->getHashedKmer(), i, it, itend));
        }

        uint64_t processedKmers = 0;
        while (nextKmer.size() > 0) {
            vector<uint32_t> colorVec;
            colorVec.clear();
            uint64_t currHash = get<0>(nextKmer.top());
            processedKmers++;
            if (processedKmers % 1000000 == 0)
                cout << processedKmers << " Kmers Processed" << endl;
            while (nextKmer.size() > 0 && get<0>(nextKmer.top()) == currHash) {
                auto colorTuple = nextKmer.top();
                nextKmer.pop();
                colorVec.push_back(get<1>(colorTuple));
                get<2>(colorTuple)->next();
                if (*get<2>(colorTuple) != *get<3>(colorTuple)) {
                    get<0>(colorTuple) = get<2>(colorTuple)->getHashedKmer();
                    nextKmer.push(colorTuple);
                } else {
                    delete get<2>(colorTuple);
                    delete get<3>(colorTuple);
                }
            }

            uint64_t prevColor = output->getCount(currHash);
            if (prevColor != 0) {
                auto res = output->getKmerDefaultColumnValue<vector<uint32_t>, insertColorColumn>(currHash);
                cout << "Error in Indexing detected at kmer " << currHash << endl;
                cout << "should be empty vector and found  ";
                for (auto a:res)
                    cout << a << " ";
                cout << endl << endl;
            }

            output->setKmerDefaultColumnValue<vector<uint32_t> &, insertColorColumn>(currHash, colorVec);

            // auto res=output->getKmerDefaultColumnValue<vector<uint32_t >, insertColorColumn>(currHash);
            // //	cout<<res.size()<<endl;
            // if(!equal(res.begin(),res.end(),colorVec.begin()))
            //   {
            //     cout<<"Error in Indexing detected at kmer "<<currHash<<endl;
            //     cout<<"expected ";
            //     for(auto a:colorVec)
            //       cout<<a<<" ";
            //     cout<<endl<<"Found ";
            //     for(auto a:res)
            //       cout<<a<<" ";
            //     cout<<endl;
            //   }

        }
        colors->populateColors();
        uint64_t noColors = colors->noColors;
        cout << noColors << " colors created" << endl;
        delete colors;

        queryColorColumn *qcolors = new queryColorColumn(input.size(), noColors, tmpFolder);
        qcolors->explainSize();
        output->changeDefaultColumnType(qcolors);
    }

    void mergeIndexes(vector<kDataFrame *> &input, kDataFrame *output) {

        vector<uint32_t> idsOffset(input.size());
        idsOffset[0] = 0;
        for (unsigned int i = 1; i < input.size(); i++) {
            idsOffset[i] = idsOffset[i - 1];
            idsOffset[i] += ((insertColorColumn *) input[i - 1]->getDefaultColumn())->noSamples;
        }
        insertColorColumn *colors = new insertColorColumn();
        output->changeDefaultColumnType(colors);


        auto compare = [](tuple<uint64_t, uint64_t, kDataFrameIterator *, kDataFrameIterator *> lhs,
                          tuple<uint64_t, uint64_t, kDataFrameIterator *, kDataFrameIterator *> rhs) {
            if (get<0>(lhs) == get<0>(rhs))
                return get<1>(lhs) > get<1>(rhs);
            return get<0>(lhs) > get<0>(rhs);
        };

        priority_queue<tuple<uint64_t, uint64_t, kDataFrameIterator *, kDataFrameIterator *>, vector<tuple<uint64_t, uint64_t, kDataFrameIterator *, kDataFrameIterator *> >, decltype(compare)> nextKmer(
                compare);

        for (unsigned int i = 0; i < input.size(); i++) {
            kDataFrameIterator *it = new kDataFrameIterator(input[i]->begin());
            kDataFrameIterator *itend = new kDataFrameIterator(input[i]->end());
            nextKmer.push(make_tuple(it->getHashedKmer(), i, it, itend));
        }

        uint64_t processedKmers = 0;
        while (nextKmer.size() > 0) {
            vector<uint32_t> colorVec;
            colorVec.clear();
            uint64_t currHash = get<0>(nextKmer.top());
            processedKmers++;
            while (nextKmer.size() > 0 && get<0>(nextKmer.top()) == currHash) {
                auto colorTuple = nextKmer.top();
                nextKmer.pop();

                uint32_t i = get<1>(colorTuple);
                vector<uint32_t> tmp = input[i]->getKmerDefaultColumnValue<vector<uint32_t>, insertColorColumn>(
                        get<0>(colorTuple));
                for (auto c:tmp)
                    colorVec.push_back(c + idsOffset[i]);

                sort(colorVec.begin(), colorVec.end());
                get<2>(colorTuple)->next();
                if (*get<2>(colorTuple) != *get<3>(colorTuple)) {
                    get<0>(colorTuple) = get<2>(colorTuple)->getHashedKmer();
                    nextKmer.push(colorTuple);
                } else {
                    delete get<2>(colorTuple);
                    delete get<3>(colorTuple);
                }
            }
            output->setKmerDefaultColumnValue<vector<uint32_t>, insertColorColumn>(currHash, colorVec);

        }
        colors->populateColors();


    }

    vector<uint64_t> estimateKmersHistogram(string fileName, int kSize, int threads) {
        std::string tmpFile = "tmp." + to_string(rand() % 1000000);
        main_ntCard(fileName, kSize, 1000, threads, tmpFile);
        string output = tmpFile + "_k" + to_string(kSize) + ".hist";
        string tag;
        uint64_t count;
        vector<uint64_t> res(1000);
        ifstream resultFile(output);
        resultFile >> tag >> count;
        resultFile >> tag >> count;
        for (int i = 0; i < 1000; i++) {
            resultFile >> tag >> res[i];
        }
        resultFile.close();
        remove(output.c_str());
        return res;
    }

    void loadFromKMC(kDataFrame *kframe, std::string KMC_DB_filename) {
        uint32 _kmer_length;
        uint32 _mode;
        uint32 _counter_size;
        uint32 _lut_prefix_length;
        uint32 _signature_len;
        uint32 _min_count;
        uint64 _max_count;
        uint64 _total_kmers;
        CKMCFile kmer_data_base;
        if (!kmer_data_base.OpenForListing(KMC_DB_filename)) {
            throw std::logic_error("Cant open KMC DB");
            return;
        }
        kmer_data_base.Info(_kmer_length, _mode, _counter_size, _lut_prefix_length, _signature_len, _min_count,
                            _max_count, _total_kmers);
        CKmerAPI kmer_object(_kmer_length);
        uint32 counter;
        std::string str;

        kframe->setkSize(_kmer_length);
        kframe->reserve(_total_kmers);
        while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
            kmer_object.to_string(str);
            kframe->insert(str, counter);
        }

    }

    void createCountColumn(kDataFrame* frame){
        frame->addColumn("count",new vectorColumn<uint32_t>(frame->size()));
        for(auto k:*frame)
        {
            frame->setKmerColumnValue<uint32_t,vectorColumn<uint32_t> >("count",k.kmer,k.count);
        }
    }
    void createColorColumn(kDataFrame* frame){
        deduplicatedColumn<vector<uint32_t>,StringColorColumn>* newColumn=new deduplicatedColumn<vector<uint32_t>,StringColorColumn>(frame->size());
        newColumn->values=(StringColorColumn*)frame->getDefaultColumn();
        uint32_t i=0;
        for(auto k:*frame)
        {
            newColumn->index[i++]=k.count;
        }
        frame->addColumn("color",newColumn);

    }

} // End of namespace kProcessor
