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

                for (int j = 0; j < reads.size(); j++) {
                    read = reads[j].first;
                    start_read:
                    if (read.size() < ksize) {
                        continue;
                    }

                    uint64_t first = 0;
                    uint64_t first_rev = 0;
                    uint64_t item = 0;

                    for (int i = 0; i < ksize; i++) {
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
        uint64_t noDistinctKmers = 0, totalNumKmers = 0;
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

    kDataFrame *transform(kDataFrame *input, kmerRow (*fn)(kmerRow it)) {
        kDataFrame *res = input->getTwin();
        kDataFrameIterator it = input->begin();
        while (it != input->end()) {
            kmerRow newkmer = fn(*it);
            res->insert(newkmer);
            it++;
        }
        return res;

    }
    kDataFrame *filter(kDataFrame *input, bool (*fn)(kmerRow it)) {
      kDataFrame *res = input->getTwin();
      kDataFrameIterator it = input->begin();
      while (it != input->end()) {
        if(fn(*it))
          res->insert(*it);
        it++;
      }
      return res;

    }
    any aggregate(kDataFrame *input, any initial ,any (*fn)(kmerRow it,any v)) {
      kDataFrameIterator it = input->begin();
      while (it != input->end()) {
        initial=fn(*it,initial);
        it++;
      }
      return initial;
    }

    void parseSequences(kmerDecoder * KD, kDataFrame* output){
        if (KD->get_kSize() != (int)output->getkSize()){
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


    void countKmersFromFile(kDataFrame *kframe, string filename, int chunk_size) {

        // kframe->reserve(100000);
        // Make a new kmerDecoder to pass the filename, but get all the params from the kDataFrame KD.
        kmerDecoder * KD = new Kmers(filename, chunk_size, kframe->ksize(), kframe->KD->hash_mode);

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

    void countKmersFromString(string sequence, kDataFrame *output) {

        kmerDecoder * KD = kmerDecoder::getInstance(output->KD->slicing_mode, output->KD->hash_mode, output->KD->string_to_params(output->KD->params_to_string()));

        std::vector<kmer_row> kmers;
        KD->seq_to_kmers(sequence, kmers);

        for (const auto &kmer : kmers) {
            output->insert(kmer.hash);
        }

        delete KD;

    }


    struct CustomKmerRow {
        bool operator()(const pair<kmerRow, int> &lhs, const pair<kmerRow, int> &rhs) {
            return lhs.first.hashedKmer > rhs.first.hashedKmer;
        }
    };

    inline void terminate_if_kDataFrameMAP(const vector<kDataFrame *> &input){
        for(const auto & kFrame : input){
            if (kFrame->get_class_name() == "MAP"){
                throw logic_error("can't apply set functions kDataFramMAP on kDataFrameMAP.");
            }
        }
    }

    void merge(const vector<kDataFrame *> &input, kDataFrame *res, kmerRow (*fn)(vector<kmerRow> &i)) {
        terminate_if_kDataFrameMAP(input);
        priority_queue<pair<kmerRow, int>, vector<pair<kmerRow, int> >, CustomKmerRow> Q;
        vector<kDataFrameIterator> iterators(input.size());
        for (int i = 0; i < input.size(); i++) {
            iterators[i] = input[i]->begin();
            if (iterators[i] != input[i]->end()) {
                //  cout<<i<<" "<<(*iterators[i]).hashedKmer<<endl;
                Q.push(make_pair(*iterators[i], i));
            }
        }
        vector<kmerRow> current(input.size());
        while (Q.size() > 0) {
            for (int i = 0; i < current.size(); i++)
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
            for (int i = 0; i < input.size(); i++) {
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
            for (int i = 1; i < input.size(); i++) {
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
            for (int i = 1; i < input.size(); i++) {
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

    
    void kmerDecoder_setHashing(kmerDecoder * KD, hashingModes hash_mode){
        KD->setHashingMode(hash_mode, KD->get_kSize());
    }

    void kmerDecoder_setHashing(kDataFrame * KF, hashingModes hash_mode){
        KF->KD->setHashingMode(hash_mode, KF->ksize());
    }


    kmerDecoder* initialize_kmerDecoder(std::string filename, int chunkSize, std::string mode, std::map<std::string, int> parse_params){

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
        } else if(mode == "protein"){
            if (parse_params.find("k_size") != parse_params.end()) {
                return new aaKmers(filename, chunkSize, parse_params["k_size"]);
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
                if(check_orf) return new Skipmers(filename, chunkSize, parse_params["m"], parse_params["n"], parse_params["k_size"], parse_params["orf"]);
                return new Skipmers(filename, chunkSize, parse_params["m"], parse_params["n"], parse_params["k_size"]);
            } else {
                std::cerr << func_name << "kmerDecoder Skipmers parameters {k_size, m, n} validation failed" << std::endl;
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

        }else{
            std::cerr << func_name << "supported kmerDecoder modes: {kmers, skipmers, minimizers}" << std::endl;
            exit(1);
        }
    }

    kmerDecoder* initialize_kmerDecoder(std::string mode, std::map<std::string, int> parse_params){

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
                if(check_orf) return new Skipmers(parse_params["m"], parse_params["n"], parse_params["k_size"], parse_params["orf"]);
                return new Skipmers(parse_params["m"], parse_params["n"], parse_params["k_size"]);
            } else {
                std::cerr << func_name << "kmerDecoder Skipmers parameters {k_size, m, n} validation failed" << std::endl;
                exit(1);
            }
        } else if (mode == "minimizers") {
            bool check_k = (parse_params.find("k_size") != parse_params.end());
            bool check_w = (parse_params.find("w") != parse_params.end());

            if (check_k && check_w) {
                return new Minimizers(parse_params["k_size"], parse_params["w"]);
            } else {
                std::cerr  << func_name << "kmerDecoder Skipmers parameters {k_size, w} validation failed" << std::endl;
                exit(1);
            }

        }else{
            std::cerr << func_name << "supported kmerDecoder modes: {kmers, skipmers, minimizers}" << std::endl;
            exit(1);
        }
    }

    kmerDecoder* initialize_kmerDecoder(int kSize, hashingModes HM = integer_hasher){
        return new Kmers(kSize, HM);
    }

    colored_kDataFrame *index(kmerDecoder *KD, string names_fileName, kDataFrame *frame) {

        if (KD->get_kSize() != (int)frame->ksize() && (KD->hash_mode != protein_hasher || KD->hash_mode != proteinDayhoff_hasher)) {
            std::cerr << "WARNING: kmerDecoder kSize must be equal to kDataFrame kSize" << std::endl;
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
            while(std::getline(iss, token, '\t'))   // but we can specify a different one
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
        int currIndex = 0;
        string kmer;
        uint64_t tagBits = 0;
        uint64_t maxTagValue = (1ULL << tagBits) - 1;
       //  kDataFrame *frame;
        int kSize = KD->get_kSize();


        uint64_t lastTag = 0;
        readID = 0;

        while (!KD->end()) {
            KD->next_chunk();

            flat_hash_map<uint64_t, uint64_t> convertMap;

            for (const auto &seq : *KD->getKmers()) {
                string readName = seq.first;

                auto it = namesMap.find(readName);
                if (it == namesMap.end()) {
                    continue;
                    // cout << "read " << readName << "dont have group. Please check the group names file." << endl;
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
                        for (int k = 1; k < colors.size(); k++) {
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
                        return NULL;
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
        colorTable *colors = new intVectorsTable();
        for (auto it : *legend) {
            colors->setColor(it.first, it.second);
        }

        colored_kDataFrame *res = new colored_kDataFrame();
        res->setColorTable(colors);
        res->setkDataFrame(frame);
        for (auto iit = namesMap.begin(); iit != namesMap.end(); iit++) {
            uint32_t sampleID = groupNameMap[iit->second];
            res->namesMap[sampleID] = iit->second;
            res->namesMapInv[iit->second] = sampleID;
        }
        return res;
    }

    colored_kDataFrame * index(kDataFrame * frame, std::map<std::string, int> parse_params, string filename, int chunk_size, string names_fileName){

        // parse_params["mode"] = 1 > Default: Kmers
        // parse_params["mode"] = 2 > Skipmers
        // parse_params["mode"] = 3 > Minimizers

        // Initialize kmerDecoder
        std::string mode = "kmers";
        bool check_mode = (parse_params.find("mode") != parse_params.end());
        if (check_mode){
            if (parse_params["mode"] == 2) mode = "skipmers";
            else if (parse_params["mode"] == 3) mode = "minimizers";
            else if (parse_params["mode"] == 4) mode = "protein";
        }

        parse_params["k_size"] = frame->ksize();
        parse_params["k"] = frame->ksize();
        kmerDecoder * KD = initialize_kmerDecoder(filename, chunk_size, mode, parse_params);

        // Clone the hashing

        kmerDecoder_setHashing(KD, frame->KD->hash_mode);

        // Processing

            if (KD->get_kSize() != (int)frame->ksize()) {
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
                while(std::getline(iss, token, '\t'))   // but we can specify a different one
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
            int currIndex = 0;
            string kmer;
            uint64_t tagBits = 0;
            uint64_t maxTagValue = (1ULL << tagBits) - 1;
            //  kDataFrame *frame;
            int kSize = KD->get_kSize();


            uint64_t lastTag = 0;
            readID = 0;
            int __batch_count = 0;
            while (!KD->end()) {
                KD->next_chunk();
                cout << "Processing Chunk(" << ++__batch_count << "): " << chunk_size*__batch_count << " seqs ..." << endl;
                flat_hash_map<uint64_t, uint64_t> convertMap;

                for (const auto &seq : *KD->getKmers()) {
                    string readName = seq.first;

                    auto it = namesMap.find(readName);
                    if (it == namesMap.end()) {
                        continue;
                        // cout << "read " << readName << "dont have group. Please check the group names file." << endl;
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
                            for (int k = 1; k < colors.size(); k++) {
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
                            return NULL;
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
            colorTable *colors = new intVectorsTable();
            for (auto it : *legend) {
                colors->setColor(it.first, it.second);
            }

            colored_kDataFrame *res = new colored_kDataFrame();
            res->setColorTable(colors);
            res->setkDataFrame(frame);
            for (auto iit = namesMap.begin(); iit != namesMap.end(); iit++) {
                uint32_t sampleID = groupNameMap[iit->second];
                res->namesMap[sampleID] = iit->second;
                res->namesMapInv[iit->second] = sampleID;
            }
            return res;
        }

    colored_kDataFrame * index(kDataFrame * frame, string filename, int chunk_size, string names_fileName){

        // parse_params["mode"] = 1 > Default: Kmers
        // parse_params["mode"] = 2 > Skipmers
        // parse_params["mode"] = 3 > Minimizers

        // Initialize kmerDecoder
        map<string, int> params = kmerDecoder::string_to_params(frame->KD->params_to_string());

        kmerDecoder * KD = kmerDecoder::getInstance(filename, chunk_size, frame->KD->slicing_mode, frame->KD->hash_mode, params);

        // Processing

            if (KD->get_kSize() != (int)frame->ksize()) {
                std::cerr << "kmerDecoder kSize must be equal to kDataFrame kSize" << std::endl;
                exit(1);
            }

            flat_hash_map<string, string> namesMap;
            flat_hash_map<string, uint64_t> tagsMap;
            flat_hash_map<string, uint64_t> groupNameMap;
            auto *legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
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
                while(std::getline(iss, token, '\t'))   // but we can specify a different one
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
            int currIndex = 0;
            string kmer;
            uint64_t tagBits = 0;
            uint64_t maxTagValue = (1ULL << tagBits) - 1;
            //  kDataFrame *frame;
            int kSize = KD->get_kSize();


            uint64_t lastTag = 0;
            readID = 0;
            int __batch_count = 0;
            while (!KD->end()) {
                KD->next_chunk();
                cout << "Processing Chunk(" << ++__batch_count << "): " << chunk_size*__batch_count << " seqs ..." << endl;
                flat_hash_map<uint64_t, uint64_t> convertMap;

                for (const auto &seq : *KD->getKmers()) {
                    string readName = seq.first;

                    auto it = namesMap.find(readName);
                    if (it == namesMap.end()) {
                        continue;
                        // cout << "read " << readName << "dont have group. Please check the group names file." << endl;
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
                            for (int k = 1; k < colors.size(); k++) {
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
                            return NULL;
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
            colorTable *colors = new intVectorsTable();
            for (auto it : *legend) {
                colors->setColor(it.first, it.second);
            }

            colored_kDataFrame *res = new colored_kDataFrame();
            res->setColorTable(colors);
            res->setkDataFrame(frame);
            for (auto iit = namesMap.begin(); iit != namesMap.end(); iit++) {
                uint32_t sampleID = groupNameMap[iit->second];
                res->namesMap[sampleID] = iit->second;
                res->namesMapInv[iit->second] = sampleID;
            }
            return res;
        }


        vector<uint64_t> estimateKmersHistogram(string fileName, int kSize ,int threads)
{
   std::string tmpFile = "tmp."+to_string(rand()%1000000);
   main_ntCard(fileName,kSize,1000,threads,tmpFile);
   string output=tmpFile+"_k"+to_string(kSize)+".hist";
   string tag;
   uint64_t count;
   vector<uint64_t> res(1000);
   ifstream resultFile(output);
   resultFile>>tag>>count;
   resultFile>>tag>>count;
   for(int i=0;i<1000;i++)
   {
     resultFile>>tag>>res[i];
   }
   resultFile.close();
   remove(output.c_str());
   return res;
}

} // End of namespace kProcessor
