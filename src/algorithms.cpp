#include "algorithms.hpp"
#include <iostream>
#include "Utils/kmer.h"
#include <limits>
#include <string>
#include <queue>
#include <functional>
#include "ntcard.hpp"
#include <cstdio>
#include <chrono>
#include "defaultColumn.hpp"
#include "kmc_file.h"

using namespace std::chrono;


using std::string;
using std::vector;
using std::cerr;
using std::cout;

using phmap::flat_hash_map;


namespace kProcessor {

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
        res->addCountColumn();
        kDataFrameIterator it = input->begin();
        while (it != input->end()) {
            kmerRow newkmer = fn(it.getKmerRow());
            res->insert(newkmer);
            it++;
        }
        // for(auto col: input->columns)
        // {
        //     string newColName= col.first;
        //     Column* column=col.second->getTwin();
        //     column->resize(res->size());
        //     res->addColumn(newColName, column);
        // }
        // for(auto kmer:(*res))
        // {
        //     for (auto col: input->columns) {
        //         string newColName = col.first;
        //         res->setKmerColumnValueFromOtherColumn(input,col.first, newColName,kmer.getKmer());
        //     }
        // }
        return res;

    }
    void transformInPlace(kDataFrame *input,  function<void (kDataFrameIterator& i)> fn)
    {
        kDataFrame *res = input;
        kDataFrameIterator it = input->begin();
        while (it != input->end()) {
            fn(it);
            it++;
        }
       return;

    }
    kDataFrame* innerJoin(vector<kDataFrame *> input, vector<uint32_t> kmersToKeep) {
        kDataFrame *res = input[0]->getTwin();
        uint64_t numKmers = 0;
        for (auto kframe:input) {
            numKmers += kframe->size();
        }
        res->reserve((uint64_t) ((double) numKmers * 0.75));
        merge(input, res, [&](vector<kDataFrameIterator*> &input) -> kmerRow {
            bool exists=false;
            for (auto i : kmersToKeep ) {
                if (input[i] != nullptr) {
                    return kmerRow("",input[i]->getHashedKmer(),1,0,nullptr);
                }
            }
            return kmerRow();
        });
        return res;
    }

    kDataFrame *filter(kDataFrame *input, function<bool (kmerRow& i)> fn) {
        kDataFrame *res = input->getTwin();
        kDataFrameIterator it = input->begin();
        while (it != input->end()) {
            kmerRow k=it.getKmerRow();
            if (fn(k))
                res->insert(k.hashedKmer);
            it++;
        }
        for(auto col: input->columns)
        {
            string newColName= col.first;
            Column* column=col.second->getTwin();
            column->resize(res->size());
            res->addColumn(newColName, column);
        }
        for(auto kmer:*res)
        {
            for (auto col: input->columns) {
                string newColName = col.first;
                res->setKmerColumnValueFromOtherColumn(input,col.first, newColName,kmer.getHashedKmer());
            }
        }

        return res;

    }
    kDataFrame *filter(kDataFrame *input, function<bool (kDataFrameIterator& i)> fn) {
        kDataFrame *res = input->getTwin();
        unordered_map<string,Column*> columns;
        for(auto col: input->columns)
        {
            columns[col.first]=col.second->getTwin();
            columns[col.first]->resize(input->size());
        }
        kDataFrameIterator it = input->begin();
        int index=0;
        while (it != input->end()) {
            if (fn(it)) {
                uint64_t hash=it.getHashedKmer();
                res->insert(it.getHashedKmer());
                for(auto col: input->columns)
                {
                    columns[col.first]->setValueFromColumn(col.second,it.getOrder(),res->getkmerOrder(hash));
                }
                index++;
            }
            it++;
        }
        for(auto col: columns)
        {
            col.second->resize(res->size());
            res->addColumn(col.first,col.second);
        }
        return res;
        


    }

    any aggregate(kDataFrame *input, any initial, function<any (kmerRow i, any v)> fn) {
        kDataFrameIterator it = input->begin();
        while (it != input->end()) {
            initial = fn(it.getKmerRow(), initial);
            it++;
        }
        return initial;
    }
    any aggregate(kDataFrame *input, any initial, function<any (kDataFrameIterator& i, any v)> fn) {
        kDataFrameIterator it = input->begin();
        while (it != input->end()) {
            initial = fn(it, initial);
            it++;
        }
        return initial;
    }

    void parseSequences(kmerDecoder *KD, kDataFrame *output) {
        output->addCountColumn();
        if (KD->get_kSize() != (int) output->getkSize()) {
            std::cerr << "kmerDecoder kSize must be equal to kDataFrame kSize" << std::endl;
            exit(1);
        }

        while (!KD->end()) {
            KD->next_chunk();
            for (const auto &seq : *KD->getKmers()) {
                for (const auto &kmer : seq.second) {
                    output->incrementCount(kmer.hash);
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
        kframe->addCountColumn();
        std::string mode = "kmers";
        bool check_mode = (parse_params.find("mode") != parse_params.end());
        if (check_mode) {
            if (parse_params["mode"] == 2) mode = "skipmers";
            else if (parse_params["mode"] == 3) mode = "minimizers";
        }
        kmerDecoder *KD = kmerDecoder::getInstance(filename, chunk_size, kframe->KD->slicing_mode, kframe->KD->hash_mode, {{"kSize", kframe->ksize()}});

        // Clone the hashing

        // kmerDecoder_setHashing(KD, kframe->KD->hash_mode);

        // Processing

        while (!KD->end()) {
            KD->next_chunk();
            for (const auto &seq : *KD->getKmers()) {
                for (const auto &kmer : seq.second) {
                    kframe->incrementCount(kmer.hash);                    
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
        output->addCountColumn();
        std::vector<kmer_row> kmers;
        KD->seq_to_kmers(sequence, kmers);

        for (const auto &kmer : kmers) {
                output->incrementCount(kmer.hash);
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

        kmerDecoder_setHashing(KD, frame->KD->hash_mode);

        if (KD->get_kSize() != (int) frame->getkSize()) {
            std::cerr << "kmerDecoder kSize must be equal to kDataFrame kSize" << std::endl;
            exit(1);
        }
        frame->addCountColumn();

        std::vector<kmer_row> kmers;
        KD->seq_to_kmers(sequence, kmers);

        for (const auto &kmer : kmers) {
                frame->incrementCount(kmer.hash);
        }

    }


    struct CustomKmerRow {
        bool operator()(const pair<kDataFrameIterator*, int> &lhs, const pair<kDataFrameIterator*, int> &rhs) {
            return lhs.first->getHashedKmer() > rhs.first->getHashedKmer();
        }
    };

    inline void terminate_if_kDataFrameMAP(const vector<kDataFrame *> &input) {
        for (const auto &kFrame : input) {
            if (kFrame->get_class_name() == "MAP") {
                throw logic_error("can't apply set functions kDataFramMAP on kDataFrameMAP.");
            }
        }
    }

    void merge(const vector<kDataFrame *> &input, kDataFrame *res,function<kmerRow (vector<kDataFrameIterator*> &i)> fn) {
        terminate_if_kDataFrameMAP(input);
        unordered_map<string,Column*> columns;
        priority_queue<pair<kDataFrameIterator*, int>, vector<pair<kDataFrameIterator*, int> >, CustomKmerRow> Q;
        vector<kDataFrameIterator> iterators(input.size());
        for (unsigned int i = 0; i < input.size(); i++) {
            for(auto col: input[i]->columns)
            {
                string newColumnName=col.first+to_string(i);
                columns[newColumnName]=col.second->getTwin();
                columns[newColumnName]->resize(res->size());
            }
            iterators[i] = input[i]->begin();
            if (iterators[i] != input[i]->end()) {
                //  cout<<i<<" "<<(*iterators[i]).hashedKmer<<endl;
                //auto tmp=iterators[i].getKmerRow();
                Q.push(make_pair(&iterators[i], i));
            }
        }
        vector<kDataFrameIterator* > current(input.size());
        uint64_t index=0;
        while (Q.size() > 0) {
            for (unsigned int i = 0; i < current.size(); i++)
                current[i] = nullptr;
            pair<kDataFrameIterator*, int> top = Q.top();
            uint64_t topHash=top.first->getHashedKmer();
            Q.pop();
//            iterators[top.second]++;
            current[top.second] = top.first;
            //  cout<<top.first.hashedKmer<<" "<<top.second<<endl;
//            if (iterators[top.second] != input[top.second]->end()) {
//                Q.push(make_pair(&iterators[top.second], top.second));
//            }

            while (!Q.empty() && topHash == Q.top().first->getHashedKmer()) {
                top = Q.top();
                Q.pop();
                current[top.second] = top.first;
                //iterators[top.second]++;
//                if (iterators[top.second] != input[top.second]->end()) {
//                    Q.push(make_pair(&iterators[top.second], top.second));
//                }
            }

            kmerRow newRow = fn(current);

            if(newRow.count > 0) {
                res->insert(newRow.hashedKmer);
                for(uint32_t i=0;i<current.size();i++)
                {
                    if(current[i]!= nullptr)
                    {
                        for(auto col: input[i]->columns)
                        {
                            string newColumnName=col.first+to_string(i);
                            if(index>=columns[newColumnName]->size())
                            {
                                columns[newColumnName]->resize(res->size());
                            }
                            columns[newColumnName]->setValueFromColumn( col.second,iterators[i].getOrder(),index);
                        }
                    }
                }
                index++;
            }
            for (unsigned int i = 0; i < current.size(); i++)
            {
                if(current[i] != nullptr) {
                    (*current[i])++;
                    if ((*current[i]) != input[i]->end())
                        Q.push(make_pair(current[i], i));

                }
            }


        }
        for(auto col: columns)
        {
            col.second->resize(res->size());
            res->addColumn(col.first,col.second);
        }

    }

    kDataFrame *kFrameUnion(const vector<kDataFrame *> &input) {
        vector<uint32_t> indexes(input.size());
        for(int i=0;i<input.size();i++)
            indexes[i]=i;
        return innerJoin(input,indexes);
    }

    kDataFrame *kFrameIntersect(const vector<kDataFrame *> &input) {
        kDataFrame *res = input[0]->getTwin();

        merge(input, res, [&](vector<kDataFrameIterator*> &input) -> kmerRow {
            bool existInAll=true;
            for(int i=0;i<input.size();i++) {
                if (input[i] == nullptr) {
                    existInAll = false;
                    break;
                }
            }
            if(existInAll)
                return kmerRow("",input[0]->getHashedKmer(),1,0,nullptr);
            return kmerRow();
        });
        return res;

    }

    kDataFrame *kFrameDiff(const vector<kDataFrame *> &input) {
        kDataFrame *res = input[0]->getTwin();

        merge(input, res, [&](vector<kDataFrameIterator*> &input) -> kmerRow {
            if(input[0] != nullptr){
                bool exist=false;
                for(int i=1;i<input.size();i++) {
                    if(input[i]!= nullptr) {
                        exist=true;
                        break;
                    }
                }
                if(!exist)
                    return kmerRow("",input[0]->getHashedKmer(),1,0,nullptr);
            }
            return kmerRow();
        });
        return res;

    }


    void kmerDecoder_setHashing(kmerDecoder * KD, hashingModes hash_mode){
        KD->setHashingMode(hash_mode, KD->get_kSize());
    }

    void kmerDecoder_setHashing(kDataFrame * KF, hashingModes hash_mode){
        KF->KD->setHashingMode(hash_mode, KF->ksize());
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

    kmerDecoder* initialize_kmerDecoder(int kSize, hashingModes HM = integer_hasher){
        return new Kmers(kSize, HM);
    }

    void index(kmerDecoder *KD, string names_fileName, kDataFrame *frame) {

        if (KD->get_kSize() != (int) frame->ksize()) {
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
                    uint64_t currentTag = frame->getkmerOrder(kmer.hash);
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
                            if (freeColors.empty()) {
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
                    if (frame->getkmerOrder(kmer.hash) != itc->second) {
                        //frame->setC(kmer,itc->second);
                        cout << "Error Founded " << kmer.str << " from sequence " << readName << " expected "
                             << itc->second << " found " << frame->getkmerOrder(kmer.hash) << endl;
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
        auto *col = new StringColorColumn();
        frame->addColumn("color",(Column *) col);
        col->colors = vector<vector<uint32_t> >(legend->size());
        col->colors.push_back(vector<uint32_t>());
        for (auto it : *legend) {
            col->colors[it.first] = it.second;
        }
        delete legend;

        for (auto & iit : namesMap) {
            uint32_t sampleID = groupNameMap[iit.second];
            col->namesMap[sampleID] = iit.second;
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

        kmerDecoder_setHashing(KD, frame->KD->hash_mode);

        // Processing

        if (KD->get_kSize() != (int) frame->ksize()) {
            std::cerr << "kmerDecoder kSize must be equal to kDataFrame kSize" << std::endl;
            exit(1);
        }
        index(KD, names_fileName, frame);


    }

    void indexPriorityQueue(vector<kDataFrame *> &input, string tmpFolder, kDataFrame *output,uint32_t num_vectors,uint32_t vector_size) {
        auto *colors =new deduplicatedColumn<vector<uint32_t>, insertColorColumn>(); 
        colors->values= new insertColorColumn(input.size(), tmpFolder,num_vectors,vector_size);
        
        output->addColumn("i",colors);
        for (unsigned int i = 0; i < input.size(); i++) {
            vector<uint32_t> tmp = {i};
            colors->values->insertAndGetIndex(tmp);
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
            auto *it = new kDataFrameIterator(input[i]->begin());
            auto *itend = new kDataFrameIterator(input[i]->end());
            nextKmer.push(make_tuple(it->getHashedKmer(), i, it, itend));
        }

        uint64_t processedKmers = 0;
        while (!nextKmer.empty()) {
            vector<uint32_t> colorVec;
            colorVec.clear();
            uint64_t currHash = get<0>(nextKmer.top());
            processedKmers++;
            if (processedKmers % 1000000 == 0)
                cout << processedKmers << " Kmers Processed" << endl;
            while (!nextKmer.empty() && get<0>(nextKmer.top()) == currHash) {
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

            uint64_t prevColor = output->getkmerOrder(currHash);
            if (prevColor != 0) {
                cout << "Error in Indexing detected at kmer " << currHash << endl;
                cout << "should be empty vector and found  " << endl;
            }
            
            output->insert(currHash);
            output->setKmerColumnValue<vector<uint32_t> , deduplicatedColumn<vector<uint32_t>, insertColorColumn>>("i",currHash, colorVec);

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
        colors->values->populateColors();
        uint64_t noColors = colors->values->noColors;
        cout << noColors << " colors created" << endl;




        auto colorColumn= new deduplicatedColumn<vector<uint32_t>, mixVectors>();
        colorColumn->values=new mixVectors(colors->values);
        colorColumn->values->explainSize();
        colorColumn->index=colors->index;
        output->removeColumn("color");
        output->addColumn("color",colorColumn);
    }

    void mergeIndexes(vector<kDataFrame *> &input, string tmpFolder, kDataFrame *output) {

        // typedef deduplicatedColumn<vector<uint32_t>, mixVectors> colorColumnType;
        // vector<uint32_t> idsOffset(input.size());
        // idsOffset[0] = 0;
        // uint32_t noSamples=0;
        // for (unsigned int i = 1; i < input.size(); i++) {
        //     idsOffset[i] = idsOffset[i - 1];
        //     idsOffset[i] += ((colorColumnType *) input[i - 1]->columns["color"])->values->noSamples;
        // }
        // noSamples+=((colorColumnType *) input[input.size() - 1]->columns["color"])->values->noSamples;

        // auto *colors = new insertColorColumn(noSamples,tmpFolder);
        // output->addColumn("i",colors);


        // auto compare = [](tuple<uint64_t, uint64_t, kDataFrameIterator *, kDataFrameIterator *> lhs,
        //                   tuple<uint64_t, uint64_t, kDataFrameIterator *, kDataFrameIterator *> rhs) {
        //     if (get<0>(lhs) == get<0>(rhs))
        //         return get<1>(lhs) > get<1>(rhs);
        //     return get<0>(lhs) > get<0>(rhs);
        // };

        // priority_queue<tuple<uint64_t, uint64_t, kDataFrameIterator *, kDataFrameIterator *>, vector<tuple<uint64_t, uint64_t, kDataFrameIterator *, kDataFrameIterator *> >, decltype(compare)> nextKmer(
        //         compare);

        // for (unsigned int i = 0; i < input.size(); i++) {
        //     auto *it = new kDataFrameIterator(input[i]->begin());
        //     auto *itend = new kDataFrameIterator(input[i]->end());
        //     nextKmer.push(make_tuple(it->getHashedKmer(), i, it, itend));
        // }

        // uint64_t processedKmers = 0;
        // while (!nextKmer.empty()) {
        //     vector<uint32_t> colorVec;
        //     colorVec.clear();
        //     uint64_t currHash = get<0>(nextKmer.top());
        //     processedKmers++;
        //     while (!nextKmer.empty() && get<0>(nextKmer.top()) == currHash) {
        //         auto colorTuple = nextKmer.top();
        //         nextKmer.pop();

        //         uint32_t i = get<1>(colorTuple);
        //         auto tmp = input[i]->getKmerColumnValue<vector<uint32_t >, deduplicatedColumn<vector<uint32_t>, mixVectors> >("color",get<0>(colorTuple));
        //         for (auto c:tmp)
        //             colorVec.push_back(c + idsOffset[i]);

        //         sort(colorVec.begin(), colorVec.end());
        //         get<2>(colorTuple)->next();
        //         if (*get<2>(colorTuple) != *get<3>(colorTuple)) {
        //             get<0>(colorTuple) = get<2>(colorTuple)->getHashedKmer();
        //             nextKmer.push(colorTuple);
        //         } else {
        //             delete get<2>(colorTuple);
        //             delete get<3>(colorTuple);
        //         }
        //     }
        //     output->setKmerColumnValue<vector<uint32_t>, insertColorColumn>("i",currHash, colorVec);

        // }
        // colors->populateColors();
        // uint64_t noColors = colors->noColors;
        // cout << noColors << " colors created!" << endl;


        // auto colorColumn= new deduplicatedColumn<vector<uint32_t>, mixVectors>();
        // colorColumn->values=new mixVectors(colors);
        // output->addColumn("color",colorColumn);
        // colorColumn->values->explainSize();
        // colorColumn->index=vector<uint32_t>(output->size());
        // for(auto k:(*output))
        // {
        //     colorColumn->index[k.getOrder()]=k.getCount();
        // }
        
        // delete colors;

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
        kframe->addCountColumn();
        while (kmer_data_base.ReadNextKmer(kmer_object, counter)) {
            kmer_object.to_string(str);
            kframe->setCount(str, counter);
        }

    }

    

} // End of namespace kProcessor
