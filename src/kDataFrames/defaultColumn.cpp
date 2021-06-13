#include <defaultColumn.hpp>
#include<iostream>
#include<fstream>
#include <cereal/types/vector.hpp>
#include <cereal/types/memory.hpp>
#include "parallel_hashmap/phmap_dump.h"
#include <cereal/archives/binary.hpp>
#include <stack>
#include <queue>
#include <iterator>
#include <regex>
#include <sdsl/util.hpp>
#include "mum.h"
#include <unistd.h>

template
class vectorColumn<int>;

template
class vectorColumn<bool>;

template
class vectorColumn<double>;


template
class vectorColumn<uint32_t>;

template
class deduplicatedColumn<vector<uint32_t>, StringColorColumn>;

bool is_file_exist(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

Column *Column::getContainerByName(std::size_t hash) {
    if (hash == typeid(vectorColumn<int>).hash_code()) {
        return new vectorColumn<int>();
    } else if (hash == typeid(vectorColumn<double>).hash_code()) {
        return new vectorColumn<double>();
    } else if (hash == typeid(vectorColumn<bool>).hash_code()) {
        return new vectorColumn<bool>();
    } else if (hash == typeid(insertColorColumn).hash_code()) {
        return new insertColorColumn();
    } else if (hash == typeid(StringColorColumn).hash_code()) {
        return new StringColorColumn();
    } else if (hash == typeid(mixVectors).hash_code()) {
        return new mixVectors();
    } else if (hash == typeid(prefixTrie).hash_code()) {
        return new prefixTrie();
    } else {
        throw logic_error("Failed to load Unknown Column " + hash);
    }
}

template<typename T>
uint32_t vectorColumn<T>::insertAndGetIndex(T item) {
    dataV.push_back(item);
    return dataV.size() - 1;
}


template<typename T>
T vectorColumn<T>::getWithIndex(uint32_t index) {
    return dataV[index];
}

template<typename T>
void vectorColumn<T>::serialize(string filename) {
//    ofstream wf(filename, ios::out | ios::binary);
//    uint32_t  size=dataV.size();
//    wf.write((char*)&size, sizeof(uint32_t));
//    for(auto t:dataV)
//        wf.write((char*)&t, sizeof(T));
//    wf.close();
    std::ofstream os(filename, std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(dataV);
    os.close();
}


template<typename T>
void vectorColumn<T>::deserialize(string filename) {
//    ifstream rf(filename, ios::out | ios::binary);
//    uint32_t  size;
//    rf.read((char*)&size, sizeof(uint32_t));
//    dataV=vector<T>(size);
//    T item;
//    for(int i=0;i<size;i++) {
//        rf.read((char *) &(item), sizeof(T));
//        dataV[i] = item;
//    }
//    rf.close();
    std::ifstream os(filename, std::ios::binary);
    cereal::BinaryInputArchive iarchive(os);
    iarchive(dataV);
    os.close();
}

template<typename T>
void vectorColumn<T>::insert(T item, uint32_t index) {
    dataV[index] = item;
}

template<typename T>
T vectorColumn<T>::get(uint32_t index) {
    return dataV[index];
}

template<typename T>
Column *vectorColumn<T>::getTwin() {
    return new vectorColumn<T>();
}



template<typename T>
void vectorColumn<T>::setValueFromColumn(Column *Container, uint32_t inputOrder, uint32_t outputOrder) {
    T val = ((vectorColumn<T> *) Container)->get(inputOrder);
    insert(val, outputOrder);
}


uint32_t insertColorColumn::insertAndGetIndex(vector<uint32_t> &item) {
    //return colorInv.getColorID(item);
    if (colorInv.hasColorID(item))
        return colorInv.getColorID(item);

    noColors++;
    uint32_t c = colorInv.getColorID(item);
    uint32_t colorSize = item.size();
    uint32_t maxSizeForNextColor;
    if (colorSize < NUM_VECTORS - 1) {
        // array of arrays fixed size
        colors[colorSize][colorsTop[colorSize]++] = c;
        for (auto s:item)
            colors[colorSize][colorsTop[colorSize]++] = s;

        maxSizeForNextColor = colorsTop[colorSize] + 1 + colorSize;

    } else {
        // array of arrays mixed sizes
        colorSize = NUM_VECTORS - 1;
        colors[colorSize][colorsTop[colorSize]++] = c;
        colors[colorSize][colorsTop[colorSize]++] = item.size();
        for (auto s:item)
            colors[colorSize][colorsTop[colorSize]++] = s;

        maxSizeForNextColor = colorsTop[colorSize] + 2 + noSamples;
    }
    if (maxSizeForNextColor > VECTOR_SIZE) {
        if (colorsTop[colorSize] != VECTOR_SIZE)
            colors[colorSize].resize(colorsTop[colorSize]);

        string colorsFileName = tmpFolder + "insertOnlyColumn." + to_string(colorSize) +
                                "." + to_string(vecCount[colorSize]++);
        sdsl::store_to_file(colors[colorSize], colorsFileName.c_str());
        if (colors[colorSize].size() != VECTOR_SIZE)
            colors[colorSize].resize(VECTOR_SIZE);
        colorsTop[colorSize] = 0;

    }
    return c;
}

/*
 * * */
vector<uint32_t> insertColorColumn::getWithIndex(uint32_t index) {
    throw logic_error("it is insert only color column");
//    return colors[index];
}


void insertColorColumn::serialize(string filename) {
    colorInv.serialize(filename);
}

void insertColorColumn::deserialize(string filename) {
    colorInv.deserialize(filename);
    populateColors();
    noSamples = colorInv.noSamples;
}

void insertColorColumn::cleanFiles() {
    for (uint32_t colorSize = 1; colorSize < NUM_VECTORS; colorSize++) {

        string colorsFileName = tmpFolder + "insertOnlyColumn." + to_string(colorSize) + "." +
                                to_string(vecCount[colorSize]);

        unlink(colorsFileName.c_str());

    }
//    colorInv.populateColors(colors);
}
void insertColorColumn::populateColors() {
    for (uint32_t colorSize = 1; colorSize < NUM_VECTORS; colorSize++) {
        if (colorsTop[colorSize] != VECTOR_SIZE)
            colors[colorSize].resize(colorsTop[colorSize]);

        string colorsFileName = tmpFolder + "insertOnlyColumn." + to_string(colorSize) + "." +
                                to_string(vecCount[colorSize]);

        sdsl::store_to_file(colors[colorSize], colorsFileName.c_str());

    }
//    colorInv.populateColors(colors);
}

Column *insertColorColumn::getTwin() {
    return new insertColorColumn();
}


void insertColorColumn::resize(uint32_t size) {

}

vector<string> StringColorColumn::getWithIndex(uint32_t index) {
    vector<string> res(colors[index].size());
    for (unsigned int i = 0; i < colors[index].size(); i++)
        res[i] = namesMap[colors[index][i]];
    return res;

}

vector<uint32_t> StringColorColumn::get(uint32_t index) {
    return colors[index];
}

void StringColorColumn::serialize(string filename) {
    std::ofstream os(filename + ".colors", std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(colors);
    os.close();

    ofstream namesMapOut(filename + ".namesMap");
    namesMapOut << namesMap.size() << endl;
    for (auto it:namesMap) {
        namesMapOut << it.first << " " << it.second << endl;
    }
    namesMapOut.close();
}


void StringColorColumn::deserialize(string filename) {
    std::ifstream os(filename + ".colors", std::ios::binary);
    cereal::BinaryInputArchive iarchive(os);
    iarchive(colors);

    ifstream namesMapIn(filename + ".namesMap");
    uint64_t size;
    namesMapIn >> size;
    for (unsigned int i = 0; i < size; i++) {
        uint32_t color;
        string name;
        namesMapIn >> color >> name;
        namesMap[color] = name;

    }
    namesMapIn.close();

}

Column *StringColorColumn::getTwin() {
    new StringColorColumn();
}


void StringColorColumn::resize(uint32_t size) {

}

bool inExactColorIndex::hasColorID(vector<uint32_t> &v) {
    uint64_t hash = mum_hash(v.data(), v.size() * 4, 4495203915657755407);
    if (colors.find(hash) == colors.end()) {
        return false;
    }
    return true;
}

uint32_t inExactColorIndex::getColorID(vector<uint32_t> &v) {
    uint64_t hash = mum_hash(v.data(), v.size() * 4, 4495203915657755407);
    auto it = colors.find(hash);
    if (it != colors.end()) {
        return it->second;
    }
    colors[hash] = ++lastColor;
    return colors[hash];


}


void mixVectors::insert(vector<uint32_t> &item, uint32_t index) {
    auto it = lower_bound(colors.begin(), colors.end(), index + 1,
                          [](vectorBase *lhs, uint32_t rhs) -> bool { return lhs->beginID < rhs; });
    // if(it==colors.end())
    it--;
    (*it)->set(index - (*it)->beginID, item);


}

uint32_t queryColorColumn::insertAndGetIndex(vector<uint32_t> &item) {
    throw std::logic_error("insertAndGetIndex is not supported in mixVectors");
    return 0;

}

vector<uint32_t> mixVectors::getWithIndex(uint32_t index) {
    if (index == 0)
        return vector<uint32_t>();
    index = idsMap[index];
    vector<uint32_t> res;
    stack<uint32_t> q;
    auto it = lower_bound(colors.begin(), colors.end(), index + 1,
                          [](vectorBase *lhs, uint32_t rhs) -> bool { return lhs->beginID < rhs; });
    //  if(it==colors.end())
    it--;
    return (*it)->get(index - (*it)->beginID);
    vector<uint32_t> compressedColor = (*it)->get(index - (*it)->beginID);
    return compressedColor;
    for (unsigned int i = 0; i < compressedColor.size(); i++) {
        if (compressedColor[i] >= noSamples) {
            q.push(compressedColor[i]);
        } else {
            res.push_back(compressedColor[i]);
        }
    }
    while (q.size() > 0) {
        index = q.top();
        // cout<<index<<endl;
        it = lower_bound(colors.begin(), colors.end(), index + 1,
                         [](vectorBase *lhs, uint32_t rhs) -> bool { return lhs->beginID < rhs; });
        //      if(it==colors.end())
        it--;
        vector<uint32_t> compressedColor = (*it)->get(index - (*it)->beginID);
        q.pop();
        for (unsigned int i = 0; i < compressedColor.size(); i++) {
            if (compressedColor[i] >= noSamples) {
                q.push(compressedColor[i]);
            } else {
                res.push_back(compressedColor[i]);
            }
        }

    }
    return res;
}

mixVectors::mixVectors(insertColorColumn *col) {
    col->populateColors();
    colors.push_back(new vectorOfVectors(0, 1));
    this->noSamples = col->noSamples;
    uint32_t colorId = 1;
    numColors = col->noColors;
    idsMap = sdsl::int_vector<>(numColors + 1);
    for (int colorSize = 1; colorSize < NUM_VECTORS - 1; colorSize++) {
        int chunkNum = 0;
        string colorsFileName = col->tmpFolder + "insertOnlyColumn." + to_string(colorSize) + "." +
                                to_string(chunkNum++);
        while (is_file_exist(colorsFileName.c_str())) {
            fixedSizeVector *f = new fixedSizeVector(colorId, colorSize);
            f->loadFromInsertOnly(colorsFileName, idsMap);
            colors.push_back(f);
            colorId += f->size();
            colorsFileName = col->tmpFolder + "insertOnlyColumn." + to_string(colorSize) + "." +
                             to_string(chunkNum++);
        }
    }

    int chunkNum = 0;
    string colorsFileName = col->tmpFolder + "insertOnlyColumn." + to_string(NUM_VECTORS - 1) + "." +
                            to_string(chunkNum++);
    while (is_file_exist(colorsFileName.c_str())) {
        vectorOfVectors *f = new vectorOfVectors(colorId);
        f->loadFromInsertOnly(colorsFileName, idsMap);
        colorId += f->size();
        colors.push_back(f);
        colorsFileName = col->tmpFolder + "insertOnlyColumn." + to_string(NUM_VECTORS - 1) + "." +
                         to_string(chunkNum++);
    }


}

void mixVectors::sortColors() {
#pragma omp parallel for
    for (unsigned int i = 0; i < colors.size(); i++)
        colors[i]->sort(idsMap);
}


Column *mixVectors::getTwin() {
    return new mixVectors();
}



void mixVectors::resize(uint32_t size) {

}


void fixedSizeVector::loadFromInsertOnly(string path, sdsl::int_vector<> &idsMap) {
    sdsl::int_vector<> curr;
    sdsl::load_from_file(curr, path);
    uint32_t noColors = curr.size() / (colorsize + 1);
    sdsl::int_vector<> tmpVec(noColors * colorsize);
    uint32_t top = 0;
    uint32_t colorId = beginID;
    for (unsigned int i = 0; i < curr.size(); i += colorsize + 1) {
        idsMap[curr[i]] = colorId++;
        for (int j = 0; j < colorsize; j++)
            tmpVec[top++] = curr[i + 1 + j];
    }
    //    vec=sdsl::enc_vector<>(tmpVec);
    vec = sdsl::int_vector<>(tmpVec);
    ((fixedSizeVectorIterator *) (endIterator->iterator))->it = vec.end();

}

void vectorOfVectors::loadFromInsertOnly(string path, sdsl::int_vector<> &idsMap) {
    sdsl::int_vector<> curr;
    sdsl::load_from_file(curr, path);
    vector<vector<uint32_t> > bigColors;
    vector<uint32_t> bigColorsIds;
    uint32_t i = 0;
    uint32_t colorId = beginID;
    uint32_t total = 0;
    while (i < curr.size()) {
        bigColorsIds.push_back(curr[i++]);
        bigColors.push_back(vector<uint32_t>(curr[i++]));
        total += bigColors.back().size();
        for (unsigned int j = 0; j < bigColors.back().size(); j++) {
            bigColors.back()[j] = curr[i++];
        }
    }

    sdsl::int_vector<> tmpStarts(bigColors.size());
    sdsl::int_vector<> tmpvecs(total);
    uint32_t curri = 0;
    for (unsigned int i = 0; i < bigColors.size(); i++) {
        tmpStarts[i] = curri;
        idsMap[bigColorsIds[i]] = colorId++;
        for (unsigned int j = 0; j < bigColors[i].size(); j++) {
            tmpvecs[curri++] = bigColors[i][j];
        }
    }
    vecs = vectype(tmpvecs);
    starts = vectype(tmpStarts);
}

void mixVectors::serialize(string filename) {
    ofstream out(filename.c_str());
    out.write((char *) (&(noSamples)), sizeof(uint32_t));
    idsMap.serialize(out);

    uint32_t tmp = colors.size();
    out.write((char *) (&(tmp)), sizeof(uint32_t));
    for (auto v:colors) {
        if (dynamic_cast<fixedSizeVector *>(v))
            tmp = 0;
        else if (dynamic_cast<vectorOfVectors *>(v))
            tmp = 1;
        else
            throw logic_error("Not supported vector");

        out.write((char *) (&(tmp)), sizeof(uint32_t));
        v->serialize(out);
        //  cout<<v->beginID<<" "<<tmp<<endl;
    }
    out.close();

}
vector<uint32_t > mixVectors::get(uint32_t index) {
    throw std::logic_error("get is not implemented use getWithIndex instead");

}
void mixVectors::deserialize(string filename) {
    ifstream input(filename.c_str());
    input.read((char *) (&(noSamples)), sizeof(uint32_t));
    idsMap.load(input);

    uint32_t colorsSize;
    input.read((char *) (&(colorsSize)), sizeof(uint32_t));
    uint32_t currColor = 0;
    for (uint32_t i = 0; i < colorsSize; i++) {
        uint32_t vecType;
        input.read((char *) (&(vecType)), sizeof(uint32_t));
        vectorBase *vec;
        if (vecType == 0)
            vec = new fixedSizeVector();
        else if (vecType == 1)
            vec = new vectorOfVectors();
        else
            throw logic_error("Not supported vector");
        vec->deserialize(input);
        vec->beginID = currColor;;
        // cout<<vec->beginID<<" "<<vecType<<endl;
        currColor += vec->size();
        colors.push_back(vec);
    }
    input.close();
    numColors = size();

}

uint64_t mixVectors::sizeInBytes() {
    uint64_t res = 0;
    for (auto vec:colors) {
        res += vec->sizeInBytes();
    }
    res += sdsl::size_in_bytes(idsMap);
    // cout<<"Ids Size = "<<sdsl::size_in_bytes(idsMap)/(1024.0*1024.0)<<"MB"<<endl;
    return res;
}

void mixVectors::explainSize() {
    cout << "Query Column" << endl;
    uint64_t numIntegers = 0;
    cout << "Ids Size = " << sdsl::size_in_bytes(idsMap) / (1024.0 * 1024.0) << "MB" << endl;
    double vMBBytes = 0.0;
    for (auto vec:colors) {
        vMBBytes += vec->sizeInMB();
        vec->explainSize();
        numIntegers += vec->numIntegers();
    }
    cout << "Arrays sizes = " << vMBBytes << "MB" << endl;
    cout << "Total = " << sizeInBytes() / (1024.0 * 1024.0) << "MB" << endl;
    cout << "Num Integers = " << numIntegers << endl;
    // cout<<"Ids Size = "<<sdsl::size_in_bytes(idsMap)/(1024.0*1024.0)<<"MB"<<endl;
}


void mixVectors::optimizeRLE() {
    explainSize();
    for (unsigned int i = 0; i < colors.size(); i++) {
        if (dynamic_cast<fixedSizeVector *>(colors[i])) {
            fixedSizeVector *tmp = (fixedSizeVector *) colors[i];

            RLEfixedSizeVector *tmp2 = new RLEfixedSizeVector(tmp, idsMap);
            colors[i] = tmp2;
//            if(tmp2->sizeInBytes() < tmp->sizeInBytes()) {
//                colors[i]=tmp2;
//                delete tmp;
//            }
//            else{
//                delete tmp2;
//            }
        }
    }
    explainSize();

}

fixedSizeVector::fixedSizeVector() {
    auto it = new fixedSizeVectorIterator(this);
    it->it = vec.end();
    endIterator = new vectorBaseIterator(it);
}

fixedSizeVector::fixedSizeVector(uint32_t beginId, uint32_t colorsize)
        : vectorBase(beginId) {
    auto it = new fixedSizeVectorIterator(this);
    it->it = vec.end();
    endIterator = new vectorBaseIterator(it);
    this->colorsize = colorsize;
    // vec.resize(noColors*size);
}

vectorBaseIterator fixedSizeVector::begin() {
    return vectorBaseIterator(new fixedSizeVectorIterator(this));
}

vectorBaseIterator fixedSizeVector::end() {
    ((fixedSizeVectorIterator *) endIterator->iterator)->it = vec.end();
    return *endIterator;
}

void fixedSizeVector::serialize(ofstream &f) {
    f.write((char *) (&(colorsize)), sizeof(uint32_t));
    vec.serialize(f);
}

void fixedSizeVector::deserialize(ifstream &f) {
    f.read((char *) (&(colorsize)), sizeof(uint32_t));
    vec.load(f);

//    sdsl::int_vector<> tmp;
//    tmp.load(f);
//    vec=sdsl::enc_vector<>(tmp);
}

void fixedSizeVector::sort(sdsl::int_vector<> &idsMap) {
    uint32_t numColors = size();
    unordered_map<uint32_t, uint32_t> idsINV;
    for (unsigned int i = 0; i < idsMap.size(); i++) {
        if (idsMap[i] >= beginID && idsMap[i] < beginID + numColors) {
            idsINV[idsMap[i]] = i;
        }
    }
    vector<pair<vector<uint32_t>, uint32_t> > aux(numColors);
    auto it = vec.begin();
    for (unsigned int i = 0; i < numColors; i++) {
        aux[i] = make_pair(vector<uint32_t>(colorsize), idsINV[beginID + i]);
        for (unsigned int j = 0; j < colorsize; j++) {
            aux[i].first[j] = *it;
            it++;
        }
    }
    std::sort(aux.begin(), aux.end(), []
            (const pair<vector<uint32_t>, uint32_t> &lhs, pair<vector<uint32_t>, uint32_t> &rhs) -> bool {
        if (lhs.first[0] > rhs.first[0])
            return true;
        else if (lhs.first[0] < rhs.first[0])
            return false;

        for (unsigned int i = 1; i < lhs.first.size() && i < rhs.first.size(); i++)
            if (lhs.first[i] < rhs.first[i])
                return true;
            else if (lhs.first[i] > rhs.first[i])
                return lhs.first.size() < rhs.first.size();
        return lhs.first.size() < rhs.first.size();
    });
    for (unsigned int i = 0; i < numColors; i++) {
        idsMap[aux[i].second] = beginID + i;
    }

    vector<uint32_t> tmpVec(numColors * colorsize);
    for (unsigned int i = 0; i < numColors; i++) {
        for (unsigned int j = 0; j < colorsize; j++)
            tmpVec[i * colorsize + j] = aux[i].first[j];

    }
    vec = vectype(tmpVec);
}


vectorOfVectors::vectorOfVectors() {
    endIterator = new vectorBaseIterator(new vectorOfVectorsIterator(this));
}

vectorOfVectors::vectorOfVectors(uint32_t beginId)
        : vectorBase(beginId) {
    endIterator = new vectorBaseIterator(new vectorOfVectorsIterator(this));
}

vectorOfVectors::vectorOfVectors(uint32_t beginId, uint32_t noColors)
        : vectorBase(beginId) {
    endIterator = new vectorBaseIterator(new vectorOfVectorsIterator(this));
    starts = vectype(sdsl::int_vector<>(noColors));
}

vectorBaseIterator vectorOfVectors::begin() {
    return vectorBaseIterator(new vectorOfVectorsIterator(this));
}

vectorBaseIterator vectorOfVectors::end() {
    ((vectorOfVectorsIterator *) endIterator->iterator)->vecsIt = vecs.end();
    ((vectorOfVectorsIterator *) endIterator->iterator)->startsIt = starts.end();
    return *endIterator;
}

void vectorOfVectors::serialize(ofstream &f) {
    vecs.serialize(f);
    starts.serialize(f);
}

void vectorOfVectors::deserialize(ifstream &f) {
    vecs.load(f);
    starts.load(f);
}


void vectorOfVectors::sort(sdsl::int_vector<> &idsMap) {
    uint32_t numColors = size();
    unordered_map<uint32_t, uint32_t> idsINV;
    for (unsigned int i = 0; i < idsMap.size(); i++) {
        if (idsMap[i] >= beginID && idsMap[i] < beginID + numColors) {
            idsINV[idsMap[i]] = i;
        }
    }
    vector<pair<vector<uint32_t>, uint32_t> > aux(numColors);
    auto it = vecs.begin();
    auto itStart = starts.begin();
    unsigned int i = 0;
    while (itStart != starts.end()) {
        uint32_t start = *itStart;
        uint32_t end = vecs.size();
        if ((itStart + 1) != starts.end()) {
            end = *(itStart + 1);
        }
        aux[i] = make_pair(vector<uint32_t>(end - start), idsINV[beginID + i]);
        for (unsigned int j = 0; j < end - start; j++) {
            aux[i].first[j] = *it;
            it++;
        }
        i++;
        itStart++;
    }
    std::sort(aux.begin(), aux.end(), []
            (const pair<vector<uint32_t>, uint32_t> &lhs, pair<vector<uint32_t>, uint32_t> &rhs) -> bool {
        if (lhs.first[0] > rhs.first[0])
            return true;
        else if (lhs.first[0] < rhs.first[0])
            return false;

        for (unsigned int i = 1; i < lhs.first.size() && i < rhs.first.size(); i++)
            if (lhs.first[i] < rhs.first[i])
                return true;
            else if (lhs.first[i] > rhs.first[i])
                return false;
        return lhs.first.size() < rhs.first.size();
    });
    for (unsigned int i = 0; i < numColors; i++) {
        idsMap[aux[i].second] = beginID + i;
    }

    vector<uint32_t> tmpVec(vecs.size());
    vector<uint32_t> tmpStarts(starts.size());
    uint32_t curr = 0;
    for (unsigned int i = 0; i < numColors; i++) {
        tmpStarts[i] = curr;
        for (unsigned int j = 0; j < aux[i].first.size(); j++)
            tmpVec[curr++] = aux[i].first[j];

    }
    vecs = vectype(tmpVec);
    starts = vectype(tmpStarts);
}

RLEfixedSizeVector::RLEfixedSizeVector(fixedSizeVector *fv, sdsl::int_vector<> &idsMap)
        : vectorBase(fv->beginID) {
    colorsize = fv->colorsize;
    numColors = fv->size();
    unordered_map<uint32_t, uint32_t> idsINV;
    for (unsigned int i = 0; i < idsMap.size(); i++) {
        if (idsMap[i] >= beginID && idsMap[i] < beginID + numColors) {
            idsINV[idsMap[i]] = i;
        }
    }
    vector<pair<vector<uint32_t>, uint32_t> > aux(numColors);
    auto it = fv->vec.begin();
    for (unsigned int i = 0; i < numColors; i++) {
        aux[i] = make_pair(vector<uint32_t>(colorsize), idsINV[beginID + i]);
        for (unsigned int j = 0; j < colorsize; j++) {
            aux[i].first[j] = *it;
            it++;
        }
    }
    std::sort(aux.begin(), aux.end(), []
            (const pair<vector<uint32_t>, uint32_t> &lhs, pair<vector<uint32_t>, uint32_t> &rhs) -> bool {
        for (unsigned int i = 0; i < lhs.first.size() && i < rhs.first.size(); i++)
            if (lhs.first[i] < rhs.first[i])
                return true;
            else if (lhs.first[i] > rhs.first[i])
                return false;
        return false;
    });
    for (unsigned int i = 0; i < numColors; i++) {
        idsMap[aux[i].second] = beginID + i;
    }

    deque<uint32_t> tmpVec;
    deque<uint32_t> tmpStart;
    uint32_t last = -1;
    for (unsigned int j = 0; j < colorsize; j++) {
        for (unsigned int i = 0; i < numColors; i++)
            if (aux[i].first[j] != last) {
                last = aux[i].first[j];
                tmpVec.push_back(last);
                tmpStart.push_back(i + j * numColors);
            }
    }
    vec = vectype(tmpVec);
    starts = vectype(tmpStart);
}

vector<uint32_t> RLEfixedSizeVector::get(uint32_t index) {
    vector<uint32_t> result(colorsize);
    for (unsigned int i = 0; i < colorsize; i++) {
        uint32_t tindex = i * numColors + index;
        auto it = lower_bound(starts.begin(), starts.end(), tindex + 1);
        it--;
        result[i] = vec[it - starts.begin()];
    }
    return result;
}


void RLEfixedSizeVector::serialize(ofstream &of) {}

void RLEfixedSizeVector::deserialize(ifstream &iif) {}

prefixTrie::prefixTrie(insertColorColumn *col)
{
    mixVectors* qq= new mixVectors(col);
    loadFromQueryColorColumn(qq);
    delete qq;
}
prefixTrie::prefixTrie(mixVectors *col) {
    loadFromQueryColorColumn(col);

}
void prefixTrie::loadFromQueryColorColumn(mixVectors *col) {
    noSamples = col->noSamples;
    numColors = col->size();
    col->sortColors();
    cerr << "Colors Sorted" << endl;
    idsMap = sdsl::int_vector<64>(col->idsMap.size());
    sdsl::int_vector<> invIdsMap(col->idsMap.size());
#pragma omp parallel for
    for (unsigned int i = 0; i < col->idsMap.size(); i++) {
        invIdsMap[col->idsMap[i]] = i;
    }
    cerr << "Inverted Ids is calculated" << endl;


    auto compare = [](tuple<vector<uint32_t>, uint32_t, vectorBaseIterator *, vectorBaseIterator *> lhs,
                      tuple<vector<uint32_t>, uint32_t, vectorBaseIterator *, vectorBaseIterator *> rhs) {
        if (std::get<0>(lhs)[0] < std::get<0>(rhs)[0])
            return true;
        else if (std::get<0>(lhs)[0] > std::get<0>(rhs)[0])
            return false;


        for (unsigned int i = 1; i < std::get<0>(lhs).size() && i < std::get<0>(rhs).size(); i++)
            if (std::get<0>(lhs)[i] > std::get<0>(rhs)[i])
                return true;
            else if (std::get<0>(lhs)[i] < std::get<0>(rhs)[i])
                return false;
        return std::get<0>(lhs).size() > std::get<0>(rhs).size();
    };
    priority_queue<tuple<vector<uint32_t>, uint32_t, vectorBaseIterator *, vectorBaseIterator *>, vector<tuple<vector<uint32_t>, uint32_t, vectorBaseIterator *, vectorBaseIterator *> >, decltype(compare)> nextColor(
            compare);
    for (auto c:col->colors) {
        vectorBaseIterator *it = new vectorBaseIterator(c->begin());
        vectorBaseIterator *itEnd = new vectorBaseIterator(c->end());
        if (*it != *itEnd) {
            vector<uint32_t> arr = **it;
            if (!arr.empty())
                nextColor.push(make_tuple(arr, it->getID(), it, itEnd));
        }
    }

    uint64_t tmpSize = col->numIntegers() / 10;


    uint64_t tmpEdgesTop = 0;
    uint64_t tmpTreeTop = 0;
    uint64_t currTree = 0;
    starts = sdsl::int_vector<64>(noSamples);
    sdsl::int_vector<> tmp_edges(tmpSize);

    starts[currTree] = 0;



    unordered_map<int, uint64_t> addedEdgesHisto;
    for (int i = 0; i <= noSamples; i++)
        addedEdgesHisto[i] = 0;

    uint64_t processedColors = 0;
    uint64_t printChunk = numColors / 20;
    uint64_t rank = 0;
    while (!nextColor.empty()) {
        starts[currTree] = rank;
        // generate the colors to be compressed and put them in colorsToProcess
        deque<tuple<vector<uint32_t>, uint32_t, vectorBaseIterator *, vectorBaseIterator *> > colorsToProcess;
        uint32_t firstNode = *std::get<0>(nextColor.top()).begin();
        while (!nextColor.empty() && firstNode == *std::get<0>(nextColor.top()).begin()) {
            auto colorTuple=nextColor.top();
            colorsToProcess.push_back(colorTuple);
            nextColor.pop();
            std::get<2>(colorTuple)->next();
            if (*std::get<2>(colorTuple) != *std::get<3>(colorTuple)) {
                std::get<0>(colorTuple) = *(*(std::get<2>(colorTuple)));
                std::get<1>(colorTuple) = std::get<2>(colorTuple)->getID();
                nextColor.push(colorTuple);
            } else {
                delete std::get<2>(colorTuple);
                delete std::get<3>(colorTuple);
            }

        }


        //process the colors
        deque<uint64_t> pastNodes;
        tmpEdgesTop = 0;
        tmpTreeTop = 0;
        tmpSize = tmp_edges.size() * 2;
        //tmp_edges.resize(tmpSize);
        tree.push_back(new sdsl::bit_vector(tmpSize * 2));
        for (unsigned c=0;c<colorsToProcess.size();c++) {
            auto colorTuple=colorsToProcess[c];
            vector<uint32_t> prevColor;
            if(c>0)
                prevColor=std::get<0>(colorsToProcess[c-1]);
            vector<uint32_t> currColor(std::get<0>(colorTuple).begin(), std::get<0>(colorTuple).end());
            unsigned int i = 0;
            for (; i < prevColor.size() && i < std::get<0>(colorTuple).size(); i++) {
                if (prevColor[i] != std::get<0>(colorTuple)[i])
                    break;

            }
            unordered_set<uint32_t> unneededNodes(prevColor.begin() + i, prevColor.end());
            unordered_set<uint32_t> neededNodes(currColor.begin(), currColor.end());
            deque<uint32_t> toBAdded;
            toBAdded.clear();
            bool hasUnneeded = false;
            unsigned int j = 0;
            for (; j < pastNodes.size(); j++) {

                for (auto c:nodesCache[pastNodes[j]])
                    if (unneededNodes.find(c) != unneededNodes.end()) {
                        hasUnneeded = true;
                        break;
                    }

                if (hasUnneeded)
                    break;

            }
            for (unsigned int k = j; k < pastNodes.size(); k++) {
                for (auto t:nodesCache[pastNodes[k]]) {
                    if (neededNodes.find(t) != neededNodes.end())
                        toBAdded.push_back(t);
                }
                nodesCache.erase(pastNodes[k]);
                rank++;
                (*tree.back())[tmpTreeTop++] = 0;
                if (tmpTreeTop == tree.back()->size()) {
                    cerr << "Tmp bp_tree of size " << tree.back()->size() << "("
                         << sdsl::size_in_mega_bytes(*tree.back())
                         << "MB) is full! size will doubled" << endl;
                    tree.back()->resize(tree.back()->size() * 2);
                    tmpSize *= 2;
                }
            }
            pastNodes.erase(pastNodes.begin() + j, pastNodes.end());

            // vector<uint32_t> tobeAdded(std::get<0>(colorTuple).size()-i);
            //   std::copy(std::get<0>(colorTuple).begin(),std::get<0>(colorTuple).end(),tobeAdded.begin());
//        for(auto it=currColor.begin() + i;it !=currColor.end();it++)
//            toBAdded.
            for (auto it = currColor.begin() + i; it != currColor.end(); it++) {
                toBAdded.push_back(*it);
            }
            std::sort(toBAdded.begin(), toBAdded.end());
            auto last = std::unique(toBAdded.begin(), toBAdded.end());
            toBAdded.erase(last, toBAdded.end());
            addedEdgesHisto[toBAdded.size()] += 1;
            uint32_t inputSize = toBAdded.size();
            deque<uint32_t> shortened;
            shortened.clear();
            shorten(toBAdded, shortened);
            uint32_t outputSize = 0;
            for (auto s:shortened)
                outputSize += nodesCache[s].size();
            if (outputSize != inputSize) {
                cerr << "Build error in rank " << rank << endl;
            }
            for (auto sample:shortened) {
                rank++;
                pastNodes.push_back(sample);
                (*tree.back())[tmpTreeTop++] = 1;
                if (tmpTreeTop == tree.back()->size()) {
                    cerr << "Tmp bp_tree of size " << tree.back()->size() << "("
                         << sdsl::size_in_mega_bytes(*tree.back())
                         << "MB) is full! size will doubled" << endl;
                    tree.back()->resize(tree.back()->size() * 2);
                }

                tmp_edges[tmpEdgesTop++] = sample;
                if (tmpEdgesTop == tmp_edges.size()) {
                    cerr << "Tmp edges of size (" << sdsl::size_in_mega_bytes(tmp_edges)
                         << "MB) is full! size will doubled"
                         << endl;
                    tmp_edges.resize(tmp_edges.size() * 2);
                }
            }

            idsMap[invIdsMap[std::get<1>(colorTuple)]] = rank - 1;

            processedColors++;
        }
        cout << "Processed " << processedColors << " / " << numColors << endl;

        // 
        for (unsigned int k = 0; k < pastNodes.size(); k++) {
            nodesCache.erase(pastNodes[k]);
            rank++;
            (*tree.back())[tmpTreeTop++] = false;
            if (tmpTreeTop == tree.back()->size()) {
                cerr << "Tmp bp_tree of size " << tree.back()->size() << "(" << sdsl::size_in_mega_bytes(*tree.back())
                     << "MB) is full! size will doubled" << endl;
                tree.back()->resize(tree.back()->size() * 2);
                tmpSize *= 2;
            }
        }
        //tmp_edges.resize(tmpEdgesTop);
        //vector<uint64_t> tmpVec(tmpEdgesTop);
        unCompressedEdges.push_back(new sdsl::int_vector<>(tmpEdgesTop));
        std::copy(tmp_edges.begin(), tmp_edges.begin() + tmpEdgesTop, unCompressedEdges.back()->begin());
        sdsl::util::bit_compress(*unCompressedEdges.back());
        //edges.push_back(new vectype(tmpVec));

        tree.back()->resize(tmpTreeTop);
        bp_tree.push_back(new sdsl::bp_support_sada<>(tree.back()));


        currTree++;


    }

    double unCompressedSize = 0.0;
    for (auto e:unCompressedEdges) {
        unCompressedSize += sdsl::size_in_mega_bytes(*e);
        edges.push_back(new vectype(*e));
        delete e;
    }
    unCompressedEdges.clear();
    cout << "Uncompressed edges size = " << unCompressedSize << endl;


    uint64_t edgesSum = 0;
    for (auto a:addedEdgesHisto) {
        edgesSum += (a.first - 1) * (a.second);
    }
    cout << "Possible saving " << edgesSum << endl;

}


uint32_t prefixTrie::insertAndGetIndex(vector<uint32_t> &item) {
    throw std::logic_error("insertAndGetIndex is not supported in mixVectors");

}

vector<uint32_t > prefixTrie::get(uint32_t index) {
    throw std::logic_error("get is not implemented use getWithIndex instead");

}

vector<uint32_t> prefixTrie::getWithIndex(uint32_t index) {
    deque<uint32_t> tmp;
    queue<uint64_t> Q;
    Q.push(idsMap[index]);
    // cout<<idsMap[index]<<" -> ";
    while (!Q.empty()) {
        uint64_t bigIndex = Q.front();
        Q.pop();
        auto it = lower_bound(starts.begin(), starts.end(), bigIndex + 1);
        it--;
        uint32_t tIndex = it - starts.begin();
        bigIndex -= *it;
        while (bigIndex != bp_tree[tIndex]->size()) {
            uint64_t edgeIndex = bp_tree[tIndex]->rank(bigIndex) - 1;
            uint64_t node = (*edges[tIndex])[edgeIndex];
            //        cout<<node<<" ";
            if (node < noSamples)
                tmp.push_back(node);
            else
                Q.push(node - noSamples);
            bigIndex = bp_tree[tIndex]->enclose(bigIndex);
        }
    }
    //cout<<endl;
    sort(tmp.begin(), tmp.end());
    vector<uint32_t> res(tmp.size());
    for (unsigned int i = 0; i < res.size(); i++)
        res[i] = tmp[i];
    return res;
}

void prefixTrie::insert(vector<uint32_t> &item, uint32_t index) {
    throw std::logic_error("insertAndGetIndex is not supported in mixVectors");

}
//vector<uint32_t > prefixTrie::get(uint32_t index){
//    return vector<uint32_t >();
//}

void prefixTrie::serialize(string filename) {
    ofstream out(filename.c_str());
    out.write((char *) (&(noSamples)), sizeof(uint32_t));
    out.write((char *) (&(numColors)), sizeof(uint32_t));
    idsMap.serialize(out);
    starts.serialize(out);
    for (uint32_t i = 0; i < tree.size(); i++) {
        tree[i]->serialize(out);
        bp_tree[i]->serialize(out);
        edges[i]->serialize(out);
    }

    out.close();
}

void prefixTrie::deserialize(string filename) {
    ifstream input(filename.c_str());
    input.read((char *) (&(noSamples)), sizeof(uint32_t));
    input.read((char *) (&(numColors)), sizeof(uint32_t));
    idsMap.load(input);
    starts.load(input);
    for (uint32_t i = 0; i < starts.size(); i++) {
        tree.push_back(new sdsl::bit_vector());
        tree.back()->load(input);
        bp_tree.push_back(new sdsl::bp_support_sada<>());
        bp_tree.back()->load(input, tree.back());
        edges.push_back(new vectype());
        edges.back()->load(input);
    }

    input.close();

}


uint32_t prefixTrie::getNumColors() {
    return numColors;
}

uint64_t prefixTrie::sizeInBytes() {
    uint64_t res = 0;
    for (auto t:tree)
        res += sdsl::size_in_bytes(*t);
    for (auto b:bp_tree)
        res += sdsl::size_in_bytes(*b);
    res += sdsl::size_in_bytes(idsMap);
    for (auto e:edges)
        res += sdsl::size_in_bytes(*e);
    return res;
}

void prefixTrie::explainSize() {
    double treeSize = 0;
    double treeCompressedSize = 0;
    double bpSize = 0;
    double eSize = 0;
    uint64_t numE = 0;
    for (auto t:tree) {
        treeSize += sdsl::size_in_mega_bytes(*t);
        sdsl::rrr_vector<> cvector(*t);
        treeCompressedSize += sdsl::size_in_mega_bytes(cvector);
    }
    for (auto b:bp_tree)
        bpSize += sdsl::size_in_mega_bytes(*b);
    for (auto e:edges) {
        eSize += sdsl::size_in_mega_bytes(*e);
        numE += e->size();
    }


    cout << "Prefix Trie index" << endl;
    cout << "Bit Tree = " << treeSize << "MB\n";
    cout << "Bit Tree Compressed = " << treeCompressedSize << "MB\n";
    cout << "BpSupport = " << bpSize << "MB\n";
    cout << "Ids Map = " << sdsl::size_in_mega_bytes(idsMap) << "MB\n";

    cout << "edges = " << eSize << "MB\n";
    cout << "edges # Integers= " << numE << "\n";
    double total = eSize + sdsl::size_in_mega_bytes(idsMap)
                   + bpSize +
                   treeSize;
    cout << "Total = " << total << "MB" << endl;


}


void prefixTrie::shorten(deque<uint32_t> &input, deque<uint32_t> &output) {
    if (input.size() == 1) {
        output.push_back(input[0]);
        nodesCache[input[0]] = {input[0]};
        return;
    }
    uint32_t treeIndex = noSamples - input[0] - 1;
    if (treeIndex >= unCompressedEdges.size()) {
        output.push_back(input[0]);
        nodesCache[input[0]] = {input[0]};
        input.erase(input.begin());
        if (!input.empty())
            shorten(input, output);
        return;
    }
    deque<uint32_t> remaining;
    vector<uint32_t> chosen;
    if ((*unCompressedEdges[treeIndex])[0] != input[0]) {
        cerr << "Wrong tree " << (*unCompressedEdges[treeIndex])[0] << endl;
        return;
    }
    auto i = input.begin();
    uint64_t treePos = 0;
    uint64_t result = tree.size();
    while (i != input.end() && treePos < tree[treeIndex]->size() && (*tree[treeIndex])[treePos] == 1) {
        uint64_t edgeIndex = bp_tree[treeIndex]->rank(treePos) - 1;
        uint32_t currNode = (*unCompressedEdges[treeIndex])[edgeIndex];
        auto it = lower_bound(i, input.end(), currNode);
        if (it == input.end() || *it != currNode) {
            treePos = bp_tree[treeIndex]->find_close(treePos) + 1;
        } else {
            for (; i < it; i++) {
                remaining.push_back(*i);
            }
            chosen.push_back(currNode);
            result = treePos;
            treePos++;
            i++;
        }
    }
    if (result == 0) {
        output.push_back(input[0]);
        nodesCache[input[0]] = {input[0]};
    } else {
        uint64_t ptr = result + starts[treeIndex] + noSamples;
        output.push_back(ptr);
        nodesCache[ptr] = chosen;
    }

    for (; i < input.end(); i++) {
        remaining.push_back(*i);
    }
    if (!remaining.empty()) {
        shorten(remaining, output);
    }

}

void prefixTrie::exportTree(string prefix, int treeIndex) {
    string outFilename = prefix + to_string(treeIndex);
    ofstream out(outFilename.c_str());
    int tabs = 0;
    out << "graph \"\"" << endl;
    tabs++;
    for (int i = 0; i < tabs; i++)
        out << "\t";
    out << "{" << endl;
    tabs++;
    for (int i = 0; i < tabs; i++)
        out << "\t";
    string bp = "";
    for (uint32_t i = 0; i < tree[treeIndex]->size(); i++) {
        if ((*tree[treeIndex])[i] == 0)
            bp += ')';
        else
            bp += '(';
    }
    out << "label =\"" << bp << "\"" << endl;
    //out<<"label =\"Tree "<<treeIndex<<"\""<<endl;
    uint32_t pos = 0;
    stack<uint32_t> parents;
    parents.push(0);
    for (int i = 0; i < tabs; i++)
        out << "\t";
    out << "n" << pos << " ;" << endl;
    for (int i = 0; i < tabs; i++)
        out << "\t";
    out << "n" << pos << " [label=\"" << (*edges[treeIndex])[pos] << "\"] ;" << endl;
    pos++;
    while (pos < tree[treeIndex]->size()) {
        for (int i = 0; i < tabs; i++)
            out << "\t";
        out << "n" << parents.top() << " -- n" << pos << " ;" << endl;
        for (int i = 0; i < tabs; i++)
            out << "\t";
        uint64_t edgeIndex = bp_tree[treeIndex]->rank(pos) - 1;
        out << "n" << pos << " [label=\"" << (*edges[treeIndex])[edgeIndex] << "\"] ;" << endl;
        parents.push(pos);
        pos++;
        while (pos < tree[treeIndex]->size() && (*tree[treeIndex])[pos] == 0) {
            pos++;
            parents.pop();
        }


    }
    tabs--;
    for (int i = 0; i < tabs; i++)
        out << "\t";
    out << "}" << endl;


}


Column *prefixTrie::getTwin() {
    return new prefixTrie();
}


void prefixTrie::resize(uint32_t size) {

}


template<typename T, typename ColumnType>
void deduplicatedColumn<T, ColumnType>::serialize(string filename) {
    string indexFilename = filename + ".index";
    string containerFilename = filename + ".container";
    std::ofstream os(indexFilename, std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(index);
    os.close();
    values->serialize(containerFilename);
}


template<typename T, typename ColumnType>
void deduplicatedColumn<T, ColumnType>::deserialize(string filename) {
    string indexFilename = filename + ".index";
    string containerFilename = filename + ".container";
    std::ifstream os(indexFilename, std::ios::binary);
    cereal::BinaryInputArchive iarchive(os);
    iarchive(index);
    os.close();
    values = new ColumnType();
    values->deserialize(containerFilename);
}


template<typename T, typename ColumnType>
T deduplicatedColumn<T, ColumnType>::get(uint32_t order) {
    return values->get(index[order]);
}

template<typename T, typename ColumnType>
Column *deduplicatedColumn<T, ColumnType>::getTwin() {
    return new deduplicatedColumn<T, ColumnType>();
}



template<typename T, typename ColumnType>
void deduplicatedColumn<T, ColumnType>::resize(uint32_t size) {
    index.resize(size);
}

template<typename T, typename ColumnType>
void
deduplicatedColumn<T, ColumnType>::setValueFromColumn(Column *Container, uint32_t inputOrder, uint32_t outputOrder) {
    deduplicatedColumn<T, ColumnType> *other = ((deduplicatedColumn<T, ColumnType> *) Container);
    values = other->values;
    index[outputOrder] = other->index[inputOrder];

}
