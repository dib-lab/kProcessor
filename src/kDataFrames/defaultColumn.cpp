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
template class vectorColumn<int>;
template class vectorColumn<bool>;
template class vectorColumn<double>;


Column* Column::getContainerByName(std::size_t hash)
{
    if(hash== typeid(vectorColumn<int>).hash_code())
    {
        return new vectorColumn<int>();
    }
    else if(hash== typeid(vectorColumn<double>).hash_code())
    {
        return new vectorColumn<double>();
    }
    else if(hash== typeid(vectorColumn<bool>).hash_code())
    {
        return new vectorColumn<bool>();
    }
    else if(hash== typeid(colorColumn).hash_code())
    {
        return new colorColumn();
    }
    else if(hash== typeid(StringColorColumn).hash_code())
    {
        return new StringColorColumn();
    }
    else
    {
        throw logic_error("Failed to load Unknown Column "+hash);
    }
}

template <typename  T>
uint32_t  vectorColumn<T>::insertAndGetIndex(T item){
    dataV.push_back(item);
    return dataV.size()-1;
}



template <typename  T>
T vectorColumn<T>::getWithIndex(uint32_t index){
    return dataV[index];
}

template <typename  T>
void vectorColumn<T>::serialize(string filename)
{
//    ofstream wf(filename, ios::out | ios::binary);
//    uint32_t  size=dataV.size();
//    wf.write((char*)&size, sizeof(uint32_t));
//    for(auto t:dataV)
//        wf.write((char*)&t, sizeof(T));
//    wf.close();
    std::ofstream os(filename , std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(dataV);
    os.close();
}


template <typename  T>
void vectorColumn<T>::deserialize(string filename)
{
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

template <typename  T>
void vectorColumn<T>::insert(T item,uint32_t index){
    dataV[index]=item;
}
template <typename  T>
T vectorColumn<T>::get(uint32_t index){
    return dataV[index];
}


uint32_t  colorColumn::insertAndGetIndex(vector<uint32_t > item){
  return colorInv.getColorID(item);
    // if(colorInv.hasColorID(item))
    //   return colorInv.getColorID(item);

    // uint32_t i= colorInv.getColorID(item);
    // colors.push_back(item);
    // if(i!=colors.size()-1)
    //   cout<<"error in insert and get index "<<i<<" "<<colors.size()<<endl;
    // return i;
  
}
vector<uint32_t > colorColumn::getWithIndex(uint32_t index){
    return colors[index];
}



void colorColumn::serialize(string filename)
{
    colorInv.serialize(filename);
}
void colorColumn::deserialize(string filename)
{
    colorInv.deserialize(filename);
    populateColors();
    noSamples=colorInv.noSamples;
}

void colorColumn::populateColors(){
    colorInv.populateColors(colors);
}


vector<string > StringColorColumn::getWithIndex(uint32_t index){
    vector<string > res(colors[index].size());
    for(unsigned  int i=0;i<colors[index].size();i++)
        res[i]=namesMap[colors[index][i]];
    return res;

}


void StringColorColumn::serialize(string filename)
{
    std::ofstream os(filename+".colors", std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(colors);
    os.close();

    ofstream namesMapOut(filename+".namesMap");
    namesMapOut<<namesMap.size()<<endl;
    for(auto it:namesMap)
    {
        namesMapOut<<it.first<<" "<<it.second<<endl;
    }
    namesMapOut.close();
}





void StringColorColumn::deserialize(string filename)
{
    std::ifstream os(filename+".colors", std::ios::binary);
    cereal::BinaryInputArchive iarchive(os);
    iarchive(colors);

    ifstream namesMapIn(filename+".namesMap");
    uint64_t size;
    namesMapIn>>size;
    for(unsigned  int i=0;i<size;i++)
    {
        uint32_t color;
        string name;
        namesMapIn>>color>>name;
        namesMap[color]=name;

    }
    namesMapIn.close();

}


colorIndex::~colorIndex()
{
    stack<tuple<colorNode*,bool> > S;
    S.push(make_tuple(root,false));

    while(S.size()>0)
    {
        if(!std::get<1>(S.top()))
        {
            std::get<1>(S.top())=true;
            for(auto it:std::get<0>(S.top())->edges)
            {
                S.push(make_tuple(it.second,false));
            }
        } else{
            delete std::get<0>(S.top());
            S.pop();
        }


    }
}
bool colorIndex::hasColorID(vector<uint32_t>& v)
{
    colorNode* currNode=root;
    unsigned  int i=0;
    while(i<v.size())
    {
        auto it=currNode->edges.find(v[i]);
        if(it==currNode->edges.end())
            break;
        currNode=it->second;
        i++;
    }
    if(i!=v.size())
        return false;
    return currNode->currColor!=0;
}
uint32_t colorIndex::getColorID(vector<uint32_t>& v)
{
    colorNode* currNode=root;
    unsigned int i=0;
    while(i<v.size())
    {
        auto it=currNode->edges.find(v[i]);
        if(it==currNode->edges.end())
            break;
        currNode=it->second;
        i++;
    }
    for(;i<v.size();i++)
    {
        currNode->edges[v[i]]=new colorNode();
        currNode=currNode->edges[v[i]];
    }
    if(currNode->currColor==0)
        currNode->currColor=++lastColor;
    return currNode->currColor;
}
void colorIndex::serialize(string fileName)
{
    ofstream output(fileName.c_str(), ios::out | ios::binary);
    output.write((char*)&lastColor, sizeof(uint32_t));
    stack<tuple<colorNode*,bool> > S;
    S.push(make_tuple(root,false));
    uint32_t tmp;
    while(S.size()>0)
    {
        if(!std::get<1>(S.top()))
        {
            std::get<1>(S.top())=true;

            tmp=get<0>(S.top())->currColor;
            output.write((char*)&tmp, sizeof(uint32_t));
            tmp=get<0>(S.top())->edges.size();
            output.write((char*)&tmp, sizeof(uint32_t));
            for(auto it:std::get<0>(S.top())->edges)
            {
                tmp=it.first;
                output.write((char*)&tmp, sizeof(uint32_t));
            }
            for(auto it:std::get<0>(S.top())->edges)
            {
                S.push(make_tuple(it.second,false));
            }
        }else{
            S.pop();
        }


    }

}


void colorIndex::deserialize(string fileName){
    ifstream rf(fileName.c_str(), ios::out | ios::binary);
    stack<tuple<colorNode*,bool> > S;
    rf.read((char *) &(lastColor), sizeof(uint32_t));
    S.push(make_tuple(root,false));
    uint32_t tmp;

    while(S.size()>0)
    {
        if(!std::get<1>(S.top()))
        {
            colorNode* curr=get<0>(S.top());
            std::get<1>(S.top())=true;


            rf.read((char *) &(tmp), sizeof(uint32_t));
            curr->currColor=tmp;

            uint32_t  noEdges;
            rf.read((char *) &(noEdges), sizeof(uint32_t));

            for(uint32_t i=0;i<noEdges;i++)
            {
                rf.read((char *) &(tmp), sizeof(uint32_t));
                colorNode* child=new colorNode();
                curr->edges[tmp]=child;
                S.push(make_tuple(child,false));
            }
        } else{
            S.pop();
        }


    }
}

void colorIndex::populateColors(vector<vector<uint32_t> > &colors) {
    colors=vector<vector<uint32_t > >(lastColor+1);
    vector<uint32_t> color;
    stack<tuple<colorNode*,uint32_t ,bool> > S;
    for(auto it:root->edges)
    {
        S.push(make_tuple(it.second,it.first,false));
    }
    S.push(make_tuple(root->edges[27],27,false));
    bool debug=false;
    while(S.size()>0)
    {
        if(!std::get<2>(S.top()))
        {
            std::get<2>(S.top())=true;
            if(debug)
            cout<<"Visiting c"<<std::get<0>(S.top())->currColor<< " With edge "<<std::get<1>(S.top())<<endl;
            color.push_back(std::get<1>(S.top()));
            noSamples=max(noSamples,std::get<1>(S.top()));
            if(std::get<0>(S.top())->currColor!=0)
            {
                colors[std::get<0>(S.top())->currColor]=color;
                if(debug) {
                    cout << std::get<0>(S.top())->currColor << " : ";
                    for (auto c :colors[std::get<0>(S.top())->currColor])
                        cout << c << " ";
                    cout << endl;
                }
//
//                if(std::get<0>(S.top())->currColor) {
//                    cout << std::get<0>(S.top())->currColor << " : ";
//                    for (auto c :colors[std::get<0>(S.top())->currColor])
//                        cout << c << " ";
//                    cout << endl;
//                }
            }

            for(auto it:std::get<0>(S.top())->edges)
            {
                if(debug)
                    cout<<"Pushing c:"<<it.second->currColor<<" withe edge "<<it.first<<endl;
                S.push(make_tuple(it.second,it.first,false));
            }
        } else{
            color.pop_back();
            S.pop();
        }


    }
    noSamples++;
}



void colorIndex::optimize() {
    stats();
    vector<vector<uint32_t > > colors=vector<vector<uint32_t > >(lastColor+1);
    vector<uint32_t> color;
    stack<tuple<colorNode*,uint32_t ,bool> > S;
    flat_hash_map<uint32_t ,uint32_t > freqs;
    for(auto it:root->edges)
    {
        S.push(make_tuple(it.second,it.first,false));
    }
    while(S.size()>0)
    {
        if(!std::get<2>(S.top()))
        {
            uint32_t currSample=std::get<1>(S.top());
            std::get<2>(S.top())=true;
            color.push_back(currSample);
            if(freqs.find(currSample)==freqs.end())
                freqs[currSample]=0;
            freqs[currSample]++;
            if(std::get<0>(S.top())->currColor!=0)
            {

                colors[std::get<0>(S.top())->currColor]=color;
            }

            for(auto it:std::get<0>(S.top())->edges)
            {
                S.push(make_tuple(it.second,it.first,false));
            }
        } else{
            color.pop_back();
            S.pop();
        }


    }
    vector<uint32_t > newColor(lastColor);
    for(unsigned int i=0;i<lastColor;i++)
        newColor[i]=i;
    sort(newColor.begin(), newColor.end(),
         [&](uint32_t & a,uint32_t & b) -> bool
         {
             uint32_t aa=a;
             uint32_t bb=b;
             return freqs[aa] > freqs[bb];
         });

    colorIndex newColorIndex;
    for(auto c:colors)
    {
        if(c.size()==0)
            continue;
        vector<uint32_t> newc(c.size());
        for(unsigned int i=0;i<c.size();i++) {
            newc[i] = newColor[c[i]];
        }
        sort(newc.begin(),newc.end());
        newColorIndex.getColorID(newc);
    }
    cout<<"Optimization Done"<<endl;
    newColorIndex.stats();
   // newColorIndex.optimize();
}

void colorIndex::stats() {
    flat_hash_map<uint32_t ,uint32_t > freqs;
    stack<tuple<colorNode*,uint32_t ,bool> > S;
    uint32_t  nNodes=0;
    for(auto it:root->edges)
    {
        S.push(make_tuple(it.second,it.first,false));
    }
    while(S.size()>0)
    {
        if(!std::get<2>(S.top()))
        {
            std::get<2>(S.top())=true;
            uint32_t  currSample=std::get<1>(S.top());
            if(freqs.find(currSample)==freqs.end())
                freqs[currSample]=0;
            freqs[currSample]++;
            nNodes++;
            for(auto it:std::get<0>(S.top())->edges)
            {
                S.push(make_tuple(it.second,it.first,false));
            }
        } else{

            S.pop();
        }


    }
    cout<<"Total Number of Nodes = "<<nNodes<<endl;
    for(auto it:freqs)
    {
        cout<<it.first<<" : "<<it.second<<"\n";
    }
}


bool stringColorIndex::hasColorID(vector<uint32_t>& v){
    auto it=colors.find(toString(v));
    return it!=colors.end();
}
uint32_t stringColorIndex::getColorID(vector<uint32_t>& v)
{
    string col=toString(v);
    auto it=colors.find(col);
    if(it==colors.end()){
        colors[col]=++lastColor;
        it=colors.find(col);
    }
    return it->second;
}
void stringColorIndex::serialize(string fileName){
    vector<vector<uint32_t > > outColors(colors.size()+1);
    for(auto c:colors)
    {
        vector<uint32_t > tmp;
        tmp.clear();
        regex r(";");
        sregex_token_iterator it(c.first.begin(), c.first.end(),r, -1);

        sregex_token_iterator reg_end;
        for (; it != reg_end; ++it) {
            if(it->str()!="")
                tmp.push_back(atoi(it->str().c_str()));
        }
        outColors[c.second]=tmp;
    }
    std::ofstream os(fileName+".colors", std::ios::binary);
    cereal::BinaryOutputArchive archive(os);
    archive(outColors);
    os.close();

}
void stringColorIndex::deserialize(string filename)
{
    std::ifstream os(filename+".colors", std::ios::binary);
    cereal::BinaryInputArchive iarchive(os);
    vector<vector<uint32_t > > inColors;
    iarchive(inColors);

    for(unsigned int i=0;i<inColors.size();i++)
    {
        colors[toString(inColors[i])]=i;
    }
}

void stringColorIndex::populateColors(vector<vector<uint32_t > >& outColors)
{
    outColors=vector<vector<uint32_t > >(colors.size()+1);
    for(auto c:colors)
    {
        vector<uint32_t > tmp;
        tmp.clear();
        regex r(";");
        sregex_token_iterator it(c.first.begin(), c.first.end(),r, -1);
        sregex_token_iterator reg_end;
        for (; it != reg_end; ++it) {
            if(it->str()!="") {
                tmp.push_back(atoi(it->str().c_str()));
            }
        }
        outColors[c.second]=tmp;
    }
}


compressedColorColumn::compressedColorColumn(colorColumn* col){
    noSamples=col->noSamples;
    colors.push_back(new constantVector(noSamples));
    colors.push_back(new vectorOfVectors(noSamples+1,col->getNumColors()-noSamples+1));
    optimize3(col);
    optimize2();
}

void compressedColorColumn::optimize(colorColumn* col)
{
    numColors=col->getNumColors();
    vector<pair<uint32_t,uint32_t > > sizeAndIndex(numColors);
    uint64_t  oldColorsSum=0,newColorsSum=0;
    idsMap.resize(numColors+1);
   // vector<sdsl::bit_vector> colorsBitvectors(colors.size());
    for(uint32_t i=0;i<numColors;i++)
    {
        vector<uint32_t> tmp=col->getWithIndex(i);
        sizeAndIndex[i]=make_pair((uint32_t)col->colors[i].size(),i);
        oldColorsSum+=col->colors[i].size();
//        colorsBitvectors[i]=sdsl::bit_vector(noSamples);
//        for(auto s:tmp) {
//            colorsBitvectors[i][s] = 1;
//        }
    }
    cout<<"Old Colors sum "<<oldColorsSum<<endl;
    sort(sizeAndIndex.begin(),sizeAndIndex.end());

    for(unsigned int ii=1; ii<numColors; ii++)
    {
        unsigned int i=sizeAndIndex[ii].second;
        idsMap[sizeAndIndex[ii].second]=ii;
        vector<uint32_t > newCombination;
        newCombination.clear();
        //vector<uint32_t > tmp= col->getWithIndex(i);
        deque<uint32_t > currV;
        copy(col->colors[i].begin(),col->colors[i].end(),back_inserter(currV));
      //  sdsl::bit_vector curr(colorsBitvectors[i]);
      //  uint32_t currOneBits = sdsl::util::cnt_one_bits( curr );
        for(unsigned  int jj=ii-1; jj>0 ; jj--)
        {
            if(currV.size()<sizeAndIndex[jj].first){
                auto it=lower_bound(sizeAndIndex.begin(),sizeAndIndex.begin()+jj,make_pair((uint32_t)currV.size(),(uint32_t)colors.size()));
                jj=it-sizeAndIndex.begin();
            }
//            if(currOneBits<sizeAndIndex[jj].first){
//                auto it=lower_bound(sizeAndIndex.begin(),sizeAndIndex.begin()+jj,make_pair(currOneBits,(uint32_t)colors.size()));
//                jj=it-sizeAndIndex.begin();
//            }

            unsigned j=sizeAndIndex[jj].second;
            if(j<noSamples)
                continue;

         //   vector<uint32_t > currJ=col->colors[j];
            bool isContain=true;
            auto it=currV.begin();
            for(auto c : col->colors[j])
            {
                it=lower_bound(it,currV.end(),c);
                if(it==currV.end()||*it!=c)
                {
                    isContain=false;
                    break;
                }
            }
            if(isContain)
            {
                newCombination.push_back(idsMap[j]);
                it=currV.begin();
                for(auto c : col->colors[j])
                {
                    it=lower_bound(it,currV.end(),c);
                    currV.erase(it);
                }
                if(currV.size()==0)
                {
                    break;
                }
            }

//            sdsl::bit_vector tmp(curr);
//            tmp|=colorsBitvectors[j];
//            if(tmp == curr)
//            {
//                newCombination.push_back(j);
//                sdsl::bit_vector mask(colorsBitvectors[j]);
//                mask.flip();
//                curr&=mask;
//                currOneBits = sdsl::util::cnt_one_bits( curr );
//                if(currOneBits==0)
//                    break;
//            }

        }
        for(auto c:currV)
            newCombination.push_back(c);

//        if(currOneBits!=0) {
//            for (uint32_t k = 0; k < noSamples; k++) {
//                if (curr[k])
//                    newCombination.push_back(k);
//            }
//        }
        newColorsSum+=newCombination.size();
        sort(newCombination.begin(),newCombination.end());
        insert(newCombination,ii);

    }
    sdsl::util::bit_compress(idsMap);
    cout<<"New Colors sum "<<newColorsSum<<endl;
}

void compressedColorColumn::optimize2()
{
    vector<pair<pair<uint32_t,uint32_t> ,uint32_t > > sizeAndMaxAndIndex(numColors);

    for(uint32_t i=1;i<=colors[0]->size();i++)
    {
        sizeAndMaxAndIndex[i]=make_pair(make_pair((uint32_t)1,i-1),i);

    }
    for(uint32_t i=0;i<colors[1]->size();i++) {
        vector<uint32_t> v=colors[1]->get(i);
        uint32_t maxColor=0;
        for(auto c:v)
            maxColor=max(maxColor,c);
        uint32_t index=i+colors[1]->beginID;
        sizeAndMaxAndIndex[index]=make_pair(make_pair((uint32_t)v.size(),maxColor),index);
    }

    sort(sizeAndMaxAndIndex.begin(),sizeAndMaxAndIndex.end());
    vector<uint32_t> newMap(numColors+1);
    for(unsigned int i=0;i<numColors;i++)
    {
        newMap[sizeAndMaxAndIndex[i].second]=i;
    }
    for(unsigned int i=0;i<numColors;i++)
    {
        idsMap[i]=newMap[idsMap[i]];
    }
    uint32_t i=noSamples+1;//skip the first trivial colors
    vectorBase* vec=colors[1];
    colors.erase(colors.begin()+1);
    uint32_t start=i,end=i;
    for(uint32_t currentSize=2;currentSize<10;currentSize++)
    {

        while(end<numColors && sizeAndMaxAndIndex[end].first.first==currentSize && sizeAndMaxAndIndex[end].first.second<65536)
        {
            end++;
        }
       // cout<<start<<" "<<end<<endl;
        if(start!=end) {
            fixedSizeVector *curr = new fixedSizeVector(end - start, currentSize);
            curr->beginID = start;
            colors.push_back(curr);
            for (; i < end; i++) {

                auto tmp = vec->get(sizeAndMaxAndIndex[i].second - vec->beginID);

                for(unsigned int j=0;j<tmp.size();j++)
                    tmp[j]=newMap[tmp[j]];

                insert(tmp, i);
            }
            sdsl::util::bit_compress(curr->vec);
        }
        start=i;
        while(end<numColors && sizeAndMaxAndIndex[end].first.first==currentSize)
        {
            end++;
        }
     //   cout<<start<<" "<<end<<endl;
        if(start!=end) {
            fixedSizeVector *curr = new fixedSizeVector(end - start, currentSize);
            curr->beginID = start;
            colors.push_back(curr);
            for (; i < end; i++) {
                auto tmp = vec->get(sizeAndMaxAndIndex[i].second - vec->beginID);
                for(unsigned int j=0;j<tmp.size();j++)
                    tmp[j]=newMap[tmp[j]];
                insert(tmp, i);
            }
            sdsl::util::bit_compress(curr->vec);
        }
    }
    start=end;
    end=numColors;
    if(start!=end) {
        vectorOfVectors *curr = new vectorOfVectors(start, end - start);
        colors.push_back(curr);
        for (; i < end; i++) {
            auto tmp = vec->get(sizeAndMaxAndIndex[i].second - vec->beginID);
            for (unsigned int j = 0; j < tmp.size(); j++)
                tmp[j] = newMap[tmp[j]];
            insert(tmp, i);
        }
    }




}
///todo return the longest seq in node colors
uint32_t getLongestSubsetColor(colorColumn* col,deque<uint32_t> & color,uint32_t colorId){
    unordered_map<uint32_t ,uint32_t > nodeColors(color.size());
    stack<tuple<colorNode*,uint32_t ,bool> > S;
    S.push(make_tuple(col->colorInv.root,0,false));

    uint32_t result;
    uint32_t resultsize=0;
    while(S.size()>0)
    {
        colorNode* currNode=std::get<0>(S.top());
        uint32_t sample=std::get<1>(S.top());
        uint32_t currColor=currNode->currColor;

        bool visited=std::get<2>(S.top());
    //    cout<<sample<<"-"<<visited<<endl;
        if(!visited)
        {

            std::get<2>(S.top())=true;
          //  cout<<"Visiting "<<std::get<1>(S.top())<<endl;
            if(currColor!=0)
            {
//                cout<<currNode->currColor<<" : ";
//                for(auto c :col->colors[currNode->currColor])
//                    cout<<c<<" ";
//                cout<<endl;
 //               nodeColors[std::get<1>(S.top())]=std::get<0>(S.top())->currColor;
                if(currColor > col->noSamples && currColor!=colorId)
                    nodeColors[sample]=currNode->currColor;
                else
                    nodeColors[sample]=sample;
            }

//            for(auto it:std::get<0>(S.top())->edges)
//            {
//                if(find(color.begin(),color.end(),it.first)!=color.end())

            for(auto c:color)
            {
                auto it=currNode->edges.find(c);
                if(it!=currNode->edges.end()) {
               //     cout<<"Pushing c:"<<it->second->currColor<<" withe edge "<<it->first<<endl;

                    S.push(make_tuple(it->second, it->first, false));
                }
            }
        } else{
            if(currColor!=0) {
                uint32_t currColor = nodeColors[sample];
                uint32_t currColorLength = col->colors[currColor].size();

//            for(auto it:std::get<0>(S.top())->edges)
//            {
//                if(find(color.begin(),color.end(),it.first)!=color.end())
//                {
                for (auto c:color) {
                    auto it = currNode->edges.find(c);
                    if (it != currNode->edges.end()) {
                        uint32_t tmpColorLength = col->colors[nodeColors[it->first]].size();
                        if (tmpColorLength > currColorLength) {
                            currColorLength = tmpColorLength;
                            currColor = nodeColors[it->first];
                            nodeColors[sample] = currColor;
                        }
                    }

                }
                if(currColorLength > resultsize)
                {
                    resultsize=currColorLength;
                    result=currColor;
                }
                else if (currColorLength == resultsize && currColor<result)
                {
                    result=currColor;
                }

            }

            S.pop();
        }


    }

    return result;

}
void compressedColorColumn::optimize3(colorColumn* col)
{
    numColors=col->getNumColors();
    uint64_t  oldColorsSum=0,newColorsSum=0;
    idsMap.resize(numColors+1);

    for(unsigned int i=1; i<numColors+1; i++)
    {
        idsMap[i]=i;
        deque<uint32_t > currV;
        copy(col->colors[i].begin(),col->colors[i].end(),back_inserter(currV));

        oldColorsSum+=currV.size();
        vector<uint32_t> newColor;
        while(currV.size()>0) {
//            cout<<"input color "<<i<<":";
//            for(auto c :currV)
//                cout<<c<<" ";
//            cout<<endl;
            uint32_t subsetColor = getLongestSubsetColor(col, currV,i);
            newColor.push_back(subsetColor);
            unordered_set<uint32_t> removeList;
            if(subsetColor < col->noSamples)
                subsetColor++;

//            cout<<"subset  ";
//            for(auto c :col->colors[subsetColor])
//                cout<<c<<" ";
//            cout<<endl;

            for(auto c: col->colors[subsetColor])
                removeList.insert(c);
            deque<uint32_t > tmp;
            auto it=currV.begin();
            while(it!=currV.end())
            {
                if(removeList.find(*it)==removeList.end())
                {
                    tmp.push_back(*it);
                }
                it++;
            }
            currV=tmp;
        }

        if(i%1000==0)
            cout<<i<<endl;

//        cout<<"result  ";
//        for(auto c :newColor)
//            cout<<c<<" ";
//        cout<<endl;
        insert(newColor,i);
        newColorsSum+=newColor.size();
    }
    sdsl::util::bit_compress(idsMap);
    cout<<"old Colors sum "<<oldColorsSum<<endl;
    cout<<"New Colors sum "<<newColorsSum<<endl;
}

void  compressedColorColumn::insert(vector<uint32_t >& item,uint32_t index){
    auto it=lower_bound(colors.begin(),colors.end(),index+1,
            [](vectorBase* lhs, uint32_t rhs) -> bool { return lhs->beginID < rhs; });
   // if(it==colors.end())
        it--;
    (*it)->set(index-(*it)->beginID,item);


}
uint32_t  compressedColorColumn::insertAndGetIndex(vector<uint32_t > item){
    throw std::logic_error("insertAndGetIndex is not supported in compressedColorColumn");
    return 0;

}
vector<uint32_t > compressedColorColumn::getWithIndex(uint32_t index){
    index=idsMap[index];
    vector<uint32_t > res;
    stack<uint32_t> q;
    auto it=lower_bound(colors.begin(),colors.end(),index+1,
                        [](vectorBase* lhs, uint32_t rhs) -> bool { return lhs->beginID < rhs; });
  //  if(it==colors.end())
        it--;
    vector<uint32_t > compressedColor=(*it)->get(index-(*it)->beginID);
    for(unsigned int i=0;i<compressedColor.size();i++) {
        if(compressedColor[i]>=noSamples)
        {
            q.push(compressedColor[i]);
        }
        else{
            res.push_back(compressedColor[i]);
        }
    }
    while(q.size()>0){
        index=q.top();
       // cout<<index<<endl;
        it=lower_bound(colors.begin(),colors.end(),index+1,
                       [](vectorBase* lhs, uint32_t rhs) -> bool { return lhs->beginID < rhs; });
  //      if(it==colors.end())
            it--;
        vector<uint32_t > compressedColor=(*it)->get(index-(*it)->beginID);
        q.pop();
        for(unsigned int i=0;i<compressedColor.size();i++) {
            if(compressedColor[i]>=noSamples)
            {
                q.push(compressedColor[i]);
            }
            else{
                res.push_back(compressedColor[i]);
            }
        }

    }
    return res;
}



void compressedColorColumn::serialize(string filename)
{
    ofstream out(filename.c_str());

//    for(unsigned int i=0;i<colors.size();i++)
//        colors[i].serialize(out);
//    out.close();
//    colorInv.serialize(filename);
}
void compressedColorColumn::deserialize(string filename)
{
  //  colorInv.deserialize(filename);
   // populateColors();
}

uint64_t compressedColorColumn::sizeInBytes()
{
    uint64_t res=0;
    for(auto vec:colors)
    {
        res+=vec->sizeInBytes();
    }
    res+=sdsl::size_in_bytes(idsMap);
    return res;
}


