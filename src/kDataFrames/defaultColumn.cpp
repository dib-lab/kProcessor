#include <defaultColumn.hpp>
#include<iostream>
#include<fstream>
#include <cereal/types/vector.hpp>
#include <cereal/types/memory.hpp>
#include "parallel_hashmap/phmap_dump.h"
#include <cereal/archives/binary.hpp>
#include <stack>
#include <iterator>
#include <regex>
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
}

void colorColumn::populateColors(){
    colorInv.populateColors(colors);
}


vector<string > StringColorColumn::getWithIndex(uint32_t index){
    vector<string > res(colors[index].size());
    for(int i=0;i<colors[index].size();i++)
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
    for(int i=0;i<size;i++)
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
    int i=0;
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
    while(S.size()>0)
    {
        if(!std::get<2>(S.top()))
        {
            std::get<2>(S.top())=true;
            color.push_back(std::get<1>(S.top()));
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