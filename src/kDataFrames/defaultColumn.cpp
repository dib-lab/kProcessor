#include <defaultColumn.hpp>
#include<iostream>
#include<fstream>
#include <cereal/types/vector.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/archives/binary.hpp>
#include <stack>
template class vectorColumn<int>;
template class vectorColumn<bool>;
template class vectorColumn<double>;


Column* Column::getContainerByName(std::size_t hash)
{
    if(hash== typeid(vectorColumn<int>).hash_code())
    {
        return new vectorColumn<int>();
    }
    if(hash== typeid(vectorColumn<double>).hash_code())
    {
        return new vectorColumn<double>();
    }
    if(hash== typeid(vectorColumn<bool>).hash_code())
    {
        return new vectorColumn<bool>();
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
    throw logic_error("not implemented yet");
}
void colorColumn::deserialize(string filename)
{
    throw logic_error("not implemented yet");
}

void colorColumn::populateColors(){
    colors=vector<vector<uint32_t > >(colorInv.lastColor+1);
    vector<uint32_t> color;
    stack<tuple<colorNode*,uint32_t ,bool> > S;
    //S.push(make_tuple(colorInv.root,(uint32_t)0,false));
    for(auto it:colorInv.root->edges)
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


