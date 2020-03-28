#include <defaultColumn.hpp>


template class vectorColumn<double>;

template <typename  T>
uint32_t  vectorColumn<T>::insertAndGetIndex(T item){
    data.push_back(item);
    return data.size()-1;
}



template <typename  T>
T vectorColumn<T>::getWithIndex(uint32_t index){
    return data[index];
}