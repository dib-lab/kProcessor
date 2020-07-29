//
// Created by mostafa on 4/6/20.
//

#include <iostream>
#include <string>
#include "kDataFrame.hpp"
#include <vector>
#include "algorithms.hpp"
using namespace std;

int main(int argc, char *argv[])
{
    string framePath=argv[1];

    queryColorColumn* col=new queryColorColumn();
    col->deserialize(framePath);
    cout<<"Num of Colors "<<col->getNumColors()<<endl;
    col->serialize(framePath+".new");
    return 0;
}
