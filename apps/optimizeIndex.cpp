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


    kDataFrame* indexFrame=kDataFrame::load(framePath);
    colorColumn* colors=(colorColumn*)indexFrame->getDefaultColumn();
    colors->optimize();
    return 0;
}
