#include "gtest/gtest.h"
#include <vector>
#include "kDataFrame.hpp"
#include "testkDataframe.h"
using namespace std;




int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);



    int ret = RUN_ALL_TESTS();
    return ret;
}
