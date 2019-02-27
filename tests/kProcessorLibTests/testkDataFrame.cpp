#include "testkDataframe.h"
//add the new kDataframes here to be tested
using MyTypes = ::testing::Types<kDataFrameMQF,kDataFrameMAP>;
TYPED_TEST_SUITE(kDataFrameTest, MyTypes);

TYPED_TEST(kDataFrameTest,emptySize)
{

    kDataFrame* kframe=&this->kFrame;
    bool tmp=kframe->empty();
    EXPECT_EQ(tmp, true);
}

TYPED_TEST(kDataFrameTest,bracketOperator)
{

    kDataFrame* kframe=&this->kFrame;
    kframe.incrementCounter();
    EXPECT_EQ(tmp, true);
}
