#include "testkDataframe.h"
using MyTypes = ::testing::Types<kDataFrameMQF,kDataFrameMAP>;
TYPED_TEST_SUITE(kDataFrameTest, MyTypes);

TYPED_TEST(kDataFrameTest,emptySize)
{

    kDataFrame* kframe=&this->kFrame;
    bool tmp=kframe->empty();
    EXPECT_EQ(tmp, true);
}
