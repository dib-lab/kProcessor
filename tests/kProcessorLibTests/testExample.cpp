#include "testExample.h"

//using ::testing::Return;

FooTest::FooTest() {

}

FooTest::~FooTest() {}

void FooTest::SetUp() {}

void FooTest::TearDown() {}

TEST_F(FooTest, onePlusOnEqualsTwo) {
    EXPECT_EQ(1+1, 2);
}

TEST_F(FooTest, twoMultiplyFiveEqualsTen) {
    EXPECT_EQ(2*5, 10);
}
