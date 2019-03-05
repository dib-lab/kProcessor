#include "gtest/gtest.h"
#include "kDataFrame.hpp"
#include <unordered_map>



template <typename T_kDataFrame>
class kDataFrameTest : public ::testing::Test{
public:
  static unordered_map<string,int> kmers;
  T_kDataFrame kFrame;
protected:
  virtual void SetUp();

};
