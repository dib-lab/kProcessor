#include <iostream>
#include <string>
#include "CLI11.hpp"
#include <vector>
#include <stdint.h>
#include <gqf.h>
#include "algorithms.hpp"
#include "Utils/kmer.h"
#include "Utils/utils.hpp"
#include <cmath>
#include <functional>
#include <parallel/algorithm>

using namespace std;


int main(int argc, char *argv[]){
  kDataFrame* kframe=new kDataFrameMQF(32);
  string fileName="test.fa";
  kProcessor::countKmersFromFile(kframe, {{"mode", 1}}, fileName, 1000); // Mode 1 : kmers, KmerSize will be cloned from the kFrame
  for(auto it:*kframe)
    cout<<it.kmer<<" "<<it.count<<endl;

  kDataFrame* kframe2=new kDataFrameMQF(32);
  transform(kframe->begin(), kframe->end(), inserter(*kframe2,kframe2->begin()),
              [](kmerRow c) -> kmerRow {  c.count+=5; return c ; });

  kDataFrame* kframe3=new kDataFrameMQF(32);
  copy_if(kframe2->begin(), kframe2->end(), inserter(*kframe3,kframe3->begin()),
            [](kmerRow c) -> bool { return c.kmer == "ACATGCATGACGATGCTAGCGTGATGCTAGCT";  });

  for(auto it:*kframe3)
        cout<<it.kmer<<" "<<it.count<<endl;
  vector<int> a={1,2,3,4,5};
  vector<int> b(5);
  __gnu_parallel::transform(b.begin(), b.end(), b.begin(),
                 [](int c) -> int { return c+5; });


  copy(a.begin(),a.end(),b.begin());
  std::transform(b.begin(), b.end(), b.begin(),
                 [](int c) -> int { return c+5; });
  for(auto i:b)
    cout<<i<<" ";
  cout<<endl;


  unordered_map<char,int> aa={{'a',1},{'b',2},{'c',3},{'d',4}};
  unordered_map<char,int> bb;
  //copy(aa.begin(),aa.end(),inserter(bb,bb.begin()));

  __gnu_parallel::transform(aa.begin(), aa.end(), inserter(bb,bb.begin()),
                 [](pair<char,int> c) -> pair<char,int> {  c.second+=5; return c ; });
  unordered_map<char,int> cc;
  //unordered_map<char,int> cc=std::filter(aa,[](pair<char,int> i){return i.second>=3;});
  copy_if(bb.begin(),bb.end(),inserter(cc,cc.begin()), [](pair<char,int> i){return i.second>=7;});
  for(auto i:cc)
    cout<<i.first<<" "<<i.second<<endl;




  return 0;

}
