#include <iostream>
#include <string>
#include "ThirdParty/CLI11.hpp"
#include <vector>
#include <stdint.h>
#include <gqf.h>
#include "KmerCounter/KmerCounter.hpp"
#include "KmerCounter/kmer.h"
#include "fstream"
#include <cmath>
#include <map>
#include <seqan/seq_io.h>
#include "kDataFrame.hpp"
#include<algorithm>
using namespace std;

class statistics{
public:
  uint64_t n_unmatched;
  uint64_t n_same;
  uint64_t n_amb_same;
  uint64_t n_clear_fusion;
  uint64_t  n_ambig_fusion;
  uint64_t n_mutli_fusion;
  uint64_t n_paired_fusion;
  uint64_t n_fragments;
  statistics(){
    n_fragments=0;
    n_unmatched = 0;
    n_same = 0;
    n_amb_same = 0;
    n_clear_fusion = 0;
    n_ambig_fusion = 0;
    n_mutli_fusion = 0;
    n_paired_fusion= 0;
  }
  void print()
  {
    cout<<"No of input fragments: "<< n_fragments<<endl;
    cout<<"unmatched:"<< n_unmatched<<endl;
    cout<<"Unique:"<< n_same<<endl;
    cout<<"Ambiguous:"<< n_amb_same<<endl;
    cout<<"Single read clear fusion:"<< n_clear_fusion<<endl;
    cout<<"Single read ambiguous fusion:"<< n_ambig_fusion<<endl;
    cout<<"Single read multi fusion:"<< n_mutli_fusion<<endl;
    cout<<"paired read fusion:"<< n_paired_fusion<<endl;

  }
};

inline string vecToString(std::vector<int> v)
{
  if(v.size()==0){
    return "[]";
  }
  string res="["+to_string(v[0]);
  for(int i=1;i<v.size();i++)
  {
    res+=" ,"+to_string(v[i]);
  }
  res+="]";
  return res;
}
inline string setsToString(std::vector<int> v)
{
  if(v.size()==0){
    return "{}";
  }
  string res="{"+to_string(v[0]);
  for(int i=1;i<v.size();i++)
  {
    res+=" ,"+to_string(v[i]);
  }
  res+="}";
  return res;
}
inline string familiesToString(vector<vector<int>  > v)
{
  if(v.size()==0){
    return "[]";
  }
  string res="["+setsToString(v[0]);
  for(int i=1;i<v.size();i++)
  {
    res+=" ,"+setsToString(v[i]);
  }
  res+="]";
  return res;
}


static string readFusion(string read,kDataFrame* DB,vector<vector<int> >* families,
  vector<int>* shared_kmers,vector<int>*gaps,statistics *stats)
{
  shared_kmers->clear();
  families->clear();
  gaps->clear();
  string flag;
  vector<int> lf_ids;
  vector<int> rt_ids;
  uint8_t kSize=DB->getkSize();
  if(read.size()<kSize){
    stats->n_unmatched++;
    return "unmatched";
  }
  //# find a matching k-mer at the beginning of the read
  lf_ids=DB->getColors(read.substr(0,kSize));
  int idx = 1;

  while (idx < read.size()-kSize+1 && lf_ids.size() == 0){
      lf_ids=DB->getColors(read.substr(idx,kSize));
      idx++;
    }
  if(lf_ids.size() == 0){
      //#print('no single match')
      stats->n_unmatched++;
      flag = "unmatched";
  }
  else if(idx == read.size()-kSize+1){
      //#print('same, only last kmer matched')
      families->push_back(lf_ids);
      if(lf_ids.size() == 1){
          stats->n_same += 1;
          flag = "unique";
        }
      else{
          stats->n_amb_same += 1;
          flag = "ambiguous";
        }
  }
  else{
    // # len(lf_ids) > 0 & idx < len(hashvals)
      //# find a matching k-mer at the end of the read
      vector<int> rt_ids=DB->getColors(read.substr(read.size()-kSize,kSize));

      int idy = read.size() - 1;
      while(idy-kSize >= idx-1 && rt_ids.size() == 0){
          rt_ids=DB->getColors(read.substr(idy-kSize,kSize));
          idy --;
        }

      if(rt_ids.size() == 0){
      //    #print('same, only one non-last kmer matched ')
          families->push_back(lf_ids);
          if(lf_ids.size() == 1){
              stats->n_same += 1;
              flag = "unique";
            }
          else{
              stats->n_amb_same += 1;
              flag = "ambiguous";
            }
      }
    else{
          vector<int> intersect_ids;
          intersect_ids.clear();
          auto it=set_intersection(lf_ids.begin(),lf_ids.end(),rt_ids.begin(),rt_ids.end(),back_inserter(intersect_ids));
          if(intersect_ids.size() > 0){
              families->push_back(intersect_ids);
              if(intersect_ids.size() == 1){
                  stats->n_same += 1;
                  flag = "unique";
                }
              else{
                  stats->n_amb_same += 1;
                  flag = "ambiguous";
                }
          }
          else{
          //# fusion to be resolved
              uint64_t shared_kmer = 1;
              uint64_t gap_size = 0;
              bool gap = false;
              while (idx <= (idy+1-kSize)){
                  vector<int> temp_ids=DB->getColors(read.substr(idx,kSize));
                  if(read.substr(idx,kSize).size()!=kSize)
                  {
                    cout<<read.substr(idx,kSize)<<" "<<idx<<" < "<<read.size()<<endl;
                  }
                  if (temp_ids.size() > 0){
                      intersect_ids.clear();
                      it=set_intersection(lf_ids.begin(),lf_ids.end(),temp_ids.begin(),temp_ids.end(),back_inserter(intersect_ids));

                      if (intersect_ids.size()> 0){
                          lf_ids = intersect_ids;
                          shared_kmer += 1;
                          gap_size = 0;
                        }
                      else{
                        // # len(intersect_ids) == 0

                          families->push_back(lf_ids);
                          shared_kmers->push_back(shared_kmer);
                          lf_ids = temp_ids;
                          shared_kmer = 1;
                          gaps->push_back(gap_size);
                          gap_size = 0;
                        }
                  }
                  else{
                      gap_size += 1;
                  }
                  idx += 1;
              }
              families->push_back(lf_ids);

              shared_kmers->push_back(shared_kmer);

              if(families->size()<=1)
              {
                cerr<<"Error"<<endl;
                return "Error";
              }
              if(families->size() == 2){
                  if((*families)[0].size() == 1 && (*families)[1].size() == 1){
                      stats->n_clear_fusion += 1;
                      flag = "clear_fusion";
                  }
                  else{
                      stats->n_ambig_fusion += 1;
                      flag = "ambig_fusion";
                    }
              }
              else{
               //# len(families) > 2
                  stats->n_mutli_fusion += 1;
                  flag = "multi_fusion";
                }

          }
      }
   }

  //#if len(families) == 0:
  //#    families = "-"

  //#if len(shared_kmers) == 0:
  //#    shared_kmers = "-"

  return flag;

}
inline bool isFusion(string flag){
  return flag=="clear_fusion" || flag=="ambig_fusion" || flag == "multi_fusion";
}
inline bool isSameRef(string flag){
  return flag=="unique" || flag=="ambiguous";
}
int detectFusionGenes_main(int argc, char *argv[]){
  CLI::App app;
  string dbFile;
  int numThreads=1;
  int chunkSize=1000;
  vector<string> singleEndReads;
  vector<string> pairedEndReads;
  string outputPrefix;

  app.add_option("-d,--database", dbFile,
   "kDataFrame containing cDGB of the transpictome.")
   ->required();

  app.add_option("-s,--single-end",singleEndReads,"Single end reads file list")
  ->check(CLI::ExistingFile);

  app.add_option("-p,--paired-end",pairedEndReads,"Paired end reads file list. File for each samplibng containing both ends.")
  ->check(CLI::ExistingFile);


  app.add_option("-o,--output-prefix", outputPrefix,
   "Output Prefix.")
   ->required();

  app.add_option("-t,--threads", numThreads,
   "Number of threads used in kmer counting. Default 1");

  CLI11_PARSE(app, argc, argv);

  if(pairedEndReads.size()==0 &&singleEndReads.size()==0)
  {
    cerr<<"No reads are supplied to -s or -p"<<endl;
    return -1;
  }
  ofstream outFusionInfo(outputPrefix+"_fusion.info");
  outFusionInfo<<"fileName"<<"\t"<<"recordIndex"<<"\t"<<"whichInPair"<<
  "\t"<<"align_class"<<"\t"<<"gene_families"<<"\t"<<"shared_kmers"<<"\t"<<"gaps"<<endl;
  
  seqan::SeqFileOut fusion_fp((outputPrefix+"_fusion.fq").c_str());

  ofstream outFusionCalc(outputPrefix+"_fusion.calc");
  outFusionCalc<<"fileName"<<"\t"<<"recordIndex"<<"\t"<<"whichInPair"<<
  "\t"<<"align_class"<<"\t"<<"familiy_A"<<"\t"<< "familiy_B"<<"\t"<<"no_families"
  <<"\t"<< "len_families"<<"\t"<< "shared_kmers"<<"\t"<< "gaps"<<"\t"<< "sorted_keys"<<endl;

  ofstream outFusionPairInfo(outputPrefix+"_fusionPair.info");
  outFusionPairInfo<<"fileName"<<"\t"<<"recordIndex"<<"\t"<<"fusion_class"<<"\t"
  <<"R1_family"<<"\t"<<"R2_family"<<endl;

  seqan::SeqFileOut fusionPair_fp((outputPrefix+"_fusionPair.fq").c_str());

  ofstream outFusionPairCalc(outputPrefix+"_fusionPair.calc");
  outFusionPairCalc<<"fileName"<<"\t"<<"recordIndex"<<"\t"<<"fusion_class"<<"\t"
  <<"familiy_A"<<"\t"<<"familiy_B"<<"\t"<<"len_families"<<"\t"<<"sorted_keys"<<endl;

  cout<<"Start Loading the Database."<<endl;
  kDataFrame* DB=kDataFrame::load(dbFile+".mqf");
  ifstream inputNamesMap(dbFile+".namesMap");
  string name;
  uint64_t nameTag;
  map<uint64_t,string> namesMap;
  while(inputNamesMap>>name>>nameTag){
    namesMap.insert(make_pair(nameTag,name));
  }
  cout<<"Finished Loading."<<endl;
  if(chunkSize%2==1)
  {
    cerr<<"chunkSize should be an even number"<<endl;
    return -1;
  }

  statistics stats;
  vector<vector<int> > families0,families1,families;
  vector<int> shared_kmers0,gaps0,shared_kmers1,gaps1,shared_kmers,gaps;


  for(auto readsPath : singleEndReads){
      seqan::SeqFileIn seqIn(readsPath.c_str());
      seqan::StringSet<seqan::CharString> ids;
      seqan::StringSet<seqan::CharString> reads;
      seqan::StringSet<seqan::CharString> quals;
      int r_index=0;
      while(!atEnd(seqIn)){
        clear(reads);
        clear(ids);
        clear(quals);
        seqan::readRecords(ids, reads,quals,seqIn,chunkSize);
        stats.n_fragments+=length(reads);
        for(int i=0;i<seqan::length(reads);i+=1){
          string read=string((char*)seqan::toCString(reads[i]));
          families.clear();
          gaps.clear();
          shared_kmers.clear();
          string flag0=readFusion(read,DB,&families,&shared_kmers,&gaps, &stats);
          // if(flag=="ambig_fusion")
          //   cout<<ids[i]<<endl;
          if(isFusion(flag0)){
            outFusionInfo<<readsPath<<"\t"<< r_index<<"\t"<< "single"<<"\t"<< flag0<<"\t"<<
            familiesToString(families)<<"\t"<<vecToString(shared_kmers)<<"\t"
            << vecToString(gaps)<<"\n";

            writeRecord(fusion_fp,ids[i],reads[i],quals[i]);
            int i=families.size()-1;
            for(auto g1:families[0]){
              string g1_name=namesMap.find(g1)->second;
              for(auto g2:families[i]){
                string g2_name=namesMap.find(g2)->second;

                vector<int> familiesSizes;
                familiesSizes.clear();
                for(auto a :families)
                  familiesSizes.push_back(a.size());

                vector<int> gs;
                gs.clear();
                gs.push_back(g1);
                gs.push_back(g2);
                sort(gs.begin(),gs.end());

                outFusionCalc<<readsPath<<"\t"<< r_index<<"\t"<< "single"<<"\t"<<
                 flag0<<"\t"<< g1<<":"<<g1_name<<"\t"<< g2<<":"<<g2_name<<"\t"<<
                 families.size()<<"\t"<<vecToString(familiesSizes)<<"\t"<<
                 vecToString(shared_kmers)<<"\t"<< vecToString(gaps)<<"\t"<<
                 vecToString(gs)<<"\n";
              }
            }

          }

          r_index++;
        }
      }


  }

  for(auto readsPath : pairedEndReads){
      seqan::SeqFileIn seqIn(readsPath.c_str());
      seqan::StringSet<seqan::CharString> ids;
      seqan::StringSet<seqan::CharString> reads;
      seqan::StringSet<seqan::CharString> quals;
      int r_index=0;
      while(!atEnd(seqIn)){
        clear(reads);
        clear(ids);
        clear(quals);
        seqan::readRecords(ids, reads,quals,seqIn,chunkSize);
        stats.n_fragments+=length(reads);
        for(int i=0;i<seqan::length(reads);i+=2){
          string read1=string((char*)seqan::toCString(reads[i]));
          string read2=string((char*)seqan::toCString(reads[i+1]));
          families0.clear();
          gaps0.clear();
          shared_kmers0.clear();
          families1.clear();
          gaps1.clear();
          shared_kmers1.clear();

          string flag0=readFusion(read1,DB,&families0,&shared_kmers0,&gaps0, &stats);
          string flag1=readFusion(read2,DB,&families1,&shared_kmers1,&gaps1, &stats);
          // if(flag=="ambig_fusion")
          //   cout<<ids[i]<<endl;

          if(isFusion(flag0) || isFusion(flag1)){

            outFusionInfo<<readsPath<<"\t"<< r_index<<"\t"<< "Read_1"<<"\t"<< flag0<<"\t"<<
            familiesToString(families0)<<"\t"<<vecToString(shared_kmers0)<<"\t"
            << vecToString(gaps0)<<"\n";

            writeRecord(fusion_fp,ids[i],read1,quals[i]);

            outFusionInfo<<readsPath<<"\t"<< r_index+1<<"\t"<< "Read_2"<<"\t"<< flag1<<"\t"<<
            familiesToString(families1)<<"\t"<<vecToString(shared_kmers1)<<"\t"
            << vecToString(gaps1)<<"\n";

            writeRecord(fusion_fp,ids[i+1],read2,quals[i+1]);

            if(isFusion(flag0)){

              int i=families0.size()-1;
              for(auto g1:families0[0]){
                string g1_name=namesMap.find(g1)->second;
                for(auto g2:families0[i]){
                  string g2_name=namesMap.find(g2)->second;

                  vector<int> familiesSizes;
                  familiesSizes.clear();
                  for(auto a :families0)
                    familiesSizes.push_back(a.size());

                  vector<int> gs;
                  gs.clear();
                  gs.push_back(g1);
                  gs.push_back(g2);
                  sort(gs.begin(),gs.end());

                  outFusionCalc<<readsPath<<"\t"<< r_index<<"\t"<< "Read_1"<<"\t"<<
                   flag0<<"\t"<< g1<<":"<<g1_name<<"\t"<< g2<<":"<<g2_name<<"\t"<<
                   families0.size()<<"\t"<<vecToString(familiesSizes)<<"\t"<<
                   vecToString(shared_kmers0)<<"\t"<< vecToString(gaps0)<<"\t"<<
                   vecToString(gs)<<"\n";
                }
              }
            }
            if(isFusion(flag1)){

              int i=families1.size()-1;
              for(auto g1:families1[0]){
                string g1_name=namesMap.find(g1)->second;
                for(auto g2:families1[i]){
                  string g2_name=namesMap.find(g2)->second;

                  vector<int> familiesSizes;
                  familiesSizes.clear();
                  for(auto a :families1)
                    familiesSizes.push_back(a.size());

                  vector<int> gs;
                  gs.clear();
                  gs.push_back(g1);
                  gs.push_back(g2);
                  sort(gs.begin(),gs.end());

                  outFusionCalc<<readsPath<<"\t"<< r_index+1<<"\t"<< "Read_2"<<"\t"<<
                   flag1<<"\t"<< g1<<":"<<g1_name<<"\t"<< g2<<":"<<g2_name<<"\t"<<
                   families1.size()<<"\t"<<vecToString(familiesSizes)<<"\t"<<
                   vecToString(shared_kmers1)<<"\t"<< vecToString(gaps1)<<"\t"<<
                   vecToString(gs)<<"\n";
                }
              }
            }

          }
          else if(isSameRef(flag0) && isSameRef(flag1)){
            vector<int> intersect_ids;
            intersect_ids.clear();
            auto it=set_intersection(families0[0].begin(),families0[0].end(),
            families1[0].begin(),families1[0].end(),back_inserter(intersect_ids));
            if(intersect_ids.size()==0)
            {
              stats.n_paired_fusion++;
              string fusion_class="ambig_fusion";
              if(flag0 == "unique" && flag1 == "unique")
                  fusion_class = "clear_fusion";


              outFusionPairInfo<<readsPath<<"\t"<< r_index<<"\t"<< fusion_class
              <<"\t"<< familiesToString(families0)<<"\t"<< familiesToString(families1)<<"\n";

              writeRecord(fusionPair_fp,ids[i],read1,quals[i]);
              writeRecord(fusionPair_fp,ids[i+1],read2,quals[i+1]);

              for(auto g1:families0[0]){
                string g1_name=namesMap.find(g1)->second;
                for(auto g2:families1[0]){
                  string g2_name=namesMap.find(g2)->second;

                  vector<int> familiesSizes;
                  familiesSizes.clear();
                  familiesSizes.push_back(families0[0].size());
                  familiesSizes.push_back(families1[0].size());

                  vector<int> gs;
                  gs.clear();
                  gs.push_back(g1);
                  gs.push_back(g2);
                  sort(gs.begin(),gs.end());

                  outFusionPairCalc<<readsPath<<"\t"<< r_index<<"\t"<< fusion_class<<
                  "\t"<< g1<<":"<<g1_name<<"\t"<< g2<<":"<<g2_name<<"\t"<<
                  vecToString(familiesSizes)<<"\t"<<vecToString(gs)<<"\n";

                }
              }


            }

          }

          r_index+=2;
        }
      }




  }




  stats.print();



  seqan::close(fusion_fp);
  outFusionInfo.close();

  return 0;
}
