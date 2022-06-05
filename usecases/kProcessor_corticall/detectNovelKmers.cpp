#include <iostream>
#include <string>
#include "CLI11.hpp"
#include "kDataFrame.hpp"
#include <vector>
#include <algorithm>
#include "algorithms.hpp"
#include <any>
#include "Utils/utils.hpp"

using namespace std;



int main(int argc, char *argv[]){
    CLI::App app;
    string kframePath;
    string inputLst;
    string sample;
    vector<string> parents;
    string progenyTreeFile;
    string outPath;
    app.add_option("-k,--kframe", kframePath,
                   "Kdataframe path")->required();

    app.add_option("-l,--list", inputLst,
                   " input list of names of kdataframes used for joining")->required();

    app.add_option("-t,--progenyTree", progenyTreeFile,
                   "tab delimited file describing the progeny tree:<sample>\t<parent cross>")->required();

    app.add_option("-s,--sample", sample,
                   "sample to be studied")->required();

    app.add_option("-p,--parent", parents,
                   "parents of the sample")->required();

    app.add_option("-o,--out", outPath,
                  "output for kDataframe")->required();


    CLI11_PARSE(app, argc, argv);

    ifstream inputLstFile(inputLst);
    string line;
    vector<string> ids;
    while(inputLstFile>>line)
    {
        ids.push_back(line);
    }
    uint32_t sampleID=std::find(ids.begin(),ids.end(),sample)-ids.begin();
    vector<uint32_t> parentsID;
    for(auto p: parents)
    {
        uint32_t pID=std::find(ids.begin(),ids.end(),p)-ids.begin();
        parentsID.push_back(pID);
    }
    unordered_map<string,string> progenyTree;
    ifstream inputProgeny(progenyTreeFile);
    string line2;
    while(inputProgeny >> line>> line2){
        progenyTree[line]=line2;
    }

    string currParent=progenyTree[sample];
    unordered_set<uint32_t> sibilingsID;
    unordered_set<uint32_t> allParentsID;
    for(auto s: progenyTree)
    {
        if(s.second == currParent)
        {
            uint32_t sid=std::find(ids.begin(),ids.end(),s.first)-ids.begin();
            sibilingsID.insert(sid);
        }
        if(s.second== "none")
        {
            uint32_t sid=std::find(ids.begin(),ids.end(),s.first)-ids.begin();
            allParentsID.insert(sid);
        }
    }


    kDataFrame* kframe= kDataFrame::load(kframePath);
    uint32_t kSize=kframe->ksize();
    unordered_map<string,Column*> newColumns;

    // renaming columns to be sample_count and index column to be color
    for(auto c: kframe->columns)
    {
        if(c.first.size()>5 && c.first.substr(0,5)=="count")
        {
            uint32_t i=atoi(c.first.substr(5,c.first.size()-5).c_str());
            string newName=ids[i]+"_count";
            newColumns[newName]=c.second;
        }
        else if(c.first.size()>5 && c.first.substr(0,5)=="color")
        {
            string newName="color";
            newColumns[newName]=c.second;
        }
        else{
            newColumns[c.first]=c.second;
        }
    }
    kframe->columns=newColumns;
    string novelColumnName= sample+"_novel";
    string countColumnName=sample+"_count";

    kframe->addColumn(novelColumnName,new vectorColumn<bool>(kframe->size()+1));
    kProcessor::transformInPlace(kframe,[&](kDataFrameIterator& it)
    {
        //detect novel Kmers
        vector<uint32_t> color;
        it.getColumnValue<deduplicatedColumn<prefixTrie> >("color",color);
        bool hasSample=std::binary_search(color.begin(),color.end(),sampleID);
        bool hasParent=false;
        for(auto p: parentsID)
            hasParent |= std::binary_search(color.begin(),color.end(),p);
        bool isNovel=hasSample && !hasParent;
        it.setColumnValue<vectorColumn<bool>>(novelColumnName,isNovel );
        if(!isNovel) return;

        // filter low cov
        uint32_t count;
        it.getColumnValue<vectorColumn<uint32_t >>(novelColumnName,count);
        bool highCov =  count>=6;
        isNovel&= highCov;
        it.setColumnValue<vectorColumn<bool>>(novelColumnName,isNovel );
        if(!isNovel) return;


        //filter shared kmers
        bool sharedKmer=false;
        for(auto c : color)
        {
            bool isParent=allParentsID.find(c) != allParentsID.end();
            bool isSibling= sibilingsID.find(c) != sibilingsID.end();
            sharedKmer|= !(isParent || isSibling);
        }
        isNovel &= !sharedKmer;
        it.setColumnValue<vectorColumn<bool>>(novelColumnName,isNovel );
        if(!isNovel) return;



    });


    auto kframe2=kProcessor::filter(kframe,[&](kDataFrameIterator& k)
    {
        bool isNovel;
        k.getColumnValue<vectorColumn<bool>>(novelColumnName,isNovel);
        return isNovel;
    });
    cout<<"Number of Novel Kmers Found = "<<kframe2->size() <<endl;



    kframe2->save(outPath);





    return 0;
}
