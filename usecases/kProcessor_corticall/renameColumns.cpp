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

    app.add_option("-k,--kframe", kframePath,
                   "Kdataframe path")->required();

    app.add_option("-l,--list", inputLst,
                   " input list of names of kdataframes used for joining")->required();



    CLI11_PARSE(app, argc, argv);

    ifstream inputLstFile(inputLst);
    string line;
    vector<string> ids;
    while(inputLstFile>>line)
    {
        ids.push_back(line);
    }

    kDataFrame* kframe= kDataFrame::load(kframePath);
    uint32_t kSize=kframe->ksize();
    unordered_map<string,Column*> newColumns;
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

    for(auto c: kframe->columns)
    {
        cout<<c.first<<" "<<c.second->size()<<endl;
    }

    kframe->save(kframePath);



    return 0;
}
