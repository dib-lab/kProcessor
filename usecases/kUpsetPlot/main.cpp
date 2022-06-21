#include <iostream>
#include <string>
#include "CLI11.hpp"
#include "kDataFrame.hpp"
#include <vector>
#include <algorithm>
#include "algorithms.hpp"
#include <any>
using namespace std;




void upsetPlot(vector<string>& genomeFileNames, string outputFileName)
{

}

int main(int argc, char *argv[]){
    CLI::App app;
    vector<string> genomesFileNames;
    string outputFilename;


    app.add_option("-g,--genome", genomesFileNames,
                   "Genome File fasta")->required();
    app.add_option("-o,--output", outputFilename,
                   "Output text to be ploted by upsetplot")->required();


    CLI11_PARSE(app, argc, argv);


    //differntialExpression(genes_file,samples,control,outputFilename);
    return 0;
}
