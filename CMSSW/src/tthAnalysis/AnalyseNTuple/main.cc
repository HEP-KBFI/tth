#include "AnalyseNtuple.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include <iostream>
#include <string>  

using namespace std;

int main ( int argc, char * argv[]) {

    string fn=argv[1]; cout << "Processing file: " << fn << endl;

    TFile f(fn.c_str()); //if (f==0) { cout << "Error: cannot open " << fn << endl;}

    TTree* tree = (TTree*) f.GetObjectChecked("tree", "TTree");

    AnalyseNtuple an(tree);

    vector<string> tags {"tth125", "dyjets", "ttjets", "ttwjets", "ttzjets", "t_t", "t_s",
	                 "t_tw", "tbar_t", "tbar_s", "tbar_tw", "ww", "wz", "zz", "singleMu", "singleEl"};

    for ( int i=0; i<tags.size(); i++){ if (fn.find(tags[i]) != std::string::npos) an.flag(i); }

    an.Loop();

    return 0;
}
