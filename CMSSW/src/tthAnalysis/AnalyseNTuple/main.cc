#include "AnalyseNtuple.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "stdlib.h"
#include <iostream>
#include <string>  

using namespace std;

int main ( int argc, char * argv[]) {
	
	//cout << "Number of args: " << argc << endl; for(int i=0; i < argc; i++) cout << i << ": `" << argv[i] << "`" << endl;
	 
	TFile f(argv[1]); //load root file

	string fn=argv[1]; cout << "Processing file: " << fn << endl;

	TTree* tree = (TTree*) gROOT->FindObject("tree");

	AnalyseNtuple an(tree);

        string tag[]={"tth125", "dyjets", "ttjets", "ttwjets", "ttzjets", "t_t", "t_s",
                      "t_tw", "tbar_t", "tbar_s", "tbar_tw", "ww", "wz", "zz", "singleMu", "singleEl"};

	for (int i=0; i<16; i++){ if (fn.find(tag[i]) != std::string::npos) an.flag(i); }

	an.Loop();

	return EXIT_SUCCESS;
}
