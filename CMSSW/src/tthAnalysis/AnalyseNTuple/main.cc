#include "AnalyseNtuple.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char * argv[])
{
	string filename=argv[1];
	cout << "Processing file: " << filename << endl;

	TFile f(filename.c_str());
	TTree* tree = (TTree*) f.GetObjectChecked("tree", "TTree");

	AnalyseNtuple an(tree);

	vector<string> tags {
		"tth125", "dyjets", "ttjets", "ttwjets", "ttzjets", "t_t", "t_s",
		"t_tw", "tbar_t", "tbar_s", "tbar_tw", "ww", "wz", "zz", "singleMu",
		"singleEl"
	};

	for(int i=0; i<tags.size(); i++){
		if(filename.find(tags[i]) != std::string::npos) {
			an.flag(i);
		}
	}

	an.Loop();

	return EXIT_SUCCESS;
}
