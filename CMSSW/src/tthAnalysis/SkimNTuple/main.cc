#include "SkimNtuple.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "stdlib.h"
#include <fstream>
#include <iostream>
#include <string.h>

using namespace std;

void MergeTrees(TChain& tree, std::ifstream& fileFist);

int main (int argc, char * argv[]) {
        // for a single root file 
        //TFile f("/hdfs/cms/store/user/calpas/tth/ntuple/v1/ww/output_110_1_QOY.root");
        //f.cd("demo");
        //TTree* tree = (TTree*) gROOT->FindObject("tree");
	//SkimNtuple an(tree);

	//cout << "argc: " << argc << endl;
	//for(int i=0; i < argc; i++) { cout << "argv "<< i <<": "<< argv[i] << endl;}

        // for mutiple root file
	TChain tree("demo/tree");

	std::ifstream ifs (argv[1]);

	MergeTrees(tree, ifs);

	SkimNtuple an(&tree);

	an.Loop();

	return 0;
}


void MergeTrees(TChain& tree, std::ifstream& ifs){

	if (not (ifs.is_open())) {
		std::cerr << "failed in open." << std::endl; exit(EXIT_FAILURE);
	}

	char line[BUFSIZ]; int counter = 0; 

	while ( !ifs.eof() ) {
	//while ( counter<5) { //eof flag for end of file
		  ifs.getline(line, sizeof(line));
		  if (strcmp(line, "") == 0) continue;
		  TFile a(line);
		  if (!a.IsOpen()) {
			  std::cerr << "failed in open :" << line << std::endl; exit(EXIT_FAILURE);
		  }
		  counter++;
		  tree.AddFile(line);
		  //std::cout<<"file number: " << counter << std::endl;
	}
	cout << "nfiles: "<< counter << endl; // to remove last empty line

	return;
}


