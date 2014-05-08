#include "SkimNtuple.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "stdlib.h"
#include <fstream>
#include <iostream>
#include <string.h>

void MergeTrees(TChain& tree, std::ifstream& filelist);

int main () {

  //several files
  TChain tree("demo/tree");
  //orig
  std::ifstream filelist("FileToSkim/TTH_HToTauTau_M-125.txt"); 
  //std::ifstream filelist("FileToSkim/TTJets.txt");

  MergeTrees(tree, filelist);

  SkimNtuple an(&tree);
  an.Loop();


/*
//With one file
TFile f("FileToSkim/TTH_HToTauTau_M-125.root");
//TFile f("/hdfs/cms/store/user/calpas/TTH/slim/v3/TTH_HToTauTau_M-125_8TeV_pythia6/TTH_HToTauTau_M-125_27_1_I7s.root");
f.cd("demo");
//f.ls();
//gDirectory->pwd();
//gFile->pwd();
 TTree* tree = (TTree*) gROOT->FindObject("tree");
 SkimNtuple an(tree);
 an.Loop();
*/

 return EXIT_SUCCESS;
}


void
MergeTrees(TChain& tree, std::ifstream& filelist)
{
  if (not (filelist.is_open())) {
    std::cerr << "failed in open." << std::endl;
    exit(EXIT_FAILURE);
  }

  char line[BUFSIZ];
  int counter = 0; 

 //while (not (filelist.eof())) { //eof flag for end of file
  while ( counter<1) { //eof flag for end of file
    filelist.getline(line, sizeof(line));

    if (strcmp(line, "") == 0) continue;

    TFile a(line);
    if (!a.IsOpen()) {
      std::cerr << "failed in open :" << line << std::endl;
      exit(EXIT_FAILURE);
    }
    
    //ERROR MESSAGE FOR NOTHING!!
    //   if (gROOT->FindObject(tree.GetName()) == 0) {
    //   //if (gROOT->FindObject("tree") == 0) {
    //     std::cerr << tree.GetName() << " is not found in :" << line << std::endl;
    //     exit(EXIT_FAILURE);
    //   }

    counter ++ ;
    tree.AddFile(line); 
    std::cout << "file number: "<< counter << std::endl;
}
  std::cout << counter << " files are loaded." << std::endl;
 
 return;
}
