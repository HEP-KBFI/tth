//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr  2 17:28:56 2014 by ROOT version 5.32/00
// from TTree tree/NTuple
// found on file: ntuple.root
//////////////////////////////////////////////////////////

#ifndef SkimNtuple_h
#define SkimNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <iostream>
#include <iomanip>
 
#include <TH1.h>
#include <TLorentzVector.h>

const int val=27;
const int gval=77;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class SkimNtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           ev;
   Int_t           np;
   Int_t           lumi;
   Int_t           run;
   Float_t         disc[val];   //[np]
   Double_t        e[val];   //[np]
   Double_t        evx;
   Double_t        evy;
   Double_t        evz;
   Int_t           gnp;
   Int_t           gpid[gval];   //[gnp]
   Double_t        ge[gval];   //[gnp]
   Int_t           gmother[gval];   //[gnp]
   Double_t        gpx[gval];   //[gnp]
   Double_t        gpy[gval];   //[gnp]
   Double_t        gpz[gval];   //[gnp]
   Int_t           gstatus[gval];   //[gnp]
   Int_t           id[val];   //[np]
   Float_t         iso[val];   //[np]
   Double_t        metcov[4];
   Double_t        metx;
   Double_t        mety;
   Double_t        npu;
   Int_t           nvtx;
   Int_t           pid[val];   //[np]
   Double_t        px[val];   //[np]
   Double_t        py[val];   //[np]
   Double_t        pz[val];   //[np]
   Double_t        rho;
   Double_t        vx[val];   //[np]
   Double_t        vy[val];   //[np]
   Double_t        vz[val];   //[np]

   // List of branches
   TBranch        *b_ev;   //!
   TBranch        *b_np;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_run;   //!
   TBranch        *b_disc;   //!
   TBranch        *b_e;   //!
   TBranch        *b_evx;   //!
   TBranch        *b_evy;   //!
   TBranch        *b_evz;   //!
   TBranch        *b_gnp;   //!
   TBranch        *b_gpid;   //!
   TBranch        *b_ge;   //!
   TBranch        *b_gmother;   //!
   TBranch        *b_gpx;   //!
   TBranch        *b_gpy;   //!
   TBranch        *b_gpz;   //!
   TBranch        *b_gstatus;   //!
   TBranch        *b_id;   //!
   TBranch        *b_iso;   //!
   TBranch        *b_metcov;   //!
   TBranch        *b_metx;   //!
   TBranch        *b_mety;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_pid;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!

   SkimNtuple(TTree *tree=0);
   virtual ~SkimNtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   //
   int nbe, nbmu, nbtau;
   int nev, ntmu, nte, ntau, nj, nbj;
   int idtmu, idte, idltau, idlj;

   std::vector<int> vgtau, vrtau;

   std::vector<int>::iterator igtau;

   TLorentzVector Lrtau, Lgtau;

   bool matchGen;

   Long64_t nentries;

   TH1F  *hpt, *h0pt, *h1pt, *h2pt, *h3pt, *h4pt, *h8pt;

   TFile *skimfile;
   TTree *skimtree;

};

#endif

#ifdef SkimNtuple_cxx
SkimNtuple::SkimNtuple(TTree *tree) : fChain(0) 
{
	Init(tree);

        nev=0; ntmu=0; nte=0; ntau=0, nj=0; nbj=0; 
	idtmu  = 1<<1;
	idte   = 1<<1;
	idltau = (1<<0) + (1<<1) + (1<<4) + (1<<8); //dlll
	idlj   = 1<<0;

	nentries =0;

	skimfile =new TFile("skim_TTH_HToTauTau_M-125.root","RECREATE");

	hpt  = new TH1F("htaupt",   "taupt",       100, 0, 300);
	h0pt = new TH1F("htaupt0",  "taupt 1<<0",  100, 0, 300);
	h1pt = new TH1F("htaupt1",  "taupt 1<<1",  100, 0, 300);
	h2pt = new TH1F("htaupt2",  "taupt 1<<2",  100, 0, 300);
	h3pt = new TH1F("htaupt3",  "taupt 1<<3",  100, 0, 300);
	h4pt = new TH1F("htaupt4",  "taupt 1<<4",  100, 0, 300);
	h8pt = new TH1F("htaupt8",  "taupt 1<<8",  100, 0, 300);

	skimtree =new TTree("tree","skim tree");
	skimtree->Branch("ev",&ev,"ev/I");
	skimtree->Branch("np",&np,"np/I");
	skimtree->Branch("lumi",&lumi,"lumi/I");
	skimtree->Branch("run",&run,"run/I");
	skimtree->Branch("disc",disc,"disc[np]/F");
	skimtree->Branch("e",e,"e[np]/D");
	skimtree->Branch("evx",&evx,"evx/D");
	skimtree->Branch("evy",&evy,"evy/D");
	skimtree->Branch("evz",&evz,"evz/D");
	skimtree->Branch("gnp",&gnp,"gnp/I");
	skimtree->Branch("gpid",gpid,"gpid[gnp]/I");
	skimtree->Branch("ge",ge,"ge[gnp]/D");
	skimtree->Branch("gmother",gmother,"gmother[gnp]/I");
	skimtree->Branch("gpx",gpx,"gpx[gnp]/D");
	skimtree->Branch("gpy",gpy,"gpy[gnp]/D");
	skimtree->Branch("gpz",gpz,"gpz[gnp]/D");
	skimtree->Branch("gstatus",gstatus,"gstatus[gnp]/I");
	skimtree->Branch("id",id,"id[np]/I");
	skimtree->Branch("iso",iso,"iso[np]/F");
	skimtree->Branch("metcov",metcov,"metcov[4]/D");
	skimtree->Branch("metx",&metx,"metx/D");
	skimtree->Branch("mety",&mety,"mety/D");
	skimtree->Branch("npu",&npu,"npu/D");
	skimtree->Branch("nvtx",&nvtx,"nvtx/I");
	skimtree->Branch("pid",pid,"pid[np]/I");
	skimtree->Branch("px",px,"px[np]/D");
	skimtree->Branch("py",py,"py[np]/D");
	skimtree->Branch("pz",pz,"pz[np]/D");
	skimtree->Branch("rho",&rho,"rho/D");
	skimtree->Branch("vx",vx,"vx[np]/D");
	skimtree->Branch("vy",vy,"vy[np]/D");
	skimtree->Branch("vz",vz,"vz[np]/D");
}

SkimNtuple::~SkimNtuple()
{
	skimtree->Write();
        skimfile->Write(); 
        skimfile->Close(); 

	//delete skimfile;
	//if (!fChain) return;
	//delete fChain->GetCurrentFile(); //crash

	//print out
	std::cout<<std::setw(20)<<"Selection"       <<std::setw(12)<<"Evt"       <<std::setw(12)<<"C.Eff.%"<<std::endl;
	std::cout<<"=========================================================="  <<std::endl;
	std::cout<<std::setw(20)<<"Initial:"        <<std::setw(12)<<nev         <<std::setw(12)<<"100"<<std::endl;
	std::cout<<std::setw(20)<<"tight muon:"     <<std::setw(12)<<ntmu        <<std::setw(12)<< 100*ntmu/nentries <<std::endl;
	std::cout<<std::setw(20)<<"tight ele:"      <<std::setw(12)<<nte         <<std::setw(12)<< 100*nte/nentries  <<std::endl;
	std::cout<<std::setw(20)<<"loose tau:"      <<std::setw(12)<<ntau        <<std::setw(12)<< 100*ntau/nentries <<std::endl;
}

Int_t SkimNtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SkimNtuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void SkimNtuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev, &b_ev);
   fChain->SetBranchAddress("np", &np, &b_np);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("disc", disc, &b_disc);
   fChain->SetBranchAddress("e", e, &b_e);
   fChain->SetBranchAddress("evx", &evx, &b_evx);
   fChain->SetBranchAddress("evy", &evy, &b_evy);
   fChain->SetBranchAddress("evz", &evz, &b_evz);
   fChain->SetBranchAddress("gnp", &gnp, &b_gnp);
   fChain->SetBranchAddress("gpid", gpid, &b_gpid);
   fChain->SetBranchAddress("ge", ge, &b_ge);
   fChain->SetBranchAddress("gmother", gmother, &b_gmother);
   fChain->SetBranchAddress("gpx", gpx, &b_gpx);
   fChain->SetBranchAddress("gpy", gpy, &b_gpy);
   fChain->SetBranchAddress("gpz", gpz, &b_gpz);
   fChain->SetBranchAddress("gstatus", gstatus, &b_gstatus);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("iso", iso, &b_iso);
   fChain->SetBranchAddress("metcov", metcov, &b_metcov);
   fChain->SetBranchAddress("metx", &metx, &b_metx);
   fChain->SetBranchAddress("mety", &mety, &b_mety);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("pid", pid, &b_pid);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   Notify();
}

Bool_t SkimNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SkimNtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SkimNtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SkimNtuple_cxx
