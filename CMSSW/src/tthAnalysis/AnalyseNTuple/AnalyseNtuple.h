//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr  2 17:28:56 2014 by ROOT version 5.32/00
// from TTree tree/NTuple
// found on file: ntuple.root
//////////////////////////////////////////////////////////

#ifndef AnalyseNtuple_h
#define AnalyseNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TH2.h>
#include <TLorentzVector.h>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip> 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class AnalyseNtuple {
    public :
	TTree          *fChain;   //!pointer to the analyzed TTree or TChain
	Int_t           fCurrent; //!current Tree number in a TChain

	// Declaration of leaf types
	Int_t           ev;
	Int_t           np;
	Int_t           lumi;
	Int_t           run;
	Float_t         disc[27];   //[np]
	Double_t        e[27];   //[np]
	Double_t        evx;
	Double_t        evy;
	Double_t        evz;
	Int_t           gnp;
	Int_t           gpid[77];   //[gnp]
	Double_t        ge[77];   //[gnp]
	Int_t           gmother[77];   //[gnp]
	Double_t        gpx[77];   //[gnp]
	Double_t        gpy[77];   //[gnp]
	Double_t        gpz[77];   //[gnp]
	Int_t           gstatus[77];   //[gnp]
	Int_t           id[27];   //[np]
	Float_t         iso[27];   //[np]
	Double_t        metcov[4];
	Double_t        metx;
	Double_t        mety;
	Double_t        npu; //was Double_t
	Int_t           nvtx;
	Int_t           pid[27];   //[np]
	Double_t        px[27];   //[np]
	Double_t        py[27];   //[np]
	Double_t        pz[27];   //[np]
	Double_t        rho;
	Int_t           trig;
	Double_t        vx[27];   //[np]
	Double_t        vy[27];   //[np]
	Double_t        vz[27];   //[np]

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
	TBranch        *b_trig;   //!
	TBranch        *b_vx;   //!
	TBranch        *b_vy;   //!
	TBranch        *b_vz;   //!

	AnalyseNtuple(TTree *tree=0);
	virtual ~AnalyseNtuple();
	virtual Int_t    Cut(Long64_t entry);
	virtual Int_t    GetEntry(Long64_t entry);
	virtual Long64_t LoadTree(Long64_t entry);
	virtual void     Init(TTree *tree);
	virtual void     Loop();
	virtual Bool_t   Notify();
	virtual void     Show(Long64_t entry = -1);

	Long64_t nentries;
	int nev, nlmu, nvmu, nte, nle, nve, nttau, nltau;
	int nbtmu, nblmu, nbvmu, nbte, nble, nbve, nbbj, nbttau, nbltau, nbtau, cj, cbj;

	double metpt;

	// efficiency
	int nsmutrig, nseltrig, nmuORe, nltmue1, nltmue2, n1tmu, nolmu, noel, ntmu, nmet, njet, n2b, ntau, n2tau, n2mu, n2e, n2lep, n2tmu, n2te; //all channel
	int ntaujet, n2j, ngj, ngb, nalljet, nsumltmue, ntrash;

	bool findtau;
	bool ra, rb, rc, rd;
	std::vector<bool> region;

	int idtmu, idlmu, idvmu, idte, idle, idve, idtj, idlj, dlll, dmll, matchTauGen, matchGen; 

	bool istmu, islmu, isvmu, iste, isle, isve, isj, islb, istb, isttau, isltau, istau;
	bool is2tmu, is2te, istmue, istlmu, istle, istmule, istelmu;
	bool isej, ismuj;

	std::vector<int> vtau;
	std::vector<int> vtmu, vlmu, vvmu;
	std::vector<int> vte, vle, vve;
	std::vector<int> vj, vjclean, vbclean;

	TLorentzVector Ltmu, Lte, Ltau,  Ltau0, Ltau1, Lj, Lmet;
	TLorentzVector Lb0, Lb1, Lj0, Lj1;
	TLorentzVector m2_Ltmu0, m2_Ltmu1, m2_Lj0, m2_Lj1, m2_Ltau0;

	// taus
	TH1F *htau0pt_sm[4], *htau0eta_sm[4]; TH1F *htau1pt_sm[4], *htau1eta_sm[4], *hdr2tau_sm[4];
	TH1F *htau0pt_se[4], *htau0eta_se[4]; TH1F *htau1pt_se[4], *htau1eta_se[4], *hdr2tau_se[4];
	// jet
	TH1F *hzeta_sm[4], *hzeta_se[4]; 
	TH1F *hjmult_sm[4], *hjpt_sm[4], *hjeta_sm[4], *hj0pt_sm[4], *hj0eta_sm[4], *hj0e_sm[4], *hj1pt_sm[4], *hj1eta_sm[4], *hj1e_sm[4], *h2jm_sm[4], *hjdr01_sm[4];
	TH1F *hbmult_sm[4], *hbpt_sm[4], *hb0pt_sm[4], *hb0eta_sm[4], *hb0e_sm[4], *hb1pt_sm[4], *hb1eta_sm[4], *hb1e_sm[4], *h2bm_sm[4], *hbdr01_sm[4];
	TH1F *hjmult_se[4], *hjpt_se[4], *hjeta_se[4], *hj0pt_se[4], *hj0eta_se[4], *hj0e_se[4], *hj1pt_se[4], *hj1eta_se[4], *hj1e_se[4], *h2jm_se[4], *hjdr01_se[4];
	TH1F *hbmult_se[4], *hbpt_se[4], *hb0pt_se[4], *hb0eta_se[4], *hb0e_se[4], *hb1pt_se[4], *hb1eta_se[4], *hb1e_se[4], *h2bm_se[4], *hbdr01_se[4];
	// lepton
	TH1F *hmupt_sm[4], *hmueta_sm[4], *hmuphi_sm[4], *hept_se[4], *heeta_se[4], *hephi_se[4];
	// met
	TH1F *hmetpt_sm[4], *hmetpt_se[4]; 
	// vtx
	TH1F *hnvtx_sm[4], *hnvtx_se[4];
	// DR
	TH1F *hdrtaumu_sm[4], *hdrtaumu_se[4];   

	TFile *file;

	// sample name
	int sID;
	int flag (int i) { sID = 1<<i; return sID; } 

	// pileup reweight
	edm::LumiReWeighting LumiWeights_;
	float pu;
	double w, wpu;
};

#endif

#ifdef AnalyseNtuple_cxx
AnalyseNtuple::AnalyseNtuple(TTree *tree) : fChain(0) 
{
    Init(tree);

    file = new TFile("output.root", "RECREATE");
    file->cd();

    std::string nvtx[8]={"nvtx_sma", "nvtx_smb", "nvtx_smc", "nvtx_smd", "nvtx_sea", "nvtx_seb", "nvtx_sec", "nvtx_sed"};
    std::string met[8]={"metpt_sma", "metpt_smb", "metpt_smc", "metpt_smd", "metpt_sea", "metpt_seb", "metpt_sec", "metpt_sed"};
    std::string tau0pt[8]={"tau0pt_sma", "tau0pt_smb", "tau0pt_smc", "tau0pt_smd", "tau0pt_sea", "tau0pt_seb", "tau0pt_sec", "tau0pt_sed" };
    std::string tau0eta[8]={"tau0eta_sma", "tau0eta_smb", "tau0eta_smc", "tau0eta_smd", "tau0eta_sea", "tau0eta_seb", "tau0eta_sec", "tau0eta_sed"};
    std::string tau1pt[8]={"tau1pt_sma", "tau1pt_smb", "tau1pt_smc", "tau1pt_smd", "tau1pt_sea", "tau1pt_seb", "tau1pt_sec", "tau1pt_sed"};
    std::string tau1eta[8]={"tau1eta_sma", "tau1eta_smb", "tau1eta_smc", "tau1eta_smd", "tau1eta_sea", "tau1eta_seb", "tau1eta_sec", "tau1eta_sed"};
    std::string dr2tau[8]={"dr2tau_sma", "dr2tau_smb", "dr2tau_smc", "dr2tau_smd", "dr2tau_sea", "dr2tau_seb", "dr2tau_sec", "dr2tau_sed"};
    std::string drtaumu[8]={"drtaumu_sma", "drtaumu_smb", "drtaumu_smc", "drtaumu_smd", "drtaumu_sea", "drtaumu_seb", "drtaumu_sec", "drtaumu_sed"};
    std::string jmult[8]={"jmult_sma", "jmult_smb", "jmult_smc", "jmult_smd", "jmult_sea", "jmult_seb", "jmult_sec", "jmult_sed"};
    std::string jpt[8]={"jpt_sma", "jpt_smb", "jpt_smc", "jpt_smd", "jpt_sea", "jpt_seb", "jpt_sec", "jpt_sed"};
    std::string jeta[8]={"jeta_sma", "jeta_smb", "jeta_smc", "jeta_smd", "jeta_sea", "jeta_seb", "jeta_sec", "jeta_sed"};
    std::string j0pt[8]={"j0pt_sma", "j0pt_smb", "j0pt_smc", "j0pt_smd", "j0pt_sea", "j0pt_seb", "j0pt_sec", "j0pt_sed"};
    std::string j0eta[8]={"j0eta_sma", "j0eta_smb", "j0eta_smc", "j0eta_smd", "j0eta_sea", "j0eta_seb", "j0eta_sec", "j0eta_sed"};
    std::string j0e[8]={"j0e_sma", "j0e_smb", "j0e_smc", "j0e_smd", "j0e_sea", "j0e_seb", "j0e_sec", "j0e_sed"};
    std::string j1pt[8]={"j1pt_sma", "j1pt_smb", "j1pt_smc", "j1pt_smd", "j1pt_sea", "j1pt_seb", "j1pt_sec", "j1pt_sed"};
    std::string j1eta[8]={"j1eta_sma", "j1eta_smb", "j1eta_smc", "j1eta_smd", "j1eta_sea", "j1eta_seb", "j1eta_sec", "j1eta_sed"};
    std::string j1e[8]={"j1e_sma", "j1e_smb", "j1e_smc", "j1e_smd", "j1e_sea", "j1e_seb", "j1e_sec", "j1e_sed"};
    std::string zeta[8]={"zeta_sma", "zeta_smb", "zeta_smc", "zeta_smd", "zeta_sea", "zeta_seb", "zeta_sec", "zeta_sed"};
    std::string jm[8]={"2jm_sma", "2jm_smb", "2jm_smc", "2jm_smd", "2jm_sea", "2jm_seb", "2jm_sec", "2jm_sed"};
    std::string jdr01[8]={"jdr01_sma", "jdr01_smb", "jdr01_smc", "jdr01_smd", "jdr01_sea", "jdr01_seb", "jdr01_sec", "jdr01_sed"};
    std::string bmult[8]={"bmult_sma", "bmult_smb", "bmult_smc", "bmult_smd", "bmult_sea", "bmult_seb", "bmult_sec", "bmult_sed"};
    std::string bpt[8]={"bpt_sma", "bpt_smb", "bpt_smc", "bpt_smd", "bpt_sea", "bpt_seb", "bpt_sec", "bpt_sed"};
    std::string b0pt[8]={"b0pt_sma", "b0pt_smb", "b0pt_smc", "b0pt_smd", "b0pt_sea", "b0pt_seb", "b0pt_sec", "b0pt_sed"};
    std::string b0eta[8]={"b0eta_sma", "b0eta_smb", "b0eta_smc", "b0eta_smd", "b0eta_sea", "b0eta_seb", "b0eta_sec", "b0eta_sed"};
    std::string b0e[8]={"b0e_sma", "b0e_smb", "b0e_smc", "b0e_smd", "b0e_sea", "b0e_seb", "b0e_sec", "b0e_sed"};
    std::string b1pt[8]={"b1pt_sma", "b1pt_smb", "b1pt_smc", "b1pt_smd", "b1pt_sea", "b1pt_seb", "b1pt_sec", "b1pt_sed"};
    std::string b1eta[8]={"b1eta_sma", "b1eta_smb", "b1eta_smc", "b1eta_smd", "b1eta_sea", "b1eta_seb", "b1eta_sec", "b1eta_sed"};
    std::string b1e[8]={"b1e_sma", "b1e_smb", "b1e_smc", "b1e_smd", "b1e_sea", "b1e_seb", "b1e_sec", "b1e_sed"};
    std::string bm[8]={"2bm_sma", "2bm_smb", "2bm_smc", "2bm_smd", "2bm_sea", "2bm_seb", "2bm_sec", "2bm_sed"};
    std::string bdr01[8]={"bdr01_sma", "bdr01_smb", "bdr01_smc", "bdr01_smd", "bdr01_sea", "bdr01_seb", "bdr01_sec", "bdr01_sed"};
    std::string pt[8]={"mupt_sma", "mupt_smb", "mupt_smc", "mupt_smd", "ept_sea", "ept_seb", "ept_sec", "ept_sed"};
    std::string eta[8]={"mueta_sma", "mueta_smb", "mueta_smc", "mueta_smd", "eeta_sea", "eeta_seb", "eeta_sec", "eeta_sed"};
    std::string phi[8]={"muphi_sma", "muphi_smb", "muphi_smc", "muphi_smd", "ephi_sea", "ephi_seb", "ephi_sec", "ephi_sed"};

    for ( int i=0; i<4; i++){
	//single muon
	hnvtx_sm[i]        = new TH1F(nvtx[i].c_str(),"primary vertices", 50, 0, 50);
	hmetpt_sm[i]       = new TH1F(met[i].c_str(),"met pt", 30, 0, 200);
	htau0pt_sm[i]      = new TH1F(tau0pt[i].c_str(),"leading taus pt", 30, 0, 400);
	htau0eta_sm[i]     = new TH1F(tau0eta[i].c_str(),"leading taus eta", 30, -3, 3);
	htau1pt_sm[i]      = new TH1F(tau1pt[i].c_str(),"trealing taus pt", 40, 0, 400);
	htau1eta_sm[i]     = new TH1F(tau1eta[i].c_str(),"trealing taus eta", 40, -3, 3);
	hdr2tau_sm[i]      = new TH1F(dr2tau[i].c_str(),"dr(2tau)", 30, 0, 10);
	hdrtaumu_sm[i]     = new TH1F(drtaumu[i].c_str(),"dr(tau,muon)", 30, 0, 10);
	hjmult_sm[i]       = new TH1F(jmult[i].c_str(),"jet multiplicity", 10, 0, 10);
	hzeta_sm[i]        = new TH1F(zeta[i].c_str(),"all jets disc", 50, 0, 1);
	hjpt_sm[i]         = new TH1F(jpt[i].c_str(),"jets pt", 50, 20, 400);
	hjeta_sm[i]        = new TH1F(jeta[i].c_str(),"jets eta", 40, -3, 3);
	hj0pt_sm[i]        = new TH1F(j0pt[i].c_str(),"leading jet pt", 50, 20, 400);
	hj0eta_sm[i]       = new TH1F(j0eta[i].c_str(),"leading jet eta", 50, -5, 5);
	hj0e_sm[i]         = new TH1F(j0e[i].c_str(),"leading jet energy", 50, 0, 500);
	hj1pt_sm[i]        = new TH1F(j1pt[i].c_str(),"trealing jet pt", 50, 20, 400);
	hj1eta_sm[i]       = new TH1F(j1eta[i].c_str(),"trealing jet eta", 50, -5, 5);
	hj1e_sm[i]         = new TH1F(j1e[i].c_str(),"trealing jet energy", 50, 0, 500);
	h2jm_sm[i]         = new TH1F(jm[i].c_str(),"2 jet mass", 50, 0, 600);
	hjdr01_sm[i]       = new TH1F(jdr01[i].c_str(),"dr(j0j1)", 50, 0, 10);
	hbmult_sm[i]       = new TH1F(bmult[i].c_str(),"b multiplicity", 10, 0, 10);
	hbpt_sm[i]         = new TH1F(bpt[i].c_str(),"b pt", 100, 20, 400);
	hb0pt_sm[i]        = new TH1F(b0pt[i].c_str(),"leading b pt", 100, 20, 400);
	hb0eta_sm[i]       = new TH1F(b0eta[i].c_str(),"leading b eta", 100, -5, 5);
	hb0e_sm[i]         = new TH1F(b0e[i].c_str(),"leading b energy", 100, 0, 500);
	hb1pt_sm[i]        = new TH1F(b1pt[i].c_str(),"trealing b pt", 100, 20, 400);
	hb1eta_sm[i]       = new TH1F(b1eta[i].c_str(),"trealing b eta", 50, -5, 5);
	hb1e_sm[i]         = new TH1F(b1e[i].c_str(),"trealing b energy", 50, 0, 500);
	h2bm_sm[i]         = new TH1F(bm[i].c_str(),"2 b mass", 50, 0, 450);
	hbdr01_sm[i]       = new TH1F(bdr01[i].c_str(),"dr(b0b1)", 50, 0, 10);
	hmupt_sm[i]        = new TH1F(pt[i].c_str(),"muon pt", 30, 0, 300);
	hmueta_sm[i]       = new TH1F(eta[i].c_str(),"muon eta", 30, -2.5, 2.5);
	hmuphi_sm[i]       = new TH1F(phi[i].c_str(),"muon phi", 30, -3.14, 3.14);
	//single el
	hnvtx_se[i]        = new TH1F(nvtx[i+4].c_str(),"primary vertices", 50, 0, 50);
	hmetpt_se[i]       = new TH1F(met[i+4].c_str(),"met pt", 30, 0, 200);
	htau0pt_se[i]      = new TH1F(tau0pt[i+4].c_str(),"leading taus pt", 30, 0, 400);
	htau0eta_se[i]     = new TH1F(tau0eta[i+4].c_str(),"leading taus eta", 30, -3, 3);
	htau1pt_se[i]      = new TH1F(tau1pt[i+4].c_str(),"trealing taus pt", 40, 0, 400);
	htau1eta_se[i]     = new TH1F(tau1eta[i+4].c_str(),"trealing taus eta", 40, -3, 3);
	hdr2tau_se[i]      = new TH1F(dr2tau[i+4].c_str(),"dr(2tau)", 30, 0, 10);
	hdrtaumu_se[i]     = new TH1F(drtaumu[i+4].c_str(),"dr(tau,muon)", 30, 0, 10);
	hjmult_se[i]       = new TH1F(jmult[i+4].c_str(),"jet multiplicity", 10, 0, 10);
	hzeta_se[i]        = new TH1F(zeta[i+4].c_str(),"all jets disc", 50, 0, 1);
	hjpt_se[i]         = new TH1F(jpt[i+4].c_str(),"jets pt", 50, 20, 400);
	hjeta_se[i]        = new TH1F(jeta[i+4].c_str(),"jets eta", 40, -3, 3);
	hj0pt_se[i]        = new TH1F(j0pt[i+4].c_str(),"leading jet pt", 50, 20, 400);
	hj0eta_se[i]       = new TH1F(j0eta[i+4].c_str(),"leading jet eta", 50, -5, 5);
	hj0e_se[i]         = new TH1F(j0e[i+4].c_str(),"leading jet energy", 50, 0, 500);
	hj1pt_se[i]        = new TH1F(j1pt[i+4].c_str(),"trealing jet pt", 50, 20, 400);
	hj1eta_se[i]       = new TH1F(j1eta[i+4].c_str(),"trealing jet eta", 50, -5, 5);
	hj1e_se[i]         = new TH1F(j1e[i+4].c_str(),"trealing jet energy", 50, 0, 500);
	h2jm_se[i]         = new TH1F(jm[i+4].c_str(),"2 jet mass", 50, 0, 600);
	hjdr01_se[i]       = new TH1F(jdr01[i+4].c_str(),"dr(j0j1)", 50, 0, 10);
	hbmult_se[i]       = new TH1F(bmult[i+4].c_str(),"b multiplicity", 10, 0, 10);
	hbpt_se[i]         = new TH1F(bpt[i+4].c_str(),"b pt", 100, 20, 400);
	hb0pt_se[i]        = new TH1F(b0pt[i+4].c_str(),"leading b pt", 100, 20, 400);
	hb0eta_se[i]       = new TH1F(b0eta[i+4].c_str(),"leading b eta", 100, -5, 5);
	hb0e_se[i]         = new TH1F(b0e[i+4].c_str(),"leading b energy", 100, 0, 500);
	hb1pt_se[i]        = new TH1F(b1pt[i+4].c_str(),"trealing b pt", 100, 20, 400);
	hb1eta_se[i]       = new TH1F(b1eta[i+4].c_str(),"trealing b eta", 50, -5, 5);
	hb1e_se[i]         = new TH1F(b1e[i+4].c_str(),"trealing b energy", 50, 0, 500);
	h2bm_se[i]         = new TH1F(bm[i+4].c_str(),"2 b mass", 50, 0, 450);
	hbdr01_se[i]       = new TH1F(bdr01[i+4].c_str(),"dr(b0b1)", 50, 0, 10);
	hept_se[i]         = new TH1F(pt[i+4].c_str(),"electron pt", 30, 0, 300);
	heeta_se[i]        = new TH1F(eta[i+4].c_str(),"electron eta", 30, -2.5, 2.5);
	hephi_se[i]        = new TH1F(phi[i+4].c_str(),"electron phi", 30, -3.14, 3.14);

    }




    idtmu= 1<<1; idlmu=1<<0; idvmu=1<<2; 
    idte=1<<1; idle=1<<0; idve=1<<2; 
    idlj=1<<0; 
    dlll= (1<<0) + (1<<1) + (1<<4) + (1<<8);  //decay, loose iso, loose e rej, loose mu rej
    dmll= (1<<0) + (1<<2) + (1<<4) + (1<<8);  //decay, loose iso, loose e rej, loose mu rej
    matchTauGen=1<<10;
    matchGen=1<<3;

    nentries=0; nev=0;

    nlmu=0; nte=0; nle=0; nttau=0; nltau=0; 

    nsmutrig=0; nseltrig=0; nmuORe=0; nltmue1=0; nltmue2=0; n1tmu=0; nolmu=0; noel=0; ntmu=0; 
    nmet=0; ntaujet=0; njet=0; n2j=0; n2b=0; ntau=0; n2tau=0; //all channel

    nsumltmue=0; ntrash=0;

    n2mu=0; n2e=0; n2lep=0; n2tmu=0; n2te=0; 

    //
    metpt=0.;       

    nbtmu=0; nblmu=0; nbvmu=0; nbte=0; nble=0; nbve=0; nbttau=0; nbltau=0; nbtau=0; cj=0; cbj=0;

    nvmu=0; nve=0;

    // initialize 1-D reweighting
    LumiWeights_ = edm::LumiReWeighting(
	    "pileupReweighting/pileup_mc.root",
	    "pileupReweighting/pileup_data.root",
	    "pileup",
	    "pileup");

}

AnalyseNtuple::~AnalyseNtuple()
{

    file->Write();
    file->Close();

    //if (!fChain) return;
    //delete fChain->GetCurrentFile(); //crash

    /*
    //print out
    std::cout<< std::endl;

    std::cout<<"======================================================================"<<std::endl;
    std::cout<<std::setw(34)<<"From skim (tight lepton >=1):" <<std::endl;  
    std::cout<<"======================================================================"<<std::endl;
    std::cout<<std::setw(34)<<"selection"       <<std::setw(12)<<"Evt"                 <<std::setw(12)<<"C.Eff.%" <<std::endl;  
    std::cout<<"======================================================================"<<std::endl;
    //std::cout<<std::setw(34)<<"tight muon>=1:"<<std::setw(12)<< std::setprecision(8)<<nmuORe*weight  <<std::setw(12)<< std::setprecision(4) << 100.*nmuORe/nev <<std::endl;
    // std::cout<<"======================================================================"<<std::endl;
    //  std::cout<<std::setw(34)<<"jet>=2 before tau removal:" <<std::setw(12)<< std::setprecision(8)<<ntaujet*weight  <<std::setw(12)<< std::setprecision(4) << 100.*ntaujet/nev<<std::endl; 
    //  std::cout<<"======================================================================"<<std::endl;
    //  std::cout<<std::setw(34)<<"jet>=2 after tau removal:" <<std::setw(12)<< std::setprecision(8)<<njet*weight  <<std::setw(12)<< std::setprecision(4) << 100.*njet/nev<<std::endl; 
    //  std::cout<<"======================================================================"<<std::endl;
    //  std::cout<<std::setw(34)<<"jet>=2 SCV<0.679:" <<std::setw(12)<< std::setprecision(8)<<nj*weight  <<std::setw(12)<< std::setprecision(4) << 100.*nj/nev<<std::endl; 
    //  std::cout<<"======================================================================"<<std::endl;
    std::cout<<std::setw(34)<<"b-jet>=2 SCV>0.679:" <<std::setw(12)<< std::setprecision(8)<<n2b*weight  <<std::setw(13)<< std::setprecision(4) << 100.*n2b/nev <<std::endl;
    std::cout<<"======================================================================"<<std::endl;
    std::cout<<std::setw(34)<<"met pt>40 GeV:"  <<std::setw(12)<< std::setprecision(8)<<nmet*weight  <<std::setw(12)<< std::setprecision(4) << 100.*nmet/nev <<std::endl;  
    std::cout<<"======================================================================"<<std::endl;
    std::cout<<std::setw(34)<<"single muon trigger:"<<std::setw(12)<< std::setprecision(8)<<nsmutrig*weight  <<std::setw(12)<< std::setprecision(4) << 100.*nsmutrig/nev <<std::endl;
    std::cout<<"======================================================================"<<std::endl;
    std::cout<<std::setw(34)<<"tight muon=1:"  <<std::setw(12)<< std::setprecision(8)<<n1tmu*weight  <<std::setw(12)<< std::setprecision(4) << 100.*n1tmu/nev <<std::endl;  
    std::cout<<"======================================================================"<<std::endl;
    std::cout<<std::setw(34)<<"loose muon=0:"  <<std::setw(12)<< std::setprecision(8)<<nolmu*weight  <<std::setw(12)<< std::setprecision(4) << 100.*nolmu/nev <<std::endl;  
    std::cout<<"======================================================================"<<std::endl;
    std::cout<<std::setw(34)<<"electron=0:"  <<std::setw(12)<< std::setprecision(8)<<noel*weight  <<std::setw(12)<< std::setprecision(4) << 100.*noel/nev <<std::endl;  
    std::cout<<"======================================================================"<<std::endl;
    std::cout<<std::setw(34)<<"jet>=2:"  <<std::setw(12)<< std::setprecision(8)<<n2j*weight  <<std::setw(12)<< std::setprecision(4) << 100.*n2j/nev <<std::endl;  
    std::cout<<"======================================================================"<<std::endl;
    std::cout<<std::setw(34)<<"muon (pt>30, |eta|<2.1, iso<0.1):"  <<std::setw(12)<< std::setprecision(8)<<n2j*weight  <<std::setw(12)<< std::setprecision(4) << 100.*n2j/nev <<std::endl;  
    std::cout<<"======================================================================"<<std::endl;
    std::cout<<std::setw(34)<<"tau>=2:"  <<std::setw(12)<< std::setprecision(8)<<n2tau*weight  <<std::setw(12)<< std::setprecision(4) << 100.*n2tau/nev <<std::endl;  
    std::cout<<"======================================================================"<<std::endl;


    //std::cout<<"======================================================================"<<std::endl;
    //  std::cout<<std::setw(34)<<"lepton+jet analysis" <<std::endl;  
    //  std::cout<<"======================================================================"<<std::endl;
    //  std::cout<<std::setw(34)<<"1 tight mu/ele + loose jet >=2:"<<std::setw(12)<< std::setprecision(5)<<(ntmu+nte)*weight  <<std::setw(12)<< std::setprecision(4) << 100.*(ntmu+nte)/nev<<std::endl; 
    //  std::cout<<"======================================================================"<<std::endl<< std::endl;
     */
    /*
       std::cout<<std::setw(34)<<"di-lepton analysis" <<std::endl;  
       std::cout<<"======================================================================"<<std::endl;
       std::cout<<std::setw(34)<<"2 tight muon:"<<std::setw(12)<< std::setprecision(5)<<n2tmu*weight  <<std::setw(12)<< std::setprecision(4) << 100.*n2tmu/nev<<std::endl; 
       std::cout<<"======================================================================"<<std::endl;
       std::cout<<std::setw(34)<<"2 tight elec:"<<std::setw(12)<< std::setprecision(5)<<n2te*weight  <<std::setw(12)<< std::setprecision(4) << 100.*n2te/nev <<std::endl;
       std::cout<<"======================================================================"<<std::endl;
       std::cout<<std::setw(34)<<"2 tight (1 muon, 1 elec):"<<std::setw(12)<< std::setprecision(5)<<n2lep*weight  <<std::setw(12)<< std::setprecision(4) << 100.*n2lep/nev<<std::endl; 
       std::cout<<"======================================================================"<<std::endl;
       std::cout<<std::setw(34)<<"(1 tight, 1 loose) muon:"<<std::setw(12)<< std::setprecision(5)<<n2mu*weight <<std::setw(12)<< std::setprecision(4) << 100.*n2mu/nev<<std::endl; 
       std::cout<<"======================================================================"<<std::endl;
       std::cout<<std::setw(34)<<"(1 tight, 1 loose) ele:"<<std::setw(12)<< std::setprecision(5)<<n2e*weight <<std::setw(12)<< std::setprecision(4) << 100.*n2e/nev<<std::endl; 


       std::cout<<"======================================================================"<<std::endl;
       std::cout<<std::setw(34)<<"1 tight (e/mu) + 1 loose (mu/e):"<<std::setw(12)<< std::setprecision(5)<<nltmue1*weight <<std::setw(12)<< std::setprecision(4) << 100.*nltmue1/nev<<std::endl; 
       std::cout<<"======================================================================"<<std::endl;
       std::cout<<std::setw(34)<<"1 tight (e/mu) + 1 loose (mu/e):"<<std::setw(12)<< std::setprecision(5)<<nltmue2*weight <<std::setw(12)<< std::setprecision(4) << 100.*nltmue2/nev<<std::endl; 
       std::cout<<"======================================================================"<<std::endl;

       std::cout << std::endl;
       std::cout<<std::setw(34)<<"other:"<<std::setw(12)<< std::setprecision(5)<<ntrash*weight  <<std::setw(12)<< std::setprecision(4) << 100.*ntrash/nev<<std::endl; 
       std::cout<<"======================================================================"<<std::endl;
     */

}

Int_t AnalyseNtuple::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t AnalyseNtuple::LoadTree(Long64_t entry)
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

void AnalyseNtuple::Init(TTree *tree)
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
    fChain->SetBranchAddress("trig", &trig, &b_trig);
    fChain->SetBranchAddress("vx", vx, &b_vx);
    fChain->SetBranchAddress("vy", vy, &b_vy);
    fChain->SetBranchAddress("vz", vz, &b_vz);
    Notify();
}

Bool_t AnalyseNtuple::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void AnalyseNtuple::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t AnalyseNtuple::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif // #ifdef AnalyseNtuple_cxx
