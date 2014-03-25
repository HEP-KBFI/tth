// -*- C++ -*-

//
// Package:    CreateNTuple
// Class:      CreateNTuple
// 
/**\class CreateNTuple CreateNTuple.cc TTH/CreateNTuple/src/CreateNTuple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Betty Calpas
//         Created:  Thu Mar 13 15:08:40 EET 2014
// $Id$
//
//


// system include files
#include <memory>
// user include files
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h" // should before Frameworkfwd.h 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/Common/interface/Handle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TTree.h"
#include <TLorentzVector.h>
#include <algorithm> 
#include <string>
//
// class declaration
//
#define MAXPART 1000

class CreateNTuple : public edm::EDAnalyzer {
public:
  explicit CreateNTuple(const edm::ParameterSet&);
  ~CreateNTuple();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // ----------member data ---------------------------
  edm::InputTag genLabel;
  edm::InputTag tauLabel;
  edm::InputTag muonLabel;
  edm::InputTag elecLabel;
  edm::InputTag jetLabel;
  edm::InputTag metLabel;
  edm::InputTag vtxLabel;

  //Evt info
  int run, lumi, ev;
  //Count
  int np, nd, nj;
  //Common 
  double Pt[MAXPART], Eta[MAXPART], Phi[MAXPART], Px[MAXPART], Py[MAXPART], Pz[MAXPART], E[MAXPART], M[MAXPART], Q[MAXPART];
  //Gen
  int pid[MAXPART], status[MAXPART], daughterId[MAXPART], ndaughter[MAXPART];
  //Vertex
  double Vx[MAXPART], Vy[MAXPART], Vz[MAXPART]; //??
  //PileUp
  double npu;


  //Tau
  double tauMinPt, tauMinEta, tauMass;
  int tauId[MAXPART], tau_byDecayModeFinding[MAXPART];
  int tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[MAXPART];
  int tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[MAXPART];
  int tau_byTightCombinedIsolationDeltaBetaCorr3Hits[MAXPART];
  int tau_againstElectronLooseMVA3[MAXPART], tau_againstElectronMediumMVA3[MAXPART];
  int tau_againstElectronTightMVA3[MAXPART], tau_againstElectronVTightMVA3[MAXPART];
  int tau_againstMuonLoose3[MAXPART], tau_againstMuonMedium3[MAXPART],tau_againstMuonTight3[MAXPART];
  //Muon
  double muonMinPt, muonMinEta;
  double sumChHadPt[MAXPART],sumNeutHadEt[MAXPART],sumPhotonEt[MAXPART],sumPUPt[MAXPART];
  double sumChPartPt[MAXPART],muonIsoPflow[MAXPART],muonIsoPflowPUcorr[MAXPART], PFMuonIso[MAXPART];
  bool isLooseMuon[MAXPART], isSoftMuon[MAXPART], isGoodMuon[MAXPART], isTightMuon[MAXPART]; 
  //Electron
  double elecMinPt, elecMinEta;
  bool isVetoEle[MAXPART], isLooseEle[MAXPART], isMediumEle[MAXPART], isTightEle[MAXPART];
  bool passCutBasedPresel[MAXPART], ele_passMVAPresel[MAXPART], passPflowPresel[MAXPART];
  double ele_MVA[MAXPART];
  double DEtaSupClusTrkAtVtx[MAXPART], DPhiSupClusTrkAtVtx[MAXPART], Sigma2Ieta[MAXPART], HadOverEm[MAXPART], ele_PfIso[MAXPART];
  double TrkDxy[MAXPART], TrkDz[MAXPART], ConvRej[MAXPART], MissConRej[MAXPART], PfIso[MAXPART], InvEminusInvPin[MAXPART];
  //Jet
  double jetMinPt, jetMinEta;
  double jetPt[MAXPART], jetEta[MAXPART], jetPhi[MAXPART], jetPx[MAXPART], jetPy[MAXPART], jetPz[MAXPART], jetE[MAXPART], jetM[MAXPART];
  double neutHadEnFrac[MAXPART], neutEmEnFrac[MAXPART], chHadEnFrac[MAXPART],chMultiplicity[MAXPART], chEmEnFrac[MAXPART];
  int nDaughter[MAXPART];
  bool isjetLoose[MAXPART], isjetTight[MAXPART];
  float jetB[MAXPART];
  //MET  
  double metMinPt;
  double metPt, metEta, metPhi, metPx, metPy;
  double metCov00, metCov10, metCov01, metCov11;
  
  TH1D *h;
  TTree *t;  
  
  JetIDSelectionFunctor jetIDLoose, jetIDTight;

};
//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
CreateNTuple::CreateNTuple(const edm::ParameterSet& iConfig):
  //Label 
  genLabel     (iConfig.getParameter<edm::InputTag>("genLabel")),
  tauLabel     (iConfig.getParameter<edm::InputTag>("tauLabel")),
  muonLabel    (iConfig.getParameter<edm::InputTag>("muonLabel")),
  elecLabel    (iConfig.getParameter<edm::InputTag>("elecLabel")),
  jetLabel     (iConfig.getParameter<edm::InputTag>("jetLabel")),
  metLabel     (iConfig.getParameter<edm::InputTag>("metLabel")),
  vtxLabel     (iConfig.getParameter<edm::InputTag>("vtxLabel")),
  //Cut
  tauMinPt     (iConfig.getParameter<double>("tauMinPt")),
  tauMinEta    (iConfig.getParameter<double>("tauMinEta")),
  muonMinPt    (iConfig.getParameter<double>("muonMinPt")),
  muonMinEta   (iConfig.getParameter<double>("muonMinEta")),
  elecMinPt    (iConfig.getParameter<double>("elecMinPt")),
  elecMinEta   (iConfig.getParameter<double>("elecMinEta")),
  jetMinPt     (iConfig.getParameter<double>("jetMinPt")),
  jetMinEta    (iConfig.getParameter<double>("jetMinEta")),
  metMinPt     (iConfig.getParameter<double>("metMinPt")),
  jetIDLoose(JetIDSelectionFunctor::PURE09,JetIDSelectionFunctor::LOOSE),
  jetIDTight(JetIDSelectionFunctor::PURE09,JetIDSelectionFunctor::TIGHT)
  
{
  //now do what ever initialization is needed
}


CreateNTuple::~CreateNTuple()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void
CreateNTuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace reco;
  using namespace edm;
  using namespace pat;
  using namespace std;

  run=0; lumi=0; ev=0; np=0; nd=0; nj=0;

  for (int i=0; i<MAXPART; i++) {
    //Common
    Pt[i]=0; Eta[i]=0; Phi[i]=0; Px[i]=0; Py[i]=0; Pz[i]=0; E[i]=0; M[i]=0; Q[i]=0;
    //Gen
    pid[i]=0; status[i]=0; ndaughter[i]=0; daughterId[i]=0; 
    //Vtx
    Vx[i]=0; Vy[i]=0; Vz[i]=0; 
    //tau
    tauId[i]=1;
 
    tau_byDecayModeFinding[i]=0;
    tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[i]=0; 
    tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[i]=0;
    tau_byTightCombinedIsolationDeltaBetaCorr3Hits[i]=0;
    tau_againstElectronLooseMVA3[i]=0;
    tau_againstElectronMediumMVA3[i]=0;
    tau_againstElectronTightMVA3[i]=0;
    tau_againstElectronVTightMVA3[i]=0;
    tau_againstMuonLoose3[i]=0;
    tau_againstMuonMedium3[i]=0;
    tau_againstMuonTight3[i]=0;

    //muon   
    sumChHadPt[i]=0; 
    sumNeutHadEt[i]=0; 
    sumPhotonEt[i]=0;
    sumPUPt[i]=0;
    sumChPartPt[i]=0; 
    muonIsoPflow[i]=0; 
    muonIsoPflowPUcorr[i]=0;
    PFMuonIso[i]=0;
    isLooseMuon[i]={false};
    isSoftMuon[i]={false};
    isGoodMuon[i]={false};
    isTightMuon[i]={false};
    //elec
    ele_MVA[i]=0;
    ele_passMVAPresel[i]=false;
    ele_PfIso[i]=0;
    isVetoEle[i]=false;    
    isLooseEle[i]=false;    
    isMediumEle[i]=false;    
    isTightEle[i]=false;    
    DEtaSupClusTrkAtVtx[i]=0;
    DPhiSupClusTrkAtVtx[i]=0;
    Sigma2Ieta[i]=0;
    HadOverEm[i]=0;
    TrkDxy[i]=0;
    TrkDz[i]=0;
    ConvRej[i]=0;
    MissConRej[i]=0;
    PfIso[i]=0;
    InvEminusInvPin[i]=0;
    
    jetPt[i]=0;
    jetEta[i]=0;
    jetPhi[i]=0;
    jetPx[i]=0;
    jetPy[i]=0;
    jetPz[i]=0;
    jetE[i]=0;
    jetM[i]=0;
    isjetLoose[i]=false;
    isjetTight[i]=false;
    neutHadEnFrac[i]=0;
    neutEmEnFrac[i]=0;
    nDaughter[i]=0;
    chHadEnFrac[i]=0;
    chMultiplicity[i]=0;
    chEmEnFrac[i]=0;
    jetB[i]=0;
  }

//Evt
  run = iEvent.id().run();
  lumi = iEvent.id().luminosityBlock();
  ev = iEvent.id().event();

//Gen
  Handle<vector<reco::GenParticle> > genReco; 	
  iEvent.getByLabel(genLabel,genReco);
  for (vector<reco::GenParticle>::const_iterator it = genReco->begin(); it != genReco->end(); it++) {
     status[np]     = it->status();  
     ndaughter[np]  = it->numberOfDaughters();
     //
     pid[np] = it->pdgId();
     Pt[np]  = it->pt();
     Eta[np] = it->eta();
     Phi[np] = it->phi();
     M[np]   = it->mass();
     E[np]   = it->energy();
     Px[np]  = it->px();
     Py[np]  = it->py();
     Pz[np]  = it->pz();
     Vx[np]  = it->vx();
     Vy[np]  = it->vy();
     Vz[np]  = it->vz();
     Q[np]   = it->charge();
     np++ ;
//mother
     for(unsigned int im = 0; im < it->numberOfMothers(); ++im) {
       pid[np] = it->motherRef(im)->pdgId();
       Pt[np]  = it->motherRef(im)->pt();
       Eta[np] = it->motherRef(im)->eta();
       Phi[np] = it->motherRef(im)->phi();
       M[np]   = it->motherRef(im)->mass();
       E[np]   = it->motherRef(im)->energy(); 
       Px[np]  = it->motherRef(im)->px();
       Py[np]  = it->motherRef(im)->py();
       Pz[np]  = it->motherRef(im)->pz();
       Vx[np]  = it->motherRef(im)->vx();
       Vy[np]  = it->motherRef(im)->vy();
       Vz[np]  = it->motherRef(im)->vz();
       Q[np]   = it->motherRef(im)->charge();
       np++; 
}
//daughter
     for(unsigned int id = 0; id < it->numberOfDaughters(); ++id) {
       pid[np] = it->daughterRef(id)->pdgId();
       Pt[np]  = it->daughterRef(id)->pt();
       Eta[np] = it->daughterRef(id)->eta();
       Phi[np] = it->daughterRef(id)->phi();
       M[np]   = it->daughterRef(id)->mass();
       E[np]   = it->daughterRef(id)->energy(); 
       Px[np]  = it->daughterRef(id)->px();
       Py[np]  = it->daughterRef(id)->py();
       Pz[np]  = it->daughterRef(id)->pz();
       Vx[np]  = it->daughterRef(id)->vx();
       Vy[np]  = it->daughterRef(id)->vy();
       Vz[np]  = it->daughterRef(id)->vz();
       Q[np]   = it->daughterRef(id)->charge();
       np++; 
     }
   }

  //tau
  Handle<pat::TauCollection> tauPat;
  iEvent.getByLabel(tauLabel,tauPat);

  cout<< "new evt: "<< endl <<endl;

  for (pat::TauCollection::const_iterator it = tauPat->begin(); it != tauPat->end(); it++) {

     pid[np] = it->pdgId();
     Pt[np]  = it->pt();
     Eta[np] = it->eta();
     Phi[np] = it->phi();
     M[np]   = it->mass();
     E[np]   = it->energy();
     Px[np]  = it->px();
     Py[np]  = it->py();
     Pz[np]  = it->pz();
     Vx[np]  = it->vx();
     Vy[np]  = it->vy();
     Vz[np]  = it->vz();
     Q[np]   = it->charge();

    tau_byDecayModeFinding[np]                          =  it->tauID("decayModeFinding"); 
    tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[np]  =  it->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
    tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[np] =  it->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
    tau_byTightCombinedIsolationDeltaBetaCorr3Hits[np]  =  it->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
    tau_againstElectronLooseMVA3[np]                    =  it->tauID("againstElectronLooseMVA3");
    tau_againstElectronMediumMVA3[np]                   =  it->tauID("againstElectronMediumMVA3");
    tau_againstElectronTightMVA3[np]                    =  it->tauID("againstElectronTightMVA3");
    tau_againstElectronVTightMVA3[np]                   =  it->tauID("againstElectronVTightMVA3");
    tau_againstMuonLoose3[np]                           =  it->tauID("againstMuonLoose3");
    //tau_againstMuonMedium3[np]                        =  it->tauID("againstMuonMedium3"); //!!missing!!
    tau_againstMuonTight3[np]                           =  it->tauID("againstMuonTight3");

cout<< "New tau:  tauid: "<<tauId[np] << endl;

    if(tau_byDecayModeFinding[np]==1) tauId[np] += 1<<0;
cout<< "decay dis: "<< tau_byDecayModeFinding[np] <<"  tauid decay: "<<tauId[np] << endl;

    if(tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[np]==1) tauId[np]+= 1<<1;
cout<< "Loose dis: "<< tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[np] <<"   tauid loose: "<<tauId[np] << endl;

    if(tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[np]==1) tauId[np]+= 1<<2;
cout<< "Medium dis: "<< tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[np] <<"   tauid medium: "<<tauId[np] << endl;

   if(tau_byTightCombinedIsolationDeltaBetaCorr3Hits[np]==1) tauId[np]+= 1<<3;
cout<< "Tight dis: "<< tau_byTightCombinedIsolationDeltaBetaCorr3Hits[np] <<"   tauid Tight: "<<tauId[np] << endl;

    if(tau_againstElectronLooseMVA3[np]==1) tauId[np]+= 1<<4;
cout<<"against Loose ele: "<<tau_againstElectronLooseMVA3[np]<< "  tauId: "<<tauId[np] << endl;

    if(tau_againstElectronMediumMVA3[np]==1) tauId[np]+= 1<<5;
cout<<"against Medium ele: "<<tau_againstElectronMediumMVA3[np]<< "  tauId: "<<tauId[np] << endl;

    if(tau_againstElectronTightMVA3[np]==1) tauId[np]+= 1<<6;
cout<<"against Tight ele: "<<tau_againstElectronTightMVA3[np]<< "  tauId: "<<tauId[np] << endl;

    if(tau_againstElectronVTightMVA3[np]==1) tauId[np]+= 1<<7;
cout<<"against VTight ele: "<<tau_againstElectronVTightMVA3[np]<< "  tauId: "<<tauId[np] << endl;

    if(tau_againstMuonLoose3[np]==1) tauId[np]+= 1<<8;
cout<<"against Loose muon: "<<tau_againstMuonLoose3[np]<<"  tauId: "<<tauId[np] << endl;

    if(tau_againstMuonMedium3[np]==1) tauId[np]+= 1<<9;
cout<<"against Medium muon: "<<tau_againstMuonMedium3[np]<<"  tauId: "<<tauId[np] << endl;

    if(tau_againstMuonTight3[np]==1) tauId[np]+= 1<<10;
cout<<"against Tight muon: "<<tau_againstMuonTight3[np]<<"  tauId: "<<tauId[np] << endl;
    np++;
  }//tau
  

  //muon
  Handle<pat::MuonCollection> muonPat;
  iEvent.getByLabel(muonLabel,muonPat);
  Handle<reco::VertexCollection> vtxRec;
  iEvent.getByLabel(vtxLabel,vtxRec);

  for (pat::MuonCollection::const_iterator it = muonPat->begin(); it != muonPat->end(); it++) {
  isLooseMuon[np]    =  it->isLooseMuon();
  isSoftMuon[np]     =  it->isSoftMuon(*vtxRec->begin());
  isGoodMuon[np]     =  it->isGood("GlobalMuonPromptTight"); 
  isTightMuon[np]    =  it->isTightMuon(*vtxRec->begin());
  sumChHadPt[np]     =  it->pfIsolationR04().sumChargedHadronPt;
  sumNeutHadEt[np]   =  it->pfIsolationR04().sumNeutralHadronEt;
  sumPhotonEt[np]    =  it->pfIsolationR04().sumPhotonEt;
  sumPUPt[np]        =  it->pfIsolationR04().sumPUPt;
  sumChPartPt[np]    =  it->pfIsolationR04().sumChargedParticlePt;
  muonIsoPflow[np]        =  (sumChHadPt[np]+sumNeutHadEt[np]+sumPhotonEt[np])/sumChPartPt[np];
  muonIsoPflowPUcorr[np]  =  (sumChHadPt[np]+std::max(0., (sumNeutHadEt[np]+sumPhotonEt[np]-0.5*sumPUPt[np])))/sumChPartPt[np];

     pid[np] = it->pdgId();
     Pt[np]  = it->pt();
     Eta[np] = it->eta();
     Phi[np] = it->phi();
     M[np]   = it->mass();
     E[np]   = it->energy();
     Px[np]  = it->px();
     Py[np]  = it->py();
     Pz[np]  = it->pz();
     Vx[np]  = it->vx();
     Vy[np]  = it->vy();
     Vz[np]  = it->vz();
     Q[np]   = it->charge();

    np++;
  
  }//muon


  //electron
  //Handle<vector<pat::Electron> > elecPat;

  Handle<ElectronCollection> elecPat;
  //Handle<View<pat::Electron>> elecPat;  
  //edm::Handle<edm::View<pat::Electron> > elecPat;

cout<< "toto1"<< elecLabel << "," << elecPat << endl;
  iEvent.getByLabel(elecLabel,elecPat); //bug
cout<< "toto2"<<endl;
  np=0;
 
  for (pat::ElectronCollection::const_iterator it = elecPat->begin(); it != elecPat->end(); it++) {

     pid[np]        = it->pdgId();
     Pt[np]  = it->pt();
     Eta[np] = it->eta();
     Phi[np] = it->phi();
     M[np]   = it->mass();
     E[np]   = it->energy();
     Px[np]  = it->px();
     Py[np]  = it->py();
     Pz[np]  = it->pz();
     Vx[np]  = it->vx();
     Vy[np]  = it->vy();
     Vz[np]  = it->vz();
     Q[np]   = it->charge();
     ele_passMVAPresel[np]       = it->passingMvaPreselection();
     ele_MVA[np] = it->mva();
    //ele_PfIso[np] = it->pfIsolationVariables();
    np++;
  }

  Handle<pat::JetCollection> jetPat;  
  iEvent.getByLabel(jetLabel,jetPat);

  for (pat::JetCollection::const_iterator it = jetPat->begin(); it != jetPat->end(); it++) {
    neutHadEnFrac[nj]  = it->neutralHadronEnergyFraction();
    neutEmEnFrac[nj]   = it->neutralEmEnergyFraction();
    nDaughter[nj]      = it->numberOfDaughters();
    chHadEnFrac[nj]    = it->chargedHadronEnergyFraction();
    chMultiplicity[nj] = it->chargedMultiplicity();
    chEmEnFrac[nj]     = it->chargedEmEnergyFraction();
    jetPt[nj]          =  it->pt();
    jetEta[nj]         =  it->eta();
    jetPhi[nj]         =  it->phi();
    jetPx[nj]          = it->px();
    jetPy[nj]          = it->py();
    jetPz[nj]          = it->pz();
    jetE[nj]           = it->energy();
    jetM[nj]           = it->mass();
    //btag
    jetB[nj]           = it->bDiscriminator("CombinedSecondaryVertex");
    nj++;
  }
  
  //MET
  Handle<pat::METCollection> met;
  iEvent.getByLabel(metLabel,met);
  metPt=met->begin()->pt();
  metEta=met->begin()->eta();
  metPhi=met->begin()->phi();
  metPx=met->begin()->px();
  metPy=met->begin()->py();
  metCov00= (met->front()).getSignificanceMatrix()(0,0);
  metCov10= (met->front()).getSignificanceMatrix()(1,0);
  metCov01= (met->front()).getSignificanceMatrix()(0,1);
  metCov11= (met->front()).getSignificanceMatrix()(1,1);

  //cout<< "np tot: " << np << endl;

 h->Fill(np);
 t->Fill(); 
}


// ------------ method called once each job just before starting event loop  ------------
void 
CreateNTuple::beginJob()
{
  edm::Service<TFileService> fs;
  h = fs->make<TH1D>("nlep", "Number of leptons", MAXPART, 0, MAXPART);
  t = new TTree("cands","Flat-tuple");
  //Evt info
  t->Branch("run",&run,"run/I");
  t->Branch("lumi",&lumi,"lumi/I");
  t->Branch("ev",&ev,"ev/I");
  t->Branch("np",&np,"np/I");
/*
  //electron
  //t->Branch("passCutBasedPresel",passCutBasedPresel,"passCutBasedPresel[np]/O");
  t->Branch("ele_passMVAPresel",ele_passMVAPresel,"ele_passMVAPresel[np]/O");
  t->Branch("ele_MVA",ele_MVA,"ele_MVA[np]/D");
  //t->Branch("ele_PfIso",ele_PfIso,"ele_PfIso[np]/D");
  //t->Branch("passPflowPresel",passPflowPresel,"passPflowPresel[np]/O");
  //t->Branch("isVetoEle",isVetoEle,"isVetoEle[np]/O");
  //t->Branch("isLooseEle",isLooseEle,"isLooseEle[np]/O");
  //t->Branch("isMediumEle",isMediumEle,"isMediumEle[np]/O");
  //t->Branch("isTightEle",isMediumEle,"isTightEle[np]/O");
*/
  //jet
  t->Branch("njet",&nj,"nj/I");
  t->Branch("isjetLoose",isjetLoose,"isjetLoose[nj]/O");
  t->Branch("isjetTight",isjetTight,"isjetTight[nj]/O");
  t->Branch("neutHadEnFrac",neutHadEnFrac,"neutHadEnFrac[nj]/D");
  t->Branch("neutEmEnFrac",neutEmEnFrac,"neutEmEnFrac[nj]/D");
  t->Branch("nDaughter",nDaughter,"nDaughter[nj]/I");
  t->Branch("chHadEnFrac",chHadEnFrac,"chHadEnFrac[nj]/D");
  t->Branch("chMultiplicity",chMultiplicity,"chMultiplicity[nj]/D");
  t->Branch("chEmEnFrac",chEmEnFrac,"chEmEnFrac[nj]/D");
  t->Branch("jetPt",jetPt,"jetPt[nj]/D");
  t->Branch("jetEta",jetEta,"jetEta[nj]/D");
  t->Branch("jetPhi",jetPhi,"jetPhi[nj]/D");
  t->Branch("jetPx",jetPx,"jetPx[nj]/D");
  t->Branch("jetPy",jetPy,"jetPy[nj]/D");
  t->Branch("jetPz",jetPz,"jetPz[nj]/D");
  t->Branch("jetE",jetE,"jetE[nj]/D");
  t->Branch("jetM",jetM,"jetM[nj]/D");
  t->Branch("jetB",jetB,"jetB[nj]/D");
  //met
  t->Branch("metPt",&metPt,"met_pt/D");
  t->Branch("metEta",&metEta,"met_eta/D");
  t->Branch("metPhi",&metPhi,"met_phi/D");
  t->Branch("metPx",&metPx,"met_px/D");
  t->Branch("metPy",&metPy,"met_py/D");
  t->Branch("metCov00",&metCov00,"metCov/D");
  t->Branch("metCov10",&metCov10,"metCov/D");
  t->Branch("metCov01",&metCov01,"metCov/D");
  t->Branch("metCov11",&metCov11,"metCov/D");
  //muon
  t->Branch("isLooseMuon",isLooseMuon,"isLooseMuon[np]/O");
  t->Branch("isSoftMuon",isSoftMuon,"isSoftMuon[np]/O");
  t->Branch("isGoodMuon",isGoodMuon,"isGoodMuon[np]/O");
  t->Branch("isTightMuon",isTightMuon,"isTightMuon[np]/O");
  //t->Branch("sumChHadPt",sumChHadPt,"sumChHadPt[np]/O");
  //t->Branch("sumNeutHadEt",sumNeutHadEt,"sumNeutHadEt[np]/O");
  //t->Branch("sumPhotonEt",sumPhotonEt,"sumPhotonEt[np]/O");
  //t->Branch("sumPUPt",sumPUPt,"sumPUPt[np]/O");
  //t->Branch("sumChPartPt",sumChPartPt,"sumChPartPt[np]/O");
  //t->Branch("muonIsoPflow",muonIsoPflow,"muonIsoPflow[np]/D");
  t->Branch("muonIsoPflowPUcorr",muonIsoPflowPUcorr,"muonIsoPflowPUcorr[np]/D");
  //tau
  t->Branch("pdgId",pid,"pid[np]/I");
  t->Branch("tauId",tauId,"tauId[np]/I");
  t->Branch("Pt",Pt,"Pt[np]/D");
  t->Branch("Eta",Eta,"Eta[np]/D");
  t->Branch("Phi",Phi,"Phi[np]/D");
  t->Branch("Px",Px,"Px[np]/D");
  t->Branch("Py",Py,"Py[np]/D");
  t->Branch("Pz",Pz,"Pz[np]/D");
  t->Branch("E",E,"E[np]/D");
  t->Branch("M",M,"M[np]/D");
  t->Branch("Q",Q,"Q[np]/D");
  t->Branch("Vx",Vx,"Vx[np]/D");
  t->Branch("Vy",Vy,"Vy[np]/D");
  t->Branch("Vz",Vz,"Vz[np]/D");
/*
  t->Branch("tau_byDecayModeFinding",tau_byDecayModeFinding,"tau_byDecayModeFinding[np]/I");
  t->Branch("tau_byLooseCombinedIsolationDeltaBetaCorr3Hits",tau_byLooseCombinedIsolationDeltaBetaCorr3Hits,"tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[np]/I");
  t->Branch("tau_byMediumCombinedIsolationDeltaBetaCorr3Hits",tau_byMediumCombinedIsolationDeltaBetaCorr3Hits,"tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[np]/I");
  t->Branch("tau_byTightCombinedIsolationDeltaBetaCorr3Hits",tau_byTightCombinedIsolationDeltaBetaCorr3Hits,"tau_byTightCombinedIsolationDeltaBetaCorr3Hits[np]/I");
  t->Branch("tau_againstElectronLooseMVA3",tau_againstElectronLooseMVA3,"tau_againstElectronLooseMVA3[np]/I");
  t->Branch("tau_againstElectronMediumMVA3",tau_againstElectronMediumMVA3,"tau_againstElectronMediumMVA3[np]/I");
  t->Branch("tau_againstElectronTightMVA3",tau_againstElectronTightMVA3,"tau_againstElectronTightMVA3[np]/I");
  t->Branch("tau_againstElectronVTightMVA3",tau_againstElectronVTightMVA3,"tau_againstElectronVTightMVA3[np]/I");
  t->Branch("tau_againstMuonLoose3",tau_againstMuonLoose3,"tau_againstMuonLoose3[np]/I");
  //t->Branch("tau_againstMuonMedium3",tau_againstMuonMedium3,"tau_againstMuonMedium3[np]/I");
  t->Branch("tau_againstMuonTight3",tau_againstMuonTight3,"tau_againstMuonTight3[np]/I"); 
*/

}

// ------------ method called once each job just after ending the event loop  ------------
void 
CreateNTuple::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
CreateNTuple::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CreateNTuple::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CreateNTuple::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CreateNTuple::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CreateNTuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CreateNTuple);



  
  //Handle<reco::TrackCollection> tracks;
  //iEvent.getByLabel("generalTracks", tracks); 
  //LogInfo("Demo") << "number of tracks "<<tracks->size();
  //hnTrack->Fill(tracks->size());
  //if( minTracks_ <= tracks->size() ) {
  //   LogInfo("Demo") << "number of tracks "<<tracks->size();
  //}


/*
    passCutBasedPresel[nelec]  = it->passingCutBasedPreselection();
    passPflowPresel[nelec]     = it->passingPflowPreselection();
    
    //Do not work: indefined ref to EgammaCutBasedEleId...
    //isVetoEle[nelec]  = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::VETO, *it, conversions_h, beamSpot, vtx_h, it->chargedHadronIso(), it->photonIso(), it->neutralHadronIso(), rhoIso);
    //isLooseEle[nelec]  = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::LOOSE, *it, conversions_h, beamSpot, vtx_h, it->chargedHadronIso(), it->photonIso(), it->neutralHadronIso(), rhoIso);
    //isMediumEle[nelec]  = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::MEDIUM, *it, conversions_h, beamSpot, vtx_h, it->chargedHadronIso(), it->photonIso(), it->neutralHadronIso(), rhoIso);
    //isTightEle[nelec]  = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, *it, conversions_h, beamSpot, vtx_h, it->chargedHadronIso(), it->photonIso(), it->neutralHadronIso(), rhoIso);
    
    DEtaSupClusTrkAtVtx[nelec] = it->deltaEtaSuperClusterTrackAtVtx();
    DPhiSupClusTrkAtVtx[nelec] = it->deltaPhiSuperClusterTrackAtVtx();
    Sigma2Ieta[nelec] = it->sigmaIetaIeta();
    HadOverEm[nelec] = it->hadronicOverEm();
    TrkDxy[nelec] = vtx_h->begin()->dxy();
    TrkDz[nelec] = vtx_h->begin()->dz();
    ConvRej[nelec] = it->conversionRejectionVariables(); 
    // number of missing hits conversion rejection 
    MissConRej[nelec] = it->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    InvEminusInvPin[nelec]= (1/it->ecalEnergy() - 1/it->trackMomentumAtVtx().p());
    */

  /*
  //jet
  Handle<reco::CaloJetCollection> jetcoll;
  iEvent.getByLabel("ak5CaloJets",jetcoll);
  Handle<reco::JetIDValueMap> jetIDmap;
  iEvent.getByLabel("ak5JetID",jetIDmap);
  int njet = 0;

  for(unsigned int ijet=0;ijet<jetcoll->size();++ijet) {
  const reco::CaloJetRef jet(jetcoll,ijet);
  if(jetIDmap->contains(jet.id())) {
  const reco::JetID & jetid = (*jetIDmap)[jet];
  pat::strbitset retLoose = jetIDLoose.getBitTemplate();
  pat::strbitset retTight = jetIDTight.getBitTemplate();
  retLoose.set(false);
  retTight.set(false);
  isjetLoose[njet] = jetIDLoose((*jetcoll)[ijet],jetid,retLoose);
  isjetTight[njet] = jetIDTight((*jetcoll)[ijet],jetid,retTight);
  jet_px[njet] = jet->px();
  jet_py[njet] = jet->py();
  jet_pz[njet] = jet->pz();
  jet_e[njet]  = jet->energy();
  jet_m[njet]  = jet->mass();
  njet++;
  }
  }
  */ 

  //bug
/*
  // conversions
  edm::Handle<reco::ConversionCollection> conversions_h;
  iEvent.getByLabel("allConversions", conversions_h);
  // beam spot
  edm::Handle<reco::BeamSpot> beamspot_h;
  iEvent.getByLabel("offlineBeamSpot", beamspot_h);
  const reco::BeamSpot &beamSpot = *(beamspot_h.product());
  // vertices
  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByLabel("offlinePrimaryVertices", vtx_h);
  // rho for isolation
  edm::Handle<double> rhoIso_h;
  iEvent.getByLabel("kt6PFJetsForIsolation", "rho", rhoIso_h);
  double rhoIso = *(rhoIso_h.product());
*/


  //SVFit bug
  /*
  // define MET
  Vector vMET(met_px, met_py, 0.);
  // define MET covariance
  TMatrixD covMET(2, 2); // PFMET significance matrix
  covMET[0][0]= (met->front() ).getSignificanceMatrix()(0,0);
  covMET[1][0]= (met->front() ).getSignificanceMatrix()(1,0);
  covMET[0][1]= (met->front() ).getSignificanceMatrix()(0,1);
  covMET[1][1]= (met->front() ).getSignificanceMatrix()(1,1);
  
  // define lepton four vectors
  NSVfitStandalone::LorentzVector l1(tau_px[0], tau_py[0], tau_pz[0], TMath::Sqrt(tau_m[0]*tau_m[0]+ tau_px[0]*tau_px[0] + tau_py[0]*tau_py[0]+ tau_pz[0]*tau_pz[0]));
  NSVfitStandalone::LorentzVector l2(tau_px[1], tau_py[1], tau_pz[1], TMath::Sqrt(tau_m[1]*tau_m[1]+ tau_px[1]*tau_px[1] + tau_py[1]*tau_py[1]+ tau_pz[1]*tau_pz[1]));

  std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
  //measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, leg1->p4()));
  //measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, leg2->p4()));
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, l1));
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, l2));
  //NSVfitStandaloneAlgorithm algo(measuredTauLeptons, met->momentum(), covMET, 0);
  
  NSVfitStandaloneAlgorithm algo(measuredTauLeptons, vMET, covMET,0);
  */
  /*  //SVFit check
  // define MET
  NSVfitStandalone::Vector MET(11.7491, -51.9172, 0.);
  // define MET covariance
  TMatrixD covMET(2, 2);
  covMET[0][0] = 787.352;
  covMET[1][0] = -178.63;
  covMET[0][1] = -178.63;
  covMET[1][1] = 179.545;
  // define lepton four vectors
  LorentzVector l3(28.9132, -17.3888, 36.6411, 49.8088); //lepton
  LorentzVector l4(-24.19, 8.77449, 16.9413, 30.8086); //tau

  std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kHadDecay, l4));
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, l3));
  // define algorithm (set the debug level to 3 for testing)
  //NSVfitStandaloneAlgorithm algo(measuredTauLeptons, MET, covMET, 0);
  //algo.addLogM(false);
  
  // integration by markov chain MC
  algo.integrateMarkovChain();

  
  double mass = algo.getMass(); // mass uncertainty not implemented yet
  if (algo.isValidSolution()) {
  std::cout << "found mass = " << mass << std::endl;
  } else {
  std::cout << "sorry -- status of NLL is not valid [" << algo.isValidSolution() << "]" << std::endl;
  }
    //return;
    */
/*  for (unsigned int i=0; i<MAXPART; i++) {
    //tau
    tauId[i]=1;
    tau_pt[i]=0;
    tau_eta[i]=0;
    tau_phi[i]=0;
    tau_px[i]=0;
    tau_py[i]=0;
    tau_pz[i]=0;
    tau_e[i]=0;
    tau_m[i]=0;
    tau_byDecayModeFinding[i]=0;
    tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[i]=0; 
    tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[i]=0;
    tau_byTightCombinedIsolationDeltaBetaCorr3Hits[i]=0;
    tau_againstElectronLooseMVA3[i]=0;
    tau_againstElectronMediumMVA3[i]=0;
    tau_againstElectronTightMVA3[i]=0;
    tau_againstElectronVTightMVA3[i]=0;
    tau_againstMuonLoose3[i]=0;
    tau_againstMuonMedium3[i]=0;
    tau_againstMuonTight3[i]=0;
    //muon   
    sumChHadPt[i]=0; 
    sumNeutHadEt[i]=0; 
    sumPhotonEt[i]=0;
    sumPUPt[i]=0;
    sumChPartPt[i]=0; 
    muonIsoPflow[i]=0; 
    muonIsoPflowPUcorr[i]=0;
    PFMuonIso[i]=0;
    isLooseMuon[i]={false};
    isSoftMuon[i]={false};
    isGoodMuon[i]={false};
    isTightMuon[i]={false};
    muon_pt[i]=0;
    muon_eta[i]=0;
    muon_phi[i]=0;
    muon_px[i]=0;
    muon_py[i]=0;
    muon_pz[i]=0;
    muon_e[i]=0;
    muon_m[i]=0;
    //elec
    elec_pt[i]=0;
    elec_eta[i]=0;
    elec_phi[i]=0;
    elec_px[i]=0;
    elec_py[i]=0;
    elec_pz[i]=0;
    elec_e[i]=0;
    elec_m[i]=0;
    ele_MVA[i]=0;
    ele_passMVAPresel[i]=false;
    ele_PfIso[i]=0;
    isVetoEle[i]=false;    
    isLooseEle[i]=false;    
    isMediumEle[i]=false;    
    isTightEle[i]=false;    
    DEtaSupClusTrkAtVtx[i]=0;
    DPhiSupClusTrkAtVtx[i]=0;
    Sigma2Ieta[i]=0;
    HadOverEm[i]=0;
    TrkDxy[i]=0;
    TrkDz[i]=0;
    ConvRej[i]=0;
    MissConRej[i]=0;
    PfIso[i]=0;
    InvEminusInvPin[i]=0;
    
    jet_pt[i]=0;
    jet_eta[i]=0;
    jet_phi[i]=0;
    jet_px[i]=0;
    jet_py[i]=0;
    jet_pz[i]=0;
    jet_e[i]=0;
    jet_m[i]=0;
    isjetLoose[i]=false;
    isjetTight[i]=false;
    neutHadEnFrac[i]=0;
    neutEmEnFrac[i]=0;
    nDaughter[i]=0;
    chHadEnFrac[i]=0;
    chMultiplicity[i]=0;
    chEmEnFrac[i]=0;
    
  }*/

