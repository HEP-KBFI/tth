// -*- C++ -*-
//
// Package:    CreateNTuple
// Class:      CreateNTuple
// 
/**\class CreateNTuple CreateNTuple.cc tthAnalysis/CreateNTuple/src/CreateNTuple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Betty Calpas
//         Created:  Tue Mar 25 14:22:19 EET 2014
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

/*
cout<< "New tau:  tauid: "<<tauId[np] << endl;

    if(tau_byDecayModeFinding[np]==1) tauId[np] += 1;
cout<< "decay dis: "<< tau_byDecayModeFinding[np] <<"  tauid decay: "<<tauId[np] << endl;

    if(tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[np]==1) tauId[np]+= 1<<1;
cout<< "Loose dis: "<< tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[np] <<"   tauid loose: "<<tauId[np] << endl;

    if(tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[np]==1) tauId[np]+= 1<<2;
cout<< "Medium dis: "<< tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[np] <<"   tauid medium: "<<tauId[np] << endl;

   if(tau_byTightCombinedIsolationDeltaBetaCorr3Hits[np]==1) tauId[np]+= 1<<3;
cout<< "Tight dis: "<< tau_byTightCombinedIsolationDeltaBetaCorr3Hits[np] <<"   tauid Tight: "<<tauId[np] << endl;

    if(tau_againstElectronLooseMVA3[np]==1) tauId[np]+= tauId[np]<<1;
cout<<"against Loose ele: "<<tau_againstElectronLooseMVA3[np]<< "  tauId: "<<tauId[np] << endl;

    if(tau_againstElectronMediumMVA3[np]==1) tauId[np]+= tauId[np]<<1;
cout<<"against Medium ele: "<<tau_againstElectronMediumMVA3[np]<< "  tauId: "<<tauId[np] << endl;

    if(tau_againstElectronTightMVA3[np]==1) tauId[np]+= tauId[np]<<1;
cout<<"against Tight ele: "<<tau_againstElectronTightMVA3[np]<< "  tauId: "<<tauId[np] << endl;

    if(tau_againstElectronVTightMVA3[np]==1) tauId[np]+= tauId[np]<<1;
cout<<"against VTight ele: "<<tau_againstElectronVTightMVA3[np]<< "  tauId: "<<tauId[np] << endl;

    if(tau_againstMuonLoose3[np]==1) tauId[np]+= tauId[np]<<1;
cout<<"against Loose muon: "<<tau_againstMuonLoose3[np]<<"  tauId: "<<tauId[np] << endl;

    if(tau_againstMuonMedium3[np]==1) tauId[np]+= tauId[np]<<1;
cout<<"against Medium muon: "<<tau_againstMuonMedium3[np]<<"  tauId: "<<tauId[np] << endl;

    if(tau_againstMuonTight3[np]==1) tauId[np]+= tauId[np]<<1;
cout<<"against Tight muon: "<<tau_againstMuonTight3[np]<<"  tauId: "<<tauId[np] << endl;
    np++;
*/
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



Handle<ElectronCollection> elecPat;
iEvent.getByLabel(elecLabel,elecPat);

  for (pat::ElectronCollection::const_iterator it = elecPat->begin(); it != elecPat->end(); it++) {
     pid[np]=it->pdgId();
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


 h->Fill(np);
 t->Fill(); 

/*
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
*/

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
