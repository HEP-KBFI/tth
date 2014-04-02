// -*- C++ -*-
//
// Package:    CreateNTuple
// Class:      CreateNTuple
// 
/**\class CreateNTuple CreateNTuple.cc tthAnalysis/CreateNTuple/src/CreateNTuple.cc

 Description: 
      Main class that writes out all the leptons and relevant details for later analysis of all kind

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
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
//
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
//
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//
#include "TTree.h"
#include "TFile.h"
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
  edm::InputTag rhoLabel, puLabel, genLabel, vtxLabel, tauLabel, muonLabel, elecLabel, jetLabel, metLabel;
  //evt 
  int run, lumi, ev, np;
  //PU
  double npu, rho;
  //gen
  int gnp, gstatus[MAXPART], gpid[MAXPART], gmother;
  double gpx[MAXPART], gpy[MAXPART], gpz[MAXPART], ge[MAXPART]; 
  //evt vtx
  int nvtx;
  double evx, evy, evz;
  reco::Vertex pv;
  //common 
  int pid[MAXPART];
  double px[MAXPART], py[MAXPART], pz[MAXPART], e[MAXPART];
  double vx[MAXPART], vy[MAXPART], vz[MAXPART];
  double iso[MAXPART], id[MAXPART], disc[MAXPART];
  //tau
  int tau_byDecayModeFinding[MAXPART];
  int tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[MAXPART];
  int tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[MAXPART];
  int tau_byTightCombinedIsolationDeltaBetaCorr3Hits[MAXPART];
  int tau_againstElectronLooseMVA3[MAXPART], tau_againstElectronMediumMVA3[MAXPART];
  int tau_againstElectronTightMVA3[MAXPART], tau_againstElectronVTightMVA3[MAXPART];
  int tau_againstMuonLoose3[MAXPART], tau_againstMuonMedium3[MAXPART],tau_againstMuonTight3[MAXPART];
  //muon
  double sumChHadPt[MAXPART],sumNeutHadEt[MAXPART],sumPhotonEt[MAXPART],sumPUPt[MAXPART];
  double sumChPartPt[MAXPART],muonIsoPflow[MAXPART],muonIsoPflowPUcorr[MAXPART];
  //elec
  double scEta[MAXPART], tIP[MAXPART], convVeto[MAXPART], mHit[MAXPART];
  double chIso03[MAXPART],nhIso03[MAXPART],phIso03[MAXPART],puChIso03[MAXPART],relIso[MAXPART],relIsodb[MAXPART],relIsorho[MAXPART];
  //jet
  int nDaughter[MAXPART];
  double jetb[MAXPART], neutHadEnFrac[MAXPART], neutEmEnFrac[MAXPART], chHadEnFrac[MAXPART],chMultiplicity[MAXPART], chEmEnFrac[MAXPART];
  //met  
  double metx, mety, metcov[4];
  //
  TTree *t;  
  TFile *f;
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
  rhoLabel     (iConfig.getParameter<edm::InputTag>("rhoLabel" )),
  puLabel      (iConfig.getParameter<edm::InputTag>("puLabel"  )),
  genLabel     (iConfig.getParameter<edm::InputTag>("genLabel" )),
  vtxLabel     (iConfig.getParameter<edm::InputTag>("vtxLabel" )),
  tauLabel     (iConfig.getParameter<edm::InputTag>("tauLabel" )),
  muonLabel    (iConfig.getParameter<edm::InputTag>("muonLabel")),
  elecLabel    (iConfig.getParameter<edm::InputTag>("elecLabel")),
  jetLabel     (iConfig.getParameter<edm::InputTag>("jetLabel" )),
  metLabel     (iConfig.getParameter<edm::InputTag>("metLabel" ))
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

  using namespace reco; using namespace edm; using namespace pat; using namespace std; using namespace isodeposit;

  // Re-initialize all variables to zero as a new event has begun
  run=0; lumi=0; ev=0; np=0; rho=0; npu=0; gnp=0; gmother=0, nvtx=0;
  evx=0; evy=0; evz=0; metx=0; mety=0; metcov[3]={0};
  //
  for (int i=0; i<MAXPART; i++) {
    //gen
    gstatus[i]=0; gpid[i]=0; gpx[i]=0; gpy[i]=0; gpz[i]=0; ge[i]=0;  
    //common
    px[i]=0; py[i]=0; pz[i]=0; e[i]=0; pid[i]=0; vx[i]=0; vy[i]=0; vz[i]=0; id[i]=0; disc[i]=0; iso[i]=0;
    //tau
    tau_byDecayModeFinding[i]=0;tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[i]=0; 
    tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[i]=0; tau_byTightCombinedIsolationDeltaBetaCorr3Hits[i]=0;
    tau_againstElectronLooseMVA3[i]=0; tau_againstElectronMediumMVA3[i]=0; 
    tau_againstElectronTightMVA3[i]=0; tau_againstElectronVTightMVA3[i]=0;
    tau_againstMuonLoose3[i]=0; tau_againstMuonMedium3[i]=0; tau_againstMuonTight3[i]=0;
    //muon   
    sumChHadPt[i]=0; sumNeutHadEt[i]=0; sumPhotonEt[i]=0; sumPUPt[i]=0; 
    sumChPartPt[i]=0; muonIsoPflow[i]=0; muonIsoPflowPUcorr[i]=0; 
    //elec
    chIso03[i]=0; nhIso03[i]=0; phIso03[i]=0; puChIso03[i]=0; relIso[i]=0; relIsodb[i]=0; 
    relIsorho[i]=0; mHit[i]=0; convVeto[i]=0; tIP[i]=0; scEta[i]=0;
    //jet
    neutHadEnFrac[i]=0; neutEmEnFrac[i]=0; nDaughter[i]=0; chHadEnFrac[i]=0; 
    chMultiplicity[i]=0; chEmEnFrac[i]=0; jetb[i]=0;
  }

  //evt
  run = iEvent.id().run(); lumi = iEvent.id().luminosityBlock(); ev = iEvent.id().event();
  //cout << "run: " << run << endl << "lumi: " << lumi << endl << "ev: " << ev << endl;   

  //PU correction
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(rhoLabel,rhoHandle);
  rho = *rhoHandle;
  //cout << "rho: " << rho << endl;

  //PU: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFastSimPileUp#Access_to_the_true_pile_up_distr
  if (!iEvent.isRealData()) {
	  Handle<std::vector< PileupSummaryInfo > >  PupInfo;
	  iEvent.getByLabel(puLabel, PupInfo);
	  npu=PupInfo->begin()->getTrueNumInteractions();
	  //cout << "npu: "<<npu<<endl;
  }

  //Gen
  Handle<vector<reco::GenParticle> > genReco; 	
  iEvent.getByLabel(genLabel,genReco);
  for (vector<reco::GenParticle>::const_iterator it = genReco->begin(); it != genReco->end(); it++) {
	  gstatus[gnp]     = it->status();  
	  gpid[gnp]        = it->pdgId();
	  ge[gnp]          = it->energy();
	  gpx[gnp]         = it->px();
	  gpy[gnp]         = it->py();
	  gpz[gnp]         = it->pz();
	  //cout << "gstatus: " << gstatus[gnp] << endl << "gpid: " << gpid[gnp] << endl << "ge: " << ge[gnp] << endl << gpx[gnp] <<", "<<gpy[gnp]<<", "<<gpz[gnp]<<endl;      
	  gnp++ ;
  }

  //evt vtx
  Handle<vector<reco::Vertex>> vtxReco;
  iEvent.getByLabel(vtxLabel,vtxReco);
  if (!vtxReco->size()) return; // We do need a vertex...
  bool foundPV = false;
  double sptSel=0;

  for (vector<reco::Vertex>::const_iterator it = vtxReco->begin(); it != vtxReco->end(); it++) {
	  nvtx++; // count number of vertex in the evt
	  double spt=0; //sum track pt of the vtx

	  for (vector<reco::TrackBaseRef>::const_iterator tr = it->tracks_begin(); tr != it->tracks_end(); tr++){
		  spt+=(*tr)->pt(); 
	  }
	  //Good vtx
	  if (!(it->isValid() && it->ndof() > 5 && fabs(it->z()) < 24 && it->position().rho() < 2)) continue; 

	  if (!foundPV) {
		  pv=*it;
		  sptSel=spt;
		  foundPV=true;
	  }
          //keep the vtx with the highest sum track pt
	  else if (spt>sptSel) {
		  pv=*it;
		  sptSel=spt;
	  }
  } //vertex

	  evx = pv.x();
	  evy = pv.y();
	  evz = pv.z();
	  //cout << "event PV vertex: " << evx << ", "<< evy << ", " << evz  << endl<<endl;

  //tau
  Handle<pat::TauCollection> tauPat;
  iEvent.getByLabel(tauLabel,tauPat);
  for (pat::TauCollection::const_iterator it = tauPat->begin(); it != tauPat->end(); it++) {
          if(!(it->isPFTau() && it->pt()>20)) continue; //tauid recommendation
	  pid[np]     = it->pdgId();
	  px[np]      = it->px();
	  py[np]      = it->py();
	  pz[np]      = it->pz();
	  e[np]       = it->energy();
	  vx[np]      = it->vx();
	  vy[np]      = it->vy();
	  vz[np]      = it->vz();
	  //cout << "pid: " <<pid[np] << endl << "e: " << e[np] << endl << px[np] <<", "<<py[np]<<", "<<pz[np]<<endl;      
	  //cout <<"tauVtx: " << vx[np] << " " << vy[np] << " " << vz[np] << endl << "taupT: "<< it->pt()<<endl;
          //tauId https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation 
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
          //
	  if(tau_byDecayModeFinding[np]==1) id[np] += 1<<0;
	  //cout<< "decay dis:  "<< tau_byDecayModeFinding[np] <<"  tauid decay: "<<id[np] << endl;
	  if(tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[np]==1)  id[np]+= 1<<1;
	  //cout<< "Loose dis:  "<< tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[np]  <<"  id loose:  "<<id[np] << endl;
	  if(tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[np]==1) id[np]+= 1<<2;
	  //cout<< "Medium dis: "<< tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[np] <<"  id medium: "<<id[np] << endl;
	  if(tau_byTightCombinedIsolationDeltaBetaCorr3Hits[np]==1)   id[np]+= 1<<3;
	  //cout<< "Tight dis:  "<< tau_byTightCombinedIsolationDeltaBetaCorr3Hits[np]  <<"  id Tight:  "<<id[np] << endl;
	  if(tau_againstElectronLooseMVA3[np]==1)  id[np]+= 1<<4;
	  //cout<<"against Loose ele:   "<<tau_againstElectronLooseMVA3[np] << "  id: "<<id[np] << endl;
	  if(tau_againstElectronMediumMVA3[np]==1) id[np]+= 1<<5;
	  //cout<<"against Medium ele:  "<<tau_againstElectronMediumMVA3[np]<< "  id: "<<id[np] << endl;
	  if(tau_againstElectronTightMVA3[np]==1)  id[np]+= 1<<6;
	  //cout<<"against Tight ele:   "<<tau_againstElectronTightMVA3[np] << "  id: "<<id[np] << endl;
	  if(tau_againstElectronVTightMVA3[np]==1) id[np]+= 1<<7;
	  //cout<<"against VTight ele:  "<<tau_againstElectronVTightMVA3[np]<< "  id: "<<id[np] << endl;
	  if(tau_againstMuonLoose3[np]==1)  id[np]+= 1<<8;
	  //cout<<"against Loose muon:  "<<tau_againstMuonLoose3[np] <<"  id: "<<id[np] << endl;
	  if(tau_againstMuonMedium3[np]==1) id[np]+= 1<<9;
	  //cout<<"against Medium muon: "<<tau_againstMuonMedium3[np]<<"  id: "<<id[np] << endl;
	  if(tau_againstMuonTight3[np]==1)  id[np]+= 1<<10;
	  //cout<<"against Tight muon:  "<<tau_againstMuonTight3[np] <<"  id: "<<id[np] << endl;
	  np++;
  }//tau

  //muon
  Handle<pat::MuonCollection> muonPat;
  iEvent.getByLabel(muonLabel,muonPat);
  for (pat::MuonCollection::const_iterator it = muonPat->begin(); it != muonPat->end(); it++) {
	  if(!(it->isPFMuon() && it->pt()>10)) continue;
	  pid[np]     = it->pdgId();
	  px[np]      = it->px();
	  py[np]      = it->py();
	  pz[np]      = it->pz();
	  e[np]       = it->energy();
	  vx[np]      = it->vx();
	  vy[np]      = it->vy();
	  vz[np]      = it->vz();
	  //cout << "pid: " <<pid[np] << endl << "e: " << e[np] << endl << px[np] <<", "<<py[np]<<", "<<pz[np]<<endl;      
	  //cout <<"vtx: " << vx[np] << " " << vy[np] << " " << vz[np] << endl;
	  //muonId: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel#Muons
          if(it->isLooseMuon() && it->pt()>20 && abs(it->eta())<2.4)
	  id[np]+= 1<<0; //ll
	  //cout<<"isLooseMuon: "<<it->isLooseMuon() << " " << it->pt() << " " << it->eta() << "id: "<<id[np] <<endl;
          if(it->isTightMuon(*vtxReco->begin()) && it->pt()>26 && abs(it->eta())<2.1)
	  id[np]+= 1<<1; //l+j
	  //cout<<"isTightMuon: "<<it->isTightMuon(*vtxReco->begin()) << " " << it->pt() << " " << it->eta()<< "id: "<<id[np] << endl;
          if(it->isLooseMuon() && it->pt()>10 && abs(it->eta())<2.5)
	  id[np]+= 1<<2; //veto l+j
	  //cout<<"isVetoLooseMuon: "<<it->isLooseMuon() << " " << it->pt() << " " << it->eta()<< "id: "<<id[np] << endl;
	  //muon combined relative isolation with Delta-Beta correction: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
	  sumChHadPt[np]          =  it->pfIsolationR04().sumChargedHadronPt;
	  sumNeutHadEt[np]        =  it->pfIsolationR04().sumNeutralHadronEt;
	  sumPhotonEt[np]         =  it->pfIsolationR04().sumPhotonEt;
	  sumPUPt[np]             =  it->pfIsolationR04().sumPUPt;
	  sumChPartPt[np]         =  it->pfIsolationR04().sumChargedParticlePt;
	  if(!(sumChPartPt[np]==0)){
		  muonIsoPflow[np]        =  (sumChHadPt[np]+sumNeutHadEt[np]+sumPhotonEt[np])/sumChPartPt[np];
		  muonIsoPflowPUcorr[np]  =  (sumChHadPt[np]+std::max(0., (sumNeutHadEt[np]+sumPhotonEt[np]-0.5*sumPUPt[np])))/sumChPartPt[np];
	  }
	  iso[np]=muonIsoPflowPUcorr[np]; //dilepton=0.2, lepton+jet=0.12, (veto 0.2 for ele/muon+jet)
          //cout << "muon CorrIso: "<< iso[np] << ", noCorrIso: " <<muonIsoPflow[np] <<endl;
	  np++;
  }//muon

  Handle<ElectronCollection> elecPat;
  iEvent.getByLabel(elecLabel,elecPat);
  for (pat::ElectronCollection::const_iterator it = elecPat->begin(); it != elecPat->end(); it++) {
	  pid[np]     = it->pdgId();
	  px[np]      = it->px();
	  py[np]      = it->py();
	  pz[np]      = it->pz();
	  e[np]       = it->energy();
	  vx[np]      = it->vx();
	  vy[np]      = it->vy();
	  vz[np]      = it->vz();
	  cout << "pid: " <<pid[np] << endl << "e: " << e[np] << endl << px[np] <<", "<<py[np]<<", "<<pz[np]<<endl;      
	  cout <<"vtx: " << vx[np] << " " << vy[np] << " " << vz[np] << endl;
          //elecId
          scEta[np]   = it->superCluster()->eta();
          tIP[np]     = it->gsfTrack()->dxy(vtxReco->begin()->position());
          convVeto[np]= it->passConversionVeto();
          mHit[np]    = it->gsfTrack()->trackerExpectedHitsInner().numberOfHits();

	  if( it->pt()>20 && abs(it->eta())<2.5 && tIP[np]<0.04 && convVeto && mHit[np]<=0 )
		  id[np] += 1<<0; //dilepton
	  //cout << "id1: "<< id[np] <<endl <<endl;
	  if( it->pt()>30 && abs(it->eta())<2.5 && (!(abs(scEta[np])>1.4442 && abs(scEta[np])<1.5660)) && tIP[np]<0.02 && convVeto && mHit[np]<=0 ) 
		  id[np] += 1<<1; //lepton+jets
	  //cout << "id2: "<< id[np] <<endl <<endl;
	  if( it->pt()>10 && abs(it->eta())<2.5)
		  id[np] += 1<<2; //veto dilepton
	  //cout << "id3: "<< id[np] <<endl <<endl;
	  if( it->pt()>20 && abs(it->eta())<2.5)
		  id[np] += 1<<3; //veto lepton+jets
	  //cout << "id4: "<< id[np] <<endl <<endl;

          //eleDisc(MVA)
          disc[np]  = it->electronID("mvaTrigV0");  //dilepton>0.5; lepton+jets>0.5(or 0.9); for veto need to be study
          //cout << "ele disc: "<< disc[np] << endl;
          //eleIso: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PfIsolation
	  float AEff03 = 0.00;
	  if(iEvent.isRealData()){
		  AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, it->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2011);
	  }else{
		  AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, it->superCluster()->eta(), ElectronEffectiveArea::kEleEAFall11MC);
	  }
	  //cout << "AEff03"<< AEff03 << endl; 
	  chIso03[np]   =  it->chargedHadronIso();
	  nhIso03[np]   =  it->neutralHadronIso();
	  phIso03[np]   =  it->photonIso();
	  puChIso03[np] =  it->puChargedHadronIso();
	  if(!(it->pt()==0)){
		  relIso[np]    =  ( chIso03[np] + nhIso03[np] + phIso03[np] ) / it->pt() ;
		  relIsodb[np]  =  ( chIso03[np] + max(0.0, nhIso03[np] + phIso03[np] - 0.5*puChIso03[np]) )/ it->pt();
		  relIsorho[np] =  ( chIso03[np] + max(0.0, nhIso03[np] + phIso03[np] - rho*AEff03) )/ it->pt();
	  }
	  iso[np] = relIsorho[np]; // ll(<0.15), l+j(<0.1), veto ll&l+j(<0.15)
	  //cout << "ele CorrIso: "<< iso[np] << ", noCorrIso: " <<relIso[np] <<endl;
	  np++;
  }

  Handle<pat::JetCollection> jetPat;  
  iEvent.getByLabel(jetLabel,jetPat);
  for (pat::JetCollection::const_iterator it = jetPat->begin(); it != jetPat->end(); it++) {
	  pid[np]     = it->pdgId();
	  px[np]      = it->px();
	  py[np]      = it->py();
	  pz[np]      = it->pz();
	  e[np]       = it->energy();
	  vx[np]      = it->vx();
	  vy[np]      = it->vy();
	  vz[np]      = it->vz();
	  // 
	  neutHadEnFrac[np]  = it->neutralHadronEnergyFraction();
	  neutEmEnFrac[np]   = it->neutralEmEnergyFraction();
	  nDaughter[np]      = it->numberOfDaughters();
	  chHadEnFrac[np]    = it->chargedHadronEnergyFraction();
	  chMultiplicity[np] = it->chargedMultiplicity();
	  chEmEnFrac[np]     = it->chargedEmEnergyFraction();
	  //id
	  if ( ( abs(it->eta())<2.4 && neutHadEnFrac[np]<0.99 && neutEmEnFrac[np]<0.99 &&
				  nDaughter[np]>1 && chHadEnFrac[np]>0 && chMultiplicity[np]>0 && chEmEnFrac[np]<0.99) || 
			  ( abs(it->eta())>2.4 && neutHadEnFrac[np]<0.99 && neutEmEnFrac[np]<0.99 && nDaughter[np]>1) );  //Loose (recommended)
	  id[np] += 1<<0;
	  //Disc(btag)
	  disc[np]  = it->bDiscriminator("combinedSecondaryVertexBJetTags"); //Loose=0.244, Medium=0.679, Tight=0.898  
	  //cout << "disc bjet: " << disc[np] << endl;
	  np++;
  }

  //MET
  Handle<pat::METCollection> met;
  iEvent.getByLabel(metLabel,met);
  metx=met->begin()->px();
  mety=met->begin()->py();
  metcov[0]= (met->front()).getSignificanceMatrix()(0,0);
  metcov[1]= (met->front()).getSignificanceMatrix()(1,0);
  metcov[2]= (met->front()).getSignificanceMatrix()(0,1);
  metcov[3]= (met->front()).getSignificanceMatrix()(1,1);
  //cout <<"met: " << metcov[0] << ", "<<metcov[1] << ", "<<metcov[2] << ", "<<metcov[3] <<endl;

  t->Fill(); 

}


// ------------ method called once each job just before starting event loop  ------------
void 
CreateNTuple::beginJob()
{
  f = new TFile("ntuple.root","RECREATE");
  t = new TTree("tree", "NTuple");
  t->Branch("ev",&ev,"ev/I");
  t->Branch("np",&np,"np/I");
  t->Branch("lumi",&lumi,"lumi/I");
  t->Branch("run",&run,"run/I");
  t->Branch("disc",disc,"disc[np]/D");
  t->Branch("e",e,"e[np]/D");
  t->Branch("evx",&evx,"evx/D");
  t->Branch("evy",&evy,"evy/D");
  t->Branch("evz",&evz,"evz/D");
  t->Branch("gnp",&gnp,"gnp/I");
  t->Branch("gpid",gpid,"pid[gnp]/I");
  t->Branch("ge",ge,"ge[gnp]/D");
  t->Branch("gpx",gpx,"gpx[gnp]/D");
  t->Branch("gpy",gpy,"gpy[gnp]/D");
  t->Branch("gpz",gpz,"gpz[gnp]/D");
  t->Branch("gstatus",gstatus,"gstatus[gnp]/I");
  t->Branch("id",id,"id[np]/I");
  t->Branch("iso",iso,"iso[np]/D");
  t->Branch("metcov",metcov,"metcov[4]/D");
  t->Branch("metx",&metx,"metx/D");
  t->Branch("mety",&mety,"mety/D");
  t->Branch("npu",&npu,"npu/D");
  t->Branch("nvtx",&nvtx,"nvtx/I");
  t->Branch("pid",pid,"pid[np]/I");
  t->Branch("px",px,"px[np]/D");
  t->Branch("py",py,"py[np]/D");
  t->Branch("pz",pz,"pz[np]/D");
  t->Branch("rho",&rho,"rho/D");
  t->Branch("vx",vx,"vx[np]/D");
  t->Branch("vy",vy,"vy[np]/D");
  t->Branch("vz",vz,"vz[np]/D");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CreateNTuple::endJob() 
{
t->Write();
f->Close();
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


