// Description: class that writes out leptons and relevant information for later analysis 
// Original Author: Betty Calpas, Morten Piibeleht
//         Created: Tue Mar 25 14:22:19 EET 2014

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

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <string>

//
// class declaration
//
const size_t maxP = 1000;

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

  int storeGenParticle(const reco::Candidate& p);

  // ----------member data ---------------------------
  edm::InputTag rhoLabel, puLabel, genLabel, vtxLabel, tauLabel, muonLabel, elecLabel, jetLabel, metLabel, trigLabel;
  std::string muLabel, elLabel;
  //evt
  int run, lumi, ev, np;
  //PU
  double npu, rho, puw, ttw;
  //gen
  int gnp, gstatus[maxP], gpid[maxP], gmother[maxP];
  const reco::Candidate * genparticle_cands[maxP];
  double gpx[maxP], gpy[maxP], gpz[maxP], ge[maxP];
  //vtx
  int nvtx;
  double evx, evy, evz;
  //common
  int pid[maxP], id[maxP];
  double px[maxP], py[maxP], pz[maxP], e[maxP];
  double vx[maxP], vy[maxP], vz[maxP];
  float iso[maxP], disc[maxP];
  //elec
  double scEta[maxP], tIP[maxP];
  bool convVeto[maxP];
  int mHit[maxP];
  float AEff03;
  //jet
  int nDaughter[maxP];
  float neutHadEnFrac[maxP], neutEmEnFrac[maxP], chHadEnFrac[maxP],chMultiplicity[maxP], chEmEnFrac[maxP];
  //met  
  double metx, mety, metcov[4];
  //other
  std::vector<std::string> muUsed, elUsed;
  int trig;
  bool muPass, elPass;
  //
  TTree *t;
  TFile *f;
  TH1D *nlep;

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
	metLabel     (iConfig.getParameter<edm::InputTag>("metLabel" )),
	trigLabel    (iConfig.getParameter<edm::InputTag>("trigLabel")),
	muLabel      (iConfig.getParameter<std::string>("muLabel")),
	elLabel      (iConfig.getParameter<std::string>("elLabel"))

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
  run=0; lumi=0; ev=0; np=0; rho=0; npu=0; puw=1, ttw=1; gnp=0; nvtx=0;
  evx=0; evy=0; evz=0; metx=0; mety=0; metcov[3]={0}; trig=0;
  
  for (size_t i=0; i<maxP; i++){ 
    gstatus[i]=0; gpid[i]=0; gpx[i]=0; gpy[i]=0; gpz[i]=0; ge[i]=0; gmother[i]=0;
    px[i]=0; py[i]=0; pz[i]=0; e[i]=0; pid[i]=-90; vx[i]=0; vy[i]=0; vz[i]=0; id[i]=0; disc[i]=0; iso[i]=0; 
    mHit[i]=0; convVeto[i]=false; tIP[i]=0; scEta[i]=0;
    neutHadEnFrac[i]=0; neutEmEnFrac[i]=0; nDaughter[i]=0; chHadEnFrac[i]=0; chMultiplicity[i]=0; chEmEnFrac[i]=0;
  }

  run= iEvent.id().run(); lumi= iEvent.id().luminosityBlock(); ev= iEvent.id().event();

  //trigger
  Handle<TriggerEvent> trigPat;
  iEvent.getByLabel(trigLabel,trigPat);
  //for (unsigned int i=0; i<trigPat->paths()->size(); i++) cout << "Path: " << trigPat->paths()->at(i).name() << endl;

  muUsed.clear(); elUsed.clear();
  for ( unsigned int j=0; j<trigPat->paths()->size(); j++) {
	  if (trigPat->paths()->at(j).name().find(muLabel) != std::string::npos) muUsed.push_back(trigPat->paths()->at(j).name());
	  if (trigPat->paths()->at(j).name().find(elLabel) != std::string::npos) elUsed.push_back(trigPat->paths()->at(j).name());
  }

  // Go over the Used list and check if any of the triggers is unprescaled and fired
  muPass=false;
  for ( unsigned int i=0; i<muUsed.size(); i++) {
          std::string t = muUsed[i];
	  if (trigPat->path(t)) { if (trigPat->path(t)->wasAccept() && trigPat->path(t)->prescale() == 1 ) muPass=true; }
  }
  if(muPass) trig+=1<<0;

  elPass=false;
  for ( unsigned int i=0; i<elUsed.size(); i++) {
	  std::string t = elUsed[i];
	  if (trigPat->path(t)) { if (trigPat->path(t)->wasAccept() && trigPat->path(t)->prescale() == 1 ) elPass=true; }
  }
  if(elPass) trig+=1<<1;


  //PU correction
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(rhoLabel,rhoHandle);
  rho = *rhoHandle;

  //PU: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFastSimPileUp#Access_to_the_true_pile_up_distr
  if (!iEvent.isRealData()) {

	  Handle<std::vector< PileupSummaryInfo > >  PupInfo; iEvent.getByLabel(puLabel, PupInfo);
	  npu=PupInfo->begin()->getTrueNumInteractions();

	  // Write out the generator level particles. Since we try to figure out the mother-daughter relationshipts also,
	  // the storeGenParticle() method has on O(n) complexity and this brings the total complexity of the loop to O(n^2). 
	  //This might become a bottleneck if the number of gen. level particles becomes large (e.g. if no filtering is done). 
	  //However, there does not seem to be a way around due to the way mother-daughter relationships are stored in the 
	  //reco::GenParticle class (it is necessary to loop over a lookup table between the object pointers and the new indices).
	  Handle<vector<reco::GenParticle> > genReco; iEvent.getByLabel(genLabel,genReco);
	  for(const reco::GenParticle& genparticle : *genReco) { storeGenParticle(genparticle); }
	  if((int)genReco->size() != gnp) {
		  LogWarning("CreateNTuple") << "genReco->size() and gnp do not match (" << genReco->size() << " vs " << gnp << ")";
	  }
  }

  Handle<vector<reco::Vertex>> vtxReco; iEvent.getByLabel(vtxLabel,vtxReco);
  if (!vtxReco.isValid()) {
	  std::cout<<"Didja hear the one about the empty vertex collection?\n";
	  return;
  }
  // require in the event that there is at least one reconstructed vertex
  if(vtxReco->size()<=0) return;
  // pick the first (i.e. highest sum pt) vertex
  const reco::Vertex* theVertex=&(vtxReco->front());
  // require that the vertex meets certain criteria
  if(theVertex->ndof()<5) return;
  if(fabs(theVertex->z())>24.0) return;
  if(fabs(theVertex->position().rho())>2.0) return;
  //get the vertex value
  evx = theVertex->x();
  evy = theVertex->y();
  evz = theVertex->z();

  std::vector<reco::Vertex>::const_iterator itv;
  int NVtx = 0;
  // now, count vertices
  for (itv = vtxReco->begin(); itv != vtxReco->end(); ++itv) {
	  // require that the vertex meets certain criteria
	  if(itv->ndof()<5) continue;
	  if(fabs(itv->z())>50.0) continue;
	  if(fabs(itv->position().rho())>2.0) continue;
	  ++NVtx;
  }
  nvtx=NVtx;


  //tau
  Handle<pat::TauCollection> tauPat; iEvent.getByLabel(tauLabel,tauPat);
  for (pat::TauCollection::const_iterator it = tauPat->begin(); it != tauPat->end(); it++) {
          if(!(it->isPFTau() && it->pt()>20 && fabs(it->pdgId())==15 )) continue; 
	  pid[np]     = it->pdgId();
	  px[np]      = it->px();
	  py[np]      = it->py();
	  pz[np]      = it->pz();
	  e[np]       = it->energy();
	  vx[np]      = it->vx();
	  vy[np]      = it->vy();
	  vz[np]      = it->vz();
	  if( it->tauID("decayModeFinding")==1) id[np] += 1<<0; 
	  if( it->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") ==1)  id[np]+= 1<<1;  
	  if( it->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")==1)  id[np]+= 1<<2; 
	  if( it->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits") ==1)  id[np]+= 1<<3; 
	  if( it->tauID("againstElectronLooseMVA3") ==1) id[np]+= 1<<4; 
	  if( it->tauID("againstElectronMediumMVA3")==1) id[np]+= 1<<5; 
	  if( it->tauID("againstElectronTightMVA3") ==1) id[np]+= 1<<6;
	  if( it->tauID("againstElectronVTightMVA3")==1) id[np]+= 1<<7;
	  if( it->tauID("againstMuonLoose3") ==1) id[np]+= 1<<8; //Medium 3 does not exist  
	  if( it->tauID("againstMuonTight3") ==1) id[np]+= 1<<9;
          if( it->genJet()) id[np]+= 1<<10;

	  np++;
  }

  //muon
  Handle<pat::MuonCollection> muonPat; iEvent.getByLabel(muonLabel,muonPat);
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
	  //muonId: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel#Muons
          if(it->isLooseMuon()   && fabs(it->eta())<2.4 /*&& it->pt()>20*/)         id[np]+= 1<<0; //ll
          if(it->isTightMuon(*theVertex) && fabs(it->eta())<2.1 /*&& it->pt()>26*/) id[np]+= 1<<1; //l+j
          if(it->isLooseMuon()   && fabs(it->eta())<2.5 /*&& it->pt()>10*/)         id[np]+= 1<<2; //veto l+j == ll
          if(it->genLepton())                                                       id[np]+= 1<<3;
	  //muon combRelIso w/ DB correction: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation
	  iso[np] = (it->pfIsolationR04().sumChargedHadronPt + std::max(0., (it->pfIsolationR04().sumNeutralHadronEt 
                     + it->pfIsolationR04().sumPhotonEt - 0.5*it->pfIsolationR04().sumPUPt )))/it->pt(); //ll=0.2, l+j=0.12, (veto 0.2 for l+j)
 
	  np++;
  }

  //ele
  Handle<ElectronCollection> elecPat; iEvent.getByLabel(elecLabel,elecPat);
  for (pat::ElectronCollection::const_iterator it = elecPat->begin(); it != elecPat->end(); it++) {
	  if(!(it->isPF() && it->pt()>10)) continue;
	  pid[np]     = it->pdgId();
	  px[np]      = it->px();
	  py[np]      = it->py();
	  pz[np]      = it->pz();
	  e[np]       = it->energy();
	  vx[np]      = it->vx();
	  vy[np]      = it->vy();
	  vz[np]      = it->vz();
          //elecId
          scEta[np]   = it->superCluster()->eta();
          tIP[np]     = it->gsfTrack()->dxy(vtxReco->begin()->position());
          convVeto[np]= it->passConversionVeto();
          mHit[np]    = it->gsfTrack()->trackerExpectedHitsInner().numberOfHits();

	  if( fabs(it->eta())<2.5 && tIP[np]<0.04 && convVeto[np]==true && mHit[np]<=0 /*&& it->pt()>20*/) id[np] += 1<<0; //ll

	  if( fabs(it->eta())<2.5 && tIP[np]<0.02 && convVeto[np]==true && mHit[np]<=0 /*&& it->pt()>30*/
              && (!(fabs(scEta[np])>1.4442 && fabs(scEta[np])<1.5660)) )                                   id[np] += 1<<1; //l+j

	  if( fabs(it->eta())<2.5 /*&& it->pt()>10||20*/) id[np] += 1<<2; //veto: ll (pt>10), l+j (pt>20)

          if( it->genLepton())                            id[np] += 1<<3;

          //eleDisc(MVA)
          disc[np]  = it->electronID("mvaTrigV0");  //ll>0.5; l+j>0.5(or 0.9); for veto need to be study
  
	  //eleIso: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PfIsolation
	  AEff03 = 0.;
	  AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, 
                   it->superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);

	  iso[np] = ( it->chargedHadronIso() + std::max(0., it->neutralHadronIso() + it->photonIso() - rho*AEff03 ))/it->pt(); // ll(<0.15), l+j(<0.1), veto ll&l+j(<0.15)

	  np++;
  }

  Handle<pat::JetCollection> jetPat; iEvent.getByLabel(jetLabel,jetPat);
  for (pat::JetCollection::const_iterator it = jetPat->begin(); it != jetPat->end(); it++) {
	  if(!(it->isPFJet() && it->pt()>10)) continue;
	  pid[np]     = 81;
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
	  if ( ( fabs(it->eta())<2.4 && neutHadEnFrac[np]<0.99 && neutEmEnFrac[np]<0.99 &&
				  nDaughter[np]>1 && chHadEnFrac[np]>0 && chMultiplicity[np]>0 && chEmEnFrac[np]<0.99) ||
			  ( fabs(it->eta())>2.4 && neutHadEnFrac[np]<0.99 && neutEmEnFrac[np]<0.99 && nDaughter[np]>1) ) id[np] += 1<<0; //L
	
	  if ( ( fabs(it->eta())<2.4 && neutHadEnFrac[np]<0.95 && neutEmEnFrac[np]<0.95 &&
				  nDaughter[np]>1 && chHadEnFrac[np]>0 && chMultiplicity[np]>0 && chEmEnFrac[np]<0.99) ||
			  ( fabs(it->eta())>2.4 && neutHadEnFrac[np]<0.95 && neutEmEnFrac[np]<0.95 && nDaughter[np]>1) ) id[np] += 1<<1; //M

	  if ( ( fabs(it->eta())<2.4 && neutHadEnFrac[np]<0.90 && neutEmEnFrac[np]<0.90 &&
				  nDaughter[np]>1 && chHadEnFrac[np]>0 && chMultiplicity[np]>0 && chEmEnFrac[np]<0.99) ||
			  ( fabs(it->eta())>2.4 && neutHadEnFrac[np]<0.90 && neutEmEnFrac[np]<0.90 && nDaughter[np]>1) )  id[np] += 1<<2; //T

          if( it->genJet()) id[np]+= 1<<3;

	  //Disc(btag)
	  disc[np]  = it->bDiscriminator("combinedSecondaryVertexBJetTags"); // > Loose=0.244, Medium=0.679, Tight=0.898  

	  np++;
  }

  //MET
  Handle<pat::METCollection> met; iEvent.getByLabel(metLabel,met);

  if( met->begin()->isPFMET()==true){
  metx=met->begin()->px();
  mety=met->begin()->py();

  metcov[0]= (met->front()).getSignificanceMatrix()(0,0);
  metcov[1]= (met->front()).getSignificanceMatrix()(1,0);
  metcov[2]= (met->front()).getSignificanceMatrix()(0,1);
  metcov[3]= (met->front()).getSignificanceMatrix()(1,1);
  //cout <<"met: " << metcov[0] << ", "<<metcov[1] << ", "<<metcov[2] << ", "<<metcov[3] <<endl;
}

  nlep->Fill(np);
  t->Fill();

}

int
CreateNTuple::storeGenParticle(const reco::Candidate& p)
{
	/*
	 * Stores a gen. particle and returns its index.
	 *
	 * If the particle is already stored, just returns the previous index.
	 * If the particle has a single mother, it first stores the mother,
	 * so that it could get the mother's index. The mother's index will be
	 * -1 if the particle does not have a mother and -2 if the particle
	 * has multiple mothers or -3 if the particle has a single mother but
	 * with issues (such as being the same as the particle or a null pointer).
	 */

	// check if it is already processed:
	for(int i=0; i<gnp; i++) {
		if(genparticle_cands[i] == &p) {
			return i;
		}
	}

	// parse the mothers
	int mother_index = -999; // -999 should never appear
	if(p.numberOfMothers() == 0) {
		mother_index = -1;
	} else if(p.numberOfMothers() == 1) {
		// first, make sure that the particle is not its own mother or did not
		// come into being from the Great Void of /dev/null (both cases are
		// apparently possible according to some scholars)
		if(p.mother() == nullptr) {
			edm::LogWarning("storeGenParticle") << "mother is a null pointer";
			mother_index = -3;
		} else if(p.mother() == &p) {
			edm::LogWarning("storeGenParticle") << "particle is its own mother";
			mother_index = -3;
		} else {
			mother_index = storeGenParticle(*p.mother());
		}
	} else {
		// otherwise, since p.numberOfMothers() is of type size_t, the value can
		// only be larger than one and we do not parse multiple mothers
		mother_index = -2;
	}

	// "reserve" an index and store necessary values
	int index = gnp++;
	genparticle_cands[index] = &p;
	gstatus[index] = p.status();
	gpid[index] = p.pdgId();
	ge[index] = p.energy();
	gpx[index] = p.px();
	gpy[index] = p.py();
	gpz[index] = p.pz();
	gmother[index] = mother_index;

	return index;
}

// ------------ method called once each job just before starting event loop  ------------
void
CreateNTuple::beginJob()
{

edm::Service<TFileService> fs;
nlep   = fs->make<TH1D>("nlep","Number of leptons",maxP,0,maxP);

t = new TTree("tree", "NTuple");

  t->Branch("ev",&ev,"ev/I");
  t->Branch("np",&np,"np/I");
  t->Branch("lumi",&lumi,"lumi/I");
  t->Branch("run",&run,"run/I");
  t->Branch("disc",disc,"disc[np]/F");
  t->Branch("e",e,"e[np]/D");
  t->Branch("evx",&evx,"evx/D");
  t->Branch("evy",&evy,"evy/D");
  t->Branch("evz",&evz,"evz/D");
  t->Branch("gnp",&gnp,"gnp/I");
  t->Branch("gpid",gpid,"gpid[gnp]/I");
  t->Branch("ge",ge,"ge[gnp]/D");
  t->Branch("gmother",gmother,"gmother[gnp]/I");
  t->Branch("gpx",gpx,"gpx[gnp]/D");
  t->Branch("gpy",gpy,"gpy[gnp]/D");
  t->Branch("gpz",gpz,"gpz[gnp]/D");
  t->Branch("gstatus",gstatus,"gstatus[gnp]/I");
  t->Branch("id",id,"id[np]/I");
  t->Branch("iso",iso,"iso[np]/F");
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
  t->Branch("trig",&trig,"trig/I");
  t->Branch("vx",vx,"vx[np]/D");
  t->Branch("vy",vy,"vy[np]/D");
  t->Branch("vz",vz,"vz[np]/D");
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

