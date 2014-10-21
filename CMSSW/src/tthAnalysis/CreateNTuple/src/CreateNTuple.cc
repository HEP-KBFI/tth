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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//
#include "TTree.h"
#include "TFile.h"
#include <string>
#include <iostream>

//
// class declaration
//
const size_t part = 10000;

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
	edm::InputTag rho_, pu_, gen_, vtx_, tau_, muon_, ele_, jet_, met_, trig_;
	// event info
	int run, lumi, ev, np;
	// pileup
	double npu, rho;
	// generator
	int gnp, gstatus[part], gpid[part], gmother[part];
	const reco::Candidate * genparticle_cands[part];
	double gpx[part], gpy[part], gpz[part], ge[part];
	// vertex
	int nvtx;
	double evx, evy, evz;
	// particle (electron, muon, tau, jet)
	int pid[part], id[part];

	double px[part], py[part], pz[part], e[part];
	double vx[part], vy[part], vz[part];
	float iso[part], disc[part];
	// missing ET  
	double metx, mety, metcov[4];
	// trigger
	std::string smtrig_, setrig_;
	std::vector<std::string> smUsed, seUsed;
	int trig;
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
    rho_     (iConfig.getParameter<edm::InputTag>("rho")),
    pu_      (iConfig.getParameter<edm::InputTag>("pu")),
    gen_     (iConfig.getParameter<edm::InputTag>("gen")),
    vtx_     (iConfig.getParameter<edm::InputTag>("vtx")),
    tau_     (iConfig.getParameter<edm::InputTag>("tau")),
    muon_    (iConfig.getParameter<edm::InputTag>("muon")),
    ele_     (iConfig.getParameter<edm::InputTag>("elec")),
    jet_     (iConfig.getParameter<edm::InputTag>("jet")),
    met_     (iConfig.getParameter<edm::InputTag>("met")),
    trig_    (iConfig.getParameter<edm::InputTag>("trig")),
    smtrig_  (iConfig.getParameter<std::string>("smtrig")),
    setrig_  (iConfig.getParameter<std::string>("setrig"))
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
 
   using namespace reco; using namespace edm; using namespace pat; using namespace std;

    // Re-initialize all variables to zero as a new event has begun
    run=0; lumi=0; ev=0; np=0; rho=0; npu=0; gnp=0; nvtx=0; evx=0; evy=0; evz=0; metx=0; mety=0; metcov[3]={0}; trig=0;
    for (size_t i=0; i<part; i++){ 
	gstatus[i]=0; gpid[i]=0; gpx[i]=0; gpy[i]=0; gpz[i]=0; ge[i]=0; gmother[i]=0;
	px[i]=0; py[i]=0; pz[i]=0; e[i]=0; pid[i]=-90; vx[i]=0; vy[i]=0; vz[i]=0; id[i]=0; disc[i]=0; iso[i]=0; 
    }

    run= iEvent.id().run(); lumi= iEvent.id().luminosityBlock(); ev= iEvent.id().event();

    /////////////
    // trigger //
    /////////////

    Handle<TriggerEvent> trigger;
    iEvent.getByLabel(trig_,trigger);
    //for (unsigned int i=0; i<trigger->paths()->size(); i++) cout << "Path: " << trigger->paths()->at(i).name() << endl;

    // look for defined triggers and store them
    smUsed.clear(); seUsed.clear();
    for ( unsigned int j=0; j<trigger->paths()->size(); j++) {
	if (trigger->paths()->at(j).name().find(smtrig_) != std::string::npos) smUsed.push_back(trigger->paths()->at(j).name());
	if (trigger->paths()->at(j).name().find(setrig_) != std::string::npos) seUsed.push_back(trigger->paths()->at(j).name());
    }

    // go over the stored trigger and select the unprescaled and fired one
    std::string name; bool smPass=false; bool sePass=false;

    for ( unsigned int i=0; i<smUsed.size(); i++) {
	name = smUsed[i];
	if (trigger->path(name)) { if(trigger->path(name)->wasAccept() && trigger->path(name)->prescale()==1) smPass=true;}
    }
    if(smPass) trig+=1<<0;

    for ( unsigned int i=0; i<seUsed.size(); i++) {
	name = seUsed[i];
	if (trigger->path(name)) { if(trigger->path(name)->wasAccept() && trigger->path(name)->prescale()==1) sePass=true;}
    }
    if(sePass) trig+=1<<1;


    ////////////
    // pileup //
    ////////////

    edm::Handle<double> rhoHandle;
    iEvent.getByLabel(rho_,rhoHandle);
    rho = *rhoHandle;

    // get true number of interaction only for MC for pileup reweighting
    if (!iEvent.isRealData()) {
	Handle<std::vector< PileupSummaryInfo > >  PupInfo; iEvent.getByLabel(pu_, PupInfo);
	npu=PupInfo->begin()->getTrueNumInteractions();

	// Write out the generator level particles. Since we try to figure out the mother-daughter relationshipts also,
	// the storeGenParticle() method has on O(n) complexity and this brings the total complexity of the loop to O(n^2). 
	// This might become a bottleneck if the number of gen. level particles becomes large (e.g. if no filtering is done). 
	// However, there does not seem to be a way around due to the way mother-daughter relationships are stored in the 
	// reco::GenParticle class (it is necessary to loop over a lookup table between the object pointers and the new indices).
	Handle<vector<reco::GenParticle> > genReco; iEvent.getByLabel(gen_,genReco);
	for(const reco::GenParticle& genparticle : *genReco) { storeGenParticle(genparticle); }
	if((int)genReco->size() != gnp) {
	    LogWarning("CreateNTuple") << "genReco->size() and gnp do not match (" << genReco->size() << " vs " << gnp << ")";
	}
    }

    ////////////
    // vertex //
    ////////////

    Handle<vector<reco::Vertex>> vtx; 
    iEvent.getByLabel(vtx_,vtx);
    if (!vtx.isValid()) return;
    // require in the event that there is at least one reconstructed vertex
    if(vtx->size()<=0) return;
    // pick the first (i.e. highest sum pt) vertex
    const reco::Vertex* theVertex=&(vtx->front());
    // require that the vertex meets certain criteria
    if(theVertex->ndof()<5) return;
    if(fabs(theVertex->z())>24.0) return;
    if(fabs(theVertex->position().rho())>2.0) return;
    // get the vertex value
    evx = theVertex->x();
    evy = theVertex->y();
    evz = theVertex->z();
    // now, count vertices
    std::vector<reco::Vertex>::const_iterator itv;
    int NVtx = 0;
    for (itv = vtx->begin(); itv != vtx->end(); ++itv) {
	// require that the vertex meets certain criteria
	if(itv->ndof()<5) continue;
	if(fabs(itv->z())>50.0) continue;
	if(fabs(itv->position().rho())>2.0) continue;
	++NVtx;
    }
    nvtx=NVtx;


    /////////
    // tau //
    /////////

    Handle<pat::TauCollection> tau; iEvent.getByLabel(tau_,tau);
    for ( auto &it : *tau ) {
	if(!( it.pt()>20 && abs(it.pdgId())==15 )) continue; 
	pid[np]     = it.pdgId();
	px[np]      = it.px();
	py[np]      = it.py();
	pz[np]      = it.pz();
	e[np]       = it.energy();
	vx[np]      = it.vx();
	vy[np]      = it.vy();
	vz[np]      = it.vz();
	if( it.tauID("decayModeFinding")==1)                             id[np]+= 1<<0; 
	if( it.tauID("byLooseCombinedIsolationDeltaBetaCorr")      ==1)  id[np]+= 1<<1;
	if( it.tauID("byMediumCombinedIsolationDeltaBetaCorr")     ==1)  id[np]+= 1<<2;
	if( it.tauID("byTightCombinedIsolationDeltaBetaCorr")      ==1)  id[np]+= 1<<3;
	if( it.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") ==1)  id[np]+= 1<<4;  
	if( it.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")==1)  id[np]+= 1<<5; 
	if( it.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits") ==1)  id[np]+= 1<<6; 
	if( it.tauID("byLooseIsolationMVA")   ==1)    id[np]+= 1<<7;
	if( it.tauID("byMediumIsolationMVA")  ==1)    id[np]+= 1<<8;
	if( it.tauID("byTightIsolationMVA")   ==1)    id[np]+= 1<<9;
	if( it.tauID("byLooseIsolationMVA2")  ==1)    id[np]+= 1<<10;
	if( it.tauID("byMediumIsolationMVA2") ==1)    id[np]+= 1<<11;
	if( it.tauID("byTightIsolationMVA2")  ==1)    id[np]+= 1<<12;
	if( it.tauID("againstElectronLooseMVA3") ==1) id[np]+= 1<<13; 
	if( it.tauID("againstElectronMediumMVA3")==1) id[np]+= 1<<14; 
	if( it.tauID("againstElectronTightMVA3") ==1) id[np]+= 1<<15;
	if( it.tauID("againstElectronVTightMVA3")==1) id[np]+= 1<<16;
	if( it.tauID("againstMuonLoose3")     ==1)    id[np]+= 1<<17;
	if( it.tauID("againstMuonTight3")     ==1)    id[np]+= 1<<18;
	if( it.genJet())                              id[np]+= 1<<19;
	np++;
    }


    ////////// 
    // muon //
    //////////
 
    Handle<pat::MuonCollection> muon; iEvent.getByLabel(muon_,muon);
    for ( auto &it : *muon) {
	if(!(it.isPFMuon() && it.pt()>10)) continue;
	pid[np]     = it.pdgId();
	px[np]      = it.px();
	py[np]      = it.py();
	pz[np]      = it.pz();
	e[np]       = it.energy();
	vx[np]      = it.vx();
	vy[np]      = it.vy();
	vz[np]      = it.vz();
	if(it.isLooseMuon() && fabs(it.eta())<2.4)           id[np]+= 1<<0; 
	if(it.isTightMuon(*theVertex) && fabs(it.eta())<2.1) id[np]+= 1<<1; 
	if(it.isLooseMuon() && fabs(it.eta())<2.5)           id[np]+= 1<<2; // veto
	if(it.genLepton())                                   id[np]+= 1<<3;
	iso[np] = (it.pfIsolationR04().sumChargedHadronPt + std::max(0., (it.pfIsolationR04().sumNeutralHadronEt 
			+ it.pfIsolationR04().sumPhotonEt - 0.5*it.pfIsolationR04().sumPUPt )))/it.pt(); 
	np++;
    }


    ///////// 
    // ele //
    /////////

    Handle<ElectronCollection> ele; iEvent.getByLabel(ele_,ele);
    for ( auto &it : *ele) {
	if(!(it.isPF() && it.pt()>10)) continue;
	pid[np]     = it.pdgId();
	px[np]      = it.px();
	py[np]      = it.py();
	pz[np]      = it.pz();
	e[np]       = it.energy();
	vx[np]      = it.vx();
	vy[np]      = it.vy();
	vz[np]      = it.vz();

	if( fabs(it.eta())<2.5 && it.gsfTrack()->dxy(vtx->begin()->position())<0.04 && 
		it.passConversionVeto()==true && it.gsfTrack()->trackerExpectedHitsInner().numberOfHits()<=0 ) id[np] += 1<<0; // L 
	if( fabs(it.eta())<2.5 && it.gsfTrack()->dxy(vtx->begin()->position())<0.02 && 
		it.passConversionVeto()==true && it.gsfTrack()->trackerExpectedHitsInner().numberOfHits()<=0 && 
		(!(fabs(it.superCluster()->eta())>1.4442 && fabs(it.superCluster()->eta())<1.5660)) )          id[np] += 1<<1; // T
	if( fabs(it.eta())<2.5 )               id[np] += 1<<2; // V
	if( it.genLepton())                    id[np] += 1<<3;

	disc[np] = it.electronID("mvaTrigV0"); 
	// eleIso: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PfIsolation
	float AEff03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, 
		it.superCluster()->eta(), ElectronEffectiveArea::kEleEAData2012);
	iso[np] = ( it.chargedHadronIso() + std::max(0., it.neutralHadronIso() + it.photonIso() - rho*AEff03 ))/it.pt();

	np++;
    }


    ////////// 
    // jets //
    //////////

    Handle<pat::JetCollection> jet; iEvent.getByLabel(jet_,jet);
    for ( auto &it : *jet ) {
	if(!(it.isPFJet() && it.pt()>20)) continue;
	pid[np]     = 81;
	px[np]      = it.px();
	py[np]      = it.py();
	pz[np]      = it.pz();
	e[np]       = it.energy();
	vx[np]      = it.vx();
	vy[np]      = it.vy();
	vz[np]      = it.vz();
	
	if ( (fabs(it.eta())<2.4 && it.neutralHadronEnergyFraction()<0.99 && it.neutralEmEnergyFraction()<0.99 &&
	      it.numberOfDaughters()>1 && it.chargedHadronEnergyFraction()>0 && it.chargedMultiplicity()>0 && 
              it.chargedEmEnergyFraction()<0.99) || ( fabs(it.eta())>2.4 && it.neutralHadronEnergyFraction()<0.99 && 
              it.neutralEmEnergyFraction()<0.99 && it.numberOfDaughters()>1) )    id[np] += 1<<0; // L

	if ( (fabs(it.eta())<2.4 && it.neutralHadronEnergyFraction()<0.95 && it.neutralEmEnergyFraction()<0.95 &&
	      it.numberOfDaughters()>1 && it.chargedHadronEnergyFraction()>0 && it.chargedMultiplicity()>0 &&
              it.chargedEmEnergyFraction()<0.99) || ( fabs(it.eta())>2.4 && it.neutralHadronEnergyFraction()<0.95 && 
              it.neutralEmEnergyFraction()<0.95 && it.numberOfDaughters()>1) )    id[np] += 1<<1; // M

	if ( (fabs(it.eta())<2.4 && it.neutralHadronEnergyFraction()<0.90 && it.neutralEmEnergyFraction()<0.90 &&
	      it.numberOfDaughters()>1 && it.chargedHadronEnergyFraction()>0 && it.chargedMultiplicity()>0 && 
              it.chargedEmEnergyFraction()<0.99) || ( fabs(it.eta())>2.4 && it.neutralHadronEnergyFraction()<0.90 && 
              it.neutralEmEnergyFraction()<0.90 && it.numberOfDaughters()>1) )    id[np] += 1<<2; // T

	if ( it.genJet()) id[np]+= 1<<3;

	disc[np]  = it.bDiscriminator("combinedSecondaryVertexBJetTags"); 

	np++;
    }
   
    /////////
    // MET //
    /////////

    Handle<pat::METCollection> met; 
    iEvent.getByLabel(met_,met);
    if( met->begin()->isPFMET()==true){
	metx=met->begin()->px();
	mety=met->begin()->py();
	metcov[0]= met->front().getSignificanceMatrix()(0,0);
	metcov[1]= met->front().getSignificanceMatrix()(1,0);
	metcov[2]= met->front().getSignificanceMatrix()(0,1);
	metcov[3]= met->front().getSignificanceMatrix()(1,1);
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
    int index = gnp++; // post incrementation!!
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
    nlep   = fs->make<TH1D>("nlep","Number of leptons",part,0,part);

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

