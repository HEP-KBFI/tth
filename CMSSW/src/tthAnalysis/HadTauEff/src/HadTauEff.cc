// -*- C++ -*-
//
// Package:    HadTauEff
// Class:      HadTauEff
// 
/**\class HadTauEff HadTauEff.cc tthAnalysis/HadTauEff/src/HadTauEff.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Betty Calpas
//         Created:  Thu May  8 18:11:48 EEST 2014
// $Id$
//
//


// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/Handle.h"

#include <TH1.h>
#include <string>
#include <TMath.h>
#include <TLorentzVector.h>

//
// class declaration
//

class HadTauEff : public edm::EDAnalyzer {
   public:
      explicit HadTauEff(const edm::ParameterSet&);
      ~HadTauEff();

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

      edm::InputTag src_;

      bool requireGenTauMatch_;
      std::string looseIsoDB3Hits_;
      std::string mediumIsoDB3Hits_;
      std::string tightIsoDB3Hits_;

      TLorentzVector LTau, LgTau;

      bool recoMatchGen, noHadTau;

      TH1 *h0, *h1, *h2, *h3;

};

//
// constants, enums and typedefs
//
const reco::GenParticle* getGenTau(const pat::Tau& patTau)
 {
   std::vector<reco::GenParticleRef> associatedGenParticles = patTau.genParticleRefs();
   for ( std::vector<reco::GenParticleRef>::const_iterator it = associatedGenParticles.begin();
         it != associatedGenParticles.end(); ++it ) {
     if ( it->isAvailable() ) {
       const reco::GenParticleRef& genParticle = (*it);
       if ( genParticle->pdgId() == -15 || genParticle->pdgId() == +15 ) return genParticle.get();
     }
   }
   return 0;
 }

//
// static data member definitions
//

//
// constructors and destructor
//
HadTauEff::HadTauEff(const edm::ParameterSet& iConfig)
{
	//now do what ever initialization is needed

	src_ = iConfig.getParameter<edm::InputTag>("src");
	requireGenTauMatch_ = iConfig.getParameter<bool>("requireGenTauMatch");
        looseIsoDB3Hits_  = iConfig.getParameter<std::string>("looseIsoDB3Hits");
        mediumIsoDB3Hits_ = iConfig.getParameter<std::string>("mediumIsoDB3Hits");
        tightIsoDB3Hits_  = iConfig.getParameter<std::string>("tightIsoDB3Hits");
}


HadTauEff::~HadTauEff()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HadTauEff::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

using namespace std;

  //tau
  edm::Handle<pat::TauCollection> PatTaus;
  iEvent.getByLabel(src_, PatTaus);
  //genPart
  edm::Handle<std::vector<reco::GenParticle> > GenPart;
  iEvent.getByLabel("savedGenParticles", GenPart);

  //recoTau
  for (pat::TauCollection::const_iterator it = PatTaus->begin(); it != PatTaus->end(); it++) {
	  if(!(it->isPFTau() && abs(it->pdgId())==15 && it->pt()>20 && abs(it->eta())<2.3  )) continue; 
	  LTau.SetPxPyPzE(it->px(), it->py(), it->pz(), it->energy());
	  //cout <<"ttau"<<endl;
	  //cout << it->pdgId()<<", "<<it->px() <<", "<<it->py()<<", "<<it->pz()<<", "<<it->energy()<<endl;

          //GenTau
	  recoMatchGen=false;
	  for (std::vector<reco::GenParticle>::const_iterator ig = GenPart->begin(); ig != GenPart->end(); ig++) {
		  if( !(abs(ig->pdgId())==15) ) continue;
		  LgTau.SetPxPyPzE(ig->px(), ig->py(), ig->pz(), ig->energy());
		  //cout <<"gtau"<<endl;
		  //cout << ig->pdgId()<<", "<<ig->px() <<", "<<ig->py()<<", "<<ig->pz()<<", "<<ig->energy()<<endl;

		  //Loop on tau daugter and remove those which have a el or mu
		  int n = ig->numberOfDaughters();
		  noHadTau=false;
		  //GenTauDauter

		  for(int j = 0; j <n; ++j) {
			  const reco::Candidate * d = ig->daughter( j );
			  if ( abs(d->pdgId())==11 || abs(d->pdgId()==13) || abs(d->pdgId())==15) noHadTau=true;
		  }
		  if (noHadTau) continue;
                  //Matching

		  if(LgTau.DeltaR(LTau)<0.3) recoMatchGen=true;

	  }
	  if(!recoMatchGen) continue;

	  h0->Fill(it->pt());
	  if( it->tauID(looseIsoDB3Hits_) ==1) h1->Fill(it->pt());
	  if( it->tauID(mediumIsoDB3Hits_)==1) h2->Fill(it->pt());
	  if( it->tauID(tightIsoDB3Hits_) ==1) h3->Fill(it->pt());

  }


}

// ------------ method called once each job just before starting event loop  ------------
void 
HadTauEff::beginJob()
{
edm::Service<TFileService> fs;
h0 = fs->make<TH1D>("h0","h0",100,0,300);
h1 = fs->make<TH1D>("h1","h1",100,0,300);
h2 = fs->make<TH1D>("h2","h2",100,0,300);
h3 = fs->make<TH1D>("h3","h3",100,0,300);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HadTauEff::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
HadTauEff::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
HadTauEff::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HadTauEff::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HadTauEff::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HadTauEff::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  //edm::ParameterSetDescription desc;
  //desc.setUnknown();
  //descriptions.addDefault(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(HadTauEff);
