// -*- C++ -*-
//
// Package:    TauAnalyserAOD
// Class:      TauAnalyserAOD
// 
/**\class TauAnalyserAOD TauAnalyserAOD.cc tthAnalysis/TauAnalyserAOD/src/TauAnalyserAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Betty Calpas
//         Created:  Mon May 12 14:22:40 EEST 2014
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
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"

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

class TauAnalyserAOD : public edm::EDAnalyzer {
   public:
      explicit TauAnalyserAOD(const edm::ParameterSet&);
      ~TauAnalyserAOD();

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
      TLorentzVector LTau, LgTau;

      TH1 *h0, *h1, *h2, *h3;

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
TauAnalyserAOD::TauAnalyserAOD(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


TauAnalyserAOD::~TauAnalyserAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TauAnalyserAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm; using namespace std;

  //tau
  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByLabel("hpsPFTauProducer", taus);
  //genPart
  edm::Handle<reco::GenParticleCollection> genp;
  iEvent.getByLabel("genParticles", genp);
  //tau disc
  edm::Handle<reco::PFTauDiscriminator> disc;
  iEvent.getByLabel("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr", disc);

  //recoTau
  //for (reco::PFTauCollection::const_iterator it = taus->begin(); it != taus->end(); it++) {
   for ( unsigned iTau = 0; iTau < taus->size(); ++iTau ) {
        reco::PFTauRef it(taus, iTau);
  
         if(!(abs(it->pdgId())==15 && it->pt()>20 && abs(it->eta())<2.3  )) continue; 
          LTau.SetPxPyPzE(it->px(), it->py(), it->pz(), it->energy());
          //cout << it->pdgId()<<", "<<it->px() <<", "<<it->py()<<", "<<it->pz()<<", "<<it->energy()<<endl;

          //GenTau
          bool recoMatchGen=false;
          for (reco::GenParticleCollection::const_iterator ig = genp->begin(); ig != genp->end(); ig++) {
                  if( !(abs(ig->pdgId())==15) ) continue;
                  LgTau.SetPxPyPzE(ig->px(), ig->py(), ig->pz(), ig->energy());
                  //cout << ig->pdgId()<<", "<<ig->px() <<", "<<ig->py()<<", "<<ig->pz()<<", "<<ig->energy()<<endl;

                  //Loop on tau daugter and remove those which have a el or mu
                  int n = ig->numberOfDaughters();
                  bool noHadTau=false;
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
          if( (*disc)[it]==1) h1->Fill(it->pt());
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
TauAnalyserAOD::beginJob()
{
edm::Service<TFileService> fs;
h0 = fs->make<TH1D>("h0","h0",100,0,300);
h1 = fs->make<TH1D>("h1","h1",100,0,300);
//h2 = fs->make<TH1D>("h2","h2",100,0,300);
//h3 = fs->make<TH1D>("h3","h3",100,0,300);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TauAnalyserAOD::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
TauAnalyserAOD::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
TauAnalyserAOD::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
TauAnalyserAOD::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TauAnalyserAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauAnalyserAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  //edm::ParameterSetDescription desc;
  //desc.setUnknown();
  //descriptions.addDefault(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
//define this as a plug-in
DEFINE_FWK_MODULE(TauAnalyserAOD);
