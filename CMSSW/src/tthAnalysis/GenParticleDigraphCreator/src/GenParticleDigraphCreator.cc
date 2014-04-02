// -*- C++ -*-
//
// Package:    GenParticleDigraphCreator
// Class:      GenParticleDigraphCreator
//
/**\class GenParticleDigraphCreator GenParticleDigraphCreator.cc tthAnalysis/GenParticleDigraphCreator/src/GenParticleDigraphCreator.cc

	Description:  Creates a dot file with graphs of gen. level particles.
*/

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <fstream>
#include <memory>

class GenParticleDigraphCreator : public edm::EDAnalyzer {
	public:
		explicit GenParticleDigraphCreator(const edm::ParameterSet&);
		~GenParticleDigraphCreator();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

		static size_t findIndex(const reco::Candidate& p, const reco::Candidate * candidates[], const size_t size);

		// ----------member data ---------------------------
		edm::InputTag collectionlabel;
		std::ofstream dotfile;
};

GenParticleDigraphCreator::GenParticleDigraphCreator(const edm::ParameterSet& iConfig)
: collectionlabel(iConfig.getParameter<edm::InputTag>("collectionlabel")),
  dotfile(iConfig.getParameter<std::string>("dotfile"))
{
}


GenParticleDigraphCreator::~GenParticleDigraphCreator()
{
	dotfile.close();
}

void
GenParticleDigraphCreator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;

	Handle< std::vector<reco::GenParticle> > genparticles;
	iEvent.getByLabel(collectionlabel, genparticles);

	const reco::Candidate * candidates[genparticles->size()];

	dotfile << "# =========================================" << std::endl;
	dotfile << "digraph event {" << std::endl;

	size_t index = 0;
	for(const reco::GenParticle& p : *genparticles) {
		LogDebug("GenParticleDigraphCreator") << ".& " << &p << " -> "  << index << std::endl;

		dotfile << '\t' << index << " [label = \"idx: " << index
		        << "\\npdg: " << p.pdgId()
		        << ", st: " << p.status()
		        << "\\n" << p.energy() << " GeV\""
		        << (p.status()==3 ? ", style=filled, fillcolor=orange" : "")
		        << "]" << std::endl;

		candidates[index++] = &p;
	}
	if(genparticles->size() != index) {
		LogWarning("GenParticleDigraphCreator") << "genparticles->size() and index do not match (" << genparticles->size() << " vs " << index << ")";
	}

	// find the mother-daughter relationships and print the edges
	const size_t ncandidates = index;
	for(const reco::GenParticle& p : *genparticles) {
		const size_t idx_p = findIndex(p, candidates, ncandidates);
		for(const reco::Candidate& d : p) {
			const size_t idx_d = findIndex(d, candidates, ncandidates);
			if(idx_d == ncandidates) {
				LogWarning("GenParticleDigraphCreator") << "no such daughter: " << &d;
			}
			dotfile << '\t' << idx_p << " -> " << idx_d << std::endl;
		}
	}

	dotfile << "}" << std::endl;
}

size_t
GenParticleDigraphCreator::findIndex(const reco::Candidate& p, const reco::Candidate * candidates[], const size_t size)
{
	for(size_t i=0; i<size; i++) {
		if(candidates[i] == &p)
			return i;
	}
	return size;
}

// ------------ method called once each job just before starting event loop  ------------
void
GenParticleDigraphCreator::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
GenParticleDigraphCreator::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
GenParticleDigraphCreator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
GenParticleDigraphCreator::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
GenParticleDigraphCreator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
GenParticleDigraphCreator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenParticleDigraphCreator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenParticleDigraphCreator);
