import FWCore.ParameterSet.Config as cms

process = cms.Process('Demo')

process.load('FWCore.MessageService.MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3) )

process.source = cms.Source('PoolSource',
	# replace 'myfile.root' with the source file you want to use
	fileNames = cms.untracked.vstring(
		'file:signal_pat.root'
	)
)

process.demo = cms.EDAnalyzer(
	'GenParticleDigraphCreator',
	collectionlabel=cms.InputTag('savedGenParticles'),
	dotfile=cms.string('output.dot')
)


process.p = cms.Path(process.demo)
