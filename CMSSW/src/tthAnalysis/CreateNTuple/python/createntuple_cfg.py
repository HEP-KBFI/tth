import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.demo = cms.EDAnalyzer('CreateNTuple',
		gen     = cms.InputTag("savedGenParticles"),
		tau     = cms.InputTag("selectedPatTaus"),
		muon    = cms.InputTag("selectedPatMuons"),
		elec    = cms.InputTag("selectedPatElectrons"),
		jet     = cms.InputTag("selectedPatJets"),
		#met     = cms.InputTag("patMETs"),
		met     = cms.InputTag("patType1CorrectedPFMet"),
		vtx     = cms.InputTag("goodOfflinePrimaryVertices"),
		rho     = cms.InputTag("kt6PFJets:rho"),
		pu      = cms.InputTag("addPileupInfo"),
		trig    = cms.InputTag("patTriggerEvent"),
	        smtrig  = cms.string("HLT_IsoMu24_eta2p1_v"), 
	        setrig  = cms.string("HLT_Ele27_WP80_v"), 
		channel = cms.string("sl")
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('output.root'))

process.p = cms.Path(process.demo)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source ("PoolSource", fileNames=cms.untracked.vstring(
'file:/hdfs/cms/store/user/calpas/TTJets_SemiLeptMGDecays_8TeV-madgraph/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1/4ff6f8bbab831f313c9802a1c3a04598/output_1110_3_nCm.root',
  ),
)

