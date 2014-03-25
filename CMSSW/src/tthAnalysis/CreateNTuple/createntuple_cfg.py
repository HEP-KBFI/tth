import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:signal_pat.root'
    )
)

process.demo = cms.EDAnalyzer('CreateNTuple',
genLabel     = cms.InputTag("savedGenParticles"),
tauLabel     = cms.InputTag("selectedPatTaus"),
muonLabel    = cms.InputTag("selectedPatMuons"),
elecLabel    = cms.InputTag("selectedPatElectrons"),
jetLabel     = cms.InputTag("selectedPatJets"),
metLabel     = cms.InputTag("patMETs"),
vtxLabel     = cms.InputTag("goodOfflinePrimaryVertices"),
minTracks    = cms.untracked.uint32(0),
tauMinPt     = cms.double(20), 
tauMinEta    = cms.double(0),  #should remove 1.460 < |eta| < 1.558 if to much e- contamination in cleanPatTau
tauMaxEta    = cms.double(6),
muonMinPt    = cms.double(0),
muonMinEta   = cms.double(0),
muonMaxEta   = cms.double(6),
elecMinPt    = cms.double(0),
elecMinEta   = cms.double(0),
elecMaxEta   = cms.double(6),
jetMinPt     = cms.double(0), 
jetMinEta    = cms.double(0),
jetMaxEta    = cms.double(6),
metMinPt     = cms.double(10)
)

process.TFileService = cms.Service("TFileService",
                                         fileName = cms.string('ntuple.root')
                                    )

process.p = cms.Path(process.demo)
