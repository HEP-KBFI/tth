import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(25) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:signal_pat.root'
        #'file:/hdfs/cms/store/user/joosep/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Mar20/87131f7abf095b2014f0aafe9c07da3d/output_100_1_i4g.root'
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
tauMinPt     = cms.double(20), #tauID recommandation
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

#Do not work
#process.load("Demo.CreateNTuple.createntuple_cfi")
#process.demo.minTracks=1000

#process.demo = cms.EDAnalyzer('MyTrackAnalyzer',
#         tracks = cms.untracked.InputTag('generalTracks')
#   )

#To see event content
#process.dump=cms.EDAnalyzer('EventContentAnalyzer')
#process.p = cms.Path(process.demo*process.dump)
#process.Tracer = cms.Service("Tracer") #track the module that worked

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('ntuple.root')
                                   )

process.p = cms.Path(process.demo)
