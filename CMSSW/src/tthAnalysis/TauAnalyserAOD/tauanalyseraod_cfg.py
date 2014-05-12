import FWCore.ParameterSet.Config as cms 
import copy

process = cms.Process("Demo")

process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] ) 
#process.GlobalTag.globaltag = cms.string('STARTUP31X_V1::All')
 
process.load('PhysicsTools.PatAlgos.patSequences_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#AOD Files
'/store/mc/Summer12_DR53X/TTH_HToTauTau_M-125_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/00918FE3-2FFC-E111-BFEB-00266CFFC980.root',
'/store/mc/Summer12_DR53X/TTH_HToTauTau_M-125_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/00B37956-1FFC-E111-AFE1-0017A4771010.root',
'/store/mc/Summer12_DR53X/TTH_HToTauTau_M-125_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/0205685C-4FFC-E111-9AF5-00266CFFCB7C.root',
'/store/mc/Summer12_DR53X/TTH_HToTauTau_M-125_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/028FFE49-B4FC-E111-B61B-002481A7329C.root'
    )
)

process.demo = cms.EDAnalyzer('TauAnalyserAOD'
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string('skim_TTH_HToTauTau_M-125.root')
)

process.p = cms.Path(process.demo)
