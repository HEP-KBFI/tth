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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#Joseep Files
"/store/user/joosep/TTH_HToTauTau_M-125_8TeV_pythia6/Mar19/87131f7abf095b2014f0aafe9c07da3d/output_100_1_fDa.root",
#"/store/user/joosep/TTH_HToTauTau_M-125_8TeV_pythia6/Mar19/87131f7abf095b2014f0aafe9c07da3d/output_101_1_ZKT.root",
#"/store/user/joosep/TTH_HToTauTau_M-125_8TeV_pythia6/Mar19/87131f7abf095b2014f0aafe9c07da3d/output_102_1_5YS.root",
#"/store/user/joosep/TTH_HToTauTau_M-125_8TeV_pythia6/Mar19/87131f7abf095b2014f0aafe9c07da3d/output_103_1_S0X.root",
#AOD Files
#'/store/mc/Summer12_DR53X/TTH_HToTauTau_M-125_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/00918FE3-2FFC-E111-BFEB-00266CFFC980.root',
#'/store/mc/Summer12_DR53X/TTH_HToTauTau_M-125_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/00B37956-1FFC-E111-AFE1-0017A4771010.root',
#'/store/mc/Summer12_DR53X/TTH_HToTauTau_M-125_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/0205685C-4FFC-E111-9AF5-00266CFFCB7C.root',
#'/store/mc/Summer12_DR53X/TTH_HToTauTau_M-125_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/028FFE49-B4FC-E111-B61B-002481A7329C.root'
    )
)

process.demo = cms.EDAnalyzer('HadTauEff',
		src = cms.InputTag('selectedPatTaus'),
		requireGenTauMatch = cms.bool(True),
                looseIsoDB3Hits  = cms.string('byLooseCombinedIsolationDeltaBetaCorr3Hits'),
                mediumIsoDB3Hits = cms.string('byMediumCombinedIsolationDeltaBetaCorr3Hits'),
                tightIsoDB3Hits  = cms.string('byTightCombinedIsolationDeltaBetaCorr3Hits'),
)

process.TFileService = cms.Service("TFileService", 
    fileName = cms.string('skim_TTH_HToTauTau_M-125.root')
)

process.p = cms.Path( process.demo )
