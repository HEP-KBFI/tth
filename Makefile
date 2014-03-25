build:
	cd CMSSW;scram b -j20

test: test_step1

test_step1: test_step1_mc_signal test_step1_mc_ttbar test_step1_mc_signal_tautau

STEP1_MC=cmsRun CMSSW/src/UserCode/TTHPAT/python/pfbreco_mc_cfg.py maxEvents=1000

test_step1_mc_signal:
	$(STEP1_MC) outputFile=step1_signal.root inputFiles=/store/mc/Summer12_DR53X/TTH_HToBB_M-125_8TeV-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/00FA9388-81FC-E111-A80D-00215E2217BE.root

test_step1_mc_signal_tautau:
	$(STEP1_MC) outputFile=step1_signal_tautau.root inputFiles=/store/mc/Summer12_DR53X/TTH_HToTauTau_M-125_8TeV_pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/00918FE3-2FFC-E111-BFEB-00266CFFC980.root

test_step1_mc_ttbar:
	$(STEP1_MC) outputFile=step1_ttbar.root inputFiles=/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v2/0000/0026AFAA-41F1-E111-80E1-00266CF2E2C8.root
