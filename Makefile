build:
	cd CMSSW;scram b -j20

test: test_step1

test_step1: test_step1_mc_signal

test_step1_mc_signal:
	cmsRun CMSSW/src/UserCode/TTHBBPAT/python/pfbreco_mc_cfg.py maxEvents=100
