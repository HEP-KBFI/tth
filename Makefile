build:
	cd CMSSW;scram b -j20

test_step1:
	cmsRun CMSSW/src/UserCode/TTHBBPAT/python/pfbreco_mc_cfg.py
