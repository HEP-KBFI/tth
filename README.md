tth
===

A concise TWiki describing this analysis is located in https://twiki.cern.ch/twiki/bin/view/Main/TTHTallinn

Setup
=====
0. be in a clean directory, with no CMSSW and no `cmsenv`
1. check out code with `git clone https://github.com/HEP-KBFI/tth.git`
2. run `./setup.sh` or if it doesn't work, run it manually line-by-line
3. `make`


Simple FWLite event loop example over PAT-tuple
=======================================
1. Code is in `CMSSW/src/UserCode/TTHPAT/bin/simple_loop.cc`
1. compile with `scram b`
2. run `$CMSSW_BASE/bin/$SCRAM_ARCH/simple_loop`, takes as input `input.root`
