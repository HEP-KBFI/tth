tth
===

A concise TWiki describing this analysis is located in https://twiki.cern.ch/twiki/bin/view/Main/TTHTallinn

Setup
=====
0. Be in a clean directory, with no CMSSW setup and no `cmsenv`. If you don't want to constantly write your password for git, add your keyfile to ssh-agent.
1. check out code with `git clone https://github.com/HEP-KBFI/tth.git`
2. run `./setup.sh` or if it doesn't work, run it manually line-by-line and see where it fails
3. `make`


Simple FWLite event loop example over PAT-tuple
=======================================
1. Code is in `CMSSW/src/UserCode/TTHPAT/bin/simple_loop.cc`
1. compile with `scram b UserCode/TTHPAT`
2. run `$CMSSW_BASE/bin/$SCRAM_ARCH/simple_loop`, takes as input `input.root`
