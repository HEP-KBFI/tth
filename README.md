tth
===

A concise TWiki describing this analysis is located in https://twiki.cern.ch/twiki/bin/viewauth/CMS/TTHTallinn.
If something fails to work, run the failing command again with `my_cmd &> log` and send the log to one of the devs (e.g. joosep.pata@cern.ch) putting as much information as you can.

Setup
=====
0. Be in a clean directory, with no CMSSW setup and no `cmsenv`. If you don't want to constantly write your password for git, add your keyfile to ssh-agent.
1. check out code with `git clone https://github.com/HEP-KBFI/tth.git`
2. run `./setup.sh` or if it doesn't work, run it manually line-by-line and see where it fails
3. make sure the last line printed out is `setup succeeded`, otherwise there were problems
4. run `make` to compile

Creating a PAT-tuple (step1)
============================
0. Make sure the code successfully compiled (see the section **Setup**)
1. create a sample tuple by running `cmsRun CMSSW/src/UserCode/TTHPAT/python/pfbreco_mc_cfg.py`

FWLite example code
=======================================
1. Code is in `CMSSW/src/UserCode/TTHPAT/bin/simple_loop.cc`
1. compile with `scram b UserCode/TTHPAT`. You should see in the output `>> Compiling  /home/joosep/test/tth/CMSSW/src/UserCode/TTHPAT/bin/simple_loop.cc`
2. run `$CMSSW_BASE/bin/$SCRAM_ARCH/simple_loop`, takes as input `input.root`. If the file does not exist compilation did not succeed!
