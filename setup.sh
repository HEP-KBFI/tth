#!/bin/bash

#abort on error
set -e

#print all lines
set -x

rm -Rf CMSSW/src
source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project -n CMSSW CMSSW CMSSW_5_3_16_patch1

cd CMSSW/src

#run cmsenv
echo "running cmsenv in  "`pwd`
eval `scramv1 runtime -sh`

echo "adding PAT"
#PAT from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATReleaseNotes52X#Add_CSCTightHaloFilter_CMSSW_5_3
git cms-addpkg PhysicsTools/PatAlgos

echo "adding ElectronTools"
#electron ID
git cms-addpkg EgammaAnalysis/ElectronTools
echo "getting electron MVAs "
cd EgammaAnalysis/ElectronTools/data/
cat download.url | xargs wget
cd $CMSSW_BASE/..

echo "running submodule init"
git checkout CMSSW/src
git submodule init
git submodule update

#submodules already added
##pull PU jet ID code
#git submodule add https://github.com/jpata/UserCode-CMG-CMGTools-External CMSSW/src/CMGTools/External
#
##pull generically useful analysis modules
#git submodule add https://github.com/jpata/AnalysisModules CMSSW/src/AnalysisModules
#
##pull common PAT code for TTH and singletop
#git submodule add https://github.com/jpata/tth-pfbreco CMSSW/src/UserCode/TTHPAT


echo "setup succeeded"
