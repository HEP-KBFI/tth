source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project -n CMSSW CMSSW CMSSW_5_3_16_patch1

#pull PU jet ID code
git submodule add https://github.com/jpata/UserCode-CMG-CMGTools-External CMSSW/src/CMGTools/External

#pull common PAT code for TTH and singletop
git submodule add https://github.com/jpata/tth-pfbreco CMSSW/src/UserCode/TTHBBPAT

#install electron ID
cd CMSSW/src
git cms-addpkg EgammaAnalysis/ElectronTools
cd EgammaAnalysis/ElectronTools/data/
cat download.url | xargs wget
cd ../..


