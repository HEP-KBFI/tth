source /cvmfs/cms.cern.ch/cmsset_default.sh
scram project -n CMSSW CMSSW CMSSW_5_3_16_patch1

#pull PU jet ID code
git submodule add https://github.com/jpata/UserCode-CMG-CMGTools-External CMSSW/src/CMGTools/External

#pull generically useful analysis modules
git submodule add https://github.com/jpata/AnalysisModules CMSSW/src/AnalysisModules

#pull common PAT code for TTH and singletop
git submodule add https://github.com/jpata/tth-pfbreco CMSSW/src/UserCode/TTHPAT

cd CMSSW/src
#PAT from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATReleaseNotes52X#Add_CSCTightHaloFilter_CMSSW_5_3
git cms-addpkg PhysicsTools/PatAlgos

#electron ID
git cms-addpkg EgammaAnalysis/ElectronTools
cd EgammaAnalysis/ElectronTools/data/
cat download.url | xargs wget
cd ../..


