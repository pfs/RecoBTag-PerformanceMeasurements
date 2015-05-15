----------------------------------------------
# RecoBTag-PerformanceMeasurements with ttbar

### Installation
```
cmsrel CMSSW_7_4_0
cd CMSSW_7_4_0/src
cmsenv
git cms-merge-topic ikrav:egm_id_74X_v0
git clone -b 7_4_X_dev git@github.com:pfs/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements
git clone -b V00-00-01 git://github.com/cms-btv-pog/cms-EventCounter.git MyAnalysis/EventCounter
git cms-merge-topic -u cms-btv-pog:TrackHistoryUpdate_from-CMSSW_7_4_0_pre7
scram b -j 8
cmsenv
```

### Producing the trees
```
cmsRun runBTagAnalyzer_cfg.py
```
Will run locally the analyzer for testing purposes
```
python submitToGrid.py -j ttbar_phys14.json -c ../runBTagAnalyzer_cfg.py -s
```
Will submit the samples described in the json to the grid.
Don't forget to init the environment for crab3
(e.g. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial)

### KIN method

### JP method

### FtM method