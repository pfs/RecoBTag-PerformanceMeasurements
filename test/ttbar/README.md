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
python submitToGrid.py -j ttbar_phys14.json -c ${CMSSW_BASE}/src/RecoBTag/PerformanceMeasurements/test/runBTagAnalyzer_cfg.py -s
```
Will submit the samples described in the json to the grid.
Partial submission can be made with -o csv_list
Don't forget to init the environment for crab3
(e.g. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial)
```
python checkProductionIntegrity.py -i /store/group/phys_btag/performance/TTbar/229a310
```
Can run as soon as ntuple production starts to end, to move from crab output directories to a more simple directory structure
which can be easily parsed by the local analysis

### Running local analysis
```
python runTTbarAnalysis.py -i /store/group/phys_btag/performance/TTbar/229a310 -j ttbar_phys14.json -n 8
```
Once grid jobs are run, and ntuples are stored in a given directory, you can run the local analysis to produce the slimmed ntuples for the efficiency measurement.
MC will be weighted by cross section. The number after -n indicates how many threads should be used.
```
python plotter.py -i analysis/ -j ttbar_phys14.json
```
Makes control plots and stores all in a ROOT file. Different options may be passed to filter plots, and show differently the plots. 

### Performance analyses

#### KIN method

#### JP method

#### FtM method