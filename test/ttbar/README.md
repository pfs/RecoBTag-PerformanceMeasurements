----------------------------------------------
# RecoBTag-PerformanceMeasurements with ttbar

### Installation
```
cmsrel CMSSW_7_4_0
cd CMSSW_7_4_0/src
cmsenv
git cms-merge-topic ikrav:egm_id_74X_v0
git cms-merge-topic -u cms-btv-pog:TrackHistoryUpdate_from-CMSSW_7_4_0_pre7
git clone -b V00-00-01 git://github.com/cms-btv-pog/cms-EventCounter.git MyAnalysis/EventCounter
git clone -b 7_4_X_dev git@github.com:pfs/RecoBTag-PerformanceMeasurements.git RecoBTag/PerformanceMeasurements
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
python checkProductionIntegrity.py -i /store/group/phys_top/psilva/BTV/d214360
```
Can run as soon as ntuple production starts to end, to move from crab output directories to a more simple directory structure
which can be easily parsed by the local analysis. 
Once production is finished you can call with --cleanup to remove original crab directories in EOS.

### Running local analysis
```
python runTTbarAnalysis.py -i /store/group/phys_top/psilva/BTV/d214360 -j ttbar_phys14.json -t taggers_phys14.json -n 8
```
Once grid jobs are run, and ntuples are stored in a given directory, you can run the local analysis to produce the slimmed ntuples for the efficiency measurement.
MC will be weighted by cross section. The number after -n indicates how many threads should be used.
NB The first time the script will produce a pickle file with the weights to be used according to the number of files found, xsec specified in the json file and 
integrated luminosity.
It prints out s.th. like "Produced normalization cache (analysis/.xsecweights.pck)"
In case you update the trees, xsec or lumi you have to remove by hand the pickle file.
```
python plotter.py -i analysis/ -j ttbar_phys14.json
```
Makes control plots and stores all in a ROOT file. Different options may be passed to filter plots, and show differently the plots. 

### Performance analyses

#### Templated methods
```
sh KIN_runClassifier.sh
```
After running the local analysis use the kin tree stored in the ttbar sample to train a kinematics discriminator for b-jets in ttbar events.
The script compiles and runs KIN_trainClassifier.C which should be modified in case different trainings are required.
```
python runTTbarAnalysis.py -i /store/group/phys_top/psilva/BTV/d214360 -j ttbar_phys14.json -t taggers_phys14.json --tmvaWgts analysis/KIN_weights/TMVAClassification_BDT.weights.xml -n 8
```
Re-run the analysis to store the KIN discriminator value per jet
```
python Templated_btagEffFitter.py -i analysis/ -o analysis/ -t taggers_phys14.json -n 8
```
Runs the fits to the templates to determine the scale factors. Valid for KIN, Mlj, JP, others one may wish to add.
The base procedure is similar for all. The first time to run will take a long time as templates need to be created.
If templates are stable and only fits need to be redone when can run with the option "--recycleTemplates"
```
python Templated_btagEffFitter.py --show analysis//close_mlj_templates/.csvivf_fits.pck,analysis//kindisc_templates/.csvivf_fits.pck,analysis//jpTagger_templates/.csvivf_fits.pck
```
Fit results are stored in pickle files which can be used to produce the final summary plots.

#### FtM method