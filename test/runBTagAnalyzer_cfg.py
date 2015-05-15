
import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
import copy

###############################
####### Parameters ############
###############################

options = VarParsing ('python')

options.register('runOnData', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('outFilename', 'JetTree',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)
options.register('reportEvery', 1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1)"
)
options.register('wantSummary', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out trigger and timing summary"
)
options.register('usePFchs', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use PFchs"
)
#options.register('mcGlobalTag', 'PHYS14_25_V1',
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "MC global tag"
#)
#options.register('dataGlobalTag', 'GR_R_70_V2',
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "Data global tag"
#)
options.register('runSubJets', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run subjets"
)
options.register('processStdAK4Jets', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Process standard AK4 jets"
)
options.register('producePtRelTemplate', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Produce PtRel template"
)
options.register('fatJetPtMin', 150.0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "Minimum pT for fat jets (default is 150 GeV)"
)
options.register('useTTbarFilter', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use TTbar filter"
)
options.register('useTopProjections', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use top projections"
)
options.register('miniAOD', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Running on miniAOD"
)
options.register('useLegacyTaggers', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use legacy taggers"
)
options.register('useExplicitJTA', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Use explicit jet-track association"
)

## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', -1)

options.parseArguments()

print "Running on data: %s"%('True' if options.runOnData else 'False')
print "Running on miniAOD: %s"%('True' if options.miniAOD else 'False')
print "Using PFchs: %s"%('True' if options.usePFchs else 'False')

## Global tag
#globalTag = options.mcGlobalTag
#if options.runOnData:
#    globalTag = options.dataGlobalTag

## Jet energy corrections
jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
jetCorrectionsAK8 = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

if not options.usePFchs:
    jetCorrectionsAK4 = ('AK4PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
    jetCorrectionsAK8 = ('AK8PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

if options.runOnData:
    jetCorrectionsAK4[1].append('L2L3Residual')
    jetCorrectionsAK8[1].append('L2L3Residual')

## b-tag infos
bTagInfosLegacy = [
'impactParameterTagInfos'
,'secondaryVertexTagInfos'
,'inclusiveSecondaryVertexFinderTagInfos'
,'softPFMuonsTagInfos'
,'softPFElectronsTagInfos'
]
bTagInfos = [
'pfImpactParameterTagInfos'
,'pfSecondaryVertexTagInfos'
,'pfInclusiveSecondaryVertexFinderTagInfos'
,'softPFMuonsTagInfos'
,'softPFElectronsTagInfos'
]
## b-tag discriminators
bTagDiscriminatorsLegacy = [
'jetBProbabilityBJetTags'
,'jetProbabilityBJetTags'
,'positiveOnlyJetBProbabilityBJetTags'
,'positiveOnlyJetProbabilityBJetTags'
,'negativeOnlyJetBProbabilityBJetTags'
,'negativeOnlyJetProbabilityBJetTags'
,'trackCountingHighPurBJetTags'
,'trackCountingHighEffBJetTags'
,'negativeTrackCountingHighEffBJetTags'
,'negativeTrackCountingHighPurBJetTags'
,'simpleSecondaryVertexHighEffBJetTags'
,'simpleSecondaryVertexHighPurBJetTags'
,'negativeSimpleSecondaryVertexHighEffBJetTags'
,'negativeSimpleSecondaryVertexHighPurBJetTags'
,'combinedSecondaryVertexV2BJetTags'
,'positiveCombinedSecondaryVertexV2BJetTags'
,'negativeCombinedSecondaryVertexV2BJetTags'
,'combinedInclusiveSecondaryVertexV2BJetTags'
,'positiveCombinedInclusiveSecondaryVertexV2BJetTags'
,'negativeCombinedInclusiveSecondaryVertexV2BJetTags'
,'softPFMuonBJetTags'
,'positiveSoftPFMuonBJetTags'
,'negativeSoftPFMuonBJetTags'
,'softPFElectronBJetTags'
,'positiveSoftPFElectronBJetTags'
,'negativeSoftPFElectronBJetTags'
]
bTagDiscriminators = [
'pfJetBProbabilityBJetTags'
,'pfJetProbabilityBJetTags'
,'pfPositiveOnlyJetBProbabilityBJetTags'
,'pfPositiveOnlyJetProbabilityBJetTags'
,'pfNegativeOnlyJetBProbabilityBJetTags'
,'pfNegativeOnlyJetProbabilityBJetTags'
,'pfTrackCountingHighPurBJetTags'
,'pfTrackCountingHighEffBJetTags'
,'pfNegativeTrackCountingHighPurBJetTags'
,'pfNegativeTrackCountingHighEffBJetTags'
,'pfSimpleSecondaryVertexHighEffBJetTags'
,'pfSimpleSecondaryVertexHighPurBJetTags'
,'pfNegativeSimpleSecondaryVertexHighEffBJetTags'
,'pfNegativeSimpleSecondaryVertexHighPurBJetTags'
,'pfCombinedSecondaryVertexV2BJetTags'
,'pfPositiveCombinedSecondaryVertexV2BJetTags'
,'pfNegativeCombinedSecondaryVertexV2BJetTags'
,'pfCombinedInclusiveSecondaryVertexV2BJetTags'
,'pfPositiveCombinedInclusiveSecondaryVertexV2BJetTags'
,'pfNegativeCombinedInclusiveSecondaryVertexV2BJetTags'
,'softPFMuonBJetTags'
,'positiveSoftPFMuonBJetTags'
,'negativeSoftPFMuonBJetTags'
,'softPFElectronBJetTags'
,'positiveSoftPFElectronBJetTags'
,'negativeSoftPFElectronBJetTags'
]

## Legacy taggers not supported with MiniAOD
if options.miniAOD and options.useLegacyTaggers:
    print "WARNING: Legacy taggers not supported with MiniAOD"
    options.useLegacyTaggers = False
## If using legacy taggers
if options.useLegacyTaggers:
    bTagInfos = bTagInfosLegacy
    bTagDiscriminators = bTagDiscriminatorsLegacy

## Postfix
postfix = "PFlow"
## Various collection names
genParticles = 'genParticles'
jetSource = 'pfJetsPFBRECO'+postfix
genJetCollection = 'ak4GenJetsNoNu'+postfix
pfCandidates = 'particleFlow'
pvSource = 'offlinePrimaryVertices'
svSource = 'inclusiveCandidateSecondaryVertices'
muSource = 'muons'
elSource = 'gedGsfElectrons'
patMuons = 'selectedPatMuons'
## If running on miniAOD
if options.miniAOD:
    genParticles = 'prunedGenParticles'
    jetSource = 'ak4PFJets'
    genJetCollection = 'ak4GenJetsNoNu'
    pfCandidates = 'packedPFCandidates'
    pvSource = 'offlineSlimmedPrimaryVertices'
    svSource = 'slimmedSecondaryVertices'
    muSource = 'slimmedMuons'
    elSource = 'slimmedElectrons'
    patMuons = 'selectedMuons'

process = cms.Process("BTagAna")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10

## Input files
process.source = cms.Source(
"PoolSource",
fileNames = cms.untracked.vstring(
  '/store/relval/CMSSW_7_4_0_pre8//RelValTTbar_13/MINIAODSIM/PU50ns_MCRUN2_74_V6-v1/00000/90A385BF-60BD-E411-85D8-0025905938A4.root',
  '/store/relval/CMSSW_7_4_0_pre8//RelValTTbar_13/MINIAODSIM/PU50ns_MCRUN2_74_V6-v1/00000/A8E0D9B1-60BD-E411-BFC1-0025905A6080.root'
#'file:/nfs/dust/cms/user/marchesi/74Xdev/testfile/22E552FD-23B7-E411-B680-002618943911.root'
)
)
if options.miniAOD:
    process.source.fileNames = [
        '/store/relval/CMSSW_7_4_0_pre8/RelValZpTT_1500_13TeV/MINIAODSIM/MCRUN2_74_V7-v1/00000/9008F5B0-54BD-E411-96FB-0025905A6110.root'
        #'file:/nfs/dust/cms/user/marchesi/74Xdev/testfile/B62A3865-39B7-E411-B76A-002618943880.root'
        ]
if options.runOnData:
    process.source.fileNames = [
        '/store/relval/CMSSW_7_4_0_pre7/SingleMu/RECO/GR_R_74_V8A_RelVal_mu2012D-v1/00000/004E151D-D8B6-E411-A889-0025905B859E.root'
        ]
    
if options.runOnData :
    if options.runSubJets :
        options.outFilename += '_data_subjets.root'
    else :
        options.outFilename += '_data.root'
else :
    if options.runSubJets :
        options.outFilename += '_mc_subjets.root'
    else :
        options.outFilename += '_mc.root'

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string(options.outFilename)
)

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Options and Output Report
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
    allowUnscheduled = cms.untracked.bool(True)
)

#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = globalTag + '::All'
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_' + ('data' if options.runOnData else 'mc'))

##############################################
# Add GenParticlePruner for Boosted b-Tagging Studies
##############################################
process.prunedGenParticlesBoost = cms.EDProducer('GenParticlePruner',
    src = cms.InputTag(genParticles),
    select = cms.vstring(
    "drop  *  ", #by default
    "keep ( status = 3 || (status>=21 && status<=29) )", #keep hard process particles
    "keep abs(pdgId) = 13 || abs(pdgId) = 15" #keep muons and taus
    )
)

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
process.load("SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackHistory.TrackHistory_cff")
process.load("SimTracker.TrackHistory.TrackClassifier_cff")
process.load("RecoBTag.Configuration.RecoBTag_cff")


#-------------------------------------
## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outFilename),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output; you need a '*' to
    # unpack the list of commands 'patEventContent'
    outputCommands = cms.untracked.vstring('drop *', *patEventContent)
)

#-------------------------------------
if not options.miniAOD:
    ## PAT Configuration
    jetAlgo="AK4"

    from PhysicsTools.PatAlgos.tools.pfTools import *
    usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=not options.runOnData, postfix=postfix,
	      jetCorrections=jetCorrectionsAK4, pvCollection=cms.InputTag(pvSource))

    ## Top projections in PF2PAT
    getattr(process,"pfPileUpJME"+postfix).checkClosestZVertex = False
    getattr(process,"pfNoPileUpJME"+postfix).enable = options.usePFchs
    if options.useTTbarFilter:
	getattr(process,"pfNoMuonJME"+postfix).enable = False
	getattr(process,"pfNoElectronJME"+postfix).enable = False
	getattr(process,"pfNoTau"+postfix).enable = False
	getattr(process,"pfNoJet"+postfix).enable = False
    else:
	getattr(process,"pfNoMuonJME"+postfix).enable = options.useTopProjections
	getattr(process,"pfNoElectronJME"+postfix).enable = options.useTopProjections
	getattr(process,"pfNoTau"+postfix).enable = False
	getattr(process,"pfNoJet"+postfix).enable = False
else:
    ## Recreate tracks and PVs for b tagging
    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
    ## Select isolated collections
    process.selectedMuons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedMuons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
       (pfIsolationR04().sumChargedHadronPt+
	max(0.,pfIsolationR04().sumNeutralHadronEt+
	pfIsolationR04().sumPhotonEt-
	0.50*pfIsolationR04().sumPUPt))/pt < 0.20 && 
	(isPFMuon && (isGlobalMuon || isTrackerMuon) )'''))
    process.selectedElectrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedElectrons"), cut = cms.string('''abs(eta)<2.5 && pt>20. &&
	gsfTrack.isAvailable() &&
	gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\') < 2 &&
	(pfIsolationVariables().sumChargedHadronPt+
	max(0.,pfIsolationVariables().sumNeutralHadronEt+
	pfIsolationVariables().sumPhotonEt-
	0.5*pfIsolationVariables().sumPUPt))/pt < 0.15'''))

    ## Do projections
    process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
    process.pfNoMuonCHS =  cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfCHS"), veto = cms.InputTag("selectedMuons"))
    process.pfNoElectronsCHS = cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfNoMuonCHS"), veto = cms.InputTag("selectedElectrons"))

    process.pfNoMuon =  cms.EDProducer("CandPtrProjector", src = cms.InputTag("packedPFCandidates"), veto = cms.InputTag("selectedMuons"))
    process.pfNoElectrons = cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfNoMuon"), veto = cms.InputTag("selectedElectrons"))

    process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
    process.ak4GenJetsNoNu = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')

    if options.useTTbarFilter:
        if options.usePFchs:
            process.ak4PFJets = ak4PFJets.clone(src = 'pfCHS', doAreaFastjet = True)
        else:
            process.ak4PFJets = ak4PFJets.clone(src = 'packedPFCandidates', doAreaFastjet = True)
    else:
        if options.usePFchs:
            if options.useTopProjections:
                process.ak4PFJets = ak4PFJets.clone(src = 'pfNoElectronsCHS', doAreaFastjet = True)
            else:
                process.ak4PFJets = ak4PFJets.clone(src = 'pfCHS', doAreaFastjet = True)
        else:
            if options.useTopProjections:
                process.ak4PFJets = ak4PFJets.clone(src = 'pfNoElectrons', doAreaFastjet = True)
            else:
                process.ak4PFJets = ak4PFJets.clone(src = 'packedPFCandidates', doAreaFastjet = True)

from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection (done in order to use the above specified b-tag infos and discriminators)
switchJetCollection(
    process,
    jetSource = cms.InputTag(jetSource),
    pfCandidates = cms.InputTag(pfCandidates),
    pvSource = cms.InputTag(pvSource),
    svSource = cms.InputTag(svSource),
    muSource = cms.InputTag(muSource),
    elSource = cms.InputTag(elSource),
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK4,
    genJetCollection = cms.InputTag(genJetCollection),
    genParticles = cms.InputTag(genParticles),
    explicitJTA = options.useExplicitJTA,
    postfix = postfix
)

## Load standard PAT objects (here we only need PAT muons but the framework will figure out what it needs to run using the unscheduled mode)
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")
#-------------------------------------

#-------------------------------------
## AK8 jets (Gen and Reco)
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak8GenJetsNoNu = ak4GenJets.clone(
    rParam = cms.double(0.8),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu") if options.miniAOD else cms.InputTag("genParticlesForJetsNoNu"+postfix))
)

from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak8PFJets = ak4PFJets.clone(
    rParam = cms.double(0.8),
    src = (getattr(process,"ak4PFJets").src if options.miniAOD else getattr(process,"pfJetsPFBRECO"+postfix).src),
    srcPVs = (getattr(process,"ak4PFJets").srcPVs if options.miniAOD else getattr(process,"pfJetsPFBRECO"+postfix).srcPVs),
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(options.fatJetPtMin)
)

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
## AK8 soft drop jets
process.ak8GenJetsNoNuSoftDrop = ak4GenJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag("packedGenParticlesForJetsNoNu"),
    useSoftDrop = cms.bool(True),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    R0   = cms.double(0.8),
    useExplicitGhosts = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)

from RecoJets.JetProducers.ak4PFJetsSoftDrop_cfi import ak4PFJetsSoftDrop
process.ak8PFJetsSoftDrop = ak4PFJetsSoftDrop.clone(
    rParam = cms.double(0.8),
    src = (getattr(process,"ak4PFJets").src if options.miniAOD else getattr(process,"pfJetsPFBRECO"+postfix).src),
    srcPVs = (getattr(process,"ak4PFJets").srcPVs if options.miniAOD else getattr(process,"pfJetsPFBRECO"+postfix).srcPVs),
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(options.fatJetPtMin)
    )

## AK8 pruned jets (Gen and Reco)
process.ak8GenJetsNoNuPruned = ak4GenJets.clone(
    SubJetParameters,
    rParam = cms.double(0.8),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu") if options.miniAOD else cms.InputTag("genParticlesForJetsNoNu"+postfix)),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak4PFJetsPruned_cfi import ak4PFJetsPruned
process.ak8PFJetsPruned = ak4PFJetsPruned.clone(
    #jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    src = (getattr(process,"ak4PFJets").src if options.miniAOD else getattr(process,"pfJetsPFBRECO"+postfix).src),
    srcPVs = (getattr(process,"ak4PFJets").srcPVs if options.miniAOD else getattr(process,"pfJetsPFBRECO"+postfix).srcPVs),
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(options.fatJetPtMin)
)

if options.runSubJets:
    ## PATify AK8 jets
    addJetCollection(
        process,
        labelName = 'AK8',
        jetSource = cms.InputTag('ak8PFJets'),
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK8,
        genJetCollection = cms.InputTag('ak8GenJetsNoNu'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = options.useExplicitJTA,
        algo = 'AK',
        rParam = 0.8,
        postfix = postfix
    )
    addJetCollection(
        process,
        labelName = 'AK8SoftDrop',
        jetSource = cms.InputTag('ak8PFJetsSoftDrop'),
        btagInfos=['None'],
        btagDiscriminators=['None'],
        jetCorrections=jetCorrectionsAK8,
        genJetCollection = cms.InputTag('ak8GenJetsNoNu'),
        genParticles = cms.InputTag(genParticles),
        getJetMCFlavour = False,
        postfix = postfix
    )
    addJetCollection(
        process,
        labelName = 'AK8SoftDropSubJets',
        jetSource = cms.InputTag('ak8PFJetsSoftDrop','SubJets'),
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK4,
        genJetCollection = cms.InputTag('ak8GenJetsNoNuSoftDrop','SubJets'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = True, # needed for subjet b tagging
        svClustering = True, # needed for subjet b tagging
        algo = 'AK',
        rParam = 0.8,
        fatJets = cms.InputTag('ak8PFJets'), # needed for subjet flavor clustering
        groomedFatJets = cms.InputTag('ak8PFJetsSoftDrop'),
        postfix = postfix
    )
    addJetCollection(
        process,
        labelName = 'AK8Pruned',
        jetSource = cms.InputTag('ak8PFJetsPruned'),
        btagInfos=['None'],
        btagDiscriminators=['None'],
        jetCorrections=jetCorrectionsAK8,
        genJetCollection = cms.InputTag('ak8GenJetsNoNu'),
        genParticles = cms.InputTag(genParticles),
        getJetMCFlavour = False,
        postfix = postfix
    )
    addJetCollection(
        process,
        labelName = 'AK8PrunedSubJets',
        jetSource = cms.InputTag('ak8PFJetsPruned','SubJets'),
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK4,
        genJetCollection = cms.InputTag('ak8GenJetsNoNuPruned','SubJets'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = True, # needed for subjet b tagging
        svClustering = True, # needed for subjet b tagging
        algo = 'AK',
        rParam = 0.8,
        fatJets = cms.InputTag('ak8PFJets'), # needed for subjet flavor clustering
        groomedFatJets = cms.InputTag('ak8PFJetsPruned'),
        postfix = postfix
    )

    ## Establish references between PATified fat jets and subjets using the BoostedJetMerger
    process.selectedPatJetsAK8SoftDropPFlowPacked = cms.EDProducer("BoostedJetMerger",
        jetSrc=cms.InputTag("selectedPatJetsAK8SoftDrop"+postfix),
        subjetSrc=cms.InputTag("selectedPatJetsAK8SoftDropSubJets"+postfix)
    )
    process.selectedPatJetsAK8PrunedPFlowPacked = cms.EDProducer("BoostedJetMerger",
        jetSrc=cms.InputTag("selectedPatJetsAK8Pruned"+postfix),
        subjetSrc=cms.InputTag("selectedPatJetsAK8PrunedSubJets"+postfix)
    )

    ## New jet flavor still requires some cfg-level adjustments for subjets until it is better integrated into PAT
    ## Adjust the jet flavor for pruned subjets
    setattr(process,'patJetFlavourAssociationAK8SoftDropSubJets'+postfix, getattr(process,'patJetFlavourAssociationAK8'+postfix).clone(
        groomedJets = cms.InputTag('ak8PFJetsSoftDrop'),
        subjets = cms.InputTag('ak8PFJetsSoftDrop','SubJets')
    ))
    getattr(process,'patJetsAK8SoftDropSubJets'+postfix).JetFlavourInfoSource = cms.InputTag('patJetFlavourAssociationAK8SoftDropSubJets'+postfix,'SubJets')
    setattr(process,'patJetFlavourAssociationAK8PrunedSubJets'+postfix, getattr(process,'patJetFlavourAssociationAK8'+postfix).clone(
        groomedJets = cms.InputTag('ak8PFJetsPruned'),
        subjets = cms.InputTag('ak8PFJetsPruned','SubJets')
    ))
    getattr(process,'patJetsAK8PrunedSubJets'+postfix).JetFlavourInfoSource = cms.InputTag('patJetFlavourAssociationAK8PrunedSubJets'+postfix,'SubJets')


#-------------------------------------

#-------------------------------------
if options.runSubJets:
  #### Add groomed masses to jet userData
  from RecoJets.JetProducers.ak8PFJetsCHS_groomingValueMaps_cfi import *
  process.ak8PFJetsSoftDropMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
      src = cms.InputTag("ak8PFJets"),
      matched = cms.InputTag("selectedPatJetsAK8SoftDropPFlowPacked"),
      distMax = cms.double(0.8),
      value = cms.string('mass')
      ) 
  process.ak8PFJetsPrunedMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
      src = cms.InputTag("ak8PFJets"),
      matched = cms.InputTag("selectedPatJetsAK8PrunedPFlowPacked"),
      distMax = cms.double(0.8),
      value = cms.string('mass')
      ) 
  getattr(process,'patJetsAK8'+postfix).userData.userFloats.src += ['ak8PFJetsSoftDropMass', 'ak8PFJetsPrunedMass']
  ## N-subjettiness
  from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness

  process.Njettiness = Njettiness.clone(
      src = cms.InputTag("ak8PFJets"),
      cone = cms.double(0.8)
  )
  getattr(process,'patJetsAK8'+postfix).userData.userFloats.src += ['Njettiness:tau1','Njettiness:tau2','Njettiness:tau3']
#-------------------------------------

#-------------------------------------
if options.runOnData and options.runSubJets:
    ## Remove MC matching when running over data
    removeMCMatching( process, ['All'] )

## Add TagInfos to PAT jets
patJets = ['patJets'+postfix]
if options.runSubJets:
    patJets += ['patJetsAK8'+postfix,'patJetsAK8SoftDropSubJets'+postfix,'patJetsAK8PrunedSubJets'+postfix]

for m in patJets:
    if hasattr(process,m):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )
        print "Switching 'addJetFlavourInfo' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addJetFlavourInfo', cms.bool(True) )
#-------------------------------------

#-------------------------------------

## Filter for removing scraping events
process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

## Filter for good primary vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag(pvSource),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)
#-------------------------------------

#-------------------------------------
if options.useTTbarFilter:
    process.load("RecoBTag.PerformanceMeasurements.TTbarSelectionFilter_cfi")
    process.load("RecoBTag.PerformanceMeasurements.TTbarSelectionProducer_cfi")

    #electron id
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
        
    if options.miniAOD:
        process.ttbarselectionproducer.electronColl = cms.InputTag('slimmedElectrons')
        process.ttbarselectionproducer.muonColl     = cms.InputTag('slimmedMuons')
        process.ttbarselectionproducer.jetColl      = cms.InputTag('selectedPatJets'+postfix)
        process.ttbarselectionproducer.metColl      = cms.InputTag('slimmedMETs')
        switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
    else:
        process.ttbarselectionproducer.electronColl = cms.InputTag('selectedPatElectrons'+postfix)
        process.ttbarselectionproducer.muonColl     = cms.InputTag('selectedPatMuons'+postfix)
        process.ttbarselectionproducer.jetColl      = cms.InputTag('selectedPatJets'+postfix)
        process.ttbarselectionproducer.metColl      = cms.InputTag('patMETs'+postfix)
        switchOnVIDElectronIdProducer(process, DataFormat.AOD)

    setupAllVIDIdsInModule(process,
                           'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
                           setupVIDElectronSelection)


    ## Conversion rejection
    ## This should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer.
    #setattr(process,'patConversions'+postfix) = cms.EDProducer("PATConversionProducer",
        #electronSource = cms.InputTag('selectedPatElectrons'+postfix)
    #)
#-------------------------------------

#-------------------------------------
from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag(pvSource))
#-------------------------------------

#-------------------------------------
from RecoBTag.PerformanceMeasurements.BTagAnalyzer_cff import *
process.btagana = bTagAnalyzer.clone()
if options.useLegacyTaggers:
    process.btagana = bTagAnalyzerLegacy.clone()
# The following combinations should be considered:
# For b-tagging performance measurements:
#   process.btagana.useSelectedTracks    = True
#   process.btagana.useTrackHistory      = False (or True for Mistag systematics with GEN-SIM-RECODEBUG samples)
#   process.btagana.produceJetTrackTree  = False
#   process.btagana.produceAllTrackTree  = False
#   process.btagana.producePtRelTemplate = False (or True for PtRel fit studies)
# or data/MC validation of jets, tracks and SVs:
#   process.btagana.useSelectedTracks    = False (or True for JP calibration)
#   process.btagana.useTrackHistory      = False
#   process.btagana.produceJetTrackTree  = True
#   process.btagana.produceAllTrackTree  = False
#   process.btagana.producePtRelTemplate = False
# or general tracks, PV and jet performance studies:
#   process.btagana.useSelectedTracks    = True
#   process.btagana.useTrackHistory      = False
#   process.btagana.produceJetTrackTree  = False
#   process.btagana.produceAllTrackTree  = True
#   process.btagana.producePtRelTemplate = False
#------------------
process.btagana.useSelectedTracks     = True  ## False if you want to run on all tracks : for commissioning studies
process.btagana.useTrackHistory       = False ## Can only be used with GEN-SIM-RECODEBUG files
process.btagana.produceJetTrackTree   = False ## True if you want to keep info for tracks associated to jets : for commissioning studies
process.btagana.produceAllTrackTree   = False ## True if you want to keep info for all tracks : for commissioning studies
process.btagana.producePtRelTemplate  = options.producePtRelTemplate  ## True for performance studies
#------------------
process.btagana.storeTagVariables     = False ## True if you want to keep TagInfo TaggingVariables
process.btagana.storeCSVTagVariables  = True  ## True if you want to keep CSV TaggingVariables
process.btagana.primaryVertexColl     = cms.InputTag(pvSource)
process.btagana.Jets                  = cms.InputTag('selectedPatJets'+postfix)
#FIXME it used to be muons, slimmedMuons from miniAOD
process.btagana.muonCollectionName    = cms.InputTag(muSource)
process.btagana.patMuonCollectionName = cms.InputTag(patMuons)
process.btagana.use_ttbar_filter      = cms.bool(options.useTTbarFilter)
process.btagana.triggerTable          = cms.InputTag('TriggerResults::HLT') # Data and MC
process.btagana.genParticles          = cms.InputTag(genParticles)
process.btagana.candidates            = cms.InputTag(pfCandidates) 

if options.runSubJets:
    process.btaganaSoftDropSubJets = process.btagana.clone(
        storeEventInfo      = cms.bool(False),
        produceJetTrackTree = cms.bool(True),
        allowJetSkipping    = cms.bool(False),
        Jets                = cms.InputTag('selectedPatJetsAK8SoftDropSubJets'+postfix),
        FatJets             = cms.InputTag('selectedPatJetsAK8'+postfix),
        GroomedFatJets      = cms.InputTag('selectedPatJetsAK8SoftDropPFlowPacked'),
        runSubJets          = options.runSubJets,
        use_ttbar_filter    = cms.bool(False)
    )
    process.btaganaPrunedSubJets = process.btagana.clone(
        storeEventInfo      = cms.bool(False),
        produceJetTrackTree = cms.bool(True),
        allowJetSkipping    = cms.bool(False),
        Jets                = cms.InputTag('selectedPatJetsAK8PrunedSubJets'+postfix),
        FatJets             = cms.InputTag('selectedPatJetsAK8'+postfix),
        GroomedFatJets      = cms.InputTag('selectedPatJetsAK8PrunedPFlowPacked'),
        runSubJets          = options.runSubJets,
        use_ttbar_filter    = cms.bool(False)
    )

#---------------------------------------

#---------------------------------------
## Trigger selection !
#import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
#process.JetHLTFilter = hlt.triggerResultsFilter.clone(
#    triggerConditions = cms.vstring(
#        "HLT_PFJet80_v*"
#    ),
#    hltResults = cms.InputTag("TriggerResults","","HLT"),
#    l1tResults = cms.InputTag( "" ),
#    throw = cms.bool( False ) #set to false to deal with missing triggers while running over different trigger menus
#)
#---------------------------------------

#---------------------------------------
## Optional MET filters:
## https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters
#process.load("RecoMET.METFilters.metFilters_cff")
#process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')
#---------------------------------------

#---------------------------------------
## Event counter
from MyAnalysis.EventCounter.eventcounter_cfi import eventCounter
process.allEvents = eventCounter.clone()
process.selectedEvents = eventCounter.clone()
#---------------------------------------

#---------------------------------------
## Define event filter sequence
process.filtSeq = cms.Sequence(
    #process.JetHLTFilter*
    #process.noscraping
    process.primaryVertexFilter
)


## Define analyzer sequence
process.analyzerSeq = cms.Sequence( )
if options.processStdAK4Jets:
    process.analyzerSeq += process.btagana
if options.runSubJets:
    process.analyzerSeq += process.btaganaSoftDropSubJets
    process.analyzerSeq += process.btaganaPrunedSubJets
if options.processStdAK4Jets and options.useTTbarFilter:
    process.analyzerSeq.replace( process.btagana, 
                                 process.egmGsfElectronIDSequence *
                                 process.ttbarselectionproducer * 
                                 process.ttbarselectionfilter * 
                                 process.btagana )
#---------------------------------------

process.p = cms.Path(
    process.allEvents
    * process.filtSeq
    * process.selectedEvents
    * process.analyzerSeq
)

# Delete predefined output module (needed for running with CRAB)
del process.out

open('pydump.py','w').write(process.dumpPython())
