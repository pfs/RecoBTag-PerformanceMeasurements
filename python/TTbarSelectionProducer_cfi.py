import FWCore.ParameterSet.Config as cms

ttbarselectionproducer = cms.EDProducer("TTbarSelectionProducer",
                                        verbose          = cms.int32(6),
                                        triggerColl      = cms.InputTag("TriggerResults","","HLT"),
                                        trigNamesToSel   = cms.vstring(),
                                        doTrigSel        = cms.bool(False),
                                        electronColl     = cms.InputTag("slimmedElectrons"),
                                        electron_cut_pt  = cms.double(20),
                                        electron_cut_eta = cms.double(2.5),
                                        electron_cut_iso = cms.double(0.15),
                                        muonColl         = cms.InputTag("slimmedMuons"),
                                        muon_cut_pt      = cms.double(20),
                                        muon_cut_eta     = cms.double(2.4),
                                        muon_cut_iso     = cms.double(0.15),
                                        jetColl          = cms.InputTag("slimmedJets"),
                                        jet_cut_pt       = cms.double(30),
                                        jet_cut_eta      = cms.double(2.5),
                                        metColl          = cms.InputTag("slimmedMETs"),
                                        met_cut          = cms.double(40)
                                        )
