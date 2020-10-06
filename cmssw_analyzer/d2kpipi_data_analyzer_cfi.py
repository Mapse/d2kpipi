import FWCore.ParameterSet.Config as cms

d2kpipi = cms.EDAnalyzer('D2Kpipi_Analyzer', 
    DoGEN           = cms.untracked.bool(False),
    DoRECO          = cms.untracked.bool(True),
    DoSinglePV      = cms.untracked.bool(True),
    trackPtMin      = cms.untracked.double(0.50),   #    0.5 - 0.25
    trackHitsMin    = cms.untracked.double(5),      #    5
    trackChi2Max    = cms.untracked.double(5.0),    #    5.0
    dcaXY           = cms.untracked.double(-0.00),  #    0.1 - 0.05
    dcaZ            = cms.untracked.double(-0.00),  #    0.2 - 0.05
    DplusMassWindow = cms.untracked.double(0.17),   #
    vtxChi2Max      = cms.untracked.double(3.0),    #    1.0 - 2.0
    vtxDispSigXYMin = cms.untracked.double(1.5),    #    2.0
    angleMax        = cms.untracked.double(0.20),   #    0.12
    rMax            = cms.untracked.double(2.0),    #    2.0
    # Diffractive analysis
    HLTPath           = cms.string("HLT_L1_BscMinBiasOR_BptxPlusORMinus"), 
    ParticleFlowTag   = cms.string("particleFlow"),
    energyThresholdHB = cms.double(1.5),    #    1.5
    energyThresholdHE = cms.double(2.0),    #    2.0
    energyThresholdHF = cms.double(4.0),    #    4.0
    comEnergy         = cms.double(7000.),  #   7000
    applyEnergyScaleHCAL  = cms.bool(False),
    EnergyScaleFactorHCAL = cms.double(1.0)    #    1.0
)
