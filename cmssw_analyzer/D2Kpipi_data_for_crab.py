import FWCore.ParameterSet.Config as cms

process = cms.Process("D2Kpipi")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('D2Kpipi')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
   limit = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
#process.load("SimTracker.TrackAssociation.TrackAssociatorByChi2_cfi")
#process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

# Data sample clean up
process.noscraping = cms.EDFilter("FilterOutScraping",
                                applyfilter = cms.untracked.bool(True),
                                debugOn = cms.untracked.bool(True),
                                numtrack = cms.untracked.uint32(10),
                                thresh = cms.untracked.double(0.25)
                                )
#process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
#                                           vertexCollection = cms.InputTag('offlinePrimaryVerticesWithBS'),
#                                           minimumNDOF = cms.uint32(4) ,
#                                           maxAbsZ = cms.double(24), 
#                                           maxd0 = cms.double(2) 
#                                           )
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')

process.GlobalTag.globaltag = 'FT_R_42_V10A::All'     # For /MinimumBias/Run2010A-Apr21ReReco-v1/AOD
#process.GlobalTag.globaltag = 'FT_R_42_V13A::All'     # For ?

# Input source
## process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
##     fileNames = cms.untracked.vstring(
         #'file:MinimumBias_Run2010A_Apr21ReReco_100evts.root'
##          'file:MinimumBias_Commissioning10_May19ReReco_100evts.root'
##    )
##)

# Use this for CRAB
# Input source
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       # '/store/user/wcarvalh/Test_CCbar_40_70_7TeV_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_MC_42_V13/Test_CCbar_40_70_7TeV_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_MC_42_V13/e44cf2cfe5a7fd2b0c9fce168dabe8a5/CCbar_40_70_7TeV_RAW_1_2_M5A.root',
] );

secFiles.extend( [
               ] )

process.load("charm.DmesonAnalyzer.d2kpipi_data_analyzer_cfi")
process.d2kpipi.PFlowThresholds = cms.PSet(
    Barrel = cms.PSet(
        hadronCharged = cms.PSet( pt = cms.double(0.0) ),
        hadronNeutral = cms.PSet( energy = cms.double(1.4) ),
        electron = cms.PSet( pt = cms.double(0.0) ),
        gamma = cms.PSet( energy = cms.double(0.9) )
    ),
    Endcap = cms.PSet(
        hadronCharged = cms.PSet( pt = cms.double(0.0) ),
        hadronNeutral = cms.PSet( energy = cms.double(2.7) ),
        electron = cms.PSet( pt = cms.double(0.0) ),
        gamma = cms.PSet( energy = cms.double(2.5) )
    ),
    Transition = cms.PSet(
        hadronCharged = cms.PSet( pt = cms.double(0.5) ),      # 0.0 before
        hadronNeutral = cms.PSet( energy = cms.double(3.8) ),
        electron = cms.PSet( pt = cms.double(0.5) ),           # 0.0 before
        gamma = cms.PSet( energy = cms.double(2.5) ),
        hadronHF = cms.PSet( energy = cms.double(4.0) ),
        emHF = cms.PSet( energy = cms.double(4.0) )            # 3.5 before
    ),
    Forward = cms.PSet(
        hadronHF = cms.PSet( energy = cms.double(4.0) ),
        emHF = cms.PSet( energy = cms.double(4.0) )            # 3.5 before
    )
)
#process.d2kpipi.minTracks = 200
#process.d2kpipi = cms.EDAnalyzer('DmesonAnalyzer',
#        minTracks = cms.untracked.uint32(150)
#    )

process.TFileService = cms.Service("TFileService",
    #fileName = cms.string('D2Kpipi_MB_Commissioning10-May19ReReco-v1.root')
    #fileName = cms.string('D2Kpipi_Run2010A-Apr21ReReco_hltTrigger.root')
    fileName = cms.string('D2Kpipi_Run2010A-Apr21ReReco.root')
)

#from HLTrigger.HLTfilters.hltHighLevel_cfi import *
#process.triggerSelection = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_L1_BscMinBiasOR_BptxPlusORMinus"])

process.p = cms.Path(#process.triggerSelection + 
                     process.noscraping       + 
                     process.HBHENoiseFilter  + 
                     process.d2kpipi)
