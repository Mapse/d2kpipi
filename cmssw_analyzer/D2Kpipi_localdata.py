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
#process.GlobalTag.globaltag = 'FT_R_42_V13A::All'
#process.GlobalTag.globaltag = 'START42_V17C::All'     # For MB Summer12 LowPU2010

# Input source
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
         #'file:MinimumBias_Run2010A_Apr21ReReco_100evts.root'
         'file:MinimumBias_Run2010A_Apr21ReReco_run136035_500evts.root'
         #'file:MinimumBias_Commissioning10_May19ReReco_100evts.root'
         #'file:MinBias_Tune4C_HFshowerLibrary_7TeV_pythia8_100evts.root'
    )
)

# Use this for CRAB
# Input source
#readFiles = cms.untracked.vstring()
#secFiles = cms.untracked.vstring() 
#process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
#readFiles.extend( [
       # '/store/user/wcarvalh/Test_CCbar_40_70_7TeV_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_MC_42_V13/Test_CCbar_40_70_7TeV_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_MC_42_V13/e44cf2cfe5a7fd2b0c9fce168dabe8a5/CCbar_40_70_7TeV_RAW_1_2_M5A.root',
       # '/store/user/wcarvalh/Test_CCbar_40_70_7TeV_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_MC_42_V13/Test_CCbar_40_70_7TeV_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_MC_42_V13/e44cf2cfe5a7fd2b0c9fce168dabe8a5/CCbar_40_70_7TeV_RAW_2_3_bN3.root' 
#] );

#secFiles.extend( [
#               ] )

process.load("charm.DmesonAnalyzer.d2kpipi_data_analyzer_cfi")
#process.load("charm.DmesonAnalyzer.d2kpipi_mc_analyzer_cfi")
#process.d2kpipi.minTracks = 200
#process.d2kpipi = cms.EDAnalyzer('DmesonAnalyzer'
#process.d2kpipi = cms.EDAnalyzer('DmesonAnalyzer',
#        minTracks = cms.untracked.uint32(150)
#    )

process.TFileService = cms.Service("TFileService",
    #fileName = cms.string('D2Kpipi_data_Commissioning10_hltTrigger.root')
    #fileName = cms.string('D2Kpipi_data_Run2010A_hltTrigger.root')
    fileName = cms.string('D2Kpipi_data_Run2010A_500evts.root')
    #fileName = cms.string('MB_pythia8_4C_7TeV_HF-SL_Summer12-LowPU2010_DR42_NoPileUp_START42_V17C-v1_100evt.root')
)

from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.triggerSelection = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_L1_BscMinBiasOR_BptxPlusORMinus"])

process.p = cms.Path(#process.triggerSelection + 
                     process.noscraping       + 
                     process.HBHENoiseFilter  + 
                     process.d2kpipi)
