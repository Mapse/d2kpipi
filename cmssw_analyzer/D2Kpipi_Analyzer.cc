// -*- C++ -*-
//
// Package:    D2Kpipi_Analyzer
// Class:      D2Kpipi_Analyzer
// 
/**\class D2Kpipi_Analyzer D2Kpipi_Analyzer.cc charm/D2Kpipi_Analyzer/src/D2Kpipi_Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Wagner De Paula Carvalho,22 1-020,+41227677686,
//         Created:  Thu Nov 29 23:54:40 CET 2012
// $Id$
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

// For MC
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

// Particle Flow
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

// JetCollection (CaloTower) - Gapside
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

//-------------------------------------------------------------
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
//--------------------------------------------------------------

// Trigger

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/RegexMatch.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
//#include "ThreeVector.h"
//
// class declaration
//

// using namespace edm;
// using namespace reco;

class D2Kpipi_Analyzer : public edm::EDAnalyzer {
   
   public:
      
      explicit D2Kpipi_Analyzer(const edm::ParameterSet&);
      ~D2Kpipi_Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
      typedef std::vector<reco::TransientTrack>  TTrackVec  ;
      typedef std::vector<TransientVertex> TVertexVec ;
      
      void analyzeGen(const edm::Event& , const edm::EventSetup&) ;
      void analyzeReco(const edm::Event& , const edm::EventSetup&) ;
      bool isD2Kpipi(const reco::GenParticle&) ;
      double findAngle(const reco::Vertex& , const TransientVertex& , const TLorentzVector& ) ;


   private:
      
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // Diffractive analysis
      
      enum calo_region_t {Barrel,Endcap,Transition,Forward};

      void initialize();
      void TriggerInfo(const edm::Event& , const edm::EventSetup&) ;
      void Gapside(const edm::Event& , const edm::EventSetup&) ;
      void resetPFThresholds(std::map<int,std::pair<double,double> >&);
      void setPFThresholds(std::map<int,std::pair<double,double> >&, edm::ParameterSet const&);
      int pflowId(std::string const& name);
      bool pflowThreshold(reco::PFCandidate const& part, std::map<int, std::map<int,std::pair<double,double> > > const& thresholdMap);
      double MassColl(reco::PFCandidateCollection const& pflowCollection, std::map<int, std::map<int,std::pair<double,double> > > const& thresholdMap);
      std::pair<double,double> xi(reco::PFCandidateCollection const& pflowCollection, double Ebeam, 
                                  std::map<int, std::map<int,std::pair<double,double> > > const& thresholdMap);
      std::pair<double,double> EPlusPz(reco::PFCandidateCollection const& pflowCollection, std::map<int, 
                                       std::map<int,std::pair<double,double> > > const& thresholdMap);
      std::pair<double,double> sumEHF(reco::PFCandidateCollection const& pflowCollection, std::map<int, 
                                       std::map<int,std::pair<double,double> > > const& thresholdMap);
      std::pair<double,double> etaMax(reco::PFCandidateCollection const& pflowCollection, std::map<int, 
                                      std::map<int,std::pair<double,double> > > const& thresholdMap);
      // void etaMaxMin(const edm::Event& iEvent, const edm::EventSetup&, std::map<int,std::map<int,std::pair<double,double> > > const& thresholdMap);
      void sumEHFTowers(const edm::Event&, const edm::EventSetup&);


   private:

      // ----------member data ---------------------------
      
      // Parameters passed through the config
      bool DoGEN_ , DoRECO_;     // Flags to process GEN and RECO information
      bool DoSinglePV_ ;         // Flag to select events with a single primary vertex
      double trackPtMin_ , trackHitsMin_ , trackChi2Max_ , dcaXY_ , dcaZ_, 
             DplusMassWindow_ , vtxChi2Max_ , vtxDispSigXYMin_ , angleMax_, rMax_ ;
      std::string hltPathName_; 

      // Set data content to be analyzed GEN and/or RECO 
      // bool DoGEN_ , DoRECO_ ;
      
      // Flags to pass information on generated Ds processes to reconstructed events
      bool hasDplus ,        //  Event has a D+ meson
           hasDplus2Kpipi ,  //  Event has a D+ -> K- pi+ pi+  (or conjugate) 
           hasD0 ,           //  Event has a D0 meson
           hasDs ,           //  Event has a Ds meson
           hasDstarplus ;    //  Event has a D*+ meson
      
      // Histograms for GEN quantities
      TH1D  *histo_nDplusGen , *histo_nD0Gen , *histo_nDsGen , 
            *histo_DplusGenVertexX , *histo_DplusGenVertexY , *histo_DplusGenVertexZ , 
            *histo_DplusDecVertexX , *histo_DplusDecVertexY , *histo_DplusDecVertexZ , 
//            *histo_DplusGenMass , *histo_DplusGenCharge ,
            *histo_DplusGenPt , *histo_DplusGenP , *histo_DplusGenPhi , *histo_DplusGenEta , 
            *histo_D2KpipiGenPt , *histo_D2KpipiGenP , *histo_D2KpipiGenPhi , *histo_D2KpipiGenEta , 
            *histo_DplusGenLength , *histo_DplusGenLengthXY , 
            *histo_KGenPt , *histo_KGenP , // *histo_ktGen , 
            *histo_piGenPt , *histo_piGenP , // *histo_piGenPerp , 
            *histo_KpiGenPtsum ; 
//            *histo_KDplusGenAngle , *histo_piDplusGenAngle , *histo_KpiGenAngle ; 
            
      // Histograms for RECO quantities
      TH1D  *histo_nTrack , *histo_trackPt , *histo_trackP , *histo_trackHits, *histo_trackNormChi2 , 
            *histo_nVertex , *histo_nTracksPerVertex , 
            *histo_PVNormChi2 , *histo_fracTracksPV , 
            *histo_PVX , *histo_PVY , *histo_PVZ ,
            *histo_PVXError , *histo_PVYError , *histo_PVZError ,
            *histo_trackDxy , *histo_trackDz , *histo_PVtrackDxy , *histo_PVtrackDz ,
            *histo_nDplus , 
            *histo_massDplus , *histo_refitMassDplus ,
            *histo_DplusVertexNormChi2 ,
            *histo_DplusVertexX , *histo_DplusVertexY , *histo_DplusVertexZ , 
            *histo_DplusVertexXError , *histo_DplusVertexYError , *histo_DplusVertexZError , 
            *histo_DplusVertexDxy , *histo_DplusVertexExy , *histo_DplusVertexSxy ,
            *histo_DplusVertexD3D , *histo_DplusVertexE3D , *histo_DplusVertexS3D ,
            *histo_angle , *histo_refittedAngle , 
            *histo_zoomedAngle , *histo_zoomedRefittedAngle , 
            *histo_matchedAngle , *histo_matchedRefittedAngle ; 
      //  For fake combinations (all same charge)
      TH1D  *histo_nFakeDplus , 
            *histo_massFakeDplus , *histo_refitMassFakeDplus ,
            *histo_FakeDplusVertexNormChi2 ,
            *histo_FakeDplusVertexX , *histo_FakeDplusVertexY , *histo_FakeDplusVertexZ , 
            *histo_FakeDplusVertexXError , *histo_FakeDplusVertexYError , *histo_FakeDplusVertexZError , 
            *histo_FakeDplusVertexDxy , *histo_FakeDplusVertexExy , *histo_FakeDplusVertexSxy ,
            *histo_FakeDplusVertexD3D , *histo_FakeDplusVertexE3D , *histo_FakeDplusVertexS3D ,
            *histo_FakeAngle , *histo_FakeRefittedAngle , 
            *histo_FakeZoomedAngle , *histo_FakeZoomedRefittedAngle , 
            *histo_FakeMatchedAngle , *histo_FakeMatchedRefittedAngle ; 
            
      // Histograms for GEN-RECO quantities
      TH1D  *histo_DplusGenRecoDR ;
      TH1D  *histo_angleGenReco , *histo_angleGenReco2 , *histo_angleGenReco3 ;
      TH1D  *histo_angleFakeGenReco , *histo_angleFakeGenReco2 , *histo_angleFakeGenReco3 ;

      TTree *genTree ;
      
      // GEN/RECO Tree variables
//      double _run ,                 //  run number
      int    _run ,                 //  run number
             _event ,               //  event number
             _procId ;              //  ID of the simulated physics proccess
      
      // GEN Tree variables
      double _D_x ,                 //  Dplus -> Kpipi  decay position X
             _D_y ,                 //  Dplus -> Kpipi  decay position Y
             _D_z ,                 //  Dplus -> Kpipi  decay position Z
             _D_xy_length ,         //  Dplus -> Kpipi  decay length in the XY plane
             _D_length ,            //  Dplus -> Kpipi  decay length 
             _D2Kpipi_p ,           //  Dplus -> Kpipi  momentum
             _D2Kpipi_pT ,          //  Dplus -> Kpipi  pT
             _D2Kpipi_eta ,         //  Dplus -> Kpipi  eta
             _D2Kpipi_phi ,         //  Dplus -> Kpipi  phi
             _D2Kpipi_pTfraction ,  //  Ratio Dplus -> Kpipi  pT over scalar pT sum of all particles
             _D2Kpipi_pT2fraction , //  Ratio Dplus -> Kpipi  pT^2 over pT^2 sum of all particles
             _mean_pT ,             //  Average particle pT
             _mean_pT2 ,            //  Average particle pT^2 
             _K_p ,                 //  p of kaon from Dplus -> Kpipi
             _K_pT ,                //  pT of kaon from Dplus -> Kpipi
             _K_phi ,               //  phi of kaon from Dplus -> Kpipi
             _K_eta ,               //  eta of kaon from Dplus -> Kpipi
             _pi1_p ,               //  p of pion from Dplus -> Kpipi
             _pi1_pT ,              //  pT of pion from Dplus -> Kpipi
             _pi1_phi ,             //  phi of pion from Dplus -> Kpipi
             _pi1_eta ,             //  eta of pion from Dplus -> Kpipi
             _pi2_p ,               //  p of pion from Dplus -> Kpipi
             _pi2_pT ,              //  pT of pion from Dplus -> Kpipi
             _pi2_phi ,             //  phi of pion from Dplus -> Kpipi
             _pi2_eta ;             //  eta of pion from Dplus -> Kpipi
//             _kt ,                  //  transverse momentum of pion/kaon with respect to Dplus line of flight
//             _K2D_angle ,           //  kaon angle with respect to Dplus line of flight
//             _pi2D_angle ;          //  pion angle with respect to Dplus line of flight
      int    _nGenPart ;          //  Number of GEN particles in the event

      TTree *theTree ;

      // RECO Tree variables
      bool   _hasdplus ,          //  Has D+ meson
             _hasdplus2kpipi ,    //  Has D-/+ -> K pi pi
             _hasd0 ,             //  Has D0 meson
             _hasds ,             //  Has Ds meson
             _hasdstarplus ;      //  Has D*+ meson
      int    _HLTPath ;           //  HLT path. 1 = selected HLT path fired
      double _massDplus ,         //  Mass of the Kpipi combination
             _refittedMassDplus , //  Mass of the Kpipi combination with refitted tracks
             _refittedMassKpi1 ,  //  Mass of the 1st K+pi combination with refitted tracks
             _refittedMassKpi2 ;  //  Mass of the 2nd K+pi combination with refitted tracks
      int    _nvtx ;              //  Number of recoed primary vertices 
      double _vtxchi2dof ,        //  Secondary vertex candidate chi2/dof
             _Dxy ,               //  Distance Lxy from secondary to the primary vertex in the XY plane
             _Sxy ,               //  Secondary to primary vertex significance ( sigma_Lxy/Lxy ) in the XY plane
             _D3D ,               //  3D distance L from secondary to the primary vertex
             _S3D ,               //  3D secondary to primary vertex significance ( sigma_L/L )
             _Rxy ,               //  Secondary vertex XY distance with respect to the nominal beam position
             _angle,              //  Matching angle between D(->Kpi) candidate momentum and displacement w/ respect to PV   
             _refittedAngle,      //  Matching angle between D(->Kpi) candidate momentum and displacement w/ respect to PV   
             _dplusp ,            //  Dplus candidate momentum
             _dpluspT ,           //  Dplus candidate pT
             _dpluseta ,          //  Dplus candidate eta
             _dplusphi ,          //  Dplus candidate phi
             _trk1dcaXY ,         //  Distance of closest approach of track 1 to PV in XY plane 
             _trk1dcaZ ,          //  Distance of closest approach of track 1 to PV in Z 
             _trk2dcaXY ,         //  Distance of closest approach of track 2 to PV in XY plane  
             _trk2dcaZ ,          //  Distance of closest approach of track 2 to PV in Z 
             _trk3dcaXY ,         //  Distance of closest approach of track 3 to PV in XY plane  
             _trk3dcaZ ,          //  Distance of closest approach of track 3 to PV in Z 
             _trk1p ,             //  Track 1 momentum 
             _trk2p ,             //  Track 2 momentum 
             _trk3p ,             //  Track 3 momentum 
             _trk1pT ,            //  Track 1 pT 
             _trk2pT ,            //  Track 2 pT 
             _trk3pT ,            //  Track 3 pT 
             _trk1chi2 ,          //  Track 1 chi2/dof 
             _trk2chi2 ,          //  Track 2 chi2/dof 
             _trk3chi2 ;          //  Track 3 chi2/dof 
//             _trk1angle ,         //  Track 1 angle with respect to the vector sum of track1 and track2 momenta
//             _trk2angle ,         //  Track 2 angle with respect to the vector sum of track1 and track2 momenta
//             _trkperp ;           //  Tracks pT with respect to the direction of the summed momentum of track1 and track2   
      int    _trk1hits ,          //  Number of valid hits of fitted track 1
             _trk2hits ,          //  Number of valid hits of fitted track 2
             _trk3hits ;          //  Number of valid hits of fitted track 3
      //Diffractive analysis
      int    _nHFPlus ,
             _nHFMinus ;
      double _HFPlusEnergyFromCaloTowers ,
             _HFMinusEnergyFromCaloTowers ,
             _HFPlusEnergyFromParticleFlow ,
             _HFMinusEnergyFromParticleFlow ,
             _EPlusPzFromPFCands ,
             _EMinusPzFromPFCands ,
             _xiPlusFromPFCands ,
             _xiMinusFromPFCands ,
             _etaMaxFromPFCands ,
             _etaMinFromPFCands ,
             _MxFromPFCands ,
             _missingMassFromXiFromPFCands ;
             // _pfsis1Eta_max ,                //  removed due redundancy with _etaMaxFromPFCands
             // _pfsis1Eta_min ,                //  removed due redundancy with _etaMinFromPFCands
             // _pfsis2Eta_max ,                //
             // _pfsis2Eta_min ,                //
             // _deltaEtapf ;                   //  redundancy (= _etaMaxFromPFCands - _etaMinFromPFCands)

      TTree *fakeTree ;

      int total_nD2KpipiGen ;
      int total_nDpluscand ;
      int total_nDpluscandFake ;

      //Diffractive analysis
      
      edm::InputTag particleFlowTag_;

      double energyThresholdHB_;
      double energyThresholdHE_;
      double energyThresholdHF_;
      double comEnergy_,Ebeam_;
      bool applyEnergyScaleHCAL_;
      double energyScaleHCAL_;

      std::map<int,std::pair<double,double> > thresholdsPFlowBarrel_;
      std::map<int,std::pair<double,double> > thresholdsPFlowEndcap_;
      std::map<int,std::pair<double,double> > thresholdsPFlowTransition_;
      std::map<int,std::pair<double,double> > thresholdsPFlowForward_;
      std::map<int,std::map<int,std::pair<double,double> > > thresholdsPFlow_;
      
};

//
// constants, enums and typedefs
//

   // Particles Masses
   
   double const MMUON  = 0.10566 ;   //  PDGID =  13
   double const MPION  = 0.139570 ;  //  PDGID = 211
   double const MKAON  = 0.49368 ;   //  PDGID = 321
   double const MPHI   = 1.01946 ;   //  PDGID = 333
   double const MD0    = 1.86484 ;   //  PDGID = 421 , 423 (D0*)
   double const MDPLUS = 1.86962 ;   //  PDGID = 411 , 413 (D+*) => M_DPLUS_STAR = 2.01027 
   double const MD_S   = 1.96850 ;   //  PDGID = 431
   double const MB0    = 5.27953 ;   //  PDGID = 511
   double const MBPLUS = 5.27927 ;   //  PDGID = 521
   
//
// static data member definitions
//

//
// constructors and destructor
//
D2Kpipi_Analyzer::D2Kpipi_Analyzer(const edm::ParameterSet& iConfig) : 
   DoGEN_(iConfig.getUntrackedParameter<bool>("DoGEN",true)),
   DoRECO_(iConfig.getUntrackedParameter<bool>("DoRECO",true)),
   DoSinglePV_(iConfig.getUntrackedParameter<bool>("DoSinglePV",false)),
   trackPtMin_(iConfig.getUntrackedParameter<double>("trackPtMin",0.5)),
   trackHitsMin_(iConfig.getUntrackedParameter<double>("trackHitsMin",5)),
   trackChi2Max_(iConfig.getUntrackedParameter<double>("trackChi2Max",5.0)),
   dcaXY_(iConfig.getUntrackedParameter<double>("dcaXY",-0.01)),
   dcaZ_(iConfig.getUntrackedParameter<double>("dcaZ",-0.01)),
   DplusMassWindow_(iConfig.getUntrackedParameter<double>("DplusMassWindow",0.20)),
   vtxChi2Max_(iConfig.getUntrackedParameter<double>("vtxChi2Max",5.0)),
   vtxDispSigXYMin_(iConfig.getUntrackedParameter<double>("vtxDispSigXYMin",1.0)),
   angleMax_(iConfig.getUntrackedParameter<double>("angleMax",3.2)),
   rMax_(iConfig.getUntrackedParameter<double>("rMax",2.0)),
   // Diffractive analysis
   hltPathName_(iConfig.getParameter<std::string>("HLTPath")),
   particleFlowTag_(iConfig.getParameter<edm::InputTag>("ParticleFlowTag")),
   energyThresholdHB_(iConfig.getParameter<double>("energyThresholdHB")),
   energyThresholdHE_(iConfig.getParameter<double>("energyThresholdHE")),
   energyThresholdHF_(iConfig.getParameter<double>("energyThresholdHF")),
   comEnergy_(iConfig.getParameter<double>("comEnergy")),
   applyEnergyScaleHCAL_(iConfig.getParameter<bool>("applyEnergyScaleHCAL")),
   energyScaleHCAL_(-1.)
{
   //now do what ever initialization is needed

   // DoGEN_ = true ;
   // DoRECO_ = true ;
   
   hasDplus = false ;
   hasDplus2Kpipi = false ;
   hasD0 = false ;
   hasDs = false ;
   hasDstarplus = false ;
   
   edm::Service<TFileService> outfile;

  if(DoRECO_)
  { 
   // Must be initialized to -1 for real data
   _procId = -1 ;
   
   //  RECO Tracks
   histo_nTrack = outfile->make<TH1D>("histo_nTrack","Track multiplicity",100,0,200);
   histo_trackPt = outfile->make<TH1D>("histo_trackPt","Track pT; (GeV/c)",100,0,25.);
   histo_trackP = outfile->make<TH1D>("histo_trackP","Track momentum; (GeV/c)",100,0,50.);
   histo_trackHits = outfile->make<TH1D>("histo_trackHits","Valid hits on track",50,0,50);
   histo_trackNormChi2 = outfile->make<TH1D>("histo_trackNormChi2","Track chi2/dof",100,0,10);
   
   histo_trackDxy = outfile->make<TH1D>("histo_trackDxy","Track XY distance to PV; (cm)",100,-0.5,0.5);
   histo_trackDz = outfile->make<TH1D>("histo_trackDz","Track Z distance to PV; (cm)",100,-1,1);
   histo_PVtrackDxy = outfile->make<TH1D>("histo_PVtrackDxy","PV Track XY distance to PV; (cm)",100,-0.5,0.5);
   histo_PVtrackDz = outfile->make<TH1D>("histo_PVtrackDz","PV Track Z distance to PV; (cm)",100,-1,1);
   
   //  RECO Vertices
   histo_nVertex = outfile->make<TH1D>("histo_nVertex","Vertex multiplicity",10,0,10);
   histo_nTracksPerVertex = outfile->make<TH1D>("histo_nTracksPerVertex","Tracks per vertex",100,0,200);
   histo_fracTracksPV = outfile->make<TH1D>("histo_fracTracksPV","Fraction of tracks from primary vertices",50,0,1);
   histo_PVNormChi2 = outfile->make<TH1D>("histo_PVNormChi2","Vertex chi2/dof",100,0,3);
   histo_PVX = outfile->make<TH1D>("histo_PVX","Vertex X position; (cm)",100,0,0.6);
   histo_PVY = outfile->make<TH1D>("histo_PVY","Vertex Y position; (cm)",100,0,0.6);
   histo_PVZ = outfile->make<TH1D>("histo_PVZ","Vertex Z position; (cm)",200,-20,20);
   histo_PVXError = outfile->make<TH1D>("histo_PVXError","Vertex X error; (cm)",100,0,0.03);
   histo_PVYError = outfile->make<TH1D>("histo_PVYError","Vertex Y error; (cm)",100,0,0.03);
   histo_PVZError = outfile->make<TH1D>("histo_PVZError","Vertex Z error; (cm)",200,0,0.2);

   histo_nDplus = outfile->make<TH1D>("histo_nDplus","Multiplicity of Dplus candidates per event",30,0,30);
//   histo_nDplus = outfile->make<TH1D>("histo_nDplus","Multiplicity of Dplus candidates per event",101,-0.5,100.5);
   histo_massDplus = outfile->make<TH1D>("histo_massDplus","Mass of all Dplus candidates; (GeV/c^{2})",50,1.75,2.00);
   histo_refitMassDplus = outfile->make<TH1D>("histo_refitMassDplus","Mass of all Dplus candidates (refitted tracks); (GeV/c^{2})",50,1.75,2.00);
   histo_DplusVertexNormChi2 = outfile->make<TH1D>("histo_DplusVertexNormChi2","Dplus candidate vertex chi2/dof",100,0,10);
   histo_DplusVertexX = outfile->make<TH1D>("histo_DplusVertexX","Dplus candidate vertex X position; (cm)",100,-1.0,1.0);
   histo_DplusVertexY = outfile->make<TH1D>("histo_DplusVertexY","Dplus candidate vertex Y position; (cm)",100,-1.0,1.0);
   histo_DplusVertexZ = outfile->make<TH1D>("histo_DplusVertexZ","Dplus candidate vertex Z position; (cm)",200,-20,20);
   histo_DplusVertexXError = outfile->make<TH1D>("histo_DplusVertexXError","Dplus candidate vertex X error; (cm)",100,0,1.0);
   histo_DplusVertexYError = outfile->make<TH1D>("histo_DplusVertexYError","Dplus candidate vertex Y error; (cm)",100,0,1.0);
   histo_DplusVertexZError = outfile->make<TH1D>("histo_DplusVertexZError","Dplus candidate vertex Z error; (cm)",100,0,1.0);
   histo_DplusVertexDxy = outfile->make<TH1D>("histo_DplusVertexDxy","Primary-secondary vertex XY displacement; (cm)",100,0,5.0);
   histo_DplusVertexExy = outfile->make<TH1D>("histo_DplusVertexExy","Primary-secondary vertex XY error; (cm)",100,0,1.0);
   histo_DplusVertexSxy = outfile->make<TH1D>("histo_DplusVertexSxy","Primary-secondary vertex XY significance ",100,0,10);
   histo_DplusVertexD3D = outfile->make<TH1D>("histo_DplusVertexD3D","Primary-secondary vertex 3D displacement; (cm)",100,0,5.0);
   histo_DplusVertexE3D = outfile->make<TH1D>("histo_DplusVertexE3D","Primary-secondary vertex 3D error; (cm)",100,0,1.0);
   histo_DplusVertexS3D = outfile->make<TH1D>("histo_DplusVertexS3D","Primary-secondary vertex 3D significance ",100,0,10);

   histo_angle = outfile->make<TH1D>("histo_angle","Momentum vs. displacement angle",330,-0.1,3.2);
   histo_refittedAngle = outfile->make<TH1D>("histo_refittedAngle","Momentum vs. displacement angle",330,-0.1,3.2);
   histo_matchedAngle = outfile->make<TH1D>("histo_matchedAngle","Momentum vs. displacement angle",100,0.,0.20);
   histo_matchedRefittedAngle = outfile->make<TH1D>("histo_matchedRefittedAngle","Momentum vs. displacement angle",100,0.,0.20);
   histo_zoomedAngle = outfile->make<TH1D>("histo_zoomedAngle","Momentum vs. displacement angle",100,0.,0.20);
   histo_zoomedRefittedAngle = outfile->make<TH1D>("histo_zoomedRefittedAngle","Momentum vs. displacement angle",100,0.,0.20);

   //  For fake combinations (all same charge)

   histo_nFakeDplus = outfile->make<TH1D>("histo_nFakeDplus","Multiplicity of fake D+ candidates per event",30,0,30);
   histo_massFakeDplus = outfile->make<TH1D>("histo_massFakeDplus","Mass of all fake D+ candidates; (GeV/c^{2})",50,1.75,2.00);
   histo_refitMassFakeDplus = outfile->make<TH1D>("histo_refitMassFakeDplus","Mass of all fake D+ candidates (refitted tracks); (GeV/c^{2})",50,1.75,2.00);
   histo_FakeDplusVertexNormChi2 = outfile->make<TH1D>("histo_FakeDplusVertexNormChi2","Fake D+ candidate vertex chi2/dof",100,0,10);
   histo_FakeDplusVertexX = outfile->make<TH1D>("histo_FakeDplusVertexX","Fake D+ candidate vertex X position; (cm)",100,-1.0,1.0);
   histo_FakeDplusVertexY = outfile->make<TH1D>("histo_FakeDplusVertexY","Fake D+ candidate vertex Y position; (cm)",100,-1.0,1.0);
   histo_FakeDplusVertexZ = outfile->make<TH1D>("histo_FakeDplusVertexZ","Fake D+ candidate vertex Z position; (cm)",200,-20,20);
   histo_FakeDplusVertexXError = outfile->make<TH1D>("histo_FakeDplusVertexXError","Fake D+ candidate vertex X error; (cm)",100,0,1.0);
   histo_FakeDplusVertexYError = outfile->make<TH1D>("histo_FakeDplusVertexYError","Fake D+ candidate vertex Y error; (cm)",100,0,1.0);
   histo_FakeDplusVertexZError = outfile->make<TH1D>("histo_FakeDplusVertexZError","Fake D+ candidate vertex Z error; (cm)",100,0,1.0);
   histo_FakeDplusVertexDxy = outfile->make<TH1D>("histo_FakeDplusVertexDxy","Primary-secondary vertex XY displacement; (cm)",100,0,5.0);
   histo_FakeDplusVertexExy = outfile->make<TH1D>("histo_FakeDplusVertexExy","Primary-secondary vertex XY error; (cm)",100,0,1.0);
   histo_FakeDplusVertexSxy = outfile->make<TH1D>("histo_FakeDplusVertexSxy","Primary-secondary vertex XY significance ",100,0,10);
   histo_FakeDplusVertexD3D = outfile->make<TH1D>("histo_FakeDplusVertexD3D","Primary-secondary vertex 3D displacement; (cm)",100,0,5.0);
   histo_FakeDplusVertexE3D = outfile->make<TH1D>("histo_FakeDplusVertexE3D","Primary-secondary vertex 3D error; (cm)",100,0,1.0);
   histo_FakeDplusVertexS3D = outfile->make<TH1D>("histo_FakeDplusVertexS3D","Primary-secondary vertex 3D significance ",100,0,10);

   histo_FakeAngle = outfile->make<TH1D>("histo_FakeAngle","Momentum vs. displacement angle",330,-0.1,3.2);
   histo_FakeRefittedAngle = outfile->make<TH1D>("histo_FakeRefittedAngle","Momentum vs. displacement angle",330,-0.1,3.2);
   histo_FakeMatchedAngle = outfile->make<TH1D>("histo_FakeMatchedAngle","Momentum vs. displacement angle",100,0.,0.20);
   histo_FakeMatchedRefittedAngle = outfile->make<TH1D>("histo_FakeMatchedRefittedAngle","Momentum vs. displacement angle",100,0.,0.20);
   histo_FakeZoomedAngle = outfile->make<TH1D>("histo_FakeZoomedAngle","Momentum vs. displacement angle",100,0.,0.20);
   histo_FakeZoomedRefittedAngle = outfile->make<TH1D>("histo_FakeZoomedRefittedAngle","Momentum vs. displacement angle",100,0.,0.20);
   
  }

   //  GEN 

  if(DoGEN_)
  { 
   
   histo_DplusGenVertexX = outfile->make<TH1D>("histo_DplusGenVertexX","Generated vertex X position; (cm)",100,-1.0,1.0);
   histo_DplusGenVertexY = outfile->make<TH1D>("histo_DplusGenVertexY","Generated vertex Y position; (cm)",100,-1.0,1.0);
   histo_DplusGenVertexZ = outfile->make<TH1D>("histo_DplusGenVertexZ","Generated vertex Z position; (cm)",200,-20,20);
   histo_DplusDecVertexX = outfile->make<TH1D>("histo_DplusDecVertexX","Generated decay vertex X position; (cm)",100,-1.0,1.0);
   histo_DplusDecVertexY = outfile->make<TH1D>("histo_DplusDecVertexY","Generated decay vertex Y position; (cm)",100,-1.0,1.0);
   histo_DplusDecVertexZ = outfile->make<TH1D>("histo_DplusDecVertexZ","Generated decay vertex Z position; (cm)",200,-20,20);

   histo_nDplusGen = outfile->make<TH1D>("histo_nDplusGen","Multiplicity of D+ generated per event",10,0,10);
   histo_nD0Gen = outfile->make<TH1D>("histo_nD0Gen","Multiplicity of D0 generated per event",10,0,10);
   histo_nDsGen = outfile->make<TH1D>("histo_nDsGen","Multiplicity of Ds generated per event",10,0,10);
   histo_DplusGenPt = outfile->make<TH1D>("histo_DplusGenPt","Generated pT; (GeV/c)",100,0,50.);
   histo_DplusGenP = outfile->make<TH1D>("histo_DplusGenP","Generated momentum; (GeV/c)",100,0,50.);
   histo_DplusGenPhi = outfile->make<TH1D>("histo_DplusGenPhi","Generated phi",100,-3.1416,3.1416);
   histo_DplusGenEta = outfile->make<TH1D>("histo_DplusGenEta","Generated eta",100,-5.0,5.0);
//   histo_DplusGenMass = outfile->make<TH1D>("histo_DplusGenMass","Generated mass",100,1.8,1.9);
//   histo_DplusGenCharge = outfile->make<TH1D>("histo_DplusGenCharge","Generated charge",4,-2,2);
   histo_DplusGenLengthXY = outfile->make<TH1D>("histo_DplusGenLengthXY","Generated track XY length; (cm)",100,0.0,1.0);
   histo_DplusGenLength = outfile->make<TH1D>("histo_DplusGenLength","Generated track length; (cm)",100,0.0,1.0);

   histo_D2KpipiGenPt = outfile->make<TH1D>("histo_D2KpipiGenPt","Generated pT of D^{+} -> K#pi#pi; (GeV/c)",100,0,50.);
   histo_D2KpipiGenP = outfile->make<TH1D>("histo_D2KpipiGenP","Generated momentum of D^{+} -> K#pi#pi; (GeV/c)",100,0,50.);
   histo_D2KpipiGenPhi = outfile->make<TH1D>("histo_D2KpipiGenPhi","Generated phi of D^{+} -> K#pi#pi",100,-3.1416,3.1416);
   histo_D2KpipiGenEta = outfile->make<TH1D>("histo_D2KpipiGenEta","Generated eta of D^{+} -> K#pi#pi",100,-5.0,5.0);
   
   histo_KGenPt = outfile->make<TH1D>("histo_KGenPt","pT of K from decay of generated Dplus; (GeV/c)",100,0,50);
   histo_KGenP = outfile->make<TH1D>("histo_KGenP","Momentum of K from decay of generated Dplus; (GeV/c)",100,0,50);
//   histo_ktGen = outfile->make<TH1D>("histo_ktGen","K/#pi pT relative to generated Dplus",50,0,1);
//   histo_KDplusGenAngle = outfile->make<TH1D>("histo_KDplusGenAngle","K/Dplus opening angle",50,0,1);
   histo_piGenPt = outfile->make<TH1D>("histo_piGenPt","pT of #pi from decay of generated Dplus; (GeV/c)",100,0,50);
   histo_piGenP = outfile->make<TH1D>("histo_piGenP","Momentum of #pi from decay of generated Dplus; (GeV/c)",100,0,50);
//   histo_piGenPerp = outfile->make<TH1D>("histo_piGenPerp","#pi pT relative to generated Dplus",50,0,1);
//   histo_piDplusGenAngle = outfile->make<TH1D>("histo_piDplusGenAngle","#pi/Dplus opening angle",50,0,1);
//   histo_KpiGenAngle = outfile->make<TH1D>("histo_KpiGenAngle","#pi/K opening angle",50,0,1);
   histo_KpiGenPtsum = outfile->make<TH1D>("histo_KpiGenPtsum","Sum of pT of K and #pi from decay of generated Dplus; (GeV/c)",100,0,50);
   
  }

   //  GEN 

  if(DoGEN_ && DoRECO_)
  { 

   histo_angleGenReco = outfile->make<TH1D>("histo_angleGenReco","Gen-Reco momentum angle",100,0.,0.10);
   histo_angleGenReco2 = outfile->make<TH1D>("histo_angleGenReco2","Gen-Reco momentum angle",100,0.,0.005);
   histo_angleGenReco3 = outfile->make<TH1D>("histo_angleGenReco3","Gen-Reco momentum angle",100,0.,3.2);
   histo_DplusGenRecoDR = outfile->make<TH1D>("histo_DplusGenRecoDR","Distance between RECO and GEN vertices; (cm)",100,0.0,1.0);
   histo_angleFakeGenReco = outfile->make<TH1D>("histo_angleFakeGenReco","Gen-Reco momentum angle",100,0.,0.10);
   histo_angleFakeGenReco2 = outfile->make<TH1D>("histo_angleFakeGenReco2","Gen-Reco momentum angle",100,0.,0.005);
   histo_angleFakeGenReco3 = outfile->make<TH1D>("histo_angleFakeGenReco3","Gen-Reco momentum angle",100,0.,3.2);
   
  }
   
   if(DoGEN_)  genTree = outfile->make<TTree>("genTree","Dplus -> Kpipi") ;
   if(DoRECO_) { 
      theTree = outfile->make<TTree>("theTree","Dplus candidates") ;
      fakeTree = outfile->make<TTree>("fakeTree","Same sign combinations") ;
   }
   
  // Diffractive analysis
  
  resetPFThresholds(thresholdsPFlowBarrel_);
  resetPFThresholds(thresholdsPFlowEndcap_);
  resetPFThresholds(thresholdsPFlowTransition_);
  resetPFThresholds(thresholdsPFlowForward_);
  
  if(iConfig.exists("PFlowThresholds")){ 
     edm::ParameterSet const& thresholdsPFPset = iConfig.getParameter<edm::ParameterSet>("PFlowThresholds");
 
     edm::ParameterSet const& thresholdsBarrel = thresholdsPFPset.getParameter<edm::ParameterSet>("Barrel");
     edm::ParameterSet const& thresholdsEndcap = thresholdsPFPset.getParameter<edm::ParameterSet>("Endcap");
     edm::ParameterSet const& thresholdsTransition = thresholdsPFPset.getParameter<edm::ParameterSet>("Transition");
     edm::ParameterSet const& thresholdsForward = thresholdsPFPset.getParameter<edm::ParameterSet>("Forward");


     setPFThresholds(thresholdsPFlowBarrel_,thresholdsBarrel);
     setPFThresholds(thresholdsPFlowEndcap_,thresholdsEndcap);
     setPFThresholds(thresholdsPFlowTransition_,thresholdsTransition);
     setPFThresholds(thresholdsPFlowForward_,thresholdsForward);
  }
  
  thresholdsPFlow_[Barrel] = thresholdsPFlowBarrel_;
  thresholdsPFlow_[Endcap] = thresholdsPFlowEndcap_; 
  thresholdsPFlow_[Transition] = thresholdsPFlowTransition_;
  thresholdsPFlow_[Forward] = thresholdsPFlowForward_;


  std::map<int,std::pair<double,double> >::const_iterator pfThreshold = thresholdsPFlowBarrel_.begin();
  std::map<int,std::pair<double,double> >::const_iterator pfThresholds_end = thresholdsPFlowBarrel_.end(); 
  
  // std::ostringstream oss;
  // oss << "Using the following PF thresholds:\n";
  edm::LogInfo("Dmeson") << "Using the following PF thresholds:\n";
  for(; pfThreshold != pfThresholds_end; ++pfThreshold){
     int key = pfThreshold->first;    
     // oss << "  " << key << ": "
     edm::LogInfo("Dmeson") << "  " << key << ": "
                 << "(" << thresholdsPFlow_[Barrel][key].first
                 << "," << thresholdsPFlow_[Barrel][key].second << ")  "
                 << "(" << thresholdsPFlow_[Endcap][key].first
                 << "," << thresholdsPFlow_[Endcap][key].second << ")  "
                 << "(" << thresholdsPFlow_[Transition][key].first
                 << "," << thresholdsPFlow_[Transition][key].second << ")  "
                 << "(" << thresholdsPFlow_[Forward][key].first
                 << "," << thresholdsPFlow_[Forward][key].second << ")\n";   
  }
  // LogDebug("Analysis") << oss.str();

  Ebeam_ = comEnergy_/2.;
  if(applyEnergyScaleHCAL_) energyScaleHCAL_ = iConfig.getParameter<double>("EnergyScaleFactorHCAL");

}


D2Kpipi_Analyzer::~D2Kpipi_Analyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
D2Kpipi_Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
  // Get primary vertices here to be able to filter events based on vertex multiplicity 
  Handle<reco::VertexCollection> vc;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",vc);
  
  //  Get run and event number
  _run = iEvent.eventAuxiliary().run() ;
  _event = iEvent.eventAuxiliary().event() ;
   
  bool processEvent = true ;
  //  Select events with a single primary vertex only 
  if(DoSinglePV_ && vc->size()!=1) processEvent = false ;
  
  if(processEvent) {
   
// --------------     GENERATOR LEVEL DISTRIBUTIONS   ----------------- //
   
      if(DoGEN_) analyzeGen(iEvent, iSetup) ;

// --------------     RECO LEVEL DISTRIBUTIONS   ----------------- //

      if(DoRECO_) 
      {
         initialize();
         TriggerInfo(iEvent,iSetup);
         Gapside(iEvent,iSetup);
         // etaMaxMin(iEvent,iSetup,thresholdsPFlow_);    // removed due to redundancy 
         sumEHFTowers(iEvent,iSetup);
         analyzeReco(iEvent, iSetup) ;
//         theTree->Fill();
      }
   
  }  // end of "if(processEvent)" call

}


// ------------ method called once each job just before starting event loop  ------------
void 
D2Kpipi_Analyzer::beginJob()
{
   // Create and open file
//   outfile = new TFile("outfile.root","RECREATE","Track Multiplicity");
//   if( !outfile->IsOpen() ) edm::LogInfo("Dmeson") << "Error opening file." ;
   
   // Create histograms
//   multihisto = new TH1D("multihisto","Track multiplicity",200,0.,200.);

  if(DoGEN_)
  { 
   genTree->Branch("run",&_run,"run/I") ;
   genTree->Branch("event",&_event,"event/I") ;
   genTree->Branch("procId",&_procId,"procId/I");
   genTree->Branch("D_x",&_D_x) ;
   genTree->Branch("D_y",&_D_y) ;
   genTree->Branch("D_z",&_D_z) ;
   genTree->Branch("D_xy_length",&_D_xy_length) ;
   genTree->Branch("D_length",&_D_length) ;
   genTree->Branch("D2Kpipi_p",&_D2Kpipi_p) ;
   genTree->Branch("D2Kpipi_pT",&_D2Kpipi_pT) ;
   genTree->Branch("D2Kpipi_phi",&_D2Kpipi_phi) ;
   genTree->Branch("D2Kpipi_eta",&_D2Kpipi_eta) ;
   genTree->Branch("D2Kpipi_pTfraction",&_D2Kpipi_pTfraction) ;
   genTree->Branch("D2Kpipi_pT2fraction",&_D2Kpipi_pT2fraction) ;
   genTree->Branch("mean_pT",&_mean_pT) ;
   genTree->Branch("mean_pT2",&_mean_pT2) ;
   genTree->Branch("nGenPart",&_nGenPart,"nGenPart/I") ;
   genTree->Branch("K_p",&_K_p) ;
   genTree->Branch("K_pT",&_K_pT) ;
   genTree->Branch("K_phi",&_K_phi) ;
   genTree->Branch("K_eta",&_K_eta) ;
   genTree->Branch("pi1_p",&_pi1_p) ;
   genTree->Branch("pi1_pT",&_pi1_pT) ;
   genTree->Branch("pi1_phi",&_pi1_phi) ;
   genTree->Branch("pi1_eta",&_pi1_eta) ;
   genTree->Branch("pi2_p",&_pi2_p) ;
   genTree->Branch("pi2_pT",&_pi2_pT) ;
   genTree->Branch("pi2_phi",&_pi2_phi) ;
   genTree->Branch("pi2_eta",&_pi2_eta) ;
//   genTree->Branch("kt",&_kt) ;
//   genTree->Branch("K2D_angle",&_K2D_angle) ;
//   genTree->Branch("pi2D_angle",&_pi2D_angle) ;
  }

   
  if(DoRECO_)
  { 
   //  Opposite sign (good) combinations
   theTree->Branch("run",&_run,"run/I") ;
   theTree->Branch("event",&_event,"event/I") ;
   theTree->Branch("hlt_bscOR_bptxOR",&_HLTPath,"hlt_bscOR_bptxOR/I");
   if(DoGEN_) theTree->Branch("procId",&_procId,"procId/I");
   if(DoGEN_) theTree->Branch("hasdplus",&_hasdplus,"hasdplus/O") ;
   if(DoGEN_) theTree->Branch("hasdplus2kpipi",&_hasdplus2kpipi,"hasdplus2kpipi/O") ;
   if(DoGEN_) theTree->Branch("hasd0",&_hasd0,"hasd0/O") ;
   if(DoGEN_) theTree->Branch("hasds",&_hasds,"hasds/O") ;
   if(DoGEN_) theTree->Branch("hasdstarplus",&_hasdstarplus,"hasdstarplus/O") ;
   // theTree->Branch("massDplus",&_massDplus) ;
   theTree->Branch("refittedMassDplus",&_refittedMassDplus) ;
   theTree->Branch("refittedMassKpi1",&_refittedMassKpi1) ;
   theTree->Branch("refittedMassKpi2",&_refittedMassKpi2) ;
   theTree->Branch("nvtx",&_nvtx,"nvtx/I") ;
   theTree->Branch("vtxchi2dof",&_vtxchi2dof) ;
   theTree->Branch("Dxy",&_Dxy) ;
   theTree->Branch("Sxy",&_Sxy) ;
   theTree->Branch("D3D",&_D3D) ;
   theTree->Branch("S3D",&_S3D) ;
   theTree->Branch("Rxy",&_Rxy) ;
   // theTree->Branch("angle",&_angle) ;
   theTree->Branch("refittedAngle",&_refittedAngle) ;
   theTree->Branch("dplusp",&_dplusp) ;
   theTree->Branch("dpluspT",&_dpluspT) ;
   theTree->Branch("dplusphi",&_dplusphi) ;
   theTree->Branch("dpluseta",&_dpluseta) ;
   theTree->Branch("trk1dcaxy",&_trk1dcaXY) ;
   theTree->Branch("trk2dcaxy",&_trk2dcaXY) ;
   theTree->Branch("trk3dcaxy",&_trk3dcaXY) ;
   theTree->Branch("trk1dcaz",&_trk1dcaZ) ;
   theTree->Branch("trk2dcaz",&_trk2dcaZ) ;
   theTree->Branch("trk3dcaz",&_trk3dcaZ) ;
   theTree->Branch("trk1p",&_trk1p) ;
   theTree->Branch("trk2p",&_trk2p) ;
   theTree->Branch("trk3p",&_trk3p) ;
   theTree->Branch("trk1pT",&_trk1pT) ;
   theTree->Branch("trk2pT",&_trk2pT) ;
   theTree->Branch("trk3pT",&_trk3pT) ;
   theTree->Branch("trk1chi2",&_trk1chi2) ;
   theTree->Branch("trk2chi2",&_trk2chi2) ;
   theTree->Branch("trk3chi2",&_trk3chi2) ;
   theTree->Branch("trk1hits",&_trk1hits,"trk1hits/I") ;
   theTree->Branch("trk2hits",&_trk2hits,"trk1hits/I") ;
   theTree->Branch("trk3hits",&_trk3hits,"trk1hits/I") ;
//   theTree->Branch("trk1angle",&_trk1angle) ;
//   theTree->Branch("trk2angle",&_trk2angle) ;
//   theTree->Branch("trkperp",&_trkperp) ;
   // Diffractive analysis
   theTree->Branch("nHFPlus",&_nHFPlus,"nHFPlus/I");
   theTree->Branch("nHFMinus",&_nHFMinus,"nHFMinus/I");
   theTree->Branch("HFPlusEnergyFromCTowers",&_HFPlusEnergyFromCaloTowers);
   theTree->Branch("HFMinusEnergyFromCTowers",&_HFMinusEnergyFromCaloTowers);
   theTree->Branch("HFPlusEnergyFromPF",&_HFPlusEnergyFromParticleFlow);
   theTree->Branch("HFMinusEnergyFromPF",&_HFMinusEnergyFromParticleFlow);
   theTree->Branch("EPlusPzFromPF",&_EPlusPzFromPFCands);
   theTree->Branch("EMinusPzFromPF",&_EMinusPzFromPFCands);
   theTree->Branch("xiPlusFromPF",&_xiPlusFromPFCands);
   theTree->Branch("xiMinusFromPF",&_xiMinusFromPFCands);
   theTree->Branch("etaMaxFromPF",&_etaMaxFromPFCands);
   theTree->Branch("etaMinFromPF",&_etaMinFromPFCands);
   theTree->Branch("MxFromPF",&_MxFromPFCands); 
   // theTree->Branch("missingMassFromXiFromPF",&_missingMassFromXiFromPFCands);     // removed, gives same value as MxFromPF
   // theTree->Branch("pfsis1EtaMax",&_pfsis1Eta_max);
   // theTree->Branch("pfsis1EtaMin",&_pfsis1Eta_min);
   // theTree->Branch("pfsis2EtaMax",&_pfsis2Eta_max);
   // theTree->Branch("pfsis2EtaMin",&_pfsis2Eta_min);
   // theTree->Branch("deltaEtaPF",&_deltaEtapf);
   // theTree->Branch("zeroEnergyHFSideFromCTowers",&_zeroEnergyHFSideFromCaloTowers);

   //  Same sign (fake) combinations
   fakeTree->Branch("run",&_run,"run/I") ;
   fakeTree->Branch("event",&_event,"event/I") ;
   fakeTree->Branch("hlt_bscOR_bptxOR",&_HLTPath,"hlt_bscOR_bptxOR/I");
   if(DoGEN_) fakeTree->Branch("procId",&_procId,"procId/I");
   if(DoGEN_) fakeTree->Branch("hasdplus",&_hasdplus,"hasdplus/O") ;
   if(DoGEN_) fakeTree->Branch("hasdplus2kpipi",&_hasdplus2kpipi,"hasdplus2kpipi/O") ;
   if(DoGEN_) fakeTree->Branch("hasd0",&_hasd0,"hasd0/O") ;
   if(DoGEN_) fakeTree->Branch("hasds",&_hasds,"hasds/O") ;
   if(DoGEN_) fakeTree->Branch("hasdstarplus",&_hasdstarplus,"hasdstarplus/O") ;
   // fakeTree->Branch("massDplus",&_massDplus) ;
   fakeTree->Branch("refittedMassDplus",&_refittedMassDplus) ;
   fakeTree->Branch("refittedMassKpi1",&_refittedMassKpi1) ;
   fakeTree->Branch("refittedMassKpi2",&_refittedMassKpi2) ;
   fakeTree->Branch("nvtx",&_nvtx,"nvtx/I") ;
   fakeTree->Branch("vtxchi2dof",&_vtxchi2dof) ;
   fakeTree->Branch("Dxy",&_Dxy) ;
   fakeTree->Branch("Sxy",&_Sxy) ;
   fakeTree->Branch("D3D",&_D3D) ;
   fakeTree->Branch("S3D",&_S3D) ;
   fakeTree->Branch("Rxy",&_Rxy) ;
   // fakeTree->Branch("angle",&_angle) ;
   fakeTree->Branch("refittedAngle",&_refittedAngle) ;
   fakeTree->Branch("dplusp",&_dplusp) ;
   fakeTree->Branch("dpluspT",&_dpluspT) ;
   fakeTree->Branch("dplusphi",&_dplusphi) ;
   fakeTree->Branch("dpluseta",&_dpluseta) ;
   fakeTree->Branch("trk1dcaxy",&_trk1dcaXY) ;
   fakeTree->Branch("trk2dcaxy",&_trk2dcaXY) ;
   fakeTree->Branch("trk3dcaxy",&_trk3dcaXY) ;
   fakeTree->Branch("trk1dcaz",&_trk1dcaZ) ;
   fakeTree->Branch("trk2dcaz",&_trk2dcaZ) ;
   fakeTree->Branch("trk3dcaz",&_trk3dcaZ) ;
   fakeTree->Branch("trk1p",&_trk1p) ;
   fakeTree->Branch("trk2p",&_trk2p) ;
   fakeTree->Branch("trk3p",&_trk3p) ;
   fakeTree->Branch("trk1pT",&_trk1pT) ;
   fakeTree->Branch("trk2pT",&_trk2pT) ;
   fakeTree->Branch("trk3pT",&_trk3pT) ;
   fakeTree->Branch("trk1chi2",&_trk1chi2) ;
   fakeTree->Branch("trk2chi2",&_trk2chi2) ;
   fakeTree->Branch("trk3chi2",&_trk3chi2) ;
   fakeTree->Branch("trk1hits",&_trk1hits,"trk1hits/I") ;
   fakeTree->Branch("trk2hits",&_trk2hits,"trk2hits/I") ;
   fakeTree->Branch("trk3hits",&_trk3hits,"trk2hits/I") ;
//   fakeTree->Branch("trk1angle",&_trk1angle) ;
//   fakeTree->Branch("trk2angle",&_trk2angle) ;
//   fakeTree->Branch("trkperp",&_trkperp) ;
   // Diffractive analysis
   fakeTree->Branch("nHFPlus",&_nHFPlus,"nHFPlus/I");
   fakeTree->Branch("nHFMinus",&_nHFMinus,"nHFMinus/I");
   fakeTree->Branch("HFPlusEnergyFromCTowers",&_HFPlusEnergyFromCaloTowers);
   fakeTree->Branch("HFMinusEnergyFromCTowers",&_HFMinusEnergyFromCaloTowers);
   fakeTree->Branch("HFPlusEnergyFromPF",&_HFPlusEnergyFromParticleFlow);
   fakeTree->Branch("HFMinusEnergyFromPF",&_HFMinusEnergyFromParticleFlow);
   fakeTree->Branch("EPlusPzFromPF",&_EPlusPzFromPFCands);
   fakeTree->Branch("EMinusPzFromPF",&_EMinusPzFromPFCands);
   fakeTree->Branch("xiPlusFromPF",&_xiPlusFromPFCands);
   fakeTree->Branch("xiMinusFromPF",&_xiMinusFromPFCands);
   fakeTree->Branch("etaMaxFromPF",&_etaMaxFromPFCands);
   fakeTree->Branch("etaMinFromPF",&_etaMinFromPFCands);
   fakeTree->Branch("MxFromPF",&_MxFromPFCands); 
  }
   
   total_nD2KpipiGen = 0;
   total_nDpluscand = 0;
   total_nDpluscandFake = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
D2Kpipi_Analyzer::endJob() 
{
//   histo_nTrack->Write();
//   outfile->Close();
   if(DoGEN_) edm::LogInfo("Dmeson") << "Total number of generated Dplus->Kpipi = " << total_nD2KpipiGen ;
   if(DoRECO_) {
              edm::LogInfo("Dmeson") << "Total number of Dpluscand = " << total_nDpluscand ;
              edm::LogInfo("Dmeson") << "Total number of DpluscandFake = " << total_nDpluscandFake ;
   }
}

// ------------ method called when starting to processes a run  ------------
void 
D2Kpipi_Analyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
D2Kpipi_Analyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
D2Kpipi_Analyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
D2Kpipi_Analyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
D2Kpipi_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// =======================================================================================
void D2Kpipi_Analyzer::analyzeGen(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// =======================================================================================
{
    // Event information
    edm::Handle<edm::HepMCProduct> hepMCProduct;
    iEvent.getByLabel("generator",hepMCProduct);

    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByLabel("generator",genEventInfo);

    int processid = -1;
    
    if(genEventInfo.isValid()){
       processid = genEventInfo->signalProcessID();
    } else {
       processid = hepMCProduct->getHepMCData().signal_process_id();
    }

    _procId = processid;

   // Gen particles
   edm::Handle<reco::GenParticleCollection> genPC;
   iEvent.getByLabel("genParticles",genPC);

   _nGenPart = genPC->size() ;
   
   double ptSum = 0 ;
   double pt2Sum = 0 ;
   for (reco::GenParticleCollection::const_iterator ip = genPC->begin() ; ip != genPC->end() ; ip++) {
      ptSum += ip->pt() ;
      pt2Sum += pow(ip->pt(),2) ;
   } 

   //  Set to zero number D0 flags
   hasDplus = false ;
   hasDplus2Kpipi = false ;
   hasD0 = false ;
   hasDs = false ;
   hasDstarplus = false ;

   int nD0Gen = 0 ;  
   int nDplusGen = 0 ;  
   int nDsGen = 0 ;
   int nDstarplusGen = 0 ;  
   // std::vector<math::XYZPoint> DplusVertexPosition ; 
   // edm::LogInfo ("Dmeson") << "\n genPC->size = " << genPC->size() ;
   for (reco::GenParticleCollection::const_iterator ip = genPC->begin() ; ip != genPC->end() ; ip++) { 
      // const reco::GenParticle &part = (*ip) ;
      // int pID = part.pdgId();
      int pID = ip->pdgId();
      int absID = abs(pID) ;
      // if(absID == 411 || absID == 421 || absID == 413 || absID == 423) {
      //    edm::LogInfo ("Dmeson") << "\n Particle ID = " << pID 
      //                     << "\n Vertex = (" << ip->vx() << "," << ip->vy() << "," << ip->vz() << ")" 
      //                     << "\n pT = " << ip->pt() << ", eta = " << ip->eta() 
      //                     << ", phi = " << ip->phi() << ", mass = " << ip->mass()  ;
      // }
      // if(ip->status() == 3) edm::LogInfo ("Dmeson") << "\n Particle having ID " << pID 
      //                                        << " and code 3 with vertex coordinates (" 
      //                                        << ip->vx() << "," << ip->vy() << "," << ip->vz() << ")" ;
      // Select Dplus's and fill histograms
      if(absID == 411) {
         nDplusGen++ ;                                 // count Dplus
         hasDplus = true ;
         const reco::Candidate *d = ip->daughter(0);
         _D_xy_length = sqrt( pow(d->vx()-ip->vx(),2) + pow(d->vy()-ip->vy(),2) );
         histo_DplusGenLengthXY->Fill(_D_xy_length) ;
         _D_length = sqrt( pow(d->vx()-ip->vx(),2) + pow(d->vy()-ip->vy(),2) + pow(d->vz()-ip->vz(),2) );
         histo_DplusGenLength->Fill(_D_length) ;
         histo_DplusGenVertexX->Fill(ip->vx()) ;
         histo_DplusGenVertexY->Fill(ip->vy()) ;
         histo_DplusGenVertexZ->Fill(ip->vz()) ;
         histo_DplusDecVertexX->Fill(d->vx()) ;
         histo_DplusDecVertexY->Fill(d->vy()) ;
         histo_DplusDecVertexZ->Fill(d->vz()) ;
         _D_x = d->vx() ;
         _D_y = d->vy() ;
         _D_z = d->vz() ;
         histo_DplusGenPt->Fill(ip->pt()) ;
         histo_DplusGenP->Fill(ip->p()) ;
         histo_DplusGenPhi->Fill(ip->phi()) ;
         histo_DplusGenEta->Fill(ip->eta()) ;
//         histo_DplusGenMass->Fill(ip->mass()) ;
//         histo_DplusGenCharge->Fill(ip->charge()) ;
         // DplusVertexPosition.push_back(d->vertex()) ;
         if(isD2Kpipi(*ip)) {
            hasDplus2Kpipi = true ;
            total_nD2KpipiGen++ ;
            _D2Kpipi_pT = ip->pt() ;
            histo_D2KpipiGenPt->Fill(ip->pt()) ;
            if(ptSum > 0) _D2Kpipi_pTfraction = _D2Kpipi_pT / ptSum ;
            _mean_pT = ptSum / genPC->size() ;
            if(pt2Sum > 0) _D2Kpipi_pT2fraction = pow(_D2Kpipi_pT,2) / pt2Sum ;
            _mean_pT2 = pt2Sum / genPC->size() ;
            _D2Kpipi_p = ip->p() ;
            histo_D2KpipiGenP->Fill(ip->p()) ;
            _D2Kpipi_phi = ip->phi() ;
            histo_D2KpipiGenPhi->Fill(ip->phi()) ;
            _D2Kpipi_eta = ip->eta() ;
            histo_D2KpipiGenEta->Fill(ip->eta()) ;
            math::XYZVector dplus_mom, k_mom, pi1_mom, pi2_mom ;
            dplus_mom = ip->momentum() ;
            if( abs(ip->daughter(0)->pdgId()) == 321 ) {
               k_mom   = ip->daughter(0)->momentum() ;
               pi1_mom = ip->daughter(1)->momentum() ;
               pi2_mom = ip->daughter(2)->momentum() ;
            } else 
            if( abs(ip->daughter(1)->pdgId()) == 321 ) {
               pi1_mom = ip->daughter(0)->momentum() ;
               k_mom   = ip->daughter(1)->momentum() ;
               pi2_mom = ip->daughter(2)->momentum() ;
            } else 
            if( abs(ip->daughter(2)->pdgId()) == 321 ) {
               pi1_mom = ip->daughter(0)->momentum() ;
               pi2_mom = ip->daughter(1)->momentum() ;
               k_mom   = ip->daughter(2)->momentum() ;
            }
            // kaon p , pT
            _K_pT = sqrt( k_mom.Perp2() ) ;
            histo_KGenPt->Fill(_K_pT);
            _K_p = sqrt( k_mom.Mag2() ) ;
            histo_KGenP->Fill(_K_p);
            _K_phi = k_mom.Phi() ;
            _K_eta = k_mom.Eta() ;
            // First pion p , pT , eta , phi
            _pi1_pT = sqrt( pi1_mom.Perp2() ) ;
            histo_piGenPt->Fill(_pi1_pT);
            _pi1_p = sqrt( pi1_mom.Mag2() ) ;
            histo_piGenP->Fill(_pi1_p);
            _pi1_phi = pi1_mom.Phi() ;
            _pi1_eta = pi1_mom.Eta() ;
            // Second pion p , pT , eta , phi
            _pi2_pT = sqrt( pi2_mom.Perp2() ) ;
            histo_piGenPt->Fill(_pi2_pT);
            _pi2_p = sqrt( pi2_mom.Mag2() ) ;
            histo_piGenP->Fill(_pi2_p);
            _pi2_phi = pi2_mom.Phi() ;
            _pi2_eta = pi2_mom.Eta() ;
            histo_KpiGenPtsum->Fill(_K_pT+_pi1_pT+_pi2_pT);
            // transverse momentum of pion/kaon with respect to Dplus line of flight
//            _kt = sqrt(  k_mom.Cross(dplus_mom).Mag2() / dplus_mom.Mag2() ) ;
//            histo_ktGen->Fill(_kt);
//            _K2D_angle = asin( _kt / sqrt(  k_mom.Mag2() ) ) ;
//            histo_KDplusGenAngle->Fill(_K2D_angle);
//            _pi2D_angle = asin( _kt / sqrt( pi1_mom.Mag2() ) ) ;
//            histo_piDplusGenAngle->Fill(_pi2D_angle);
//            double angle_Kpi = acos( k_mom.Dot(pi1_mom) / ( sqrt(k_mom.Mag2()) * sqrt(pi1_mom.Mag2()) ) ) ;
//            histo_KpiGenAngle->Fill(angle_Kpi);
            
            genTree->Fill() ;

            /*
            edm::LogInfo ("Dmeson") // << "\n Dplus mom = (" << dplus_mom.x() << "," << dplus_mom.y() << "," << dplus_mom.z() << ")" 
                             << "\n Dplus mom = (" << sqrt(dplus_mom.Perp2()) << "," << dplus_mom.z() << ")" 
                             // << "\n  K mom = (" <<  k_mom.x() << "," <<  k_mom.y() << "," <<  k_mom.z() << ")" 
                             << "\n  K mom = (" <<  sqrt(k_mom.Perp2()) << "," <<  k_mom.z() << ")" 
                             // << "\n pi mom = (" << pi1_mom.x() << "," << pi1_mom.y() << "," << pi1_mom.z() << ")" 
                             << "\n pi mom = (" << sqrt(pi1_mom.Perp2()) << "," << pi1_mom.z() << ")" 
                             << "\n  _kt with respect to Dplus mom = " <<  _kt ;  
            */
         }
      }
      if(absID == 421) { nD0Gen++        ; hasD0 = true ;        }           // count D+
      if(absID == 431) { nDsGen++        ; hasDs = true ;        }           // count Ds
      if(absID == 413) { nDstarplusGen++ ; hasDstarplus = true ; }           // count Ds
   }
   // Fill multiplicity histos for D0, D+ and Ds
   histo_nD0Gen->Fill(nD0Gen) ;
   histo_nDplusGen->Fill(nDplusGen) ;
   histo_nDsGen->Fill(nDsGen) ;
}
   

// =======================================================================================
void D2Kpipi_Analyzer::analyzeReco(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// =======================================================================================
{
   // Get tracks in event
   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel("generalTracks",tracks);
   
   // Get primary vertices
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByLabel("offlinePrimaryVerticesWithBS",vertices);

   reco::Vertex PV ;
   
   // Fill histos
   
   histo_nTrack->Fill(tracks->size());

   if(vertices->size()>1) edm::LogInfo("Dmeson") << " Event with " << vertices->size() << " vertices." ;
   histo_nVertex->Fill(vertices->size());
   
   // Loop over list of primary vertices (PV)
   unsigned int nTracksFromPrimVert = 0;
   for (reco::VertexCollection::const_iterator iv = vertices->begin() ; iv != vertices->end() ; iv++) 
   {
      if(iv->isFake()) {
        edm::LogInfo("Dmeson") << " Vertex with " << iv->tracksSize() << " tracks in event with " 
                        << vertices->size() << " vertices and " << tracks->size() << " tracks is fake." ;  
      } else {
        PV = (*iv) ;
        nTracksFromPrimVert += iv->tracksSize() ;
        histo_nTracksPerVertex->Fill(iv->tracksSize());
        histo_PVNormChi2->Fill(iv->normalizedChi2());
        histo_PVX->Fill(iv->x());
        histo_PVY->Fill(iv->y());
        histo_PVZ->Fill(iv->z());
        histo_PVXError->Fill(iv->xError());
        histo_PVYError->Fill(iv->yError());
        histo_PVZError->Fill(iv->zError());
	// edm::LogInfo("Dmeson") <<   "  Primary vertex :" 
	//                 << "\n     Normalised chi2 : " << iv->normalizedChi2() 
	//                 << "\n     Position (x,y,z) :  (" << iv->x() << "," 
	//                 << iv->y() << "," << iv->z() << ")" 
	//                 << "\n     Transversal position XY : " 
	//                 << sqrt( iv->x()*iv->x() + iv->y()*iv->y() ) ;
        for (std::vector<reco::TrackBaseRef>::const_iterator it = iv->tracks_begin() ; 
                                                             it != iv->tracks_end() ; 
                                                             it++) {
           //  Distance of closest approach of reco'ed PV track to PV position
           histo_PVtrackDxy->Fill(it->get()->dxy(PV.position())) ;
           histo_PVtrackDz->Fill(it->get()->dz(PV.position())) ;
        }

      }
   }
   
   if(tracks->size()!=0) histo_fracTracksPV->Fill( (double)nTracksFromPrimVert/tracks->size() );
   
   // Get tracks collection and look for secondary vertices made of three tracks  
   // having invariant mass consistent with Dplus decaying to Kpipi 
   
   for (reco::TrackCollection::const_iterator it = tracks->begin() ; it != tracks->end() ; it++) {
      histo_trackPt->Fill(it->pt()) ;
      histo_trackP->Fill(it->p()) ;
      histo_trackHits->Fill(it->numberOfValidHits()) ;
      histo_trackNormChi2->Fill(it->normalizedChi2()) ;
      //  Distance of closest approach of reco'ed track to PV position
      histo_trackDxy->Fill(it->dxy(PV.position())) ;
      histo_trackDz->Fill(it->dz(PV.position())) ;
   }
   
   // Get the builder:
   edm::ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
   // Do the conversion:
   TTrackVec t_tracks = (*theB).build(tracks);
   
   //  Collections to hold transient tracks which match pre-selection cuts 
   TTrackVec t_good_tracks ;
   for(TTrackVec::const_iterator itt = t_tracks.begin() ; itt != t_tracks.end() ; itt++) {
      if ( itt->track().pt() > trackPtMin_ && 
           fabs(itt->track().dxy(PV.position()))>dcaXY_ && 
           fabs(itt->track().dz(PV.position()))>dcaZ_ && 
           fabs(itt->track().numberOfValidHits())>=trackHitsMin_ && 
           fabs(itt->track().normalizedChi2())<trackChi2Max_  )     { t_good_tracks.push_back(*itt) ; 
      }
   }
   
   // edm::LogInfo("Dmeson") << "  There are " << t_good_tracks.size() << " good tracks in this event." ;
   
   
   // Lorentz Vectors
   TLorentzVector vkaon  , vpion1 , vpion2, 
                  vdplus , vkpi1  , vkpi2 ;
   reco::Track track1, track2 ;
   int nDpluscand , nFakeDplus ;
   double massDplus , refittedMassDplus, massKpi1, massKpi2, refittedMassKpi1, refittedMassKpi2 ;

   TTrackVec fitTracks ;
   TransientVertex D_cand ;
   
   nDpluscand = 0 ;
   nFakeDplus = 0 ;
   
   // Loop over tracks looking for combinations of 3 charged tracks with invariant mass around Dplus mass 
   
   if(t_good_tracks.size()>2) {
     // First track
     for(TTrackVec::const_iterator it1 = t_good_tracks.begin() ; it1 != t_good_tracks.end()-2 ; it1++) {
       // Second track
       for(TTrackVec::const_iterator it2 = it1+1 ; it2 != t_good_tracks.end()-1 ; it2++) {
         // Third track
         for(TTrackVec::const_iterator it3 = it2+1 ; it3 != t_good_tracks.end() ; it3++) {
         
           if( it1->charge()*it2->charge()>0 ) {
             //  Reset combination flags and respective masses
             massDplus = 0 ;
	 
	     // Combination 1  ( pi+ pi+ K- )
	     vpion1.SetPtEtaPhiM(it1->track().pt(), it1->track().eta(), it1->track().phi(), MPION);
	     vpion2.SetPtEtaPhiM(it2->track().pt(), it2->track().eta(), it2->track().phi(), MPION);
	      vkaon.SetPtEtaPhiM(it3->track().pt(), it3->track().eta(), it3->track().phi(), MKAON);
	     vdplus = vkaon + vpion1 + vpion2 ;
	     vkpi1 = vkaon + vpion1 ;
	     vkpi2 = vkaon + vpion2 ;
	     // edm::LogInfo("Dmeson") << " Dplus candidate 1 momentum = (" << vdplus.Px() << ","
	     //                 << vdplus.Py() << "," << vdplus.Pz() << ")" ;
	     massDplus = vdplus.M() ;
	     massKpi1 = vkpi1.M() ;
	     massKpi2 = vkpi2.M() ;
	     // Perpendicular momentum of track1 or track2 with respect to their vector sum.
	     // Same value for Combination 2 as for Combination 1. 
//	     double trkperp = vkaon.Vect().Perp(vd0_2.Vect()) ;
	     // track1 or track2 angles with respect to their vector sum.
//	     double trk1angle = vpion1.Vect().Angle(vd0_2.Vect()) ;
//	     double trk2angle = vkaon.Vect().Angle(vd0_2.Vect()) ;
	 
	     if( fabs( massDplus-MDPLUS ) < DplusMassWindow_ ) {
	       fitTracks.clear();
	       fitTracks.push_back(*it1);
	       fitTracks.push_back(*it2);
	       fitTracks.push_back(*it3);
               KalmanVertexFitter theFitter(true) ;
               D_cand = theFitter.vertex(fitTracks) ;
               if(D_cand.isValid()) {
                 double angle = findAngle(PV,D_cand,vdplus) ;
                 // if(fabs(cos(angle)) > cos(angleMax_)) {
                 if(angle < angleMax_) {
                   /////////////////////////////////////
                   //  OPPOSITE SIGN combinations
	           //        ( pi+ pi+ K- )
                   /////////////////////////////////////
                   if(it1->charge()*it3->charge()<0) {
	             double chi2 = D_cand.normalisedChiSquared();
	             histo_DplusVertexNormChi2->Fill(chi2);
	             if(chi2 < vtxChi2Max_) {
	               // edm::LogInfo("Dmeson") << " Dplus candidate 1 vertex position = (" << D_cand.position().x() << ","
	               //                 << D_cand.position().y() << "," << D_cand.position().z() << ")" ;
	               // for(std::vector<math::XYZPoint>::const_iterator j = DplusVertexPosition.begin() ;
	               // j != DplusVertexPosition.end() ; j++) {
	               //    // Calculate the distance between Dplus reco'ed vertex candidate and generated vertex
	               //    double dx = D_cand.position().x() - j->x() ;
	               //    double dy = D_cand.position().y() - j->y() ;
	               //    double dz = D_cand.position().z() - j->z() ;
	               //    double distance = sqrt(dx*dx + dy*dy + dz*dz) ;
	               //    histo_DplusGenRecoDR->Fill(distance) ;
	               // }
	               // XY primary-secondary vertices distance
	               VertexDistanceXY vdXY ;
	               double dXY = vdXY.distance(PV,D_cand).value() ;
	               histo_DplusVertexDxy->Fill(dXY) ;
	               double eXY = vdXY.distance(PV,D_cand).error() ;
	               histo_DplusVertexExy->Fill(eXY) ;
	               double sXY = dXY / eXY ;
	               histo_DplusVertexSxy->Fill(sXY) ;
	               // 3D primary-secondary vertices distance
	               VertexDistance3D vd3D ;
	               double d3D = vd3D.distance(PV,D_cand).value() ;
	               if(cos(angle)<0) d3D = -d3D ;     //  signed displacement
	               histo_DplusVertexD3D->Fill(d3D) ;
	               double e3D = vd3D.distance(PV,D_cand).error() ;
	               histo_DplusVertexE3D->Fill(e3D) ;
	               double s3D = d3D / e3D ;
	               histo_DplusVertexS3D->Fill(s3D) ;
	               // Position relative to the beam line
	               double x = D_cand.position().x();
	               double y = D_cand.position().y();
	               // double z = D_cand.position().z();
	               double r = sqrt(x*x + y*y);
	               // Calculate angle between momentum and displacement vectors
//	               double angle = findAngle(PV,D_cand,vdplus) ;
	               histo_angle->Fill(angle) ;
	               histo_zoomedAngle->Fill(angle) ;
	          //   if(angle<0.1)
	          //   {
	               histo_DplusVertexX->Fill(D_cand.position().x());
	               histo_DplusVertexY->Fill(D_cand.position().y());
	               histo_DplusVertexZ->Fill(D_cand.position().z());
	               histo_DplusVertexXError->Fill(D_cand.positionError().cxx());
	               histo_DplusVertexYError->Fill(D_cand.positionError().cyy());
	               histo_DplusVertexZError->Fill(D_cand.positionError().czz());
	          //   }
	               double refittedAngle = -0.09 ;
	               if( r < rMax_ && sXY > vtxDispSigXYMin_) {
	                 if(D_cand.hasRefittedTracks()) {
	                   //  Recalculate the invariant masses with reffited tracks
	                   TTrackVec tracksRefitted = D_cand.refittedTracks() ;
	                   // Combination 1
	                   vpion1.SetPtEtaPhiM(tracksRefitted[0].track().pt(), tracksRefitted[0].track().eta(), tracksRefitted[0].track().phi(), MPION);
	                   vpion2.SetPtEtaPhiM(tracksRefitted[1].track().pt(), tracksRefitted[1].track().eta(), tracksRefitted[1].track().phi(), MPION);
	                    vkaon.SetPtEtaPhiM(tracksRefitted[2].track().pt(), tracksRefitted[2].track().eta(), tracksRefitted[2].track().phi(), MKAON);
	                   vdplus = vkaon + vpion1 + vpion2 ;
	                   refittedMassDplus = vdplus.M() ;
	                   vkpi1 = vkaon + vpion1 ;
	                   refittedMassKpi1 = vkpi1.M() ;
	                   vkpi2 = vkaon + vpion2 ;
	                   refittedMassKpi2 = vkpi2.M() ;
	                   // Calculate angle between momentum and displacement vectors
	                   refittedAngle = findAngle(PV,D_cand,vdplus) ;
	                   histo_refittedAngle->Fill(refittedAngle) ;
	                   histo_zoomedRefittedAngle->Fill(refittedAngle) ;
	                 } else {
	                   //  Take the same mass
	                   refittedMassDplus = massDplus ;
	                   refittedMassKpi1 = massKpi1 ;
	                   refittedMassKpi2 = massKpi2 ;
	                 }
	                 // Angle between p_GEN and p_RECO
	                 if(DoGEN_) {
	                   // Get p_GEN for Dplus->Kpipi
	                   edm::Handle<reco::GenParticleCollection> gen_pc;
	                   iEvent.getByLabel("genParticles",gen_pc);
	                   TVector3 v3_D2Kpipi ;
	                   int nD2Kpipigen=0 ; 
	                   double delta = -9 ;
	                   for (reco::GenParticleCollection::const_iterator ip = gen_pc->begin() ; ip != gen_pc->end() ; ip++) {
	                     // if(abs(ip->pdgId())==421) edm::LogInfo("Dmeson") << "\n Dplus found !!! " ;
	                     if( abs(ip->pdgId())==421 && isD2Kpipi(*ip) ) { 
	                       TVector3 v3_D2Kpipi(ip->px(),ip->py(),ip->pz()) ;
	                       delta = vdplus.Vect().Angle(v3_D2Kpipi) ;
	                       edm::LogInfo("Dmeson") << "\n Angle between p_gen and p_reco for Dplus->Kpi = " << delta ;
	                       nD2Kpipigen++ ;
	                     }
	                   }
	                   if(nD2Kpipigen==1) { 
	                     histo_angleGenReco->Fill( delta ) ;
	                     histo_angleGenReco2->Fill( delta ) ;
	                     histo_angleGenReco3->Fill( delta ) ;
	                     //  Match criteria for p_gen and p_reco 
	                     if(delta<0.01) {
	                       histo_matchedAngle->Fill(angle) ;
	                       histo_matchedRefittedAngle->Fill(refittedAngle) ;
	                     }
	                   }
	                 }
	                 nDpluscand++ ;
	                 histo_massDplus->Fill(massDplus) ;
	                 histo_refitMassDplus->Fill(refittedMassDplus) ;
	                   
	                 // edm::LogInfo("Dmeson") <<   "  D0 candidate vertex :" 
	                 //                 << "\n     Normalised chi2 : " << chi2 
	                 //                 << "\n     Position (x,y,z) :  (" << x << "," 
	                 //                 << y << "," << z << ")" 
	                 //                 << "\n     Transversal position XY : " << r ;

	                 _hasdplus = hasDplus ; 
	                 _hasdplus2kpipi = hasDplus2Kpipi ; 
	                 _hasd0 = hasD0 ; 
	                 _hasds = hasDs ; 
	                 _hasdstarplus = hasDstarplus ; 
	                 _massDplus = massDplus ;
	                 _refittedMassDplus = refittedMassDplus ;
	                 _refittedMassKpi1 = refittedMassKpi1 ;
	                 _refittedMassKpi2 = refittedMassKpi2 ;
	                 _nvtx = vertices->size() ;
	                 _vtxchi2dof = chi2 ;
	                 _Dxy = dXY ;
	                 _Sxy = sXY ;
	                 _D3D = d3D ;
	                 _S3D = s3D ;
	                 _Rxy = r ;
	                 _angle = angle ;
	                 _refittedAngle = refittedAngle ;
	                 _dplusp = vdplus.P() ;
	                 _dpluspT = vdplus.Pt() ;
	                 _dplusphi = vdplus.Phi() ;
	                 _dpluseta = vdplus.Eta() ;
	                 _trk1dcaXY = it1->track().dxy(PV.position()) ;
	                 _trk1dcaZ  = it1->track().dz(PV.position()) ;
	                 _trk1pT = it1->track().pt() ;
	                 _trk1p = it1->track().p() ;
	                 _trk1chi2 = it1->track().normalizedChi2() ;
	                 _trk1hits = it1->track().numberOfValidHits() ;
//	                 _trk1angle = trk1angle ;
	                 _trk2dcaXY = it2->track().dxy(PV.position()) ;
	                 _trk2dcaZ  = it2->track().dz(PV.position()) ;
	                 _trk2pT = it2->track().pt() ;
	                 _trk2p = it2->track().p() ;
	                 _trk2chi2 = it2->track().normalizedChi2() ;
	                 _trk2hits = it2->track().numberOfValidHits() ;
//	                 _trk2angle = trk2angle ;
//	                 _trkperp = trkperp ;
	                 _trk3dcaXY = it3->track().dxy(PV.position()) ;
	                 _trk3dcaZ  = it3->track().dz(PV.position()) ;
	                 _trk3pT = it3->track().pt() ;
	                 _trk3p = it3->track().p() ;
	                 _trk3chi2 = it3->track().normalizedChi2() ;
	                 _trk3hits = it3->track().numberOfValidHits() ;
	             
	                 theTree->Fill();
	 
	               }     // end XY distance and significance loop
	             }     // end chi2 loop
	           } else 
                   /////////////////////////////////////
                   //  SAME SIGN combinations
	           //        ( pi+ pi+ K+ )
                   /////////////////////////////////////
	           if(it1->charge()*it3->charge()>0) {
	             double chi2 = D_cand.normalisedChiSquared();
	             histo_FakeDplusVertexNormChi2->Fill(chi2);
	             if(chi2 < vtxChi2Max_) {
	               // XY primary-secondary vertices distance
	               VertexDistanceXY vdXY ;
	               double dXY = vdXY.distance(PV,D_cand).value() ;
	               histo_FakeDplusVertexDxy->Fill(dXY) ;
	               double eXY = vdXY.distance(PV,D_cand).error() ;
	               histo_FakeDplusVertexExy->Fill(eXY) ;
	               double sXY = dXY / eXY ;
	               histo_FakeDplusVertexSxy->Fill(sXY) ;
	               // 3D primary-secondary vertices distance
	               VertexDistance3D vd3D ;
	               double d3D = vd3D.distance(PV,D_cand).value() ;
	               if(cos(angle)<0) d3D = -d3D ;     //  signed displacement
	               histo_FakeDplusVertexD3D->Fill(d3D) ;
	               double e3D = vd3D.distance(PV,D_cand).error() ;
	               histo_FakeDplusVertexE3D->Fill(e3D) ;
	               double s3D = d3D / e3D ;
	               histo_FakeDplusVertexS3D->Fill(s3D) ;
	               // Position relative to the beam line
	               double x = D_cand.position().x();
	               double y = D_cand.position().y();
	               // double z = D_cand.position().z();
	               double r = sqrt(x*x + y*y);
	               // Calculate angle between momentum and displacement vectors
	               double angle = findAngle(PV,D_cand,vdplus) ;
	               histo_FakeAngle->Fill(angle) ;
	               histo_FakeZoomedAngle->Fill(angle) ;
	            // if(angle<0.1)
	            // {
	               histo_FakeDplusVertexX->Fill(D_cand.position().x());
	               histo_FakeDplusVertexY->Fill(D_cand.position().y());
	               histo_FakeDplusVertexZ->Fill(D_cand.position().z());
	               histo_FakeDplusVertexXError->Fill(D_cand.positionError().cxx());
	               histo_FakeDplusVertexYError->Fill(D_cand.positionError().cyy());
	               histo_FakeDplusVertexZError->Fill(D_cand.positionError().czz());
	            // }
	               double refittedAngle = -0.09 ;
	               if( r < rMax_ && sXY > vtxDispSigXYMin_) {
	                 if(D_cand.hasRefittedTracks()) {
	                   //  Recalculate the invariant masses with reffited tracks
	                   TTrackVec tracksRefitted = D_cand.refittedTracks() ;
	                   // Combination 1
	                    vkaon.SetPtEtaPhiM(tracksRefitted[0].track().pt(), tracksRefitted[0].track().eta(), tracksRefitted[0].track().phi(), MKAON);
	                   vpion1.SetPtEtaPhiM(tracksRefitted[1].track().pt(), tracksRefitted[1].track().eta(), tracksRefitted[1].track().phi(), MPION);
	                   vpion2.SetPtEtaPhiM(tracksRefitted[2].track().pt(), tracksRefitted[2].track().eta(), tracksRefitted[2].track().phi(), MPION);
	                   vdplus = vkaon + vpion1 + vpion2 ;
	                   refittedMassDplus = vdplus.M() ;
	                   vkpi1 = vkaon + vpion1 ;
	                   refittedMassKpi1 = vkpi1.M() ;
	                   vkpi2 = vkaon + vpion2 ;
	                   refittedMassKpi2 = vkpi2.M() ;
	                   // Calculate angle between momentum and displacement vectors
	                   refittedAngle = findAngle(PV,D_cand,vdplus) ;
	                   histo_FakeRefittedAngle->Fill(refittedAngle) ;
	                   histo_FakeZoomedRefittedAngle->Fill(refittedAngle) ;
	                 } else {
	                   //  Take the same mass
	                   refittedMassDplus = massDplus ;
	                   refittedMassKpi1 = massKpi1 ;
	                   refittedMassKpi2 = massKpi2 ;
	                 }
	                 // Angle between p_GEN and p_RECO
	                 if(DoGEN_) {
	                   // Get p_GEN for Dplus->Kpipi
	                   edm::Handle<reco::GenParticleCollection> gen_pc;
	                   iEvent.getByLabel("genParticles",gen_pc);
	                   TVector3 v3_D2Kpipi ;
	                   int nD2Kpipigen=0 ; 
	                   double delta = -9 ;
	                   for (reco::GenParticleCollection::const_iterator ip = gen_pc->begin() ; ip != gen_pc->end() ; ip++) {
	                     // if(abs(ip->pdgId())==421) edm::LogInfo("Dmeson") << "\n Dplus found !!! " ;
	                     if( abs(ip->pdgId())==421 && isD2Kpipi(*ip) ) { 
	                       TVector3 v3_D2Kpipi(ip->px(),ip->py(),ip->pz()) ;
	                       delta = vdplus.Vect().Angle(v3_D2Kpipi) ;
	                       edm::LogInfo("Dmeson") << "\n Angle between p_gen and p_reco for Dplus->Kpi = " << delta ;
	                       nD2Kpipigen++ ;
	                     }
	                   }
	                   if(nD2Kpipigen==1) { 
	                     histo_angleFakeGenReco->Fill( delta ) ;
	                     histo_angleFakeGenReco2->Fill( delta ) ;
	                     histo_angleFakeGenReco3->Fill( delta ) ;
	                     //  Match criteria for p_gen and p_reco 
	                     if(delta<0.01) {
	                       histo_FakeMatchedAngle->Fill(angle) ;
	                       histo_FakeMatchedRefittedAngle->Fill(refittedAngle) ;
	                     }
	                   }
	                 }
	                 nFakeDplus++ ;
	                 histo_massFakeDplus->Fill(massDplus) ;
	                 histo_refitMassFakeDplus->Fill(refittedMassDplus) ;
	                 // edm::LogInfo("Dmeson") <<   "  D0 candidate vertex :" 
	                 //                 << "\n     Normalised chi2 : " << chi2 
	                 //                 << "\n     Position (x,y,z) :  (" << x << "," 
	                 //                 << y << "," << z << ")" 
	                 //                 << "\n     Transversal position XY : " << r ;

	                 _hasdplus = hasDplus ; 
	                 _hasdplus2kpipi = hasDplus2Kpipi ; 
	                 _hasd0 = hasD0 ; 
	                 _hasds = hasDs ; 
	                 _hasdstarplus = hasDstarplus ; 
	                 _massDplus = massDplus ;
	                 _refittedMassDplus = refittedMassDplus ;
	                 _refittedMassKpi1 = refittedMassKpi1 ;
	                 _refittedMassKpi2 = refittedMassKpi2 ;
	                 _nvtx = vertices->size() ;
	                 _vtxchi2dof = chi2 ;
	                 _Dxy = dXY ;
	                 _Sxy = sXY ;
	                 _D3D = d3D ;
	                 _S3D = s3D ;
	                 _Rxy = r ;
	                 _angle = angle ;
	                 _refittedAngle = refittedAngle ;
	                 _dplusp = vdplus.P() ;
	                 _dpluspT = vdplus.Pt() ;
	                 _dplusphi = vdplus.Phi() ;
	                 _dpluseta = vdplus.Eta() ;
	                 _trk1dcaXY = it1->track().dxy(PV.position()) ;
	                 _trk1dcaZ  = it1->track().dz(PV.position()) ;
	                 _trk1pT = it1->track().pt() ;
	                 _trk1p = it1->track().p() ;
	                 _trk1chi2 = it1->track().normalizedChi2() ;
	                 _trk1hits = it1->track().numberOfValidHits() ;
//	                 _trk1angle = trk1angle ;
	                 _trk2dcaXY = it2->track().dxy(PV.position()) ;
	                 _trk2dcaZ  = it2->track().dz(PV.position()) ;
	                 _trk2pT = it2->track().pt() ;
	                 _trk2p = it2->track().p() ;
	                 _trk2chi2 = it2->track().normalizedChi2() ;
	                 _trk2hits = it2->track().numberOfValidHits() ;
//	                 _trk2angle = trk2angle ;
//	                 _trkperp = trkperp ;
	                 _trk3dcaXY = it3->track().dxy(PV.position()) ;
	                 _trk3dcaZ  = it3->track().dz(PV.position()) ;
	                 _trk3pT = it3->track().pt() ;
	                 _trk3p = it3->track().p() ;
	                 _trk3chi2 = it3->track().normalizedChi2() ;
	                 _trk3hits = it3->track().numberOfValidHits() ;
	             
	                 fakeTree->Fill();
	 
	               }     // end XY distance and significance loop
	             }     // end chi2 loop
	           }     // end charge if
	         }     //  end angle
	       }     // end D_cand.isValid() loop
	     }     // end "DplusMassWindow" loop
           }     // end same charge pions if
           else if( it1->charge()*it2->charge()<0 ) {
             //  Reset combination flags and respective masses
             massDplus = 0 ;
	 
	     if( it1->charge()*it3->charge()>0 ) {
	       // Combination 2  ( pi+ K- pi+ )
	       vpion1.SetPtEtaPhiM(it1->track().pt(), it1->track().eta(), it1->track().phi(), MPION);
	        vkaon.SetPtEtaPhiM(it2->track().pt(), it2->track().eta(), it2->track().phi(), MKAON);
	       vpion2.SetPtEtaPhiM(it3->track().pt(), it3->track().eta(), it3->track().phi(), MPION);
	     }
	     if( it2->charge()*it3->charge()>0 ) {
	       // Combination 3  ( K- pi+ pi+ )
	        vkaon.SetPtEtaPhiM(it1->track().pt(), it1->track().eta(), it1->track().phi(), MKAON);
	       vpion1.SetPtEtaPhiM(it2->track().pt(), it2->track().eta(), it2->track().phi(), MPION);
	       vpion2.SetPtEtaPhiM(it3->track().pt(), it3->track().eta(), it3->track().phi(), MPION);
	     }
	     vdplus = vkaon + vpion1 + vpion2 ;
	     vkpi1 = vkaon + vpion1 ;
	     vkpi2 = vkaon + vpion2 ;
	     // edm::LogInfo("Dmeson") << " Dplus candidate 1 momentum = (" << vdplus.Px() << ","
	     //                 << vdplus.Py() << "," << vdplus.Pz() << ")" ;
	     massDplus = vdplus.M() ;
	     massKpi1 = vkpi1.M() ;
	     massKpi2 = vkpi2.M() ;
	     // Perpendicular momentum of track1 or track2 with respect to their vector sum.
	     // Same value for Combination 2 as for Combination 1. 
//	     double trkperp = vkaon.Vect().Perp(vd0_2.Vect()) ;
	     // track1 or track2 angles with respect to their vector sum.
//	     double trk1angle = vpion1.Vect().Angle(vd0_2.Vect()) ;
//	     double trk2angle = vkaon.Vect().Angle(vd0_2.Vect()) ;
	 
	     if( fabs( massDplus-MDPLUS ) < DplusMassWindow_ ) {
	       fitTracks.clear();
	       fitTracks.push_back(*it1);
	       fitTracks.push_back(*it2);
	       fitTracks.push_back(*it3);
               KalmanVertexFitter theFitter(true) ;
               D_cand = theFitter.vertex(fitTracks) ;
               if(D_cand.isValid()) {
                 double angle = findAngle(PV,D_cand,vdplus) ;
                 // if(fabs(cos(angle)) > cos(angleMax_)) {
                 if(angle < angleMax_) {
	             double chi2 = D_cand.normalisedChiSquared();
	             histo_DplusVertexNormChi2->Fill(chi2);
	             if(chi2 < vtxChi2Max_) {
	               // edm::LogInfo("Dmeson") << " Dplus candidate 1 vertex position = (" << D_cand.position().x() << ","
	               //                 << D_cand.position().y() << "," << D_cand.position().z() << ")" ;
	               // for(std::vector<math::XYZPoint>::const_iterator j = DplusVertexPosition.begin() ;
	               // j != DplusVertexPosition.end() ; j++) {
	               //    // Calculate the distance between Dplus reco'ed vertex candidate and generated vertex
	               //    double dx = D_cand.position().x() - j->x() ;
	               //    double dy = D_cand.position().y() - j->y() ;
	               //    double dz = D_cand.position().z() - j->z() ;
	               //    double distance = sqrt(dx*dx + dy*dy + dz*dz) ;
	               //    histo_DplusGenRecoDR->Fill(distance) ;
	               // }
	               // XY primary-secondary vertices distance
	               VertexDistanceXY vdXY ;
	               double dXY = vdXY.distance(PV,D_cand).value() ;
	               histo_DplusVertexDxy->Fill(dXY) ;
	               double eXY = vdXY.distance(PV,D_cand).error() ;
	               histo_DplusVertexExy->Fill(eXY) ;
	               double sXY = dXY / eXY ;
	               histo_DplusVertexSxy->Fill(sXY) ;
	               // 3D primary-secondary vertices distance
	               VertexDistance3D vd3D ;
	               double d3D = vd3D.distance(PV,D_cand).value() ;
	               if(cos(angle)<0) d3D = -d3D ;     //  signed displacement
	               histo_DplusVertexD3D->Fill(d3D) ;
	               double e3D = vd3D.distance(PV,D_cand).error() ;
	               histo_DplusVertexE3D->Fill(e3D) ;
	               double s3D = d3D / e3D ;
	               histo_DplusVertexS3D->Fill(s3D) ;
	               // Position relative to the beam line
	               double x = D_cand.position().x();
	               double y = D_cand.position().y();
	               // double z = D_cand.position().z();
	               double r = sqrt(x*x + y*y);
	               // Calculate angle between momentum and displacement vectors
//	               double angle = findAngle(PV,D_cand,vdplus) ;
	               histo_angle->Fill(angle) ;
	               histo_zoomedAngle->Fill(angle) ;
	          //   if(angle<0.1)
	          //   {
	               histo_DplusVertexX->Fill(D_cand.position().x());
	               histo_DplusVertexY->Fill(D_cand.position().y());
	               histo_DplusVertexZ->Fill(D_cand.position().z());
	               histo_DplusVertexXError->Fill(D_cand.positionError().cxx());
	               histo_DplusVertexYError->Fill(D_cand.positionError().cyy());
	               histo_DplusVertexZError->Fill(D_cand.positionError().czz());
	          //   }
	               double refittedAngle = -0.09 ;
	               if( r < rMax_ && sXY > vtxDispSigXYMin_) {
	                 if(D_cand.hasRefittedTracks()) {
	                   //  Recalculate the invariant masses with reffited tracks
	                   TTrackVec tracksRefitted = D_cand.refittedTracks() ;
	                   // Combination 1
	                   vpion1.SetPtEtaPhiM(tracksRefitted[0].track().pt(), tracksRefitted[0].track().eta(), tracksRefitted[0].track().phi(), MPION);
	                   vpion2.SetPtEtaPhiM(tracksRefitted[1].track().pt(), tracksRefitted[1].track().eta(), tracksRefitted[1].track().phi(), MPION);
	                    vkaon.SetPtEtaPhiM(tracksRefitted[2].track().pt(), tracksRefitted[2].track().eta(), tracksRefitted[2].track().phi(), MKAON);
	                   vdplus = vkaon + vpion1 + vpion2 ;
	                   refittedMassDplus = vdplus.M() ;
	                   vkpi1 = vkaon + vpion1 ;
	                   refittedMassKpi1 = vkpi1.M() ;
	                   vkpi2 = vkaon + vpion2 ;
	                   refittedMassKpi2 = vkpi2.M() ;
	                   // Calculate angle between momentum and displacement vectors
	                   refittedAngle = findAngle(PV,D_cand,vdplus) ;
	                   histo_refittedAngle->Fill(refittedAngle) ;
	                   histo_zoomedRefittedAngle->Fill(refittedAngle) ;
	                 } else {
	                   //  Take the same mass
	                   refittedMassDplus = massDplus ;
	                   refittedMassKpi1 = massKpi1 ;
	                   refittedMassKpi2 = massKpi2 ;
	                 }
	                 // Angle between p_GEN and p_RECO
	                 if(DoGEN_) {
	                   // Get p_GEN for Dplus->Kpipi
	                   edm::Handle<reco::GenParticleCollection> gen_pc;
	                   iEvent.getByLabel("genParticles",gen_pc);
	                   TVector3 v3_D2Kpipi ;
	                   int nD2Kpipigen=0 ; 
	                   double delta = -9 ;
	                   for (reco::GenParticleCollection::const_iterator ip = gen_pc->begin() ; ip != gen_pc->end() ; ip++) {
	                     // if(abs(ip->pdgId())==421) edm::LogInfo("Dmeson") << "\n Dplus found !!! " ;
	                     if( abs(ip->pdgId())==421 && isD2Kpipi(*ip) ) { 
	                       TVector3 v3_D2Kpipi(ip->px(),ip->py(),ip->pz()) ;
	                       delta = vdplus.Vect().Angle(v3_D2Kpipi) ;
	                       edm::LogInfo("Dmeson") << "\n Angle between p_gen and p_reco for Dplus->Kpi = " << delta ;
	                       nD2Kpipigen++ ;
	                     }
	                   }
	                   if(nD2Kpipigen==1) { 
	                     histo_angleGenReco->Fill( delta ) ;
	                     histo_angleGenReco2->Fill( delta ) ;
	                     histo_angleGenReco3->Fill( delta ) ;
	                     //  Match criteria for p_gen and p_reco 
	                     if(delta<0.01) {
	                       histo_matchedAngle->Fill(angle) ;
	                       histo_matchedRefittedAngle->Fill(refittedAngle) ;
	                     }
	                   }
	                 }
	                 nDpluscand++ ;
	                 histo_massDplus->Fill(massDplus) ;
	                 histo_refitMassDplus->Fill(refittedMassDplus) ;
	                   
	                 // edm::LogInfo("Dmeson") <<   "  D0 candidate vertex :" 
	                 //                 << "\n     Normalised chi2 : " << chi2 
	                 //                 << "\n     Position (x,y,z) :  (" << x << "," 
	                 //                 << y << "," << z << ")" 
	                 //                 << "\n     Transversal position XY : " << r ;

	                 _hasdplus = hasDplus ; 
	                 _hasdplus2kpipi = hasDplus2Kpipi ; 
	                 _hasd0 = hasD0 ; 
	                 _hasds = hasDs ; 
	                 _hasdstarplus = hasDstarplus ; 
	                 _massDplus = massDplus ;
	                 _refittedMassDplus = refittedMassDplus ;
	                 _refittedMassKpi1 = refittedMassKpi1 ;
	                 _refittedMassKpi2 = refittedMassKpi2 ;
	                 _nvtx = vertices->size() ;
	                 _vtxchi2dof = chi2 ;
	                 _Dxy = dXY ;
	                 _Sxy = sXY ;
	                 _D3D = d3D ;
	                 _S3D = s3D ;
	                 _Rxy = r ;
	                 _angle = angle ;
	                 _refittedAngle = refittedAngle ;
	                 _dplusp = vdplus.P() ;
	                 _dpluspT = vdplus.Pt() ;
	                 _dplusphi = vdplus.Phi() ;
	                 _dpluseta = vdplus.Eta() ;
	                 _trk1dcaXY = it1->track().dxy(PV.position()) ;
	                 _trk1dcaZ  = it1->track().dz(PV.position()) ;
	                 _trk1pT = it1->track().pt() ;
	                 _trk1p = it1->track().p() ;
	                 _trk1chi2 = it1->track().normalizedChi2() ;
	                 _trk1hits = it1->track().numberOfValidHits() ;
//	                 _trk1angle = trk1angle ;
	                 _trk2dcaXY = it2->track().dxy(PV.position()) ;
	                 _trk2dcaZ  = it2->track().dz(PV.position()) ;
	                 _trk2pT = it2->track().pt() ;
	                 _trk2p = it2->track().p() ;
	                 _trk2chi2 = it2->track().normalizedChi2() ;
	                 _trk2hits = it2->track().numberOfValidHits() ;
//	                 _trk2angle = trk2angle ;
//	                 _trkperp = trkperp ;
	                 _trk3dcaXY = it3->track().dxy(PV.position()) ;
	                 _trk3dcaZ  = it3->track().dz(PV.position()) ;
	                 _trk3pT = it3->track().pt() ;
	                 _trk3p = it3->track().p() ;
	                 _trk3chi2 = it3->track().normalizedChi2() ;
	                 _trk3hits = it3->track().numberOfValidHits() ;
	             
	                 theTree->Fill();
	 
	               }     // end XY distance and significance loop
	             }     // end chi2 loop
	         }     //  end angle
	       }     // end D_cand.isValid() loop
	     }     // end "DplusMassWindow" loop
           }     // end same charge pions if
         }     // end third track loop 
       }     // end second track loop 
     }     // end first track loop 
   }
   
   histo_nDplus->Fill(nDpluscand) ;
   histo_nFakeDplus->Fill(nFakeDplus) ;
   
   total_nDpluscand += nDpluscand ;
   total_nDpluscandFake += nFakeDplus ;
   
}

// =========================================================
bool D2Kpipi_Analyzer::isD2Kpipi(const reco::GenParticle& p) 
// =========================================================
{
   bool d2kpipi = false;
   if(abs(p.pdgId()) == 411) {
      // edm::LogInfo("Dmeson") << " Dplus found." ;
      if(p.numberOfDaughters() == 3) {
         int id1 = p.daughter(0)->pdgId() ;
         int id2 = p.daughter(1)->pdgId() ;
         int id3 = p.daughter(2)->pdgId() ;
         if( ( abs(id1)==321 && abs(id2)==211 && abs(id3)==211 && (id1*id2)<0 && (id1*id3)<0 ) || 
             ( abs(id1)==211 && abs(id2)==321 && abs(id3)==211 && (id2*id1)<0 && (id2*id3)<0 ) || 
             ( abs(id1)==211 && abs(id2)==211 && abs(id3)==321 && (id3*id1)<0 && (id3*id2)<0 ) ) 
         { 
            d2kpipi = true ;
            // edm::LogInfo("Dmeson") <<   " ==================================================== " 
            //                      << "\n Candidate with 3 daughters of ids " <<  id1 << ", " <<  id2 << " and " << id3 
            //                      << "\n ==================================================== " ;
         }
      }
   } else { 
      // edm::LogInfo("Dmeson") << " Not a Dplus." ; 
   }
   return d2kpipi;
}


// =========================================================
double D2Kpipi_Analyzer::findAngle(const reco::Vertex& pv , const TransientVertex& sv , const TLorentzVector& v ) 
// =========================================================
{
   CLHEP::Hep3Vector displacement( sv.position().x() - pv.position().x() , 
                                   sv.position().y() - pv.position().y() , 
                                   sv.position().z() - pv.position().z() ) ; 
   CLHEP::Hep3Vector momentum( v.Px() , v.Py() , v.Pz() ) ; 
   return momentum.angle(displacement) ;
}


//=================================
void D2Kpipi_Analyzer::initialize()
//=================================
{
   // _zeroEnergyHFSideFromCaloTowers = false ;
   _nHFPlus = 0 ; 
   _nHFMinus = 0 ;
   _HFPlusEnergyFromCaloTowers=0. ; 
   _HFMinusEnergyFromCaloTowers=0. ;
   _HFPlusEnergyFromParticleFlow=0 ; 
   _HFMinusEnergyFromParticleFlow=0 ;
   _EPlusPzFromPFCands=0 ; 
   _EMinusPzFromPFCands=0 ;
   _xiPlusFromPFCands=0 ; 
   _xiMinusFromPFCands=0 ; 
   _etaMaxFromPFCands=0 ; 
   _etaMinFromPFCands=0 ; 
   _MxFromPFCands=0 ; 
   _missingMassFromXiFromPFCands=0 ;

   // _pfsis1Eta_max=0.; _pfsis1Eta_min=0.; _pfsis2Eta_max=0.; _pfsis2Eta_min=0.; _deltaEtapf=0.;
   
   _HLTPath=0 ; 

/*
     TTBit_32 = 0;
     TTBit_33 = 0;
     TTBit_34 = 0;

     TTBit_8 = 0;
     TTBit_9 = 0;
     TTBit_10 = 0;

   xiGenPlus_.clear(); xiGenMinus_.clear(); MxGen_.clear(); MxGenRange_.clear(); sumEnergyHEPlusGen_.clear(); 
   sumEnergyHEMinusGen_.clear(); sumEnergyHFPlusGen_.clear(); MxGenMinus_.clear();
   sumEnergyHFMinusGen_.clear(); etaMaxGen_.clear(); etaMinGen_.clear(); deltaEtaGen_.clear(); etaGapLow_.clear(); etaGapHigh_.clear(); MxGenPlus_.clear();

   triggers.clear();
*/
}


// =======================================================================================
void D2Kpipi_Analyzer::TriggerInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// =======================================================================================
{
   edm::Handle<edm::TriggerResults> triggerResults;
   iEvent.getByLabel(edm::InputTag("TriggerResults","","HLT"),triggerResults);

   if(triggerResults.isValid() && hltPathName_ != "")
   {
      const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
      // In case hltPathName_ is a pattern (e.g. HLT_Jet30U*)
      std::string hltPath;
      if( edm::is_glob(hltPathName_) ){
         std::vector< std::vector<std::string>::const_iterator > matches = edm::regexMatch(triggerNames.triggerNames(), hltPathName_);  
         if( matches.empty() ) throw cms::Exception("Configuration") << "Could not find any HLT path of type " << hltPathName_ << "\n";
         else if( matches.size() > 1) throw cms::Exception("Configuration") << "HLT path type " << hltPathName_ << " not unique\n";
         else hltPath = *(matches[0]);
      } else {
         hltPath = hltPathName_; 
      } 

      unsigned int idxHLT = triggerNames.triggerIndex(hltPath);

      if (idxHLT < triggerResults->size()) 
      { 
         _HLTPath = (triggerResults->wasrun(idxHLT) && triggerResults->accept(idxHLT)) ? 1 : 0 ; 
      }
      else 
      {
         edm::LogWarning("Analysis") << " Trigger index: " << idxHLT << " Trigger Results Size: " << triggerResults->size()  
                                     << " Trigger index  must be equal/more that the Trigger Results !! " ;
         _HLTPath = -1;
      }
   } else 
   {
      _HLTPath = -1;
   }
}


// =====================================================================================
void D2Kpipi_Analyzer::Gapside(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// =====================================================================================
{
  using namespace std;
  using namespace reco;
  using namespace edm;

  double energyScale = (applyEnergyScaleHCAL_) ? energyScaleHCAL_ : -1.;

  // Leave only PF-based variables
  Handle<reco::PFCandidateCollection> particleFlowCollectionH;
  iEvent.getByLabel(particleFlowTag_,particleFlowCollectionH);

  double MxFromPFCands = MassColl(*particleFlowCollectionH,thresholdsPFlow_);

  std::pair<double,double> xiFromPFCands = xi(*particleFlowCollectionH,Ebeam_,thresholdsPFlow_);

  std::pair<double,double> EPlusPzFromPFCands = EPlusPz(*particleFlowCollectionH,thresholdsPFlow_);

  std::pair<double,double> sumEHFFromPFCands = sumEHF(*particleFlowCollectionH,thresholdsPFlow_);

  double missingMassFromXiFromPFCands = 2*Ebeam_*sqrt(xiFromPFCands.first*xiFromPFCands.second);

  std::pair<double,double> etaMaxFromPFCands = etaMax(*particleFlowCollectionH,thresholdsPFlow_);


  _MxFromPFCands = MxFromPFCands ;
  _xiPlusFromPFCands = xiFromPFCands.first ;
  _xiMinusFromPFCands = xiFromPFCands.second ;
  _EPlusPzFromPFCands = EPlusPzFromPFCands.first ;
  _EMinusPzFromPFCands = EPlusPzFromPFCands.second ;
  _HFPlusEnergyFromParticleFlow = sumEHFFromPFCands.first ;
  _HFMinusEnergyFromParticleFlow = sumEHFFromPFCands.second ;
  _etaMaxFromPFCands = etaMaxFromPFCands.first ;
  _etaMinFromPFCands = etaMaxFromPFCands.second ;
  _missingMassFromXiFromPFCands = missingMassFromXiFromPFCands ;

}


// ===============================================================================================
void D2Kpipi_Analyzer::resetPFThresholds(std::map<int,std::pair<double,double> >& thresholdsPFlow)
// ===============================================================================================
{
  thresholdsPFlow[reco::PFCandidate::X] = std::make_pair(-1.,-1.);
  thresholdsPFlow[reco::PFCandidate::h] = std::make_pair(-1.,-1.);
  thresholdsPFlow[reco::PFCandidate::e] = std::make_pair(-1.,-1.);
  thresholdsPFlow[reco::PFCandidate::mu] = std::make_pair(-1.,-1.);
  thresholdsPFlow[reco::PFCandidate::gamma] = std::make_pair(-1.,-1.);
  thresholdsPFlow[reco::PFCandidate::h0] = std::make_pair(-1.,-1.);
  thresholdsPFlow[reco::PFCandidate::h_HF] = std::make_pair(-1.,-1.);
  thresholdsPFlow[reco::PFCandidate::egamma_HF] = std::make_pair(-1.,-1.);
}


// =============================================================================================
void D2Kpipi_Analyzer::setPFThresholds(std::map<int,std::pair<double,double> >& thresholdsPFlow, 
                                   edm::ParameterSet const& thresholdsPFPset)
// =============================================================================================
{
  std::vector<std::string> pfThresholdNames = thresholdsPFPset.getParameterNames();
  std::vector<std::string>::const_iterator param = pfThresholdNames.begin();
  std::vector<std::string>::const_iterator params_end = pfThresholdNames.end();
  for(; param != params_end; ++param){
     //reco::PFCandidate::ParticleType particleType = pflowId(*param);
     int particleType = pflowId(*param);
     if(particleType == -1)
        throw cms::Exception("Configuration") << "Parameter " << *param
                                              << " does not correspond to any particle type (PF)";

     edm::ParameterSet const& typePSet = thresholdsPFPset.getParameter<edm::ParameterSet>(*param);
     double ptThreshold = -1.;
     if(typePSet.exists("pt")) ptThreshold = typePSet.getParameter<double>("pt");
     double energyThreshold = -1.;
     if(typePSet.exists("energy")) energyThreshold = typePSet.getParameter<double>("energy");
     thresholdsPFlow[particleType].first = ptThreshold;
     thresholdsPFlow[particleType].second = energyThreshold;
  }
}


// ===================================================
int D2Kpipi_Analyzer::pflowId(std::string const& name)
// ===================================================
{
   // FIXME: The labels definition could go somewhere else
   std::vector<std::string> labels_X, labels_h, labels_e, labels_mu, labels_gamma, labels_h0, labels_h_HF, labels_egamma_HF;
   labels_X.push_back("X");
   labels_X.push_back("undefined");
   labels_h.push_back("h");
   labels_h.push_back("chargedHadron");
   labels_h.push_back("hadronCharged");
   labels_e.push_back("e");
   labels_e.push_back("electron");
   labels_mu.push_back("mu");
   labels_mu.push_back("muon");
   labels_gamma.push_back("gamma");
   labels_gamma.push_back("photon");
   labels_h0.push_back("h0");
   labels_h0.push_back("neutralHadron");
   labels_h0.push_back("hadronNeutral");
   labels_h_HF.push_back("h_HF");
   labels_h_HF.push_back("hadronHF");
   labels_egamma_HF.push_back("egamma_HF");
   labels_egamma_HF.push_back("emHF");
   // Find corresponding particle type   
   if( std::find(labels_X.begin(), labels_X.end(), name) != labels_X.end() )
      return reco::PFCandidate::X;
   else if( std::find(labels_h.begin(), labels_h.end(), name) != labels_h.end() )
      return reco::PFCandidate::h;
   else if( std::find(labels_e.begin(), labels_e.end(), name) != labels_e.end() )
      return reco::PFCandidate::e;
   else if( std::find(labels_mu.begin(), labels_mu.end(), name) != labels_mu.end() )
      return reco::PFCandidate::mu;
   else if( std::find(labels_gamma.begin(), labels_gamma.end(), name) != labels_gamma.end() ) 
      return reco::PFCandidate::gamma;
   else if( std::find(labels_h0.begin(), labels_h0.end(), name) != labels_h0.end() ) 
      return reco::PFCandidate::h0;
   else if( std::find(labels_h_HF.begin(), labels_h_HF.end(), name) != labels_h_HF.end() ) 
      return reco::PFCandidate::h_HF;
   else if( std::find(labels_egamma_HF.begin(), labels_egamma_HF.end(), name) != labels_egamma_HF.end() ) 
      return reco::PFCandidate::egamma_HF;
   else return -1;
}

// =================================================================================================
bool D2Kpipi_Analyzer::pflowThreshold(reco::PFCandidate const& part, std::map<int, 
                                     std::map<int,std::pair<double,double> > > const& thresholdMap)
// =================================================================================================
{
   //FIXME
   // HF eta rings 29, 40, 41
   if( ( (fabs(part.eta()) >= 2.866) && (fabs(part.eta()) < 2.976) ) || 
         (fabs(part.eta()) >= 4.730) ) return false;

   bool accept = true;

   double eta = part.eta();
   int region = -1;
   if( (fabs(eta) >= 0.) && (fabs(eta) < 1.4) ) region = Barrel;
   else if( (fabs(eta) >= 1.4) && (fabs(eta) < 2.6) ) region = Endcap;
   else if( (fabs(eta) >= 2.6) && (fabs(eta) < 3.2) ) region = Transition;
   else if( (fabs(eta) >= 3.2) ) region = Forward;
   std::map<int,std::pair<double,double> > const& thresholds = thresholdMap.find(region)->second;
   
   double ptThreshold = -1.0;
   double eThreshold = -1.0;
   int partType = part.particleId();
   std::map<int,std::pair<double,double> >::const_iterator it_threshold = thresholds.find(partType);
   if(it_threshold != thresholds.end()) {
      ptThreshold = it_threshold->second.first;
      eThreshold = it_threshold->second.second;
   }

   if(part.pt() < ptThreshold) accept = false;
   if(part.energy() < eThreshold) accept = false;

   return accept;
}

// ==========================================================================================================
double D2Kpipi_Analyzer::MassColl(reco::PFCandidateCollection const& pflowCollection, 
                                 std::map<int, std::map<int,std::pair<double,double> > > const& thresholdMap)
// ==========================================================================================================
{
   math::XYZTLorentzVector allCands(0.,0.,0.,0.);
   reco::PFCandidateCollection::const_iterator part = pflowCollection.begin();
   reco::PFCandidateCollection::const_iterator pfCands_end = pflowCollection.end();
   for(; part != pfCands_end; ++part){
      if(pflowThreshold(*part,thresholdMap)) allCands += part->p4();
   }

   return allCands.M();
}


// ======================================================================================================================
std::pair<double,double> D2Kpipi_Analyzer::xi(reco::PFCandidateCollection const& pflowCollection, double Ebeam, 
                                             std::map<int, std::map<int,std::pair<double,double> > > const& thresholdMap)
// ======================================================================================================================
{
   double xi_towers_plus = 0.;
   double xi_towers_minus = 0.;
   reco::PFCandidateCollection::const_iterator part = pflowCollection.begin();
   reco::PFCandidateCollection::const_iterator pfCands_end = pflowCollection.end();
   for(; part != pfCands_end; ++part){
     if(!pflowThreshold(*part,thresholdMap)) continue;

     xi_towers_plus += part->et()*TMath::Exp(part->eta());
     xi_towers_minus += part->et()*TMath::Exp(-part->eta());
   }

   xi_towers_plus /= 2*Ebeam;
   xi_towers_minus /= 2*Ebeam;
   
   return std::make_pair(xi_towers_plus,xi_towers_minus);
}


// ============================================================================================================================
std::pair<double,double> D2Kpipi_Analyzer::EPlusPz(reco::PFCandidateCollection const& pflowCollection, 
                                                  std::map<int, std::map<int,std::pair<double,double> > > const& thresholdMap)
// ============================================================================================================================
{
   double e_plus_pz = 0.;
   double e_minus_pz = 0.;
   reco::PFCandidateCollection::const_iterator part = pflowCollection.begin();
   reco::PFCandidateCollection::const_iterator pfCands_end = pflowCollection.end();
   for(; part != pfCands_end; ++part){
      if(!pflowThreshold(*part,thresholdMap)) continue;

      e_plus_pz += part->energy() + part->pz(); 
      e_minus_pz += part->energy() - part->pz();
   }

   return std::make_pair(e_plus_pz,e_minus_pz);
}


// ===========================================================================================================================
std::pair<double,double> D2Kpipi_Analyzer::sumEHF(reco::PFCandidateCollection const& pflowCollection, 
                                                 std::map<int, std::map<int,std::pair<double,double> > > const& thresholdMap)
// ===========================================================================================================================
{
   double sumEHFPlus = 0.;
   double sumEHFMinus = 0.;
   reco::PFCandidateCollection::const_iterator part = pflowCollection.begin();
   reco::PFCandidateCollection::const_iterator pfCands_end = pflowCollection.end();
   for(; part != pfCands_end; ++part){
     if(!pflowThreshold(*part,thresholdMap)) continue;
     if((3.0 < part->eta()) && (part->eta() < 4.9) ){
        sumEHFPlus += part->energy();
     }
     if((-4.9 < part->eta()) && (part->eta() < -3.0) ){
        sumEHFMinus += part->energy();
     }

   }

   return std::make_pair(sumEHFPlus,sumEHFMinus);
}


// ===========================================================================================================================
std::pair<double,double> D2Kpipi_Analyzer::etaMax(reco::PFCandidateCollection const& pflowCollection, 
                                                 std::map<int, std::map<int,std::pair<double,double> > > const& thresholdMap)
// ===========================================================================================================================
{
   std::vector<double> etaCands;
   reco::PFCandidateCollection::const_iterator part = pflowCollection.begin();
   reco::PFCandidateCollection::const_iterator pfCands_end = pflowCollection.end();
   for(; part != pfCands_end; ++part){                           
      if(!pflowThreshold(*part,thresholdMap)) continue;            
      etaCands.push_back( part->eta() );
   }                                                             
   double eta_max = etaCands.size() ? *( std::max_element(etaCands.begin(), etaCands.end()) ) : -999.;
   double eta_min = etaCands.size() ? *( std::min_element(etaCands.begin(), etaCands.end()) ) : -999.;

   return std::make_pair(eta_max,eta_min);
}


/*

//  This method has been temporarily disabled as the variables which are given values on it were commented 
//  on the main code, namely: _pfsis1Eta_max, _pfsis1Eta_min, _pfsis2Eta_max, _pfsis2Eta_min, _deltaEtapf .

// ==========================================================================================================
void D2Kpipi_Analyzer::etaMaxMin(const edm::Event& iEvent, const edm::EventSetup& iSetup, 
                                std::map<int, std::map<int,std::pair<double,double> > > const& thresholdMap)
// ==========================================================================================================
{
   // Particle Flow
   edm::Handle<reco::PFCandidateCollection> pfsis;
   iEvent.getByLabel("particleFlow",pfsis);

   // Declaring Variables
   int pfsissize = pfsis->size();
   int Npfsis=0;
   int ipfsis;
   const reco::PFCandidate* pfsis1_max=NULL;
   const reco::PFCandidate* pfsis2_max=NULL;
   const reco::PFCandidate* pfsis1_min=NULL;
   const reco::PFCandidate* pfsis2_min=NULL;

   // Defense
   if (pfsissize >= 2)
   {
      for(ipfsis=0; ipfsis < pfsissize; ++ipfsis)
      {
         ++Npfsis;		
         const reco::PFCandidate* pfsis_max = &((*pfsis)[ipfsis]);
         if (!pflowThreshold(*pfsis_max,thresholdMap)) continue;
         if (pfsis_max==NULL) continue;
         if (pfsis1_max==NULL) {pfsis1_max=pfsis_max; continue;}
         if (pfsis_max->eta()>pfsis1_max->eta()) {
            pfsis2_max=pfsis1_max;
            pfsis1_max=pfsis_max;
            continue;
         }
         if (pfsis2_max==NULL) {pfsis2_max=pfsis_max; continue;}
         if (pfsis_max->eta()>pfsis2_max->eta()) pfsis2_max = pfsis_max;
      }
      if(pfsis1_max!=NULL) _pfsis1Eta_max = pfsis1_max->eta();
      if(pfsis2_max!=NULL) _pfsis2Eta_max = pfsis2_max->eta();
//   }
//
//   if (pfsissize >= 2)
//   {
      for(ipfsis=0; ipfsis < pfsissize; ++ipfsis)
      {
         ++Npfsis;		
         const reco::PFCandidate* pfsis_min = &((*pfsis)[ipfsis]);
         if (!pflowThreshold(*pfsis_min,thresholdMap)) continue;     
         if (pfsis_min==NULL) continue;
         if (pfsis1_min==NULL) {pfsis1_min=pfsis_min; continue;}
         if (pfsis_min->eta()<pfsis1_min->eta()) {
            pfsis2_min=pfsis1_min;
            pfsis1_min=pfsis_min;
            continue;
         }
         if (pfsis2_min==NULL) {pfsis2_min=pfsis_min; continue;}
         if (pfsis_min->eta()<pfsis2_min->eta()) pfsis2_min = pfsis_min;
      }
      if(pfsis1_min!=NULL) _pfsis1Eta_min = pfsis1_min->eta();
      if(pfsis2_min!=NULL) _pfsis2Eta_min = pfsis2_min->eta();
 
      // Esta certo isso ??
      if(pfsis1_min!=NULL && pfsis1_max!=NULL) _deltaEtapf = fabs(pfsis1_max->eta() - pfsis1_min->eta());
    }
}
*/

// =========================================================================================
void D2Kpipi_Analyzer::sumEHFTowers(const edm::Event& iEvent, const edm::EventSetup& iSetup)
// =========================================================================================
{
   edm::Handle<CaloTowerCollection> caloTowers;
   iEvent.getByLabel("towerMaker",caloTowers);

   for ( CaloTowerCollection::const_iterator cal = caloTowers->begin(); cal != caloTowers->end(); ++ cal )
   {
      if (fabs(cal->eta()) > 3.0 && fabs(cal->eta()) < 4.9)
      {
         if (cal->energy() >= 4.0)
         {
            if (cal->zside() < 0)
            {
               _HFMinusEnergyFromCaloTowers += cal->energy();
               _nHFMinus++;
            }
            if (cal->zside() > 0)
            {
               _HFPlusEnergyFromCaloTowers += cal->energy();
               _nHFPlus++;
            }
         }
      }
   }
   // if (_HFMinusEnergyFromCaloTowers == 0. || _HFPlusEnergyFromCaloTowers == 0.) _zeroEnergyHFSideFromCaloTowers = true ;
}



//define this as a plug-in
DEFINE_FWK_MODULE(D2Kpipi_Analyzer);
