#include <iostream>
#include <fstream>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TChain.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "load_d2kpipi_MB_MC.C"

void load_d2kpipi_MB_MC(TChain*) ;

Double_t GaussPol1 (Double_t * , Double_t *) ;
Double_t GaussPol2 (Double_t * , Double_t *) ;

double binwidth;

int D2KpipiRecoMC()
///////////////////////////////////////////////////////////////////////////////////////// 
//                                                                                     //
//  Purpose :  Macro to compare distributions from 2 root files.                       //
//                                                                                     //
//  Parameters :  Set accordingly the variables                                        //
//                                                                                     //
//    fileName1                                                                        //
//    fileName2                                                                        //
//                                                                                     //
//    treeName                                                                         //
//                                                                                     //
//  Author :  Wagner Carvalho , 11/02/2009                                             //
//                                                                                     //
///////////////////////////////////////////////////////////////////////////////////////// 
{
   // TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           hlt_bscOR_bptxOR;
   Int_t           procId;
   Bool_t          hasdplus;
   Bool_t          hasdplus2kpipi;
   Bool_t          hasd0;
   Bool_t          hasds;
   Bool_t          hasdstarplus;
   Double_t        refittedMassDplus;
   Double_t        refittedMassKpi1;
   Double_t        refittedMassKpi2;
   Int_t           nvtx;
   Double_t        vtxchi2dof;
   Double_t        Dxy;
   Double_t        Sxy;
   Double_t        D3D;
   Double_t        S3D;
   Double_t        Rxy;
   Double_t        refittedAngle;
   Double_t        dplusp;
   Double_t        dpluspT;
   Double_t        dplusphi;
   Double_t        dpluseta;
   Double_t        trk1dcaxy;
   Double_t        trk2dcaxy;
   Double_t        trk3dcaxy;
   Double_t        trk1dcaz;
   Double_t        trk2dcaz;
   Double_t        trk3dcaz;
   Double_t        trk1p;
   Double_t        trk2p;
   Double_t        trk3p;
   Double_t        trk1pT;
   Double_t        trk2pT;
   Double_t        trk3pT;
   Double_t        trk1chi2;
   Double_t        trk2chi2;
   Double_t        trk3chi2;
   Int_t           trk1hits;
   Int_t           trk2hits;
   Int_t           trk3hits;
   Int_t           nHFPlus;
   Int_t           nHFMinus;
   Double_t        HFPlusEnergyFromCTowers;
   Double_t        HFMinusEnergyFromCTowers;
   Double_t        HFPlusEnergyFromPF;
   Double_t        HFMinusEnergyFromPF;
   Double_t        EPlusPzFromPF;
   Double_t        EMinusPzFromPF;
   Double_t        xiPlusFromPF;
   Double_t        xiMinusFromPF;
   Double_t        etaMaxFromPF;
   Double_t        etaMinFromPF;
   Double_t        MxFromPF;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_hlt_bscOR_bptxOR;   //!
   TBranch        *b_procId;   //!
   TBranch        *b_hasdplus;   //!
   TBranch        *b_hasdplus2kpipi;   //!
   TBranch        *b_hasd0;   //!
   TBranch        *b_hasds;   //!
   TBranch        *b_hasdstarplus;   //!
   TBranch        *b_refittedMassDplus;   //!
   TBranch        *b_refittedMassKpi1;   //!
   TBranch        *b_refittedMassKpi2;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_vtxchi2dof;   //!
   TBranch        *b_Dxy;   //!
   TBranch        *b_Sxy;   //!
   TBranch        *b_D3D;   //!
   TBranch        *b_S3D;   //!
   TBranch        *b_Rxy;   //!
   TBranch        *b_refittedAngle;   //!
   TBranch        *b_dplusp;   //!
   TBranch        *b_dpluspT;   //!
   TBranch        *b_dplusphi;   //!
   TBranch        *b_dpluseta;   //!
   TBranch        *b_trk1dcaxy;   //!
   TBranch        *b_trk2dcaxy;   //!
   TBranch        *b_trk3dcaxy;   //!
   TBranch        *b_trk1dcaz;   //!
   TBranch        *b_trk2dcaz;   //!
   TBranch        *b_trk3dcaz;   //!
   TBranch        *b_trk1p;   //!
   TBranch        *b_trk2p;   //!
   TBranch        *b_trk3p;   //!
   TBranch        *b_trk1pT;   //!
   TBranch        *b_trk2pT;   //!
   TBranch        *b_trk3pT;   //!
   TBranch        *b_trk1chi2;   //!
   TBranch        *b_trk2chi2;   //!
   TBranch        *b_trk3chi2;   //!
   TBranch        *b_trk1hits;   //!
   TBranch        *b_trk1hits;   //!
   TBranch        *b_trk1hits;   //!
   TBranch        *b_nHFPlus;   //!
   TBranch        *b_nHFMinus;   //!
   TBranch        *b_HFPlusEnergyFromCTowers;   //!
   TBranch        *b_HFMinusEnergyFromCTowers;   //!
   TBranch        *b_HFPlusEnergyFromPF;   //!
   TBranch        *b_HFMinusEnergyFromPF;   //!
   TBranch        *b_EPlusPzFromPF;   //!
   TBranch        *b_EMinusPzFromPF;   //!
   TBranch        *b_xiPlusFromPF;   //!
   TBranch        *b_xiMinusFromPF;   //!
   TBranch        *b_etaMaxFromPF;   //!
   TBranch        *b_etaMinFromPF;   //!
   TBranch        *b_MxFromPF;   //!

   // Histograms
   
   //  Parameters for mass histograms
   const int NBINS = 30 ;
   const double LOWER_MASS = 1.70 ;
   const double HIGHER_MASS = 2.04 ;
   binwidth = (HIGHER_MASS - LOWER_MASS) / NBINS ;
   
   
   TH1F *h_refitMassDplus_0 = new TH1F("h_refitMassDplus_0","Mass of D+ candidate (refitted tracks)",NBINS,LOWER_MASS,HIGHER_MASS) ; 
   TH1F *h_refitMassDplus = new TH1F("h_refitMassDplus","Mass of D+ candidate (refitted tracks)",NBINS,LOWER_MASS,HIGHER_MASS) ; 
   TH1F *h_refitMassKpi1 = new TH1F("h_refitMassKpi1","Invariant mass of K+pi1 system (refitted tracks)",120,0.60,1.80) ; 
   TH1F *h_refitMassKpi2 = new TH1F("h_refitMassKpi2","Invariant mass of K+pi2 system (refitted tracks)",120,0.60,1.80) ; 
   //  HISTOS FOR SCANNING !!!
   TH1F *h_angle = new TH1F("h_angle","Angle of D+ candidate momentum (refitted tracks)",79,0.,3.16) ; 
   // TH1F *h_vtxchi2 = new TH1F("h_vtxchi2","Chi2/d.o.f. of D+ candidate vertex (refitted tracks)",6,0.,3.) ; 
   // TH1F *h_sxy = new TH1F("h_sxy","Sxy of D+ candidate vertex",9,1.,10.) ; 
   // TH1F *h_dplusp = new TH1F("h_dplusp","Momentum of D+ candidate",50,0.,25.) ; 
   // TH1F *h_dpluspT = new TH1F("h_dpluspT","pT of D+ candidate",20,0.,10.) ; 


   TH1F *h_fake_refitMassDplus_0 = new TH1F("h_fake_refitMassDplus_0","Mass of D+ candidate (refitted tracks)",NBINS,LOWER_MASS,HIGHER_MASS) ; 
   TH1F *h_fake_refitMassDplus = new TH1F("h_fake_refitMassDplus","Mass of D+ fake candidate (refitted tracks)",NBINS,LOWER_MASS,HIGHER_MASS) ; 
   TH1F *h_fake_refitMassKpi1 = new TH1F("h_fake_refitMassKpi1","Invariant mass of K+pi1 system (refitted tracks)",120,0.60,1.80) ; 
   TH1F *h_fake_refitMassKpi2 = new TH1F("h_fake_refitMassKpi2","Invariant mass of K+pi2 system (refitted tracks)",120,0.60,1.80) ; 
   //  HISTOS FOR SCANNING !!!
   TH1F *h_fake_angle = new TH1F("h_fake_angle","Angle of D+ fake candidate momentum (refitted tracks)",79,0.,3.16) ; 
   // TH1F *h_fake_vtxchi2 = new TH1F("h_fake_vtxchi2","Chi2/d.o.f. of D+ candidate vertex (refitted tracks)",6,0.,3.) ; 
   // TH1F *h_fake_sxy = new TH1F("h_fake_sxy","Sxy of D+ candidate vertex",9,1.,10.) ; 
   // TH1F *h_fake_dplusp = new TH1F("h_fake_dplusp","Momentum of D+ candidate",50,0.,25.) ; 
   // TH1F *h_fake_dpluspT = new TH1F("h_fake_dpluspT","pT of D+ candidate",20,0.,10.) ; 


   TH1F *h_diff_refitMassDplus = new TH1F("h_diff_refitMassDplus","D+ candidate mass with refitted tracks (fake subtracted)",NBINS,LOWER_MASS,HIGHER_MASS) ; 
   TH1F *h_mass = new TH1F("h_mass","Mass of D+ candidate (refitted tracks)",30,1.7,2.0) ; 


   //  Pass the cuts through these variables
   double trkdcaxy_cut = 0.1 ;    //  no cut = 150
   double trkdcaz_cut = 0.5 ;     //  no cut = 300
   double trkhits_cut = 3 ;       //  no cut = 3
   double trk1p_cut = 0.5 ;       //  no cut = 0.5
   double trk2p_cut = 0.5 ;       //  no cut = 0.5
   double trk3p_cut = 0.5 ;       //  no cut = 0.5
   double trk1pT_cut = 0.5 ;      //  no cut = 0.5 
   double trk2pT_cut = 0.5 ;      //  no cut = 0.5
   double trk3pT_cut = 0.5 ;      //  no cut = 0.5
   double pT_1st_cut = 0.5 ;      //  std cut = 1.7 
   double pT_2nd_cut = 0.5 ;      //  std cut = 1.5 
   double pT_3rd_cut = 0.5 ;      //  std cut = 0.9 
   double dplusp_cut = 1000 ;     //  no cut = 1000.
   double dpluspT_cut = 0.0 ;     //  no cut = 0.
   double angle_cut = 0.05 ;      //  no cut = 3.1416
   double angle_cut_min = 0.000 ;
   double angle_cut_max = 0.050 ;
   double Sxy_cut = 3. ;          //  no cut = 1
   double vtxchi2dof_cut = 15. ;  //  no cut = 10 ?
   // double trkperp_cut = 0.7 ;      //  no cut = 0.
   // vtxchi2dof < 1.5  ?!?!
   // dpluspT > 3  ?!?!

/*
   //  Ewerton cuts
   trkdcaxy_cut = 0.1 ;    //  no cut = 150
   trkdcaz_cut = 0.5 ;     //  no cut = 300
   trkhits_cut = 5 ;       //  std cut = 5
   trk1p_cut = 0.5 ;       //  no cut = 0.5
   trk2p_cut = 0.5 ;       //  no cut = 0.5
   trk3p_cut = 0.5 ;       //  no cut = 0.5
   trk1pT_cut = 0.5 ;      //  no cut = 0.5 
   trk2pT_cut = 0.5 ;      //  no cut = 0.5
   trk3pT_cut = 0.5 ;      //  no cut = 0.5
   pT_1st_cut = 1.7 ;      //  std cut = 1.7 
   pT_2nd_cut = 1.5 ;      //  std cut = 1.5 
   pT_3rd_cut = 0.9 ;      //  std cut = 0.9 
   dplusp_cut = 1500 ;     //  no cut = 1000.
   dpluspT_cut = 4.3 ;     //  std cut = 0.
   angle_cut = 0.04 ;      //  std cut = 3.1416
   angle_cut_min = 0.000 ;
   angle_cut_max = 0.040 ;
   Sxy_cut = 3.6 ;         //  std cut = 1
   vtxchi2dof_cut = 200 ;  //  no cut = 10 ?

*/   
   // FOR VARIABLE SCANNING:
   bool DoSCAN = false ;
   
   double CUT_OFFSET, CUT_MIN, CUT_MAX, cut_offset, cut_min, cut_max, d_cut, cut_lim ;
   // This is for TCanvas cscan and cscan_fake
   int rows, columns, k ;
   
   if(DoSCAN) {
      CUT_OFFSET = 0 * 0.80 ;
      CUT_MIN = 0.00 + CUT_OFFSET ;
      CUT_MAX = 0.025 + CUT_OFFSET ;
      rows = 2 ;
      columns = 5 ;
   }
   
   // For signal window cut study
   double mass_window  = 0.010 ;
   double nominal_mass = 1.870 ;
   
   TH1F *h_sig_vtxchi2 = new TH1F("h_sig_vtxchi2","Chi2/d.o.f. of D+ vertex",6,0.,3.) ;
   TH1F *h_nosig_vtxchi2 = new TH1F("h_nosig_vtxchi2","Chi2/d.o.f. of D+ vertex",6,0.,3.) ;
   TH1F *h_sig_dplusp = new TH1F("h_sig_dplusp","Momentum of D+",30,0.,30.) ;
   TH1F *h_nosig_dplusp = new TH1F("h_nosig_dplusp","Momentum of D+",30,0.,30.) ;
   TH1F *h_sig_dpluspT = new TH1F("h_sig_dpluspT","pT of D+",20,0.,10.) ;
   TH1F *h_nosig_dpluspT = new TH1F("h_nosig_dpluspT","pT of D+",20,0.,10.) ;
   TH1F *h_sig_trk1chi2 = new TH1F("h_sig_trk1chi2","track chi2dof",50,0.,5.) ;
   TH1F *h_nosig_trk1chi2 = new TH1F("h_nosig_trk1chi2","track chi2dof",50,0.,5.) ;
   TH1F *h_sig_trk1hits = new TH1F("h_sig_trk1hits","track hits",40,0.,40.) ;
   TH1F *h_nosig_trk1hits = new TH1F("h_nosig_trk1hits","track hits",40,0.,40.) ;
   TH1F *h_sig_trk1p = new TH1F("h_sig_trk1p","track p",40,0.,20.) ;
   TH1F *h_nosig_trk1p = new TH1F("h_nosig_trk1p","track p",40,0.,20.) ;
   TH1F *h_sig_trk1pT = new TH1F("h_sig_trk1pT","track pT",40,0.,4.0) ;
   TH1F *h_nosig_trk1pT = new TH1F("h_nosig_trk1pT","track pT",40,0.,4.0) ;
   
   // bool keep = false ;
   
   //   Pre-selection cuts:
   //   -  trk1hits > 4      ;  trk2hits > 4
   //   -  trk1pT > 0.5 GeV  ;  trk2pT > 0.5 GeV
   //   -  trk1chi2 < 5      ;  trk2chi2 < 5
   //   -  vtxchi2dof < 3    ;  Sxy > 1       ;  Rxy < 2 cm
   // Char_t *fileName1 = "MB_test_10000evt_new_1_1_a4D.root";
   
   // Set here the TTree name
   Char_t *goodTree = "d2kpipi/theTree";
   Char_t *fakeTree = "d2kpipi/fakeTree";
   
   // Open file
   // cout << "\n    ***  Opening file \"" << fileName1 << "\"  ***" << endl;
   // TFile *file1 = new TFile(fileName1);


   // fChain = (TTree*)file1->Get(goodTree);
   fChain = new TChain(goodTree) ;
   // fChain->Add("D2Kpipi_CCbar_Pt_5_10_7TeV_MC_42_V13_All_1_1_j5g.root") ;
   // fChain->Add("D2Kpipi_MB_pythia8_4C_7TeV_HF-SL_Summer12-LowPU2010_DR42_NoPileUp_START42_V17C-v1_1_1_Srx.root") ;
   load_d2kpipi_MB_MC(fChain) ;
   // Init(fChain);
   if (!fChain) return 1;
   
   // Set branch addresses and branch pointers
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("hlt_bscOR_bptxOR", &hlt_bscOR_bptxOR, &b_hlt_bscOR_bptxOR);
   fChain->SetBranchAddress("procId", &procId, &b_procId);
   fChain->SetBranchAddress("hasdplus", &hasdplus, &b_hasdplus);
   fChain->SetBranchAddress("hasdplus2kpipi", &hasdplus2kpipi, &b_hasdplus2kpipi);
   fChain->SetBranchAddress("hasd0", &hasd0, &b_hasd0);
   fChain->SetBranchAddress("hasds", &hasds, &b_hasds);
   fChain->SetBranchAddress("hasdstarplus", &hasdstarplus, &b_hasdstarplus);
   fChain->SetBranchAddress("refittedMassDplus", &refittedMassDplus, &b_refittedMassDplus);
   fChain->SetBranchAddress("refittedMassKpi1", &refittedMassKpi1, &b_refittedMassKpi1);
   fChain->SetBranchAddress("refittedMassKpi2", &refittedMassKpi2, &b_refittedMassKpi2);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("vtxchi2dof", &vtxchi2dof, &b_vtxchi2dof);
   fChain->SetBranchAddress("Dxy", &Dxy, &b_Dxy);
   fChain->SetBranchAddress("Sxy", &Sxy, &b_Sxy);
   fChain->SetBranchAddress("D3D", &D3D, &b_D3D);
   fChain->SetBranchAddress("S3D", &S3D, &b_S3D);
   fChain->SetBranchAddress("Rxy", &Rxy, &b_Rxy);
   fChain->SetBranchAddress("refittedAngle", &refittedAngle, &b_refittedAngle);
   fChain->SetBranchAddress("dplusp", &dplusp, &b_dplusp);
   fChain->SetBranchAddress("dpluspT", &dpluspT, &b_dpluspT);
   fChain->SetBranchAddress("dplusphi", &dplusphi, &b_dplusphi);
   fChain->SetBranchAddress("dpluseta", &dpluseta, &b_dpluseta);
   fChain->SetBranchAddress("trk1dcaxy", &trk1dcaxy, &b_trk1dcaxy);
   fChain->SetBranchAddress("trk2dcaxy", &trk2dcaxy, &b_trk2dcaxy);
   fChain->SetBranchAddress("trk3dcaxy", &trk3dcaxy, &b_trk3dcaxy);
   fChain->SetBranchAddress("trk1dcaz", &trk1dcaz, &b_trk1dcaz);
   fChain->SetBranchAddress("trk2dcaz", &trk2dcaz, &b_trk2dcaz);
   fChain->SetBranchAddress("trk3dcaz", &trk3dcaz, &b_trk3dcaz);
   fChain->SetBranchAddress("trk1p", &trk1p, &b_trk1p);
   fChain->SetBranchAddress("trk2p", &trk2p, &b_trk2p);
   fChain->SetBranchAddress("trk3p", &trk3p, &b_trk3p);
   fChain->SetBranchAddress("trk1pT", &trk1pT, &b_trk1pT);
   fChain->SetBranchAddress("trk2pT", &trk2pT, &b_trk2pT);
   fChain->SetBranchAddress("trk3pT", &trk3pT, &b_trk3pT);
   fChain->SetBranchAddress("trk1chi2", &trk1chi2, &b_trk1chi2);
   fChain->SetBranchAddress("trk2chi2", &trk2chi2, &b_trk2chi2);
   fChain->SetBranchAddress("trk3chi2", &trk3chi2, &b_trk3chi2);
   fChain->SetBranchAddress("trk1hits", &trk1hits, &b_trk1hits);
   fChain->SetBranchAddress("trk2hits", &trk2hits, &b_trk1hits);
   fChain->SetBranchAddress("trk3hits", &trk3hits, &b_trk1hits);
   fChain->SetBranchAddress("nHFPlus", &nHFPlus, &b_nHFPlus);
   fChain->SetBranchAddress("nHFMinus", &nHFMinus, &b_nHFMinus);
   fChain->SetBranchAddress("HFPlusEnergyFromCTowers", &HFPlusEnergyFromCTowers, &b_HFPlusEnergyFromCTowers);
   fChain->SetBranchAddress("HFMinusEnergyFromCTowers", &HFMinusEnergyFromCTowers, &b_HFMinusEnergyFromCTowers);
   fChain->SetBranchAddress("HFPlusEnergyFromPF", &HFPlusEnergyFromPF, &b_HFPlusEnergyFromPF);
   fChain->SetBranchAddress("HFMinusEnergyFromPF", &HFMinusEnergyFromPF, &b_HFMinusEnergyFromPF);
   fChain->SetBranchAddress("EPlusPzFromPF", &EPlusPzFromPF, &b_EPlusPzFromPF);
   fChain->SetBranchAddress("EMinusPzFromPF", &EMinusPzFromPF, &b_EMinusPzFromPF);
   fChain->SetBranchAddress("xiPlusFromPF", &xiPlusFromPF, &b_xiPlusFromPF);
   fChain->SetBranchAddress("xiMinusFromPF", &xiMinusFromPF, &b_xiMinusFromPF);
   fChain->SetBranchAddress("etaMaxFromPF", &etaMaxFromPF, &b_etaMaxFromPF);
   fChain->SetBranchAddress("etaMinFromPF", &etaMinFromPF, &b_etaMinFromPF);
   fChain->SetBranchAddress("MxFromPF", &MxFromPF, &b_MxFromPF);

   int nentries = fChain->GetEntries();
   
   cout << " Nentries = " << nentries << endl;
   
   for(int i=0; i<nentries; i++) {
      // Get entry
      fChain->GetEntry(i);
      //  Reference plots with no extra cut
      h_refitMassDplus_0->Fill(refittedMassDplus) ;
      //  Apply cuts
      bool keep = false ;
      double pT_1st , pT_2nd , pT_3rd ;
      pT_1st = pT_2nd = pT_3rd = 0. ;
      if(trk1pT>trk2pT) {
         if(trk1pT>trk3pT) {
            pT_1st = trk1pT ;
            if(trk2pT>trk3pT) {
               pT_2nd = trk2pT ;
               pT_3rd = trk3pT ;
            } else {
               pT_2nd = trk3pT ;
               pT_3rd = trk2pT ;
            }
         } else { 
            pT_1st = trk3pT ;
            pT_2nd = trk1pT ;
            pT_3rd = trk2pT ;
         }
      } else {
         if(trk2pT>trk3pT) {
            pT_1st = trk2pT ;
            if(trk1pT>trk3pT) {
               pT_2nd = trk1pT ;
               pT_3rd = trk3pT ;
            } else {
               pT_2nd = trk3pT ;
               pT_3rd = trk1pT ;
            }
         } else { 
            pT_1st = trk3pT ;
            pT_2nd = trk2pT ;
            pT_3rd = trk1pT ;
         }
      }
//      keep = trk1pT>trk1pT_cut && trk2pT>trk2pT_cut && trk3pT>trk3pT_cut &&
//              trk1p>trk1p_cut  &&  trk2p>trk2p_cut  &&  trk3p>trk3p_cut  &&
      keep = refittedMassDplus>LOWER_MASS && refittedMassDplus<HIGHER_MASS && 
             pT_1st>pT_1st_cut && pT_2nd>pT_2nd_cut && pT_3rd>pT_3rd_cut &&
             trk1hits>trkhits_cut && trk2hits>trkhits_cut && trk3hits>trkhits_cut && 
             trk1dcaxy<trkdcaxy_cut && trk2dcaxy<trkdcaxy_cut && 
             trk1dcaz<trkdcaz_cut && trk2dcaz<trkdcaz_cut &&
             // hlt_bscOR_bptxOR>=hlt_cut && 
             // hasdplus==1 &&
             hasdplus2kpipi==1 &&
             // hasdplus2kpipi!=1 &&
             // hasd0==1 &&
             // hasds==1 &&
             // hasdstarplus==1 &&
             // ( hasd0==1 || hasds==1 || hasdstarplus==1 ) &&
             vtxchi2dof < vtxchi2dof_cut &&
             refittedAngle>=angle_cut_min && refittedAngle<angle_cut_max && Sxy>Sxy_cut && 
             dpluspT>dpluspT_cut && dplusp<dplusp_cut ;
             // refittedAngle>=angle_cut_min && refittedAngle<angle_cut_max && Sxy>Sxy_cut && vtxchi2dof < 1.5 ;
             // refittedAngle<angle_cut && Sxy>Sxy_cut && vtxchi2dof < 1.5 ;
      if(keep) {
         h_refitMassDplus->Fill(refittedMassDplus) ;
         h_refitMassKpi1->Fill(refittedMassKpi1) ;
         h_refitMassKpi2->Fill(refittedMassKpi2) ;
         // Mass window cut study
         if(fabs(refittedMassDplus-nominal_mass) < mass_window) {
            h_sig_vtxchi2->Fill(vtxchi2dof) ;
            h_sig_dplusp->Fill(dplusp) ;
            h_sig_dpluspT->Fill(dpluspT) ;
            h_sig_trk1chi2->Fill(trk1chi2) ;
            h_sig_trk1hits->Fill(trk1hits) ;
            h_sig_trk1p->Fill(trk1p) ;
            h_sig_trk1pT->Fill(trk1pT) ;
         } else if(fabs(refittedMassDplus-nominal_mass) > 3*mass_window) {
            h_nosig_vtxchi2->Fill(vtxchi2dof) ;
            h_nosig_dplusp->Fill(dplusp) ;
            h_nosig_dpluspT->Fill(dpluspT) ;
            h_nosig_trk1chi2->Fill(trk1chi2) ;
            h_nosig_trk1hits->Fill(trk1hits) ;
            h_nosig_trk1p->Fill(trk1p) ;
            h_nosig_trk1pT->Fill(trk1pT) ;
         } 
      }
      if( trk1pT>trk1pT_cut && trk2pT>trk2pT_cut && trk3pT>trk3pT_cut &&
           trk1p>trk1p_cut  &&  trk2p>trk2p_cut  &&  trk3p>trk3p_cut  &&
           Sxy>Sxy_cut) h_angle->Fill(refittedAngle ) ;
   }

   if(DoSCAN) {
      cut_offset = CUT_OFFSET ;
      cut_min = CUT_MIN ;
      cut_max = CUT_MAX ;
      d_cut = cut_max - cut_min ;
      cut_lim = (rows * columns - 0.5) * d_cut + cut_offset ;
      cout << "\n     cut_lim = " << cut_lim << endl ;
      k = 1 ;
      //  Create Canvas
      TCanvas *cscan = new TCanvas("cscan","Vtx chi2/dof  scan",10,10,1500,800);
      cscan->Divide(columns,rows);
      while(cut_min < cut_lim) {
        h_mass->Reset() ;
        for(int i=0; i<nentries; i++) {
          // Get entry
          fChain->GetEntry(i);
          //  Apply cuts
          bool keep = false ;
          keep = trk1pT>trk1pT_cut && trk2pT>trk2pT_cut && trk3pT>trk3pT_cut &&
                  trk1p>trk1p_cut  &&  trk2p>trk2p_cut  &&  trk3p>trk3p_cut  &&
                 refittedAngle>=cut_min && refittedAngle<cut_max && Sxy>Sxy_cut ;
          if(keep) {
             h_mass->Fill(refittedMassDplus) ;
          }
        }
        // Draw histo
        cscan->cd(k) ;
        cout << "\n  for k = " << k << " , angle = [ " << cut_min << " : " << cut_max << " )" ; 
        h_mass->DrawCopy("e1") ;
        // Update angle and counter
        cut_min += d_cut ;
        cut_max += d_cut ;
        k++ ;
      }
   }


   // fChain = (TTree*)file1->Get(fakeTree);
   fChain = new TChain(fakeTree) ;
   // fChain->Add("D2Kpipi_CCbar_Pt_5_10_7TeV_MC_42_V13_All_1_1_j5g.root") ;
   // fChain->Add("D2Kpipi_MB_pythia8_4C_7TeV_HF-SL_Summer12-LowPU2010_DR42_NoPileUp_START42_V17C-v1_1_1_Srx.root") ;
   load_d2kpipi_MB_MC(fChain) ;
   // Init(fChain);
   if (!fChain) return 1;
   
   // Set branch addresses and branch pointers
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("hlt_bscOR_bptxOR", &hlt_bscOR_bptxOR, &b_hlt_bscOR_bptxOR);
   fChain->SetBranchAddress("procId", &procId, &b_procId);
   fChain->SetBranchAddress("hasdplus", &hasdplus, &b_hasdplus);
   fChain->SetBranchAddress("hasdplus2kpipi", &hasdplus2kpipi, &b_hasdplus2kpipi);
   fChain->SetBranchAddress("hasd0", &hasd0, &b_hasd0);
   fChain->SetBranchAddress("hasds", &hasds, &b_hasds);
   fChain->SetBranchAddress("hasdstarplus", &hasdstarplus, &b_hasdstarplus);
   fChain->SetBranchAddress("refittedMassDplus", &refittedMassDplus, &b_refittedMassDplus);
   fChain->SetBranchAddress("refittedMassKpi1", &refittedMassKpi1, &b_refittedMassKpi1);
   fChain->SetBranchAddress("refittedMassKpi2", &refittedMassKpi2, &b_refittedMassKpi2);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("vtxchi2dof", &vtxchi2dof, &b_vtxchi2dof);
   fChain->SetBranchAddress("Dxy", &Dxy, &b_Dxy);
   fChain->SetBranchAddress("Sxy", &Sxy, &b_Sxy);
   fChain->SetBranchAddress("D3D", &D3D, &b_D3D);
   fChain->SetBranchAddress("S3D", &S3D, &b_S3D);
   fChain->SetBranchAddress("Rxy", &Rxy, &b_Rxy);
   fChain->SetBranchAddress("refittedAngle", &refittedAngle, &b_refittedAngle);
   fChain->SetBranchAddress("dplusp", &dplusp, &b_dplusp);
   fChain->SetBranchAddress("dpluspT", &dpluspT, &b_dpluspT);
   fChain->SetBranchAddress("dplusphi", &dplusphi, &b_dplusphi);
   fChain->SetBranchAddress("dpluseta", &dpluseta, &b_dpluseta);
   fChain->SetBranchAddress("trk1dcaxy", &trk1dcaxy, &b_trk1dcaxy);
   fChain->SetBranchAddress("trk2dcaxy", &trk2dcaxy, &b_trk2dcaxy);
   fChain->SetBranchAddress("trk3dcaxy", &trk3dcaxy, &b_trk3dcaxy);
   fChain->SetBranchAddress("trk1dcaz", &trk1dcaz, &b_trk1dcaz);
   fChain->SetBranchAddress("trk2dcaz", &trk2dcaz, &b_trk2dcaz);
   fChain->SetBranchAddress("trk3dcaz", &trk3dcaz, &b_trk3dcaz);
   fChain->SetBranchAddress("trk1p", &trk1p, &b_trk1p);
   fChain->SetBranchAddress("trk2p", &trk2p, &b_trk2p);
   fChain->SetBranchAddress("trk3p", &trk3p, &b_trk3p);
   fChain->SetBranchAddress("trk1pT", &trk1pT, &b_trk1pT);
   fChain->SetBranchAddress("trk2pT", &trk2pT, &b_trk2pT);
   fChain->SetBranchAddress("trk3pT", &trk3pT, &b_trk3pT);
   fChain->SetBranchAddress("trk1chi2", &trk1chi2, &b_trk1chi2);
   fChain->SetBranchAddress("trk2chi2", &trk2chi2, &b_trk2chi2);
   fChain->SetBranchAddress("trk3chi2", &trk3chi2, &b_trk3chi2);
   fChain->SetBranchAddress("trk1hits", &trk1hits, &b_trk1hits);
   fChain->SetBranchAddress("trk2hits", &trk2hits, &b_trk1hits);
   fChain->SetBranchAddress("trk3hits", &trk3hits, &b_trk1hits);
   fChain->SetBranchAddress("nHFPlus", &nHFPlus, &b_nHFPlus);
   fChain->SetBranchAddress("nHFMinus", &nHFMinus, &b_nHFMinus);
   fChain->SetBranchAddress("HFPlusEnergyFromCTowers", &HFPlusEnergyFromCTowers, &b_HFPlusEnergyFromCTowers);
   fChain->SetBranchAddress("HFMinusEnergyFromCTowers", &HFMinusEnergyFromCTowers, &b_HFMinusEnergyFromCTowers);
   fChain->SetBranchAddress("HFPlusEnergyFromPF", &HFPlusEnergyFromPF, &b_HFPlusEnergyFromPF);
   fChain->SetBranchAddress("HFMinusEnergyFromPF", &HFMinusEnergyFromPF, &b_HFMinusEnergyFromPF);
   fChain->SetBranchAddress("EPlusPzFromPF", &EPlusPzFromPF, &b_EPlusPzFromPF);
   fChain->SetBranchAddress("EMinusPzFromPF", &EMinusPzFromPF, &b_EMinusPzFromPF);
   fChain->SetBranchAddress("xiPlusFromPF", &xiPlusFromPF, &b_xiPlusFromPF);
   fChain->SetBranchAddress("xiMinusFromPF", &xiMinusFromPF, &b_xiMinusFromPF);
   fChain->SetBranchAddress("etaMaxFromPF", &etaMaxFromPF, &b_etaMaxFromPF);
   fChain->SetBranchAddress("etaMinFromPF", &etaMinFromPF, &b_etaMinFromPF);
   fChain->SetBranchAddress("MxFromPF", &MxFromPF, &b_MxFromPF);

   nentries = fChain->GetEntries();
   
   cout << " Nentries = " << nentries << endl;
   
   for(int i=0; i<nentries; i++) {
      // Get entry
      fChain->GetEntry(i);
      //  Reference plots with no extra cut
      h_fake_refitMassDplus_0->Fill(refittedMassDplus) ;
      //  Apply cuts
      bool keep = false ;
      double pT_1st , pT_2nd , pT_3rd ;
      pT_1st = pT_2nd = pT_3rd = 0. ;
      if(trk1pT>trk2pT) {
         if(trk1pT>trk3pT) {
            pT_1st = trk1pT ;
            if(trk2pT>trk3pT) {
               pT_2nd = trk2pT ;
               pT_3rd = trk3pT ;
            } else {
               pT_2nd = trk3pT ;
               pT_3rd = trk2pT ;
            }
         } else { 
            pT_1st = trk3pT ;
            pT_2nd = trk1pT ;
            pT_3rd = trk2pT ;
         }
      } else {
         if(trk2pT>trk3pT) {
            pT_1st = trk2pT ;
            if(trk1pT>trk3pT) {
               pT_2nd = trk1pT ;
               pT_3rd = trk3pT ;
            } else {
               pT_2nd = trk3pT ;
               pT_3rd = trk1pT ;
            }
         } else { 
            pT_1st = trk3pT ;
            pT_2nd = trk2pT ;
            pT_3rd = trk1pT ;
         }
      }
//      keep = trk1pT>trk1pT_cut && trk2pT>trk2pT_cut && trk3pT>trk3pT_cut &&
//              trk1p>trk1p_cut  &&  trk2p>trk2p_cut  &&  trk3p>trk3p_cut  &&
      keep = refittedMassDplus>LOWER_MASS && refittedMassDplus<HIGHER_MASS && 
             pT_1st>pT_1st_cut && pT_2nd>pT_2nd_cut && pT_3rd>pT_3rd_cut &&
             trk1hits>trkhits_cut && trk2hits>trkhits_cut && trk3hits>trkhits_cut && 
             trk1dcaxy<trkdcaxy_cut && trk2dcaxy<trkdcaxy_cut && 
             trk1dcaz<trkdcaz_cut && trk2dcaz<trkdcaz_cut &&
             // hlt_bscOR_bptxOR>=hlt_cut && 
             // hasdplus==1 &&
             hasdplus2kpipi==1 &&
             // hasdplus2kpipi!=1 &&
             // hasd0==1 &&
             // hasds==1 &&
             // hasdstarplus==1 &&
             // ( hasd0==1 || hasds==1 || hasdstarplus==1 ) &&
             vtxchi2dof < vtxchi2dof_cut &&
             refittedAngle>=angle_cut_min && refittedAngle<angle_cut_max && Sxy>Sxy_cut && 
             dpluspT>dpluspT_cut && dplusp<dplusp_cut ;
             // refittedAngle>=angle_cut_min && refittedAngle<angle_cut_max && Sxy>Sxy_cut && vtxchi2dof < 1.5 ;
             // refittedAngle<angle_cut && Sxy>Sxy_cut && vtxchi2dof < 1.5 ;
      if(keep) {
         h_fake_refitMassDplus->Fill(refittedMassDplus) ;
         h_fake_refitMassKpi1->Fill(refittedMassKpi1) ;
         h_fake_refitMassKpi2->Fill(refittedMassKpi2) ;
      }
      if( trk1pT>trk1pT_cut && trk2pT>trk2pT_cut && trk3pT>trk3pT_cut &&
           trk1p>trk1p_cut  &&  trk2p>trk2p_cut  &&  trk3p>trk3p_cut  &&
           Sxy>Sxy_cut) h_fake_angle->Fill(refittedAngle) ;
   }
   
   if(DoSCAN) {
      cut_offset = CUT_OFFSET ;
      cut_min = CUT_MIN ;
      cut_max = CUT_MAX ;
      d_cut = cut_max - cut_min ;
      cut_lim = (rows * columns - 0.5) * d_cut + cut_offset ;
      cout << "\n     cut_lim = " << cut_lim << endl ;
      k = 1 ;
      //  Create Canvas
      TCanvas *cscan_fake = new TCanvas("cscan_fake","Vtx chi2/dof  scan (fake candidates)",10,10,1500,800);
      cscan_fake->Divide(columns,rows);
      while(cut_min < cut_lim) {
        h_mass->Reset() ;
        for(int i=0; i<nentries; i++) {
          // Get entry
          fChain->GetEntry(i);
          //  Apply cuts
          bool keep = false ;
          keep = trk1pT>trk1pT_cut && trk2pT>trk2pT_cut && trk3pT>trk3pT_cut &&
                  trk1p>trk1p_cut  &&  trk2p>trk2p_cut  &&  trk3p>trk3p_cut  &&
                 refittedAngle>=cut_min && refittedAngle<cut_max && Sxy>Sxy_cut ;
          if(keep) {
             h_mass->Fill(refittedMassDplus) ;
          }
        }
        // Draw histo
        cscan_fake->cd(k) ;
        cout << "\n  for k = " << k << " , angle = [ " << cut_min << " : " << cut_max << " )" ; 
        h_mass->DrawCopy("e1") ;
        // Update angle and counter
        cut_min += d_cut ;
        cut_max += d_cut ;
        k++ ;
      }
   }


   // Close file 1
   // file1->Close();
   
   //======================================  PLOT HISTOGRAMS  =================================//
   
   gStyle->SetOptFit(1);
   
   Double_t Mmin = LOWER_MASS ;
   Double_t Mmax = HIGHER_MASS ;

/*
   h_diff_massDplus->Sumw2() ;
   h_diff_massDplus->Add(h_massDplus,h_fake_massDplus,1,-1) ;
   
   //  Create Canvas
   TCanvas *c1 = new TCanvas("c1","Mass distribution of D+ -> Kpipi candidates",50,150,600,800);
   c1->Divide(1,2) ;
   c1->cd(1) ;
   h_massDplus->SetLineColor(4);            // h_massDplus->SetLineStyle(1); 
   h_massDplus->Draw("e1"); 
   h_fake_massDplus->SetLineColor(2);       // h_massDplus->SetLineStyle(1); 
   h_fake_massDplus->Draw("e1same"); 
   
   // Fit histogram

   // TF1 *g1 = new TF1("g1","gaus(0)+[3]",Mmin,Mmax);
   TF1 *g1 = new TF1("g1",GaussPol1,Mmin,Mmax,5);
   
   Double_t max1   = h_massDplus->GetMaximum();
   g1->SetParameter(0,max1);
   g1->SetParameter(1,nominal_mass);
   g1->SetParameter(2,0.010);
   // g1->SetParNames("Nevents","Mean","Sigma","Const");
   g1->SetParNames("Nevents","Mean","Sigma","Const","Slope");  //  Only with GaussPol1
   
   g1->SetParLimits(0,0,10*max1);
   g1->SetParLimits(1,1.850,1.880);
   g1->SetParLimits(2,0.002,0.020);

   h_massDplus->Fit("g1","B","SAME",Mmin,Mmax);
   
   c1->cd(2) ;
   h_diff_massDplus->SetLineColor(46);
   h_diff_massDplus->Draw("e1");

   // Fit histogram

   // TF1 *gdif1 = new TF1("gdif1","gaus(0)+[3]",Mmin,Mmax);
   TF1 *gdif1 = new TF1("gdif1",GaussPol1,Mmin,Mmax,5);
   
   Double_t difmax1   = h_diff_massDplus->GetMaximum();
   gdif1->SetParameter(0,difmax1);
   gdif1->SetParameter(1,nominal_mass);
   gdif1->SetParameter(2,0.010);
   // gdif1->SetParNames("Nevents","Mean","Sigma","Const");
   gdif1->SetParNames("Nevents","Mean","Sigma","Const","Slope");  //  Only with GaussPol1
   
   gdif1->SetParLimits(0,0,10*difmax1);
   gdif1->SetParLimits(1,1.850,1.880);
   gdif1->SetParLimits(2,0.002,0.020);
   
   h_diff_massDplus->Fit("gdif1","B","SAME",Mmin,Mmax);
*/
   

   h_diff_refitMassDplus->Sumw2() ;
   h_diff_refitMassDplus->Add(h_refitMassDplus,h_fake_refitMassDplus,1,-1) ;
   
   //  Create Canvas
   TCanvas *c2 = new TCanvas("c2","Refitted mass distribution of D+ -> Kpipi candidates",700,150,600,800);
   c2->Divide(1,2) ;
   c2->cd(1) ;
   h_refitMassDplus->SetLineColor(4);       // h_refitMassDplus->SetLineStyle(1); 
   h_refitMassDplus->Draw("e"); 
   h_fake_refitMassDplus->SetLineColor(2);  // h_refitMassDplus->SetLineStyle(1);
   h_fake_refitMassDplus->SetFillColor(2);
   h_fake_refitMassDplus->SetFillStyle(3001);
   h_fake_refitMassDplus->Draw("same"); 

   // Fit histogram

   // TF1 *g2 = new TF1("g2","gaus(0)+[3]",Mmin,Mmax);
   // TF1 *g2 = new TF1("g2",GaussPol1,Mmin,Mmax,5);
   TF1 *g2 = new TF1("g2",GaussPol2,Mmin,Mmax,6);
   
   Double_t max2   = h_refitMassDplus->GetMaximum();
   g2->SetParameter(0,max2);
   g2->SetParameter(1,nominal_mass);
   g2->SetParameter(2,0.010);
   // g2->SetParNames("Nevents","Mean","Sigma","Const");
   // g2->SetParNames("Nevents","Mean","Sigma","Const","Slope");  //  with GaussPol1
   g2->SetParNames("Nevents","Mean","Sigma","Const","Linear Coef","Quad Coef");  //  with GaussPol2
   
   g2->SetParLimits(0,0,10*max2);
   g2->SetParLimits(1,nominal_mass-0.005,nominal_mass+0.005);
   g2->SetParLimits(2,0.005,0.020);

   h_refitMassDplus->Fit("g2","B","SAME",Mmin,Mmax);
   
   c2->cd(2) ;
   h_diff_refitMassDplus->SetLineColor(46);
   h_diff_refitMassDplus->Draw("e1");

   // Fit histogram

   // TF1 *gdif2 = new TF1("gdif2","gaus(0)+[3]",Mmin,Mmax);
   // TF1 *gdif2 = new TF1("gdif2",GaussPol1,Mmin,Mmax,5);
   TF1 *gdif2 = new TF1("gdif2",GaussPol2,Mmin,Mmax,6);
   
   Double_t difmax2   = h_diff_refitMassDplus->GetMaximum();
   gdif2->SetParameter(0,difmax2);
   gdif2->SetParameter(1,nominal_mass);
   gdif2->SetParameter(2,0.010);
   // gdif2->SetParNames("Nevents","Mean","Sigma","Const");
   // gdif2->SetParNames("Nevents","Mean","Sigma","Const","Slope");  //  Only with GaussPol1
   gdif2->SetParNames("Nevents","Mean","Sigma","Const","Linear Coef","Quad Coef");  //  with GaussPol2
   // gdif2->FixParameter(4,0.);
   
   gdif2->SetParLimits(0,0,10*difmax2);
   gdif2->SetParLimits(1,nominal_mass-0.005,nominal_mass+0.005);
   gdif2->SetParLimits(2,0.005,0.020);

   h_diff_refitMassDplus->Fit("gdif2","BWW","SAME",Mmin,Mmax);
   
   TH1 *h_sig_tmp, *h_nosig_tmp ;
   
   //  Create Canvas
   TCanvas *c3 = new TCanvas("c3","Mass window cut study",700,150,1400,700);
   c3->Divide(4,2) ;

   c3->cd(1) ;
   h_sig_vtxchi2->Sumw2() ;
   h_sig_tmp = h_sig_vtxchi2->DrawNormalized("e1") ;
   h_sig_tmp->SetLineColor(4);
   h_nosig_vtxchi2->Sumw2() ;
   h_nosig_tmp = h_nosig_vtxchi2->DrawNormalized("e1same") ;
   h_nosig_tmp->SetLineColor(2);

   c3->cd(3) ;
   h_sig_dplusp->Sumw2() ;
   h_sig_tmp = h_sig_dplusp->DrawNormalized("e1") ;
   h_sig_tmp->SetLineColor(4);
   h_nosig_dplusp->Sumw2() ;
   h_nosig_tmp = h_nosig_dplusp->DrawNormalized("e1same") ;
   h_nosig_tmp->SetLineColor(2);

   c3->cd(4) ;
   h_nosig_dpluspT->Sumw2() ;
   h_nosig_tmp = h_nosig_dpluspT->DrawNormalized("e1") ;
   h_nosig_tmp->SetLineColor(2);
   h_sig_dpluspT->Sumw2() ;
   h_sig_tmp = h_sig_dpluspT->DrawNormalized("e1same") ;
   h_sig_tmp->SetLineColor(4);

   c3->cd(5) ;
   h_sig_trk1chi2->Sumw2() ;
   h_sig_tmp = h_sig_trk1chi2->DrawNormalized("e1") ;
   h_sig_tmp->SetLineColor(4);
   h_nosig_trk1chi2->Sumw2() ;
   h_nosig_tmp = h_nosig_trk1chi2->DrawNormalized("e1same") ;
   h_nosig_tmp->SetLineColor(2);

   c3->cd(6) ;
   h_sig_trk1hits->Sumw2() ;
   h_sig_tmp = h_sig_trk1hits->DrawNormalized("e1") ;
   h_sig_tmp->SetLineColor(4);
   h_nosig_trk1hits->Sumw2() ;
   h_nosig_tmp = h_nosig_trk1hits->DrawNormalized("e1same") ;
   h_nosig_tmp->SetLineColor(2);

   c3->cd(7) ;
   h_sig_trk1p->Sumw2() ;
   h_sig_tmp = h_sig_trk1p->DrawNormalized("e1") ;
   h_sig_tmp->SetLineColor(4);
   h_nosig_trk1p->Sumw2() ;
   h_nosig_tmp = h_nosig_trk1p->DrawNormalized("e1same") ;
   h_nosig_tmp->SetLineColor(2);

   c3->cd(8) ;
   h_sig_trk1pT->Sumw2() ;
   h_sig_tmp = h_sig_trk1pT->DrawNormalized("e1") ;
   h_sig_tmp->SetLineColor(4);
   h_nosig_trk1pT->Sumw2() ;
   h_nosig_tmp = h_nosig_trk1pT->DrawNormalized("e1same") ;
   h_nosig_tmp->SetLineColor(2);
   
/*   
   //  Create Canvas
   TCanvas *cang = new TCanvas("cang","Angle(refitted)",500,50,600,600);
   h_angle->SetLineColor(4);       // h_angle->SetLineStyle(1); 
   h_angle->Draw("e1"); 
   h_fake_angle->SetLineColor(2);
   h_fake_angle->Draw("e1same"); 

   //  Create Canvas
   TCanvas *cang_fake = new TCanvas("cang_fake","Angle(refitted) - Mass",500,50,600,600);
   h_fake_angle->SetLineColor(4);       // h_fake_angle->SetLineStyle(1); 
   h_fake_angle->Draw("e1"); 
*/

   //  Create Canvas
   TCanvas *c4 = new TCanvas("c4","Refitted mass distribution of D+ -> Kpipi candidates",100,150,600,800);
   
   c4->Divide(1,2) ;
   
   c4->cd(1) ;
   h_refitMassKpi1->SetLineColor(4);       // h_refitMassKpi1->SetLineStyle(1); 
   h_refitMassKpi1->Draw("e"); 
   h_fake_refitMassKpi1->SetLineColor(2);  // h_fake_refitMassKpi1->SetLineStyle(1); 
   h_fake_refitMassKpi1->SetFillColor(2);
   h_fake_refitMassKpi1->SetFillStyle(3001);
   h_fake_refitMassKpi1->Draw("same"); 
   
   c4->cd(2) ;
   h_refitMassKpi2->SetLineColor(4);       // h_refitMassKpi2->SetLineStyle(1); 
   h_refitMassKpi2->Draw("e"); 
   h_fake_refitMassKpi2->SetLineColor(2);  // h_fake_refitMassKpi2->SetLineStyle(1); 
   h_fake_refitMassKpi2->SetFillColor(2);
   h_fake_refitMassKpi2->SetFillStyle(3001);
   h_fake_refitMassKpi2->Draw("same"); 

   
   return 0;
}


/////////////////////////////////////////
//                                     //
//  F(x) = Gaussian(x) + ( a + b*x )   //
//                                     //
/////////////////////////////////////////
Double_t GaussPol1 (Double_t *v , Double_t *par)
{
   Double_t arg = 0 ;
   Double_t normGauss = 1. ;
   // Double_t binwidth = 0.005 ;
   Double_t Area = par[0] ;
   Double_t Mass = par[1] ;
   Double_t Sigma = par[2] ;
   if (Sigma != 0) { 
      arg = (v[0] - Mass)/Sigma ;
      normGauss = sqrt(2*3.14159263) * Sigma ;
   }
//   Double_t fitval = Area*TMath::Exp(-0.5*arg*arg);
   Double_t gauss = (Area*binwidth) * ( exp(-0.5*arg*arg) / normGauss ) ;
   Double_t pol1 = par[3] + par[4]*v[0] ;
   
   return (gauss + pol1) ;
   
}

/////////////////////////////////////////////////
//                                             //
//  F(x) = Gaussian(x) + ( a + b*x + c*x^2 )   //
//                                             //
/////////////////////////////////////////////////
Double_t GaussPol2 (Double_t *v , Double_t *par)
{
   Double_t arg = 0 ;
   Double_t normGauss = 1. ;
   // Double_t binwidth = 0.005 ;
   Double_t Area = par[0] ;
   Double_t Mass = par[1] ;
   Double_t Sigma = par[2] ;
   if (Sigma != 0) { 
      arg = (v[0] - Mass)/Sigma ;
      normGauss = sqrt(2*3.14159263) * Sigma ;
   }
//   Double_t fitval = Area*TMath::Exp(-0.5*arg*arg);
   Double_t gauss = (Area*binwidth) * ( exp(-0.5*arg*arg) / normGauss ) ;
   Double_t pol2 = par[3] + par[4]*v[0] + par[5]*v[0]*v[0] ;
   
   return (gauss + pol2) ;
   
}
