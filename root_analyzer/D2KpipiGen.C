#include <iostream>
#include <fstream>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TChain.h>
#include <TBranch.h>
// #include <TCanvas.h>
// #include <TStyle.h>

#include "load_d2kpipi_MB_MC.C"

void load_d2kpipi_MB_MC(TChain*) ;

int D2KpipiGen()
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
   Double_t        D_x;
   Double_t        D_y;
   Double_t        D_z;
   Double_t        D_xy_length;
   Double_t        D_length;
   Double_t        D2Kpipi_p;
   Double_t        D2Kpipi_pT;
   Double_t        D2Kpipi_phi;
   Double_t        D2Kpipi_eta;
   Double_t        D2Kpipi_pTfraction;
   Double_t        D2Kpipi_pT2fraction;
   Double_t        mean_pT;
   Double_t        mean_pT2;
   Double_t        nGenPart;
   Double_t        K_p;
   Double_t        K_pT;
   Double_t        K_phi;
   Double_t        K_eta;
   Double_t        pi1_p;
   Double_t        pi1_pT;
   Double_t        pi1_phi;
   Double_t        pi1_eta;
   Double_t        pi2_p;
   Double_t        pi2_pT;
   Double_t        pi2_phi;
   Double_t        pi2_eta;

   // List of branches
   TBranch        *b_D_x;   //!
   TBranch        *b_D_y;   //!
   TBranch        *b_D_z;   //!
   TBranch        *b_D_xy_length;   //!
   TBranch        *b_D_length;   //!
   TBranch        *b_D2Kpipi_p;   //!
   TBranch        *b_D2Kpipi_pT;   //!
   TBranch        *b_D2Kpipi_phi;   //!
   TBranch        *b_D2Kpipi_eta;   //!
   TBranch        *b_D2Kpipi_pTfraction;   //!
   TBranch        *b_D2Kpipi_pT2fraction;   //!
   TBranch        *b_mean_pT;   //!
   TBranch        *b_mean_pT2;   //!
   TBranch        *b_nGenPart;   //!
   TBranch        *b_K_p;   //!
   TBranch        *b_K_pT;   //!
   TBranch        *b_K_phi;   //!
   TBranch        *b_K_eta;   //!
   TBranch        *b_pi1_p;   //!
   TBranch        *b_pi1_pT;   //!
   TBranch        *b_pi1_phi;   //!
   TBranch        *b_pi1_eta;   //!
   TBranch        *b_pi2_p;   //!
   TBranch        *b_pi2_pT;   //!
   TBranch        *b_pi2_phi;   //!
   TBranch        *b_pi2_eta;   //!

   // Histograms
   
   //  Pass the cuts through these variables

   //  Ewerton cuts
   double trketa_cut = 2.5 ;    //  std cut = 2.4
   double trk1p_cut = 0.5 ;       //  no cut = 0.5
   double trk2p_cut = 0.5 ;       //  no cut = 0.5
   double trk3p_cut = 0.5 ;       //  no cut = 0.5
   double K_pT_cut = 0.5 ;        //  no cut = 0.5 
   double pi1_pT_cut = 0.5 ;      //  no cut = 0.5
   double pi2_pT_cut = 0.5 ;      //  no cut = 0.5
   double pT_1st_cut = 0.5 ;    //  std cut = 1.7 
   double pT_2nd_cut = 0.5 ;    //  std cut = 1.5 
   double pT_3rd_cut = 0.5 ;    //  std cut = 0.9 
   double dpluseta_cut = 12 ;    //  no cut = 12.
   double dplusp_cut = 15000 ;    //  no cut = 1000.
   double dpluspT_cut = 3.5 ;   //  std cut = 4.3
   
   // bool keep = false ;
   
   // Set here the TTree name
   Char_t *goodTree = "d2kpipi/genTree";
   
   // Open file
   // cout << "\n    ***  Opening file \"" << fileName1 << "\"  ***" << endl;
   // TFile *file1 = new TFile(fileName1);


   // fChain = (TTree*)file1->Get(goodTree);
   fChain = new TChain(goodTree) ;
   load_d2kpipi_MB_MC(fChain) ;
   // Init(fChain);
   if (!fChain) return 1;
   
   // Set branch addresses and branch pointers
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("D_x", &D_x, &b_D_x);
   fChain->SetBranchAddress("D_y", &D_y, &b_D_y);
   fChain->SetBranchAddress("D_z", &D_z, &b_D_z);
   fChain->SetBranchAddress("D_xy_length", &D_xy_length, &b_D_xy_length);
   fChain->SetBranchAddress("D_length", &D_length, &b_D_length);
   fChain->SetBranchAddress("D2Kpipi_p", &D2Kpipi_p, &b_D2Kpipi_p);
   fChain->SetBranchAddress("D2Kpipi_pT", &D2Kpipi_pT, &b_D2Kpipi_pT);
   fChain->SetBranchAddress("D2Kpipi_phi", &D2Kpipi_phi, &b_D2Kpipi_phi);
   fChain->SetBranchAddress("D2Kpipi_eta", &D2Kpipi_eta, &b_D2Kpipi_eta);
   fChain->SetBranchAddress("D2Kpipi_pTfraction", &D2Kpipi_pTfraction, &b_D2Kpipi_pTfraction);
   fChain->SetBranchAddress("D2Kpipi_pT2fraction", &D2Kpipi_pT2fraction, &b_D2Kpipi_pT2fraction);
   fChain->SetBranchAddress("mean_pT", &mean_pT, &b_mean_pT);
   fChain->SetBranchAddress("mean_pT2", &mean_pT2, &b_mean_pT2);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("K_p", &K_p, &b_K_p);
   fChain->SetBranchAddress("K_pT", &K_pT, &b_K_pT);
   fChain->SetBranchAddress("K_phi", &K_phi, &b_K_phi);
   fChain->SetBranchAddress("K_eta", &K_eta, &b_K_eta);
   fChain->SetBranchAddress("pi1_p", &pi1_p, &b_pi1_p);
   fChain->SetBranchAddress("pi1_pT", &pi1_pT, &b_pi1_pT);
   fChain->SetBranchAddress("pi1_phi", &pi1_phi, &b_pi1_phi);
   fChain->SetBranchAddress("pi1_eta", &pi1_eta, &b_pi1_eta);
   fChain->SetBranchAddress("pi2_p", &pi2_p, &b_pi2_p);
   fChain->SetBranchAddress("pi2_pT", &pi2_pT, &b_pi2_pT);
   fChain->SetBranchAddress("pi2_phi", &pi2_phi, &b_pi2_phi);
   fChain->SetBranchAddress("pi2_eta", &pi2_eta, &b_pi2_eta);

   int nentries = fChain->GetEntries();
   
   cout << " Nentries = " << nentries << endl;
   
   double N_accepted = 0 ;
   double N_eta = 0 ;
   
   bool keep = false ;
   double pT_1st , pT_2nd , pT_3rd ;
   for(int i=0; i<nentries; i++) {
      // Get entry
      fChain->GetEntry(i);
      //  Apply cuts
      keep = false ;
      pT_1st = pT_2nd = pT_3rd = 0. ;
      if(K_pT>pi1_pT) {
         if(K_pT>pi2_pT) {
            pT_1st = K_pT ;
            if(pi1_pT>pi2_pT) {
               pT_2nd = pi1_pT ;
               pT_3rd = pi2_pT ;
            } else {
               pT_2nd = pi2_pT ;
               pT_3rd = pi1_pT ;
            }
         } else { 
            pT_1st = pi2_pT ;
            pT_2nd = K_pT ;
            pT_3rd = pi1_pT ;
         }
      } else {
         if(pi1_pT>pi2_pT) {
            pT_1st = pi1_pT ;
            if(K_pT>pi2_pT) {
               pT_2nd = K_pT ;
               pT_3rd = pi2_pT ;
            } else {
               pT_2nd = pi2_pT ;
               pT_3rd = K_pT ;
            }
         } else { 
            pT_1st = pi2_pT ;
            pT_2nd = pi1_pT ;
            pT_3rd = K_pT ;
         }
      }
      
      if(fabs(K_eta)<trketa_cut && fabs(pi1_eta)<trketa_cut && fabs(pi2_eta)<trketa_cut && 
         fabs(D2Kpipi_eta)<dpluseta_cut) N_eta++ ;
      
//      keep = K_pT>K_pT_cut && pi1_pT>pi1_pT_cut && pi2_pT>pi2_pT_cut &&
//              trk1p>trk1p_cut  &&  trk2p>trk2p_cut  &&  trk3p>trk3p_cut  &&
      keep = pT_1st>pT_1st_cut       && pT_2nd>pT_2nd_cut         && pT_3rd>pT_3rd_cut         &&
             fabs(K_eta)<trketa_cut  && fabs(pi1_eta)<trketa_cut  && fabs(pi2_eta)<trketa_cut  &&
             D2Kpipi_pT>dpluspT_cut  && D2Kpipi_p<dplusp_cut      && fabs(D2Kpipi_eta)<dpluseta_cut ;
      if(keep) {
         N_accepted++ ;
      }
   }

   cout << " N_initial = " << nentries 
        << "\n N_accepted = " << N_accepted 
        << "\n Acceptance = " << ((float) N_accepted)/((float) nentries) 
        << "\n Efficiency = " << ((float) N_accepted)/((float) N_eta) 
        << endl;
   
   return 0;
}

