using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"

const float PI = TMath::Pi();

// EPOS particle ids
const int ID_pip = 120;    // pi+
const int ID_pim = -120;   // pi-
const int ID_pi0 = 110;    // pi0

// track cuts
const float etaCut = 1.0;  // |eta| < etaCut for v2(pT)
const float ptMin  = 0.05; // minimum pT

void Elliptic(int cen = 5, int job = 0){	//main_function

    float Eweight = 1;
    TChain* chain = new TChain("teposevent");

    chain->Add("z-*.root");
    //chain->Add("../auau200-6/z-*.root");

    char fname_out[200];
    sprintf(fname_out,"cen%d.v2pT_pion_job%d.root",cen,job);
    TFile fout(fname_out,"RECREATE");

    //defining histograms
    // v2(pT) = <cos(2(phi - Psi))> vs pT, w.r.t. reaction plane (phir) and participant plane (phir+psi2)
    TProfile* pV2Pt_pip = new TProfile("pV2Pt_pip","v_{2}(p_{T}) #pi^{+} (RP)",60,0,3);
    TProfile* pV2Pt_pim = new TProfile("pV2Pt_pim","v_{2}(p_{T}) #pi^{-} (RP)",60,0,3);
    TProfile* pV2Pt_pi  = new TProfile("pV2Pt_pi" ,"v_{2}(p_{T}) #pi^{#pm} (RP)",60,0,3);
    TProfile* pV2Pt_pi0 = new TProfile("pV2Pt_pi0","v_{2}(p_{T}) #pi^{0} (RP)",60,0,3);

    TProfile* pV2Pt_pi_PP = new TProfile("pV2Pt_pi_PP","v_{2}(p_{T}) #pi^{#pm} (PP)",60,0,3);

    // v2 vs eta (integrated over pT) and pT-integrated v2 per event
    TProfile* pV2Eta_pi = new TProfile("pV2Eta_pi","v_{2}(#eta) #pi^{#pm} (RP)",26,-6,6);

    // QA histograms
    TH1D* Hist_Pt_pi  = new TH1D("Hist_Pt_pi","p_{T} #pi^{#pm}",60,0,3);
    TH1D* Hist_dPhi_pi = new TH1D("Hist_dPhi_pi","#phi - #Psi_{RP} #pi^{#pm}",72,-PI,PI);
    TH1D* Hist_b = new TH1D("Hist_b","Hist_b",50,0,20);
    TH1D* Hist_Psi = new TH1D("Hist_Psi","#Psi_{RP} (phir)",72,-PI,PI);

    Int_t nentries = chain->GetEntries();
    cout << nentries << "\n";
    for(int i = 0; i < nentries; i++){

        if((i+1)%1000==0) cout << "Processing entry == "<< i+1 <<" == out of "<<nentries<<".\n";
        chain->GetEntry(i);

        TLeaf* leaf_b   = chain->GetLeaf("bim");
        float bim = leaf_b->GetValue(0);
        if(bim < 8.4 || bim > 9.2){	//30-40% centrality
            continue;
        }
        TLeaf* leaf_NoTracks = chain->GetLeaf("np");
        TLeaf* leaf_Np_p= chain->GetLeaf("npartproj");
        TLeaf* leaf_Np_t= chain->GetLeaf("nparttarg");
        int Np = leaf_Np_p->GetValue(0) + leaf_Np_t->GetValue(0);
        int NPTracks = (int)leaf_NoTracks->GetValue(0);

        // event plane: "phi" = random rotation of the event in the lab (reaction plane),
        // "phir" = fluctuation of the participant plane relative to the reaction plane
        TLeaf* leaf_phi  = chain->GetLeaf("phi");
        TLeaf* leaf_phir = chain->GetLeaf("phir");
        float PsiRP = leaf_phi->GetValue(0);
        float PsiPP = PsiRP + leaf_phir->GetValue(0);

        Hist_b->Fill(bim);
        Hist_Psi->Fill(TVector2::Phi_mpi_pi(PsiRP));

        TLeaf* leaf_PxV0       = chain->GetLeaf("px");
        TLeaf* leaf_PyV0       = chain->GetLeaf("py");
        TLeaf* leaf_PzV0       = chain->GetLeaf("pz");
        TLeaf* leaf_PIDV0    = chain->GetLeaf("id");
        TLeaf* leaf_ist        = chain->GetLeaf("ist");

        for(int trk = 0; trk < NPTracks; trk++) {

            int istAsso = (int)leaf_ist->GetValue(trk);
            int pid = (int)leaf_PIDV0->GetValue(trk);
            //charged pions: final state (ist==0); pi0: all decay, so take ist==1 (decayed) instead
            bool isChargedPi = (istAsso == 0) && (pid == ID_pip || pid == ID_pim);
            bool isPi0       = (istAsso == 1) && (pid == ID_pi0);
            if(!isChargedPi && !isPi0){
                continue;
            }

            float PxAsso    = leaf_PxV0->GetValue(trk);
            float PyAsso    = leaf_PyV0->GetValue(trk);
            float PzAsso    = leaf_PzV0->GetValue(trk);

            float PtAsso = sqrt(PxAsso*PxAsso+PyAsso*PyAsso);
            float PhiAsso = atan2(PyAsso,PxAsso);
            float ThetaAsso = atan2(PtAsso,PzAsso);
            float EtaAsso = -log(tan(ThetaAsso/2.));

            if(PtAsso < ptMin) continue;

            float cos2dPhi_RP = cos(2.*(PhiAsso - PsiRP));
            float cos2dPhi_PP = cos(2.*(PhiAsso - PsiPP));

            if(pid == ID_pip || pid == ID_pim){
                pV2Eta_pi->Fill(EtaAsso, cos2dPhi_RP, Eweight);
            }

            if(EtaAsso > etaCut || EtaAsso < -etaCut) continue;

            if(pid == ID_pip){
                pV2Pt_pip->Fill(PtAsso, cos2dPhi_RP, Eweight);
            }
            if(pid == ID_pim){
                pV2Pt_pim->Fill(PtAsso, cos2dPhi_RP, Eweight);
            }
            if(pid == ID_pip || pid == ID_pim){
                pV2Pt_pi->Fill(PtAsso, cos2dPhi_RP, Eweight);
                pV2Pt_pi_PP->Fill(PtAsso, cos2dPhi_PP, Eweight);
                Hist_Pt_pi->Fill(PtAsso, Eweight);
                Hist_dPhi_pi->Fill(TVector2::Phi_mpi_pi(PhiAsso - PsiRP), Eweight);
            }
            if(pid == ID_pi0){
                pV2Pt_pi0->Fill(PtAsso, cos2dPhi_RP, Eweight);
            }
        }
    }
    fout.Write();
    return;
}
