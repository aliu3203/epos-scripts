using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <vector>
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

// EPOS particle ids (see analyze.C / idt.dt)
const int ID_pip  = 120;    // pi+
const int ID_pim  = -120;   // pi-
const int ID_Kp   = 130;    // K+
const int ID_Km   = -130;   // K-
const int ID_p    = 1120;   // proton
const int ID_pbar = -1120;  // antiproton

// POI (particle of interest) cuts
const float etaCut = 1.0;   // |eta| < etaCut
const float ptMin  = 0.05;  // minimum pT

// EP (event plane) track cuts
const float ep_ptMin  = 0.2;   // 0.2 < pT < 2 for Q-vector tracks
const float ep_ptMax  = 2.0;
const float ep_etaMax = 1.0;   // |eta| < 1
const float ep_etaGap = 0.05;  // sub-events: A = (-1,-0.05), B = (0.05,1)

void Elliptic_EP(int cen = 5, int job = 0){	//main_function

    float Eweight = 1;
    TChain* chain = new TChain("teposevent");

    chain->Add("z-*.root");
    //chain->Add("../auau200-6/z-*.root");

    //only read the branches we use -- skips decompressing x,y,z,t,zus,... (about half the file volume)
    chain->SetBranchStatus("*",0);
    const char* activeBranches[] = {"np","bim","npartproj","nparttarg","phi","phir","px","py","pz","id","ist"};
    for(int ib = 0; ib < 11; ib++) chain->SetBranchStatus(activeBranches[ib],1);

    char fname_out[200];
    sprintf(fname_out,"cen%d.v2pT_EP_job%d.root",cen,job);
    TFile fout(fname_out,"RECREATE");

    //defining histograms
    // raw v2(pT) = <cos(2(phi - Psi2))> vs pT; resolution correction is applied in MakeFigure_Elliptic_EP.C
    TProfile* pV2Pt_pi_EP  = new TProfile("pV2Pt_pi_EP" ,"raw v_{2}(p_{T}) #pi^{#pm}, full EP",60,0,3);
    TProfile* pV2Pt_K_EP   = new TProfile("pV2Pt_K_EP"  ,"raw v_{2}(p_{T}) K^{#pm}, full EP",60,0,3);
    TProfile* pV2Pt_pr_EP  = new TProfile("pV2Pt_pr_EP" ,"raw v_{2}(p_{T}) p+#bar{p}, full EP",60,0,3);

    TProfile* pV2Pt_pi_sub = new TProfile("pV2Pt_pi_sub","raw v_{2}(p_{T}) #pi^{#pm}, #eta-sub EP",60,0,3);
    TProfile* pV2Pt_K_sub  = new TProfile("pV2Pt_K_sub" ,"raw v_{2}(p_{T}) K^{#pm}, #eta-sub EP",60,0,3);
    TProfile* pV2Pt_pr_sub = new TProfile("pV2Pt_pr_sub","raw v_{2}(p_{T}) p+#bar{p}, #eta-sub EP",60,0,3);

    // truth reference: true reaction plane from the generator ("phi" event branch)
    TProfile* pV2Pt_pi_RP  = new TProfile("pV2Pt_pi_RP" ,"v_{2}(p_{T}) #pi^{#pm} (true RP)",60,0,3);
    TProfile* pV2Pt_K_RP   = new TProfile("pV2Pt_K_RP"  ,"v_{2}(p_{T}) K^{#pm} (true RP)",60,0,3);
    TProfile* pV2Pt_pr_RP  = new TProfile("pV2Pt_pr_RP" ,"v_{2}(p_{T}) p+#bar{p} (true RP)",60,0,3);

    // event-plane resolution terms
    // bin 1: <cos2(PsiA-PsiB)>   bin 2: <cos2(PsiFull-PsiRP)>
    // bin 3: <cos2(PsiA-PsiRP)>  bin 4: <cos2(PsiB-PsiRP)>
    TProfile* pRes = new TProfile("pRes","EP resolution terms",4,0.5,4.5);

    // QA
    TH1D* Hist_PsiFull = new TH1D("Hist_PsiFull","#Psi_{2}^{full} - #Psi_{RP}",72,-PI/2,PI/2);
    TH1D* Hist_PsiAB   = new TH1D("Hist_PsiAB","#Psi_{2}^{A} - #Psi_{2}^{B}",72,-PI/2,PI/2);
    TH1D* Hist_nEP     = new TH1D("Hist_nEP","N tracks in full Q-vector",300,-0.5,2999.5);
    TH1D* Hist_b       = new TH1D("Hist_b","Hist_b",50,0,20);

    // per-event POI storage (filled in the track loop, used after the EP is built)
    std::vector<float> poi_pt, poi_eta, poi_c2, poi_s2;
    std::vector<int>   poi_pid;
    std::vector<bool>  poi_inQ;

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
        int NPTracks = (int)leaf_NoTracks->GetValue(0);

        TLeaf* leaf_phi  = chain->GetLeaf("phi");
        float PsiRP = leaf_phi->GetValue(0);	//true reaction plane (generator)

        TLeaf* leaf_PxV0  = chain->GetLeaf("px");
        TLeaf* leaf_PyV0  = chain->GetLeaf("py");
        TLeaf* leaf_PzV0  = chain->GetLeaf("pz");
        TLeaf* leaf_PIDV0 = chain->GetLeaf("id");
        TLeaf* leaf_ist   = chain->GetLeaf("ist");

        //---- pass 1: build the Q-vectors and collect the POIs ----
        float QxFull = 0, QyFull = 0;
        float QxA = 0, QyA = 0, QxB = 0, QyB = 0;
        int nFull = 0, nA = 0, nB = 0;
        poi_pt.clear(); poi_eta.clear(); poi_c2.clear(); poi_s2.clear();
        poi_pid.clear(); poi_inQ.clear();

        for(int trk = 0; trk < NPTracks; trk++) {

            int istAsso = (int)leaf_ist->GetValue(trk);
            if (istAsso != 0) continue;	//final state particles only

            int pid = (int)leaf_PIDV0->GetValue(trk);
            //charged pi/K/p: the dominant charged species (same selection as analyze.C)
            if(!(pid == ID_pip || pid == ID_pim || pid == ID_Kp || pid == ID_Km
              || pid == ID_p   || pid == ID_pbar)){
                continue;
            }

            float PxAsso = leaf_PxV0->GetValue(trk);
            float PyAsso = leaf_PyV0->GetValue(trk);
            float PzAsso = leaf_PzV0->GetValue(trk);

            float PtAsso = sqrt(PxAsso*PxAsso+PyAsso*PyAsso);
            float PhiAsso = atan2(PyAsso,PxAsso);
            float ThetaAsso = atan2(PtAsso,PzAsso);
            float EtaAsso = -log(tan(ThetaAsso/2.));

            float c2 = cos(2.*PhiAsso);
            float s2 = sin(2.*PhiAsso);

            //Q-vector tracks
            bool inQ = (PtAsso > ep_ptMin && PtAsso < ep_ptMax && fabs(EtaAsso) < ep_etaMax);
            if(inQ){
                QxFull += c2; QyFull += s2; nFull++;
                if(EtaAsso < -ep_etaGap) { QxA += c2; QyA += s2; nA++; }
                if(EtaAsso >  ep_etaGap) { QxB += c2; QyB += s2; nB++; }
            }

            //POIs: same species, |eta|<1, pT>0.05
            if(PtAsso > ptMin && fabs(EtaAsso) < etaCut){
                poi_pt.push_back(PtAsso);
                poi_eta.push_back(EtaAsso);
                poi_c2.push_back(c2);
                poi_s2.push_back(s2);
                poi_pid.push_back(pid);
                poi_inQ.push_back(inQ);
            }
        }

        if(nA < 2 || nB < 2 || nFull < 4) continue;

        float PsiA = 0.5*atan2(QyA,QxA);
        float PsiB = 0.5*atan2(QyB,QxB);
        float PsiFull = 0.5*atan2(QyFull,QxFull);

        //resolution terms
        pRes->Fill(1., cos(2.*(PsiA-PsiB)));
        pRes->Fill(2., cos(2.*(PsiFull-PsiRP)));
        pRes->Fill(3., cos(2.*(PsiA-PsiRP)));
        pRes->Fill(4., cos(2.*(PsiB-PsiRP)));

        Hist_b->Fill(bim);
        Hist_PsiFull->Fill(TVector2::Phi_mpi_pi(2.*(PsiFull-PsiRP))/2.);
        Hist_PsiAB->Fill(TVector2::Phi_mpi_pi(2.*(PsiA-PsiB))/2.);
        Hist_nEP->Fill(nFull);

        //---- pass 2: correlate the POIs with the event planes ----
        for(unsigned int ip = 0; ip < poi_pt.size(); ip++) {

            int pid = poi_pid[ip];
            float PtAsso = poi_pt[ip];
            float PhiAsso = 0.5*atan2(poi_s2[ip],poi_c2[ip]);	//phi mod pi is enough for cos(2(phi-Psi))

            //full EP, removing this particle's own contribution from Q (no autocorrelation)
            float qx = QxFull, qy = QyFull;
            if(poi_inQ[ip]) { qx -= poi_c2[ip]; qy -= poi_s2[ip]; }
            float PsiExcl = 0.5*atan2(qy,qx);
            float cos2dPhi_EP = cos(2.*(PhiAsso - PsiExcl));

            //eta-sub EP: use the sub-event in the opposite hemisphere (no autocorrelation by construction)
            float PsiSub = (poi_eta[ip] >= 0) ? PsiA : PsiB;
            float cos2dPhi_sub = cos(2.*(PhiAsso - PsiSub));

            //truth
            float cos2dPhi_RP = cos(2.*(PhiAsso - PsiRP));

            if(pid == ID_pip || pid == ID_pim){
                pV2Pt_pi_EP->Fill(PtAsso, cos2dPhi_EP, Eweight);
                pV2Pt_pi_sub->Fill(PtAsso, cos2dPhi_sub, Eweight);
                pV2Pt_pi_RP->Fill(PtAsso, cos2dPhi_RP, Eweight);
            }
            if(pid == ID_Kp || pid == ID_Km){
                pV2Pt_K_EP->Fill(PtAsso, cos2dPhi_EP, Eweight);
                pV2Pt_K_sub->Fill(PtAsso, cos2dPhi_sub, Eweight);
                pV2Pt_K_RP->Fill(PtAsso, cos2dPhi_RP, Eweight);
            }
            if(pid == ID_p || pid == ID_pbar){
                pV2Pt_pr_EP->Fill(PtAsso, cos2dPhi_EP, Eweight);
                pV2Pt_pr_sub->Fill(PtAsso, cos2dPhi_sub, Eweight);
                pV2Pt_pr_RP->Fill(PtAsso, cos2dPhi_RP, Eweight);
            }
        }
    }
    fout.Write();
    return;
}
