using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TMath.h"
#include "TProfile.h"
#include "TVector2.h"

const float PI = TMath::Pi();

// EPOS particle ids (see idt.dt): CHARGED MESONS = pi+/- and K+/- (baryons excluded),
// matching STAR's default POI (Figs 5-15 of HEPData ins2928164; Fig 16 is the exception
// that adds p/pbar). Charge separation is by sign of id (= electric charge sign for these).
const int ID_pip = 120;    // pi+
const int ID_pim = -120;   // pi-
const int ID_Kp  = 130;    // K+
const int ID_Km  = -130;   // K-
const float MPI_CHG = 0.13957;   // charged-pion mass (GeV)
const float MKA_CHG = 0.493677;  // charged-kaon mass (GeV) -- needed for the |y|<1 cut

// POI cuts, matching the STAR/ESS analysis: |y|<1, 0.2<pT<2, true reaction plane
const float yCut  = 1.0;
const float ptMin = 0.2;
const float ptMax = 2.0;

// same observables as Correlators.C (v2, gamma_OS/SS, delta_OS/SS and the OS-SS differences),
// but the POI set is charged mesons (pi+K). Charge-separated Q-vectors are built over all
// positive mesons (pi+,K+) and all negative mesons (pi-,K-); species enters only through the
// particle mass used to convert to rapidity for the |y|<1 acceptance cut.

static inline float MesonMass(int pid){
    int a = abs(pid);
    if(a == 130) return MKA_CHG;   // kaon
    return MPI_CHG;                // pion (default)
}

void Correlators_mesons(int cen = 5, int job = 0){

    TChain* chain = new TChain("teposevent");
    chain->Add("z-auau_run_*.root");

    //only read the branches we use
    chain->SetBranchStatus("*",0);
    const char* activeBranches[] = {"np","bim","phi","px","py","pz","id","ist"};
    for(int ib = 0; ib < 8; ib++) chain->SetBranchStatus(activeBranches[ib],1);

    char fname_out[200];
    sprintf(fname_out,"cen%d.correlators_mesons_job%d.root",cen,job);
    TFile fout(fname_out,"RECREATE");

    //----- differential elliptic flow (per-particle fill; shape/QA) -----
    TProfile* pV2pT  = new TProfile("pV2pT" ,"v_{2}(p_{T}) ch. mesons (RP);p_{T} (GeV/c);v_{2}",60,0,3);
    TProfile* pV2eta = new TProfile("pV2eta","v_{2}(#eta) ch. mesons (RP);#eta;v_{2}",30,-3,3);

    //----- integrated observables, filled ONCE PER EVENT (error = event-to-event spread) -----
    TProfile* pV2int   = new TProfile("pV2int"  ,"integrated v_{2} ch. mesons",1,0,1);
    TProfile* pGammaOS = new TProfile("pGammaOS","#gamma_{OS}",1,0,1);
    TProfile* pGammaSS = new TProfile("pGammaSS","#gamma_{SS}",1,0,1);
    TProfile* pDGamma  = new TProfile("pDGamma" ,"#Delta#gamma = #gamma_{OS}-#gamma_{SS}",1,0,1);
    TProfile* pDeltaOS = new TProfile("pDeltaOS","#delta_{OS}",1,0,1);
    TProfile* pDeltaSS = new TProfile("pDeltaSS","#delta_{SS}",1,0,1);
    TProfile* pDDelta  = new TProfile("pDDelta" ,"#Delta#delta = #delta_{OS}-#delta_{SS}",1,0,1);

    //pair-weighted accumulators (the paper's pair+event average convention) for cross-check
    double tNgOS=0, tNgSS=0, tNdOS=0, tNdSS=0, tDOS=0, tDSS=0;
    double sumNp=0, sumNm=0;
    Long64_t nUsed = 0;

    Int_t nentries = chain->GetEntries();
    cout << nentries << " entries in chain\n";

    for(int i = 0; i < nentries; i++){

        if((i+1)%1000==0) cout << "Processing entry == "<< i+1 <<" == out of "<<nentries<<".\n";
        chain->GetEntry(i);

        float bim = chain->GetLeaf("bim")->GetValue(0);
        if(bim < 8.4 || bim > 9.2) continue;	//30-40% centrality (float compare!)

        int NPTracks = (int)chain->GetLeaf("np")->GetValue(0);
        float PsiRP  = chain->GetLeaf("phi")->GetValue(0);	//true reaction plane
        float C2 = cos(2.*PsiRP), S2 = sin(2.*PsiRP);

        TLeaf* leaf_px  = chain->GetLeaf("px");
        TLeaf* leaf_py  = chain->GetLeaf("py");
        TLeaf* leaf_pz  = chain->GetLeaf("pz");
        TLeaf* leaf_id  = chain->GetLeaf("id");
        TLeaf* leaf_ist = chain->GetLeaf("ist");

        //charge-separated Q-vectors over POI charged mesons
        double pP=0, qP=0, p2P=0, q2P=0;  int Np=0;   // positive mesons (pi+, K+)
        double pM=0, qM=0, p2M=0, q2M=0;  int Nm=0;   // negative mesons (pi-, K-)

        for(int trk = 0; trk < NPTracks; trk++) {
            if((int)leaf_ist->GetValue(trk) != 0) continue;	//final-state only
            int pid = (int)leaf_id->GetValue(trk);
            int a = abs(pid);
            if(!(a == 120 || a == 130)) continue;	//charged mesons only (pi+/-, K+/-)

            float px = leaf_px->GetValue(trk);
            float py = leaf_py->GetValue(trk);
            float pz = leaf_pz->GetValue(trk);

            float pt = sqrt(px*px+py*py);
            float phi = atan2(py,px);
            float theta = atan2(pt,pz);
            float eta = -log(tan(theta/2.));
            float m = MesonMass(pid);
            float E = sqrt(m*m+px*px+py*py+pz*pz);
            float y = 0.5*log((E+pz)/(E-pz));

            float cos2dphi = cos(2.*(phi - PsiRP));

            //differential v2 (loosened one cut each so the shapes are informative)
            if(pt > ptMin && pt < ptMax) pV2eta->Fill(eta, cos2dphi);
            if(fabs(y) < yCut)           pV2pT->Fill(pt,  cos2dphi);

            //POI selection for the integrated observables and correlators
            if(pt < ptMin || pt > ptMax) continue;
            if(fabs(y) > yCut) continue;

            double c1 = cos(phi),   s1 = sin(phi);
            double c2 = cos(2.*phi), s2 = sin(2.*phi);
            if(pid > 0){ pP+=c1; qP+=s1; p2P+=c2; q2P+=s2; Np++; }  // positive meson
            else       { pM+=c1; qM+=s1; p2M+=c2; q2M+=s2; Nm++; }  // negative meson
        }

        if(Np < 2 || Nm < 2) continue;	//need >=2 of each charge for all pair classes

        double NpairOS = (double)Np*Nm;
        double NpairSS = 0.5*Np*(Np-1) + 0.5*Nm*(Nm-1);

        //---- gamma numerators ----
        double NumgOS = (pP*pM - qP*qM)*C2 + (qP*pM + pP*qM)*S2;
        double NumgSS = 0.5*( (pP*pP - qP*qP - p2P)*C2 + (2.*pP*qP - q2P)*S2 )
                      + 0.5*( (pM*pM - qM*qM - p2M)*C2 + (2.*pM*qM - q2M)*S2 );
        //---- delta numerators ----
        double NumdOS = pP*pM + qP*qM;
        double NumdSS = 0.5*(pP*pP + qP*qP - Np) + 0.5*(pM*pM + qM*qM - Nm);

        double gOS = NumgOS/NpairOS, gSS = NumgSS/NpairSS;
        double dOS = NumdOS/NpairOS, dSS = NumdSS/NpairSS;

        //integrated v2 = <cos2(phi-Psi)> over the POI (Sum cos2phi etc. already in Q-vectors)
        double v2ev = ((p2P+p2M)*C2 + (q2P+q2M)*S2) / (double)(Np+Nm);

        //fill per-event (single-bin) profiles -> GetBinError = event-level statistical error
        pV2int  ->Fill(0.5, v2ev);
        pGammaOS->Fill(0.5, gOS);
        pGammaSS->Fill(0.5, gSS);
        pDGamma ->Fill(0.5, gOS - gSS);
        pDeltaOS->Fill(0.5, dOS);
        pDeltaSS->Fill(0.5, dSS);
        pDDelta ->Fill(0.5, dOS - dSS);

        //pair-weighted accumulators
        tNgOS += NumgOS; tNdOS += NumdOS; tDOS += NpairOS;
        tNgSS += NumgSS; tNdSS += NumdSS; tDSS += NpairSS;
        sumNp += Np; sumNm += Nm; nUsed++;
    }

    //store metadata + pair-weighted central values for the figure/projection macro
    TH1D* hMeta = new TH1D("hMeta","meta: 1 Nev, 2 <N+>, 3 <N->, 4 gOS_pw, 5 gSS_pw, 6 dOS_pw, 7 dSS_pw",7,0.5,7.5);
    hMeta->SetBinContent(1, (double)nUsed);
    hMeta->SetBinContent(2, nUsed? sumNp/nUsed : 0);
    hMeta->SetBinContent(3, nUsed? sumNm/nUsed : 0);
    hMeta->SetBinContent(4, tDOS? tNgOS/tDOS : 0);
    hMeta->SetBinContent(5, tDSS? tNgSS/tDSS : 0);
    hMeta->SetBinContent(6, tDOS? tNdOS/tDOS : 0);
    hMeta->SetBinContent(7, tDSS? tNdSS/tDSS : 0);

    cout << "\n==== correlators (30-40% Au+Au 200 GeV, CHARGED MESONS pi+K, true RP) ====\n";
    cout << "events used            = " << nUsed << "\n";
    cout << "<N+>, <N-> per event   = " << (nUsed? sumNp/nUsed:0) << ", " << (nUsed? sumNm/nUsed:0) << "\n\n";
    printf("  %-10s %14s %14s\n","observable","per-event","pair-weighted");
    printf("  %-10s %10.4e +/- %8.1e   (%10.4e)\n","v2",      pV2int  ->GetBinContent(1), pV2int  ->GetBinError(1), 0.0);
    printf("  %-10s %10.4e +/- %8.1e   %10.4e\n","gamma_OS", pGammaOS->GetBinContent(1), pGammaOS->GetBinError(1), tNgOS/tDOS);
    printf("  %-10s %10.4e +/- %8.1e   %10.4e\n","gamma_SS", pGammaSS->GetBinContent(1), pGammaSS->GetBinError(1), tNgSS/tDSS);
    printf("  %-10s %10.4e +/- %8.1e\n",         "Dgamma",   pDGamma ->GetBinContent(1), pDGamma ->GetBinError(1));
    printf("  %-10s %10.4e +/- %8.1e   %10.4e\n","delta_OS", pDeltaOS->GetBinContent(1), pDeltaOS->GetBinError(1), tNdOS/tDOS);
    printf("  %-10s %10.4e +/- %8.1e   %10.4e\n","delta_SS", pDeltaSS->GetBinContent(1), pDeltaSS->GetBinError(1), tNdSS/tDSS);
    printf("  %-10s %10.4e +/- %8.1e\n",         "Ddelta",   pDDelta ->GetBinContent(1), pDDelta ->GetBinError(1));

    fout.Write();
    return;
}
