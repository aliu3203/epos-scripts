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

// EPOS particle ids (see idt.dt): charged pions are the POI of arXiv:2307.14997
const int ID_pip = 120;    // pi+
const int ID_pim = -120;   // pi-
const float MPI_CHG = 0.13957;  // charged-pion mass (GeV); POI are exclusively pi+/-, so no need for the 'e' branch

// POI cuts, matching arXiv:2307.14997 (Xu et al.): |y|<1, 0.2<pT<2, true reaction plane
const float yCut  = 1.0;
const float ptMin = 0.2;
const float ptMax = 2.0;

// CME / flow correlators computed here (all vs the true reaction plane Psi_RP = "phi" branch):
//   v2      = <cos(2(phi - Psi_RP))>                       single-particle elliptic flow
//   gamma   = <cos(phi_a + phi_b - 2 Psi_RP)>   (Eq. 2)    charge-separation (CME) correlator
//   delta   = <cos(phi_a - phi_b)>              (App. B)   two-particle correlator
// split into OS (opposite-sign) and SS (same-sign) charge combinations, and the differences
//   Dgamma = gamma_OS - gamma_SS   (Eq. 3),  Ddelta = delta_OS - delta_SS
//
// gamma and delta are evaluated per event from charge-separated Q-vectors (O(N), no pair loop):
//   p_c  = Sum cos(phi),  q_c  = Sum sin(phi)      over pions of charge c
//   p2_c = Sum cos(2phi), q2_c = Sum sin(2phi)
//   NpairOS = N+ N- ,   NpairSS = N+(N+-1)/2 + N-(N--1)/2
//   Num(gamma_OS) = (p+ p- - q+ q-)cos2Psi + (q+ p- + p+ q-)sin2Psi
//   Num(gamma_SS) = 1/2 Sum_c [ (p_c^2 - q_c^2 - p2_c)cos2Psi + (2 p_c q_c - q2_c)sin2Psi ]
//   Num(delta_OS) = p+ p- + q+ q-
//   Num(delta_SS) = 1/2 Sum_c ( p_c^2 + q_c^2 - N_c )
// The statistical error is taken from the event-to-event spread (events are independent, pairs are not),
// so each single-bin TProfile's GetBinError() is the correct standard error of the mean over events.

void Correlators(int cen = 5, int job = 0){

    TChain* chain = new TChain("teposevent");
    chain->Add("z-auau_run_*.root");

    //only read the branches we use
    chain->SetBranchStatus("*",0);
    const char* activeBranches[] = {"np","bim","phi","px","py","pz","id","ist"};
    for(int ib = 0; ib < 8; ib++) chain->SetBranchStatus(activeBranches[ib],1);

    char fname_out[200];
    sprintf(fname_out,"cen%d.correlators_job%d.root",cen,job);
    TFile fout(fname_out,"RECREATE");

    //----- differential elliptic flow (per-particle fill; shape/QA) -----
    TProfile* pV2pT  = new TProfile("pV2pT" ,"v_{2}(p_{T}) #pi^{#pm} (RP);p_{T} (GeV/c);v_{2}",60,0,3);
    TProfile* pV2eta = new TProfile("pV2eta","v_{2}(#eta) #pi^{#pm} (RP);#eta;v_{2}",30,-3,3);

    //----- integrated observables, filled ONCE PER EVENT (error = event-to-event spread) -----
    TProfile* pV2int   = new TProfile("pV2int"  ,"integrated v_{2} #pi^{#pm}",1,0,1);
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

    //---- subgroup (sub-sampling) accumulators for the statistical error of derived ratios ----
    // kappa112 = Dgamma/(v2*Ddelta) is a nonlinear combination of three event-averaged
    // quantities that are all measured on the SAME events, so Dgamma and Ddelta are
    // correlated and a naive independent error propagation is wrong. The sub-sampling
    // method (used by STAR) splits the events into NSUB independent subgroups, evaluates
    // the observable in each, and takes the spread of the subgroup values -> this carries
    // the full covariance automatically. We accumulate the per-event v2/Dgamma/Ddelta into
    // NSUB subgroups here; kappa (and, as a cross-check, v2/Dgamma/Ddelta) errors are
    // formed from the subgroup spread after the loop.
    const int NSUB = 20;
    double sV2[NSUB]={0}, sDg[NSUB]={0}, sDd[NSUB]={0};
    long   cSub[NSUB]={0};

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

        //charge-separated Q-vectors over POI pions
        double pP=0, qP=0, p2P=0, q2P=0;  int Np=0;   // pi+
        double pM=0, qM=0, p2M=0, q2M=0;  int Nm=0;   // pi-

        for(int trk = 0; trk < NPTracks; trk++) {
            if((int)leaf_ist->GetValue(trk) != 0) continue;	//final-state only
            int pid = (int)leaf_id->GetValue(trk);
            if(!(pid == ID_pip || pid == ID_pim)) continue;

            float px = leaf_px->GetValue(trk);
            float py = leaf_py->GetValue(trk);
            float pz = leaf_pz->GetValue(trk);

            float pt = sqrt(px*px+py*py);
            float phi = atan2(py,px);
            float theta = atan2(pt,pz);
            float eta = -log(tan(theta/2.));
            float E = sqrt(MPI_CHG*MPI_CHG+px*px+py*py+pz*pz);
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
            if(pid == ID_pip){ pP+=c1; qP+=s1; p2P+=c2; q2P+=s2; Np++; }
            else             { pM+=c1; qM+=s1; p2M+=c2; q2M+=s2; Nm++; }
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

        //subgroup accumulators (round-robin so subgroups are equal-sized & unbiased in time)
        int sub = nUsed % NSUB;
        sV2[sub] += v2ev;
        sDg[sub] += (gOS - gSS);
        sDd[sub] += (dOS - dSS);
        cSub[sub]++;

        sumNp += Np; sumNm += Nm; nUsed++;
    }

    //---- form kappa112 and its statistical error from the subgroup spread ----
    // Full-sample central values (same as the single-bin profiles):
    double v2f = pV2int->GetBinContent(1);
    double Dgf = pDGamma->GetBinContent(1);
    double Ddf = pDDelta->GetBinContent(1);
    double kappa112 = (v2f*Ddf!=0) ? Dgf/(v2f*Ddf) : 0.0;

    // Per-subgroup estimates -> standard error of the mean over subgroups.
    // err = sqrt( Sum(x_m - xbar)^2 / (M*(M-1)) ), valid for equal-sized subgroups.
    // Done for kappa112 (needs the covariance) and, as a cross-check, for v2/Dgamma/Ddelta
    // (these should agree with the profile GetBinError values).
    double kappaErr=0, v2ErrSub=0, DgErrSub=0, DdErrSub=0;
    {
        int M=0; double meanK=0, meanV=0, meanG=0, meanD=0;
        double kk[NSUB], vv[NSUB], gg[NSUB], dd[NSUB];
        for(int m=0;m<NSUB;m++){
            if(cSub[m]<=0) continue;
            double v=sV2[m]/cSub[m], g=sDg[m]/cSub[m], d=sDd[m]/cSub[m];
            double k=(v*d!=0)? g/(v*d) : 0.0;
            kk[M]=k; vv[M]=v; gg[M]=g; dd[M]=d;
            meanK+=k; meanV+=v; meanG+=g; meanD+=d; M++;
        }
        if(M>=2){
            meanK/=M; meanV/=M; meanG/=M; meanD/=M;
            double sK=0,sV=0,sG=0,sD=0;
            for(int m=0;m<M;m++){
                sK+=(kk[m]-meanK)*(kk[m]-meanK);
                sV+=(vv[m]-meanV)*(vv[m]-meanV);
                sG+=(gg[m]-meanG)*(gg[m]-meanG);
                sD+=(dd[m]-meanD)*(dd[m]-meanD);
            }
            double denom = (double)M*(M-1);
            kappaErr = sqrt(sK/denom);
            v2ErrSub = sqrt(sV/denom);
            DgErrSub = sqrt(sG/denom);
            DdErrSub = sqrt(sD/denom);
        }
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

    //kappa112 = Dgamma / (v2 * Ddelta); value + sub-sampling statistical error
    TH1D* hKappa = new TH1D("hKappa","#kappa_{112}=#Delta#gamma/(v_{2}#Delta#delta); ; value",1,0.5,1.5);
    hKappa->SetBinContent(1, kappa112);
    hKappa->SetBinError(1, kappaErr);

    cout << "\n==== correlators (30-40% Au+Au 200 GeV, charged pions, true RP) ====\n";
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

    cout << "\n---- derived: kappa112 = Dgamma / (v2 * Ddelta) ----\n";
    printf("  v2      = %+.4e +/- %.1e\n", v2f, pV2int ->GetBinError(1));
    printf("  Dgamma  = %+.4e +/- %.1e\n", Dgf, pDGamma->GetBinError(1));
    printf("  Ddelta  = %+.4e +/- %.1e\n", Ddf, pDDelta->GetBinError(1));
    printf("  kappa112= %+.4e +/- %.1e   (sub-sampling, NSUB=%d)\n", kappa112, kappaErr, NSUB);
    printf("  [subgroup cross-check errors  v2:%.1e Dgamma:%.1e Ddelta:%.1e -- should match the profile errors above]\n",
           v2ErrSub, DgErrSub, DdErrSub);

    fout.Write();
    return;
}
