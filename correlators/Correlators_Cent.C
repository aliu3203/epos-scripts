using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TMath.h"
#include "TProfile.h"
#include "EposPID.h"

// -----------------------------------------------------------------------------
// Centrality-differential CME/flow correlators vs the true reaction plane, for
// direct comparison with STAR (HEPData ins2928164, Au+Au 200 GeV):
//
//   v2        = <cos2(phi - Psi_RP)>                        (Fig. 5)
//   gamma112  = <cos(phi_a + phi_b - 2 Psi_RP)>             (Fig. 6)
//   gamma132  = <cos(phi_a - 3 phi_b + 2 Psi_RP)>           (Fig. 7)
//   delta     = <cos(phi_a - phi_b)>                        (Fig. 8)
//   kappa112  = Dgamma112 / (v2 * Ddelta)                   (Fig. 13)
//   kappa132  = Dgamma132 / (v2 * Ddelta)                   (Fig. 13)
// (the kappa definitions were verified against STAR's own published numbers:
//  their gamma/delta/v2 values reproduce their quoted kappa to 0.2%)
//
// HOW THE CORRELATORS ARE BUILT
// -----------------------------
// Every term is computed in an explicit loop over the N accepted tracks of the event
// and is Filled into a TProfile individually -- one Fill per term, never a per-event
// average.  The profile therefore holds the full distribution of terms, and
// GetBinError(1) is the standard deviation of that distribution divided by
// sqrt(number of terms).  That is the quoted error.  Nothing else is modelled.
//
//   gamma112, delta:  for(a = 0; a < N; a++) for(b = a+1; b < N; b++)
//                     -> N(N-1)/2 terms.  The pair term is symmetric in (a,b),
//                        so each unordered pair is visited once.
//
//   gamma132:         for(a = 0; a < N; a++) for(b = 0; b < N; b++) if(a==b) continue;
//                     -> N(N-1) terms.  cos(phi_a - 3 phi_b + 2Psi) is NOT symmetric
//                        in (a,b), so (a,b) and (b,a) are two distinct terms and both
//                        are kept.
//
// gamma132 therefore has twice the terms of gamma112 and its error comes out smaller
// by ~sqrt(2).  That is a direct consequence of the definitions, not a bug.
//
//   v2:               one Fill per accepted track, cos2(phi - Psi_RP).
//
// Each term goes into the OS or the SS profile according to the charge product of
// (a,b), because the physics quantities are the differences
//   Dgamma = gamma_OS - gamma_SS,  Ddelta = delta_OS - delta_SS,
// whose errors are the quadrature sum of the two profile errors.
//
// Centrality comes from refmult percentile cuts in cen_cuts.txt (see
// CentralityCalib.C). If that file is absent the macro runs in INCLUSIVE mode:
// one bin holding every event that passes, which is the correct behaviour for a
// b-restricted (single-centrality) sample like auau200-7.
// -----------------------------------------------------------------------------

const float MPI_CHG = 0.13957;
const float yCut = 1.0, ptMin = 0.2, ptMax = 2.0;

// poiMode 0 = charged pions (the earlier auau200-7 results)
//         1 = all charged hadrons, including p/pbar (STAR's Fig.16 definition,
//             used only at 11.5-19.6 GeV -- NOT the 200 GeV figures)
//         2 = charged mesons, pi+/- and K+/- (STAR's POI for the 200 GeV
//             Figs. 5,6,7,8,13 -- "v2 of charged mesons"; p/pbar excluded).
//             THIS is the mode to use when comparing with the 200 GeV data.
void Correlators_Cent(int job = 0, int poiMode = 0, const char* filepat = "z-auau_run_*.root"){

    if(!EposPID::Load("idt.dt")) return;

    //------------------- centrality cuts (optional) -------------------
    int cut[NCENT]; bool useCent = false;
    {
        ifstream in("cen_cuts.txt");
        if(in){
            int n=0; string line;
            while(getline(in,line)){
                if(line.empty() || line[0]=='#') continue;
                if(n<NCENT) cut[n++] = atoi(line.c_str());
            }
            if(n==NCENT){ useCent=true;
                printf("centrality cuts from cen_cuts.txt: ");
                for(int i=0;i<NCENT;i++) printf("%d ",cut[i]);
                printf("\n");
            }
        }
    }
    const int NB = useCent ? NCENT : 1;
    if(!useCent)
        printf("\n*** cen_cuts.txt not found -> INCLUSIVE mode (single bin, all events).\n"
               "*** This is correct for a b-restricted sample; run CentralityCalib.C on a\n"
               "*** minimum-bias sample first to enable centrality-differential output.\n\n");

    //------------------- input -------------------
    TChain* chain = new TChain("teposevent");
    chain->Add(filepat);
    chain->SetBranchStatus("*",0);
    const char* act[] = {"np","bim","phi","px","py","pz","id","ist","npartproj","nparttarg"};
    for(int i=0;i<10;i++) chain->SetBranchStatus(act[i],1);

    char fname_out[200];
    sprintf(fname_out,"cen.correlators_cent_poi%d_job%d.root",poiMode,job);
    TFile fout(fname_out,"RECREATE");

    //------------------- output profiles, one set per centrality bin -------------------
    // Every correlator profile is filled once per TERM.  The event-level metadata
    // profiles (Npart, refmult, b) are filled once per event.
    TProfile *pV2[NCENT];
    TProfile *pG112OS[NCENT], *pG112SS[NCENT];
    TProfile *pG132OS[NCENT], *pG132SS[NCENT];
    TProfile *pDOS[NCENT],    *pDSS[NCENT];
    TProfile *pNpart[NCENT], *pRefmult[NCENT], *pBim[NCENT];
    for(int ic=0; ic<NB; ic++){
        pV2[ic]    = new TProfile(Form("pV2_c%d",ic),    Form("v_{2} %s",CentLabel(ic)),1,0,1);
        pG112OS[ic]= new TProfile(Form("pG112OS_c%d",ic),Form("#gamma^{112}_{OS} %s",CentLabel(ic)),1,0,1);
        pG112SS[ic]= new TProfile(Form("pG112SS_c%d",ic),Form("#gamma^{112}_{SS} %s",CentLabel(ic)),1,0,1);
        pG132OS[ic]= new TProfile(Form("pG132OS_c%d",ic),Form("#gamma^{132}_{OS} %s",CentLabel(ic)),1,0,1);
        pG132SS[ic]= new TProfile(Form("pG132SS_c%d",ic),Form("#gamma^{132}_{SS} %s",CentLabel(ic)),1,0,1);
        pDOS[ic]   = new TProfile(Form("pDOS_c%d",ic),   Form("#delta_{OS} %s",CentLabel(ic)),1,0,1);
        pDSS[ic]   = new TProfile(Form("pDSS_c%d",ic),   Form("#delta_{SS} %s",CentLabel(ic)),1,0,1);
        pNpart[ic] = new TProfile(Form("pNpart_c%d",ic), Form("N_{part} %s",CentLabel(ic)),1,0,1);
        pRefmult[ic]=new TProfile(Form("pRefmult_c%d",ic),Form("refmult %s",CentLabel(ic)),1,0,1);
        pBim[ic]   = new TProfile(Form("pBim_c%d",ic),   Form("b %s",CentLabel(ic)),1,0,1);
    }

    // per-track trig, rebuilt each event: cos/sin of phi and of 3phi, plus the charge
    vector<double> c1, s1, c3, s3;
    vector<int>    qq;

    Long64_t nentries = chain->GetEntries();
    cout << nentries << " entries in chain\n";
    Long64_t nUsed = 0;

    for(Long64_t i=0;i<nentries;i++){
        if((i+1)%10000==0) cout<<"Processing entry == "<<i+1<<" == out of "<<nentries<<".\n";
        chain->GetEntry(i);

        int   NPTracks = (int)chain->GetLeaf("np")->GetValue(0);
        float PsiRP    = chain->GetLeaf("phi")->GetValue(0);
        float bim      = chain->GetLeaf("bim")->GetValue(0);
        int   Npart    = (int)chain->GetLeaf("npartproj")->GetValue(0)
                       + (int)chain->GetLeaf("nparttarg")->GetValue(0);
        double C2 = cos(2.*PsiRP), S2 = sin(2.*PsiRP);

        TLeaf* lpx=chain->GetLeaf("px"); TLeaf* lpy=chain->GetLeaf("py");
        TLeaf* lpz=chain->GetLeaf("pz"); TLeaf* lid=chain->GetLeaf("id");
        TLeaf* lst=chain->GetLeaf("ist");

        c1.clear(); s1.clear(); c3.clear(); s3.clear(); qq.clear();
        int refmult = 0, Np = 0, Nm = 0;

        //---- track loop: select the POI and cache their trig ----
        for(int trk=0; trk<NPTracks; trk++){
            if((int)lst->GetValue(trk) != 0) continue;              // final state
            int pid = (int)lid->GetValue(trk);
            int chg = EposPID::Charge(pid);
            if(chg == 0) continue;                                   // neutral: no use here

            float px=lpx->GetValue(trk), py=lpy->GetValue(trk), pz=lpz->GetValue(trk);
            float pt2 = px*px+py*py;

            // ---- refmult (centrality estimator) ----
            if(pt2 >= REFMULT_PTMIN*REFMULT_PTMIN && pz*pz < SINH_ETACUT_SQ*pt2) refmult++;

            // ---- POI selection ----
            int apid = abs(pid);
            bool isPOI;
            if      (poiMode==0) isPOI = (apid==120);                 // pi+/- only
            else if (poiMode==2) isPOI = (apid==120 || apid==130);    // charged mesons
            else                 isPOI = true;                        // all charged hadrons
            if(!isPOI) continue;
            float pt = sqrt(pt2);
            if(pt < ptMin || pt > ptMax) continue;
            float m  = (poiMode==0) ? MPI_CHG : EposPID::Mass(pid);
            float E  = sqrt(m*m + pt2 + pz*pz);
            float y  = 0.5*log((E+pz)/(E-pz));
            if(fabs(y) > yCut) continue;

            float phi = atan2(py,px);
            c1.push_back(cos(phi));    s1.push_back(sin(phi));
            c3.push_back(cos(3.*phi)); s3.push_back(sin(3.*phi));
            if(chg > 0){ qq.push_back(+1); Np++; }
            else       { qq.push_back(-1); Nm++; }
        }

        if(Np < 2 || Nm < 2) continue;

        int ic = useCent ? CentBinFromRefmult(refmult, cut) : 0;
        if(ic < 0) continue;                                        // more peripheral than 80%

        const int N = (int)c1.size();
        TProfile* g112[2] = {pG112OS[ic], pG112SS[ic]};
        TProfile* g132[2] = {pG132OS[ic], pG132SS[ic]};
        TProfile* del [2] = {pDOS[ic],    pDSS[ic]};

        //---- v2: one term per track ----
        // cos(2 phi - 2 Psi) = cos2phi cos2Psi + sin2phi sin2Psi
        for(int a=0;a<N;a++){
            double c2a = c1[a]*c1[a] - s1[a]*s1[a];      // cos 2phi
            double s2a = 2.*c1[a]*s1[a];                 // sin 2phi
            pV2[ic]->Fill(0.5, c2a*C2 + s2a*S2);
        }

        //---- gamma112 and delta: b = a+1 .. N-1  ->  N(N-1)/2 terms ----
        for(int a=0;a<N;a++){
            const double ca=c1[a], sa=s1[a];
            const int    qa=qq[a];
            for(int b=a+1;b<N;b++){
                const int k = (qa*qq[b] < 0) ? 0 : 1;    // 0 = opposite sign, 1 = same sign
                const double cb=c1[b], sb=s1[b];

                // cos(phi_a + phi_b - 2Psi)
                double Rp = ca*cb - sa*sb, Ip = sa*cb + ca*sb;
                g112[k]->Fill(0.5, Rp*C2 + Ip*S2);

                // cos(phi_a - phi_b)
                del[k]->Fill(0.5, ca*cb + sa*sb);
            }
        }

        //---- gamma132: b = 0 .. N-1, skipping b == a  ->  N(N-1) terms ----
        for(int a=0;a<N;a++){
            const double ca=c1[a], sa=s1[a];
            const int    qa=qq[a];
            for(int b=0;b<N;b++){
                if(a==b) continue;
                const int k = (qa*qq[b] < 0) ? 0 : 1;

                // cos(phi_a - 3 phi_b + 2Psi) = Re[ e^{i phi_a} conj(e^{3i phi_b}) e^{2iPsi} ]
                double X = ca*c3[b] + sa*s3[b];
                double Y = sa*c3[b] - ca*s3[b];
                g132[k]->Fill(0.5, X*C2 - Y*S2);
            }
        }

        pNpart[ic] ->Fill(0.5, Npart);
        pRefmult[ic]->Fill(0.5, refmult);
        pBim[ic]   ->Fill(0.5, bim);
        nUsed++;
    }

    //------------------- results -------------------
    TH1D* hK112 = new TH1D("hKappa112","#kappa^{112} vs centrality;centrality [%];#kappa^{112}",NB,0,NB);
    TH1D* hK132 = new TH1D("hKappa132","#kappa^{132} vs centrality;centrality [%];#kappa^{132}",NB,0,NB);
    TH1D* hNev  = new TH1D("hNev","events per centrality bin",NB,0,NB);

    printf("\n================ EPOS correlators%s ================\n",
           useCent? " vs centrality" : " (inclusive)");
    const char* poiName = (poiMode==0)? "charged pions (pi+/-)"
                        : (poiMode==2)? "charged mesons (pi+/-, K+/-)  [STAR 200 GeV POI]"
                                      : "all charged hadrons incl. p/pbar";
    printf("POI = %s,  |y|<1, %.1f<pT<%.1f,  true reaction plane\n", poiName, ptMin, ptMax);
    printf("gamma112, delta : b = a+1..N-1   (N(N-1)/2 terms)\n");
    printf("gamma132        : b = 0..N-1, b != a   (N(N-1) terms)\n");
    printf("errors are TProfile GetBinError = std dev of the terms / sqrt(N terms)\n");
    printf("events used = %lld\n\n", nUsed);

    for(int ic=0; ic<NB; ic++){
        double nev = pBim[ic]->GetBinEntries(1);
        if(nev <= 0) continue;

        double v2  = pV2[ic]->GetBinContent(1),  ev2 = pV2[ic]->GetBinError(1);
        double a1o = pG112OS[ic]->GetBinContent(1), e1o = pG112OS[ic]->GetBinError(1);
        double a1s = pG112SS[ic]->GetBinContent(1), e1s = pG112SS[ic]->GetBinError(1);
        double a3o = pG132OS[ic]->GetBinContent(1), e3o = pG132OS[ic]->GetBinError(1);
        double a3s = pG132SS[ic]->GetBinContent(1), e3s = pG132SS[ic]->GetBinError(1);
        double ddo = pDOS[ic]->GetBinContent(1),    edo = pDOS[ic]->GetBinError(1);
        double dds = pDSS[ic]->GetBinContent(1),    eds = pDSS[ic]->GetBinError(1);

        double Dg112 = a1o - a1s, eDg112 = sqrt(e1o*e1o + e1s*e1s);
        double Dg132 = a3o - a3s, eDg132 = sqrt(e3o*e3o + e3s*e3s);
        double Ddel  = ddo - dds, eDdel  = sqrt(edo*edo + eds*eds);

        double den   = v2*Ddel;
        double k112  = den!=0 ? Dg112/den : 0.0;
        double k132  = den!=0 ? Dg132/den : 0.0;
        // relative errors added in quadrature
        double rel   = (v2!=0 && Ddel!=0) ? sqrt(ev2*ev2/(v2*v2) + eDdel*eDdel/(Ddel*Ddel)) : 0.0;
        double ek112 = (Dg112!=0) ? fabs(k112)*sqrt(eDg112*eDg112/(Dg112*Dg112) + rel*rel) : 0.0;
        double ek132 = (Dg132!=0) ? fabs(k132)*sqrt(eDg132*eDg132/(Dg132*Dg132) + rel*rel) : 0.0;

        printf("---- %s  (Nev=%.0f, <refmult>=%.1f, <Npart>=%.1f, <b>=%.2f fm) ----\n",
               useCent? CentLabel(ic):"inclusive", nev,
               pRefmult[ic]->GetBinContent(1), pNpart[ic]->GetBinContent(1),
               pBim[ic]->GetBinContent(1));
        printf("   %-12s %14s %12s %16s %12s\n","quantity","value","error","N terms","terms/event");
        struct { const char* nm; TProfile* p; } row[7] = {
            {"v2",       pV2[ic]},
            {"g112_OS",  pG112OS[ic]}, {"g112_SS", pG112SS[ic]},
            {"g132_OS",  pG132OS[ic]}, {"g132_SS", pG132SS[ic]},
            {"delta_OS", pDOS[ic]},    {"delta_SS",pDSS[ic]} };
        for(int r=0;r<7;r++){
            double nt = row[r].p->GetBinEntries(1);
            printf("   %-12s %+14.6e %12.4e %16.0f %12.1f\n",
                   row[r].nm, row[r].p->GetBinContent(1), row[r].p->GetBinError(1), nt, nt/nev);
        }
        printf("   %-12s %+14.6e %12.4e\n","Dgamma112", Dg112, eDg112);
        printf("   %-12s %+14.6e %12.4e\n","Dgamma132", Dg132, eDg132);
        printf("   %-12s %+14.6e %12.4e\n","Ddelta",    Ddel,  eDdel);
        printf("   %-12s %+14.6f %12.4f\n","kappa112",  k112,  ek112);
        printf("   %-12s %+14.6f %12.4f\n","kappa132",  k132,  ek132);
        printf("   term ratio g132/g112 = %.4f (OS), %.4f (SS);  err ratio Dg112/Dg132 = %.4f\n\n",
               pG132OS[ic]->GetBinEntries(1)/pG112OS[ic]->GetBinEntries(1),
               pG132SS[ic]->GetBinEntries(1)/pG112SS[ic]->GetBinEntries(1),
               eDg132>0 ? eDg112/eDg132 : 0.0);

        hK112->SetBinContent(ic+1,k112); hK112->SetBinError(ic+1,ek112);
        hK132->SetBinContent(ic+1,k132); hK132->SetBinError(ic+1,ek132);
        hNev ->SetBinContent(ic+1,nev);
        if(useCent){
            hK112->GetXaxis()->SetBinLabel(ic+1,CentLabel(ic));
            hK132->GetXaxis()->SetBinLabel(ic+1,CentLabel(ic));
        }
    }

    fout.Write();
    printf("wrote %s\n", fname_out);
}
