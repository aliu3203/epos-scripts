using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <cmath>
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include "EposPID.h"
#include "STAR200.h"

// -----------------------------------------------------------------------------
// Overlay the EPOS correlators on the STAR Au+Au 200 GeV measurements
// (HEPData ins2928164) as a function of centrality.
//
// Works in both modes of Correlators_Cent.C:
//   * centrality-differential (9 bins, needs a minimum-bias sample) -> full curves
//   * inclusive (1 bin, b-restricted sample)                        -> single point,
//     drawn at `inclusiveCent` (default 35%, i.e. the 8.4-9.2 fm slice of auau200-7)
// -----------------------------------------------------------------------------

static TGraphErrors* mkg(int n, const double* x, const double* y, const double* e,
                         int col, int mstyle, double msize=1.3){
    TGraphErrors* g = new TGraphErrors(n);
    for(int i=0;i<n;i++){ g->SetPoint(i,x[i],y[i]); g->SetPointError(i,0,e?e[i]:0); }
    g->SetMarkerStyle(mstyle); g->SetMarkerSize(msize);
    g->SetMarkerColor(col); g->SetLineColor(col); g->SetLineWidth(2);
    return g;
}

// STAR total error = stat (+) sys in quadrature
static void startot(const double* st, const double* sy, double* out){
    for(int i=0;i<NSTAR;i++) out[i]=sqrt(st[i]*st[i]+sy[i]*sy[i]);
}

void MakeFigure_Cent_STAR(int poiMode = 0, int job = 0, double inclusiveCent = 35.0){

    char fname[200];
    sprintf(fname,"cen.correlators_cent_poi%d_job%d.root",poiMode,job);
    TFile* f = TFile::Open(fname);
    if(!f || f->IsZombie()){ printf("cannot open %s\n",fname); return; }

    TH1D* hNev  = (TH1D*)f->Get("hNev");
    TH1D* hK112 = (TH1D*)f->Get("hKappa112");
    TH1D* hK132 = (TH1D*)f->Get("hKappa132");
    if(!hNev){ printf("hNev missing\n"); return; }
    int NB = hNev->GetNbinsX();
    bool inclusive = (NB==1);
    printf("EPOS file %s : %d centrality bin(s)%s\n", fname, NB, inclusive?" (inclusive mode)":"");

    // ---- pull EPOS values into arrays ----
    double ex[NCENT], ev2[NCENT],ev2e[NCENT], eg1OS[NCENT],eg1OSe[NCENT], eg1SS[NCENT],eg1SSe[NCENT];
    double eg3OS[NCENT],eg3OSe[NCENT], eg3SS[NCENT],eg3SSe[NCENT];
    double edOS[NCENT],edOSe[NCENT], edSS[NCENT],edSSe[NCENT];
    double eDg1[NCENT],eDg1e[NCENT], eDg3[NCENT],eDg3e[NCENT], eDd[NCENT],eDde[NCENT];
    double ek1[NCENT],ek1e[NCENT], ek3[NCENT],ek3e[NCENT];

    auto grab=[&](const char* pre,int ic,double& v,double& e){
        TProfile* p=(TProfile*)f->Get(Form("%s_c%d",pre,ic));
        v = p? p->GetBinContent(1):0; e = p? p->GetBinError(1):0;
    };
    for(int ic=0; ic<NB; ic++){
        ex[ic] = inclusive ? inclusiveCent : CENT_MID[ic];
        grab("pV2",ic,ev2[ic],ev2e[ic]);
        grab("pG112OS",ic,eg1OS[ic],eg1OSe[ic]);  grab("pG112SS",ic,eg1SS[ic],eg1SSe[ic]);
        grab("pG132OS",ic,eg3OS[ic],eg3OSe[ic]);  grab("pG132SS",ic,eg3SS[ic],eg3SSe[ic]);
        grab("pDOS",ic,edOS[ic],edOSe[ic]);       grab("pDSS",ic,edSS[ic],edSSe[ic]);
        // Correlators_Cent.C no longer stores per-event Delta profiles: every term is
        // filled individually into the OS/SS profiles, so Delta is formed here as
        // OS - SS with the two profile errors added in quadrature.
        eDg1[ic] = eg1OS[ic]-eg1SS[ic];
        eDg1e[ic]= sqrt(eg1OSe[ic]*eg1OSe[ic] + eg1SSe[ic]*eg1SSe[ic]);
        eDg3[ic] = eg3OS[ic]-eg3SS[ic];
        eDg3e[ic]= sqrt(eg3OSe[ic]*eg3OSe[ic] + eg3SSe[ic]*eg3SSe[ic]);
        eDd[ic]  = edOS[ic]-edSS[ic];
        eDde[ic] = sqrt(edOSe[ic]*edOSe[ic] + edSSe[ic]*edSSe[ic]);
        ek1[ic]=hK112->GetBinContent(ic+1); ek1e[ic]=hK112->GetBinError(ic+1);
        ek3[ic]=hK132->GetBinContent(ic+1); ek3e[ic]=hK132->GetBinError(ic+1);
    }

    // ---- STAR derived quantities ----
    double sDg1[NSTAR],sDg1e[NSTAR], sDg3[NSTAR],sDg3e[NSTAR], sDd[NSTAR],sDde[NSTAR];
    double t_g1OS[NSTAR],t_g1SS[NSTAR],t_g3OS[NSTAR],t_g3SS[NSTAR],t_dOS[NSTAR],t_dSS[NSTAR];
    double t_v2[NSTAR],t_k1[NSTAR],t_k3[NSTAR];
    startot(STAR_g112OS_stat,STAR_g112OS_sys,t_g1OS);
    startot(STAR_g112SS_stat,STAR_g112SS_sys,t_g1SS);
    startot(STAR_g132OS_stat,STAR_g132OS_sys,t_g3OS);
    startot(STAR_g132SS_stat,STAR_g132SS_sys,t_g3SS);
    startot(STAR_dOS_stat,STAR_dOS_sys,t_dOS);
    startot(STAR_dSS_stat,STAR_dSS_sys,t_dSS);
    startot(STAR_v2_stat,STAR_v2_sys,t_v2);
    startot(STAR_k112_stat,STAR_k112_sys,t_k1);
    startot(STAR_k132_stat,STAR_k132_sys,t_k3);
    for(int i=0;i<NSTAR;i++){
        sDg1[i]=STAR_g112OS[i]-STAR_g112SS[i]; sDg1e[i]=sqrt(t_g1OS[i]*t_g1OS[i]+t_g1SS[i]*t_g1SS[i]);
        sDg3[i]=STAR_g132OS[i]-STAR_g132SS[i]; sDg3e[i]=sqrt(t_g3OS[i]*t_g3OS[i]+t_g3SS[i]*t_g3SS[i]);
        sDd[i] =STAR_dOS[i]-STAR_dSS[i];       sDde[i] =sqrt(t_dOS[i]*t_dOS[i]+t_dSS[i]*t_dSS[i]);
    }

    // ---- numeric comparison table at the EPOS centrality/centralities ----
    printf("\n============ EPOS vs STAR (Au+Au 200 GeV) ============\n");
    printf("POI = %s\n", poiMode==0? "charged pions (pi+/-)"
                       : poiMode==2? "charged mesons (pi+/-, K+/-)  [matches STAR 200 GeV]"
                                   : "all charged hadrons incl. p/pbar");
    auto near=[&](double c){ int b=0; double d=1e9;
        for(int i=0;i<NSTAR;i++){ double dd=fabs(STAR_cent[i]-c); if(dd<d){d=dd;b=i;} } return b; };
    for(int ic=0; ic<NB; ic++){
        int s = near(ex[ic]);
        printf("\n--- centrality %.0f%% (STAR bin %.0f%%) ---\n", ex[ic], STAR_cent[s]);
        printf("  %-11s %14s %10s   %14s %10s   %7s\n","obs","EPOS","+/-","STAR","+/-","n-sigma");
        auto row=[&](const char* n,double e,double ee,double s_,double se){
            double d=sqrt(ee*ee+se*se);
            printf("  %-11s %+14.5e %10.3e   %+14.5e %10.3e   %7.1f\n",n,e,ee,s_,se, d>0?fabs(e-s_)/d:0.0);
        };
        row("v2",       ev2[ic],ev2e[ic],   STAR_v2[s],      t_v2[s]);
        row("g112_OS",  eg1OS[ic],eg1OSe[ic],STAR_g112OS[s], t_g1OS[s]);
        row("g112_SS",  eg1SS[ic],eg1SSe[ic],STAR_g112SS[s], t_g1SS[s]);
        row("Dgamma112",eDg1[ic],eDg1e[ic], sDg1[s],         sDg1e[s]);
        row("g132_OS",  eg3OS[ic],eg3OSe[ic],STAR_g132OS[s], t_g3OS[s]);
        row("g132_SS",  eg3SS[ic],eg3SSe[ic],STAR_g132SS[s], t_g3SS[s]);
        row("Dgamma132",eDg3[ic],eDg3e[ic], sDg3[s],         sDg3e[s]);
        row("delta_OS", edOS[ic],edOSe[ic], STAR_dOS[s],     t_dOS[s]);
        row("delta_SS", edSS[ic],edSSe[ic], STAR_dSS[s],     t_dSS[s]);
        row("Ddelta",   eDd[ic],eDde[ic],   sDd[s],          sDde[s]);
        row("kappa112", ek1[ic],ek1e[ic],   STAR_k112[s],    t_k1[s]);
        row("kappa132", ek3[ic],ek3e[ic],   STAR_k132[s],    t_k3[s]);
    }

    // ------------------------------- figure -------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetTitleSize(0.06,"t");
    TCanvas* c = new TCanvas("cCent","EPOS vs STAR vs centrality",1500,1200);
    c->Divide(3,3);

    int ipad=1;
    auto panel=[&](const char* title,const char* ytit,
                   const double* sy1,const double* se1,const double* sy2,const double* se2,
                   const double* ey1,const double* ee1,const double* ey2,const double* ee2,
                   const char* l1,const char* l2){
        c->cd(ipad++); gPad->SetGridy(); gPad->SetLeftMargin(0.16);
        // build all graphs first to find the y-range
        double lo=1e30,hi=-1e30;
        auto upd=[&](const double* y,const double* e,int n){ if(!y)return;
            for(int i=0;i<n;i++){ double a=y[i]-(e?e[i]:0), b=y[i]+(e?e[i]:0); if(a<lo)lo=a; if(b>hi)hi=b; } };
        upd(sy1,se1,NSTAR); upd(sy2,se2,NSTAR); upd(ey1,ee1,NB); upd(ey2,ee2,NB);
        double pad=0.18*(hi-lo); if(pad<=0) pad=0.1*fabs(hi)+1e-12;

        TGraphErrors* gs1 = mkg(NSTAR,STAR_cent,sy1,se1,kBlack,20);
        gs1->SetTitle(Form("%s;centrality [%%];%s",title,ytit));
        gs1->GetXaxis()->SetLimits(0,80);
        gs1->GetYaxis()->SetRangeUser(lo-pad,hi+pad);
        gs1->GetYaxis()->SetTitleOffset(1.6);
        gs1->Draw("APL");
        TGraphErrors* gs2=0; if(sy2){ gs2=mkg(NSTAR,STAR_cent,sy2,se2,kGray+2,24); gs2->Draw("PL same"); }
        TGraphErrors* ge1 = mkg(NB,ex,ey1,ee1,kRed+1,29,2.0); ge1->Draw("PL same");
        TGraphErrors* ge2=0; if(ey2){ ge2=mkg(NB,ex,ey2,ee2,kAzure+2,30,2.0); ge2->Draw("PL same"); }

        TLegend* lg=new TLegend(0.18,0.68,0.62,0.89);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.043);
        lg->AddEntry(gs1,Form("STAR %s",l1),"lp");
        if(gs2) lg->AddEntry(gs2,Form("STAR %s",l2),"lp");
        lg->AddEntry(ge1,Form("EPOS %s",l1),"lp");
        if(ge2) lg->AddEntry(ge2,Form("EPOS %s",l2),"lp");
        lg->Draw();
    };

    panel("v_{2}","v_{2}",                       STAR_v2,t_v2,0,0,          ev2,ev2e,0,0,            "","");
    panel("#gamma^{112}","#gamma^{112}",         STAR_g112OS,t_g1OS,STAR_g112SS,t_g1SS, eg1OS,eg1OSe,eg1SS,eg1SSe,"OS","SS");
    panel("#gamma^{132}","#gamma^{132}",         STAR_g132OS,t_g3OS,STAR_g132SS,t_g3SS, eg3OS,eg3OSe,eg3SS,eg3SSe,"OS","SS");
    panel("#delta","#delta",                     STAR_dOS,t_dOS,STAR_dSS,t_dSS,         edOS,edOSe,edSS,edSSe,"OS","SS");
    panel("#Delta#gamma^{112}","#Delta#gamma^{112}", sDg1,sDg1e,0,0,        eDg1,eDg1e,0,0,          "","");
    panel("#Delta#gamma^{132}","#Delta#gamma^{132}", sDg3,sDg3e,0,0,        eDg3,eDg3e,0,0,          "","");
    panel("#Delta#delta","#Delta#delta",         sDd,sDde,0,0,              eDd,eDde,0,0,            "","");
    panel("#kappa^{112}","#kappa^{112}",         STAR_k112,t_k1,0,0,        ek1,ek1e,0,0,            "","");
    panel("#kappa^{132}","#kappa^{132}",         STAR_k132,t_k3,0,0,        ek3,ek3e,0,0,            "","");

    char png[200]; sprintf(png,"correlators_vs_centrality_poi%d.png",poiMode);
    c->SaveAs(png);
    printf("\nwrote %s\n",png);
    if(inclusive)
        printf("NOTE: EPOS shown as a single point at %.0f%% -- auau200-7 is a b-restricted\n"
               "      (30-40%%) sample. Generate minimum-bias events for the full curve.\n", inclusiveCent);
}
