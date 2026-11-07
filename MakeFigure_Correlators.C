using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraphErrors.h"

// Read cen5.correlators_job0.root, report v2 / gamma / delta with statistical
// uncertainties, and project how many events are needed to reach target precision.
// Statistical error scales as 1/sqrt(N_events), so for a current error sigma from
// Nev events, reaching target sigma_t needs Nev*(sigma/sigma_t)^2 events.

static void proj(const char* name, double val, double err, double Nev){
    printf("\n  %s = %+.4e +/- %.4e   (current significance %.2f sigma, Nev=%.0f)\n",
           name, val, err, err>0? fabs(val)/err : 0.0, Nev);
    printf("    to measure the observed central value at:\n");
    double sig[] = {3.0, 5.0};
    for(double ns : sig){
        double st = fabs(val)/ns;                 // target error for ns-sigma
        double Nn = Nev*(err/st)*(err/st);        // events needed
        printf("       %.0f sigma : sigma_t=%.2e  -> %.3e events  (%.0f files, x%.1f more)\n",
               ns, st, Nn, Nn/1000., Nn/Nev);
    }
    printf("    to reach a target absolute uncertainty:\n");
    double tgt[] = {1e-4, 5e-5, 2e-5, 1e-5};
    for(double st : tgt){
        double Nn = Nev*(err/st)*(err/st);
        printf("       sigma_t=%.0e -> %.3e events  (%.0f files, x%.1f more)\n",
               st, Nn, Nn/1000., Nn/Nev);
    }
}

void MakeFigure_Correlators(int cen = 5, int job = 0){

    char fname[200];
    sprintf(fname,"cen%d.correlators_job%d.root",cen,job);
    TFile* f = TFile::Open(fname);
    if(!f || f->IsZombie()){ printf("cannot open %s\n",fname); return; }

    TProfile* pV2pT   = (TProfile*)f->Get("pV2pT");
    TProfile* pV2eta  = (TProfile*)f->Get("pV2eta");
    TProfile* pV2int  = (TProfile*)f->Get("pV2int");
    TProfile* pGammaOS= (TProfile*)f->Get("pGammaOS");
    TProfile* pGammaSS= (TProfile*)f->Get("pGammaSS");
    TProfile* pDGamma = (TProfile*)f->Get("pDGamma");
    TProfile* pDeltaOS= (TProfile*)f->Get("pDeltaOS");
    TProfile* pDeltaSS= (TProfile*)f->Get("pDeltaSS");
    TProfile* pDDelta = (TProfile*)f->Get("pDDelta");
    TH1D*     hMeta   = (TH1D*)    f->Get("hMeta");
    TH1D*     hKappa  = (TH1D*)    f->Get("hKappa");

    double Nev = hMeta->GetBinContent(1);

    printf("\n=================  30-40%% Au+Au 200 GeV, charged pions, true RP  =================\n");
    printf("events used = %.0f   ( <N+>=%.1f  <N->=%.1f per event )\n",
           Nev, hMeta->GetBinContent(2), hMeta->GetBinContent(3));
    printf("\n  %-10s %13s %13s %10s\n","observable","value","stat.err","signif.");
    auto row = [&](const char* n, TProfile* p){
        double v=p->GetBinContent(1), e=p->GetBinError(1);
        printf("  %-10s %+13.4e %13.4e %9.2f\n", n, v, e, e>0?fabs(v)/e:0.0);
    };
    row("v2",       pV2int);
    row("gamma_OS", pGammaOS);
    row("gamma_SS", pGammaSS);
    row("Dgamma",   pDGamma);
    row("delta_OS", pDeltaOS);
    row("delta_SS", pDeltaSS);
    row("Ddelta",   pDDelta);

    // kappa112 = Dgamma / (v2 * Ddelta), with its sub-sampling statistical error
    // (computed in Correlators.C, which carries the Dgamma-Ddelta covariance correctly).
    double kap=0, kapE=0;
    if(hKappa){ kap=hKappa->GetBinContent(1); kapE=hKappa->GetBinError(1); }
    printf("  %-10s %+13.4e %13.4e %9.2f\n", "kappa112", kap, kapE, kapE>0?fabs(kap)/kapE:0.0);
    printf("            ( kappa112 = Dgamma/(v2*Ddelta); error via sub-sampling, incl. Dgamma-Ddelta corr. )\n");

    printf("\n----------------  event-count projection (stat. error ~ 1/sqrt(Nev))  ----------------\n");
    proj( "Dgamma",   pDGamma->GetBinContent(1), pDGamma->GetBinError(1), Nev);
    proj( "Ddelta",   pDDelta->GetBinContent(1), pDDelta->GetBinError(1), Nev);
    if(hKappa) proj( "kappa112", kap, kapE, Nev);

    //------------------------------------ figures ------------------------------------
    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas("c1","v2",900,400);
    c1->Divide(2,1);
    c1->cd(1); pV2pT ->SetMarkerStyle(20); pV2pT ->SetMarkerColor(kBlue+1); pV2pT ->SetLineColor(kBlue+1); pV2pT ->Draw("E");
    c1->cd(2); pV2eta->SetMarkerStyle(20); pV2eta->SetMarkerColor(kRed+1);  pV2eta->SetLineColor(kRed+1);  pV2eta->Draw("E");
    c1->SaveAs("v2_pion_correlators.png");

    //gamma and delta OS-vs-SS with error bars, on separate y-scales (gamma is ~10x smaller than delta)
    TCanvas* c2 = new TCanvas("c2","gamma_delta",1000,450);
    c2->Divide(2,1);

    auto panel = [&](int pad, const char* title, TProfile* pOS, TProfile* pSS, TProfile* pD,
                     const char* lOS, const char* lSS, const char* lD){
        c2->cd(pad);
        gPad->SetGridy();
        double x[3]  = {1,2,3};
        double y[3]  = {pOS->GetBinContent(1), pSS->GetBinContent(1), pD->GetBinContent(1)};
        double ey[3] = {pOS->GetBinError(1),   pSS->GetBinError(1),   pD->GetBinError(1)};
        double ex[3] = {0,0,0};
        TGraphErrors* g = new TGraphErrors(3,x,y,ex,ey);
        g->SetTitle(Form("%s;;correlator",title));
        g->SetMarkerStyle(21); g->SetMarkerSize(1.6); g->SetMarkerColor(kAzure+2);
        g->SetLineColor(kAzure+2); g->SetLineWidth(2);
        g->GetXaxis()->SetLimits(0.5,3.5);
        g->GetXaxis()->SetNdivisions(0);
        //pad the y-range so points+errors and the zero line are visible
        double lo=1e30, hi=-1e30;
        for(int i=0;i<3;i++){ lo=min(lo,y[i]-ey[i]); hi=max(hi,y[i]+ey[i]); }
        lo=min(lo,0.0); hi=max(hi,0.0);
        double pad_=0.25*(hi-lo); g->GetYaxis()->SetRangeUser(lo-pad_,hi-0+pad_);
        g->Draw("AP");
        TLatex tl; tl.SetTextSize(0.045); tl.SetTextAlign(22);
        const char* lab[3]={lOS,lSS,lD};
        for(int i=0;i<3;i++) tl.DrawLatex(x[i], lo-pad_+0.04*(hi-lo+2*pad_), lab[i]);
    };
    panel(1,"#gamma correlator (30-40% Au+Au 200 GeV)",pGammaOS,pGammaSS,pDGamma,
          "#gamma_{OS}","#gamma_{SS}","#Delta#gamma");
    panel(2,"#delta correlator (30-40% Au+Au 200 GeV)",pDeltaOS,pDeltaSS,pDDelta,
          "#delta_{OS}","#delta_{SS}","#Delta#delta");
    if(hKappa){
        c2->cd(1);
        TLatex kl; kl.SetNDC(); kl.SetTextSize(0.040); kl.SetTextColor(kBlack);
        kl.DrawLatex(0.14,0.86,Form("#kappa_{112}=#Delta#gamma/(v_{2}#Delta#delta) = %.2f #pm %.2f", kap, kapE));
    }
    c2->SaveAs("gamma_delta_pion.png");

    printf("\nwrote v2_pion_correlators.png and gamma_delta_pion.png\n");
}
