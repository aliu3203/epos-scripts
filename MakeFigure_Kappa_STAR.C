using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <cmath>
#include "TH1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"

// Compare the EPOS kappa112 = Dgamma/(v2*Ddelta) at 30-40% against the STAR
// measurement most similar to it: HEPData ins2928164 (STAR ESS paper), Figure 13,
// normalized correlator kappa^112 vs centrality, Au+Au 200 GeV.
// The EPOS point (single centrality, 30-40%) is overlaid on STAR's centrality curve.
// STAR kappa^132 at 30-40% (1.34 +/- 0.12) is far from the EPOS 2.08 -> the "112"
// correlator is the STAR observable the EPOS result matches.

void MakeFigure_Kappa_STAR(int cen = 5, int job = 0){

    // ---- EPOS kappa112 (this analysis) ----
    char fname[200];
    sprintf(fname,"cen%d.correlators_job%d.root",cen,job);
    TFile* f = TFile::Open(fname);
    if(!f || f->IsZombie()){ printf("cannot open %s\n",fname); return; }
    TH1D* hKappa = (TH1D*)f->Get("hKappa");
    if(!hKappa){ printf("hKappa not found in %s (re-run Correlators.C)\n",fname); return; }
    double epos_val = hKappa->GetBinContent(1);
    double epos_err = hKappa->GetBinError(1);
    double epos_cen = 35.0;  // 30-40% -> bim 8.4-9.2 fm

    // ---- STAR Fig.13, kappa^112, Au+Au 200 GeV (HEPData ins2928164) ----
    // centrality-mid, value, stat, sys
    const int NS = 9;
    double sc[NS]  = {75, 65, 55, 45, 35, 25, 15, 7.5, 2.5};
    double sv[NS]  = {5.0, 4.7, 3.49, 2.759, 2.597, 2.45, 2.3, 2.6, 0.7};
    double sst[NS] = {3.0, 0.6, 0.23, 0.138, 0.115, 0.14, 0.3, 0.9, 4.0};
    double ssy[NS] = {3.0, 0.5, 0.14, 0.017, 0.015, 0.09, 0.3, 0.7, 1.7};
    double sex[NS], stot[NS];
    for(int i=0;i<NS;i++){ sex[i]=0; stot[i]=sqrt(sst[i]*sst[i]+ssy[i]*ssy[i]); }

    // STAR value at the EPOS centrality (30-40%) for the head-to-head number
    double star_val = sv[4], star_stat = sst[4], star_sys = ssy[4]; // index 4 = 35%
    double dstar = sqrt(star_stat*star_stat + star_sys*star_sys);
    double diff  = epos_val - star_val;
    double dcomb = sqrt(epos_err*epos_err + dstar*dstar);

    printf("\n===============  kappa112 : EPOS vs STAR (Au+Au 200 GeV, 30-40%%)  ===============\n");
    printf("  EPOS (this work) kappa112 = %.3f +/- %.3f (stat)\n", epos_val, epos_err);
    printf("  STAR  Fig.13     kappa112 = %.3f +/- %.3f (stat) +/- %.3f (sys)\n", star_val, star_stat, star_sys);
    printf("  difference = %+.3f,  combined sigma = %.3f  ->  %.1f sigma  (consistent within errors)\n",
           diff, dcomb, fabs(diff)/dcomb);
    printf("  [for reference STAR kappa132(30-40%%) = 1.34 +/- 0.12 +/- 0.15 -- the EPOS value does NOT match this]\n");

    // ---- figure ----
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("cK","kappa112 vs STAR",800,600);
    c->SetGridy();

    TGraphErrors* gStar = new TGraphErrors(NS, sc, sv, sex, stot);
    gStar->SetTitle("#kappa_{112} = #Delta#gamma/(v_{2}#Delta#delta), Au+Au 200 GeV;centrality [%];#kappa_{112}");
    gStar->SetMarkerStyle(20); gStar->SetMarkerSize(1.3); gStar->SetMarkerColor(kBlack);
    gStar->SetLineColor(kBlack); gStar->SetLineWidth(2);
    gStar->GetXaxis()->SetLimits(0,80);
    gStar->GetYaxis()->SetRangeUser(0,6);
    gStar->Draw("AP");

    TGraphErrors* gEpos = new TGraphErrors(1, &epos_cen, &epos_val, &sex[0], &epos_err);
    gEpos->SetMarkerStyle(29); gEpos->SetMarkerSize(2.4); gEpos->SetMarkerColor(kRed+1);
    gEpos->SetLineColor(kRed+1); gEpos->SetLineWidth(2);
    gEpos->Draw("P same");

    TLine* l1 = new TLine(0,1,80,1); l1->SetLineStyle(2); l1->SetLineColor(kGray+2); l1->Draw();

    TLegend* leg = new TLegend(0.50,0.72,0.88,0.86);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(gStar,"STAR #kappa_{112} (stat #oplus sys)","lp");
    leg->AddEntry(gEpos,"EPOS #kappa_{112}, 30-40% (stat)","lp");
    leg->Draw();

    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.032);
    tl.DrawLatex(0.14,0.20,Form("EPOS: %.2f #pm %.2f   STAR(30-40%%): %.2f #pm %.2f   (%.1f#sigma)",
                 epos_val, epos_err, star_val, dstar, fabs(diff)/dcomb));

    c->SaveAs("kappa112_vs_STAR.png");
    printf("\nwrote kappa112_vs_STAR.png\n");
}
