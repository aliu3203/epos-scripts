#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <TChain.h>
#include "TLeaf.h"
#include "TF1.h"
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
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

// Turn a TProfile into a TGraphErrors, skipping empty/low-stat bins
TGraphErrors* ProfToGraph(TProfile* p, double minEntries = 20) {
    TGraphErrors* g = new TGraphErrors();
    int ip = 0;
    for(int i = 1; i <= p->GetNbinsX(); i++) {
        if(p->GetBinEntries(i) < minEntries) continue;
        g->SetPoint(ip, p->GetBinCenter(i), p->GetBinContent(i));
        g->SetPointError(ip, 0, p->GetBinError(i));
        ip++;
    }
    return g;
}

void DrawPanel(TVirtualPad* pad, TProfile* prof, const char* xtitle, const char* ytitle,
               const char* tag, int color, int marker, double xmax, double ymax) {
    pad->cd();
    pad->SetTopMargin(0.04);
    pad->SetRightMargin(0.03);
    pad->SetBottomMargin(0.15);
    pad->SetLeftMargin(0.16);
    pad->SetTickx();
    pad->SetTicky();

    TH1F* fr = new TH1F(Form("fr_%s",tag), "", 24, 0, xmax);
    fr->SetMinimum(0.0);
    fr->SetMaximum(ymax);
    fr->GetXaxis()->SetTitle(xtitle);
    fr->GetYaxis()->SetTitle(ytitle);
    fr->GetYaxis()->SetTitleOffset(1.05);
    fr->GetYaxis()->SetTitleSize(0.065);
    fr->GetXaxis()->SetTitleSize(0.065);
    fr->GetXaxis()->SetTitleOffset(1.00);
    fr->GetXaxis()->SetNdivisions(505);
    fr->GetYaxis()->SetNdivisions(505);
    fr->GetXaxis()->SetLabelSize(0.055);
    fr->GetYaxis()->SetLabelSize(0.055);
    fr->Draw();

    // linear fit to guide the eye
    TF1* fit = new TF1(Form("fit_%s",tag),"[0]+[1]*x",0,xmax);
    fit->SetLineColor(color);
    fit->SetLineStyle(2);
    fit->SetLineWidth(2);

    TGraphErrors* g = ProfToGraph(prof);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.3);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
    g->Fit(fit,"Q0");
    g->Draw("pe1");
    fit->Draw("same");

    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextColor(color);
    tex->SetTextSize(0.060);
    tex->DrawLatex(0.20, 0.88, tag);

    tex->SetTextColor(1);
    tex->SetTextSize(0.048);
    tex->DrawLatex(0.20, 0.80, Form("slope = %.3f #pm %.3f", fit->GetParameter(1), fit->GetParError(1)));
    tex->DrawLatex(0.20, 0.73, Form("intercept = %.3f", fit->GetParameter(0)));
}

void MakeFigure_Elliptic_ESS(int cen=5, int job=0) {

char fname[200];
sprintf(fname,"cen%d.v2q2_ESS_job%d.root",cen,job);
TFile *f = new TFile(fname);

TProfile* pSS = (TProfile*)f->Get("pV2s_q2s");	// single v2 vs single q2^2
TProfile* pSP = (TProfile*)f->Get("pV2p_q2s");	// pair   v2 vs single q2^2
TProfile* pPS = (TProfile*)f->Get("pV2s_q2p");	// single v2 vs pair   q2^2
TProfile* pPP = (TProfile*)f->Get("pV2p_q2p");	// pair   v2 vs pair   q2^2
TH1D* hConst  = (TH1D*)f->Get("hConst");

double v2_2  = hConst ? hConst->GetBinContent(1) : 0;
double v2p_2 = hConst ? hConst->GetBinContent(2) : 0;
cout << "v2{2} = " << v2_2 << ",  v2pair{2} = " << v2p_2 << endl;

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    gStyle->SetOptDate(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);

// x ranges chosen so each q2^2 panel shows its well-populated region
double xmaxS = 3.5;	// single q2^2
double xmaxP = 3.5;	// pair q2^2
double ymax  = 0.18;

TCanvas* c = new TCanvas("cESS","ESS v2 vs q2^2",1000,900);
c->Divide(2,2,0.001,0.001);

// (a) single q2^2  ->  single v2   [unmixed]
DrawPanel(c->cd(1), pSS, "single q_{2}^{2}", "single v_{2}", "(a) single q_{2}^{2}, single v_{2}",
          2, kFullCircle, xmaxS, ymax);
// (b) single q2^2  ->  pair v2     [mixed]
DrawPanel(c->cd(2), pSP, "single q_{2}^{2}", "pair v_{2}",   "(b) single q_{2}^{2}, pair v_{2}",
          4, kFullSquare, xmaxS, ymax);
// (c) pair q2^2    ->  single v2   [mixed]
DrawPanel(c->cd(3), pPS, "pair q_{2}^{2}",   "single v_{2}", "(c) pair q_{2}^{2}, single v_{2}",
          kGreen+2, kFullTriangleUp, xmaxP, ymax);
// (d) pair q2^2    ->  pair v2     [unmixed]
DrawPanel(c->cd(4), pPP, "pair q_{2}^{2}",   "pair v_{2}",   "(d) pair q_{2}^{2}, pair v_{2}",
          kMagenta+1, kFullDiamond, xmaxP, ymax);

c->cd(1);
TLatex* head = new TLatex();
head->SetNDC();
head->SetTextSize(0.050);
head->SetTextColor(1);
head->DrawLatex(0.20, 0.29, "30-40% Au+Au 200 GeV (EPOS)");
head->DrawLatex(0.20, 0.23, "#pi^{#pm}, |y|<1, 0.2<p_{T}<2");
head->DrawLatex(0.20, 0.17, Form("v_{2}{2}=%.3f, v_{2,pair}{2}=%.3f", v2_2, v2p_2));

c->SaveAs("v2_vs_q2_ESS.png");

return;
}
