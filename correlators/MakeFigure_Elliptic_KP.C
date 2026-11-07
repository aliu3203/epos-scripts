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

// Turn a TProfile into a TGraphErrors, skipping empty bins
TGraphErrors* ProfileToGraphKP(TProfile* p, float xshift = 0.) {
    TGraphErrors* g = new TGraphErrors();
    int ip = 0;
    for(int i = 1; i <= p->GetNbinsX(); i++) {
        if(p->GetBinEntries(i) < 1) continue;
        g->SetPoint(ip, p->GetBinCenter(i) + xshift, p->GetBinContent(i));
        g->SetPointError(ip, 0, p->GetBinError(i));
        ip++;
    }
    return g;
}

// Styled canvas + axis frame shared by all the v2 figures
TCanvas* MakeFrame(const char* cname, TH1F*& frame, float xlo, float xhi,
                   float ylo, float yhi, const char* xtitle, const char* ytitle) {
    TCanvas* can = new TCanvas(cname, cname, 40,40,740,580);
    can-> SetTopMargin(0.06);
    can-> SetRightMargin(0.02);
    can-> SetBottomMargin(0.17);
    can-> SetLeftMargin(0.13);
    can->Draw();

    gPad->SetGridx(0);
    gPad->SetGridy(0);
    gPad->SetTickx();
    gPad->SetTicky();
    frame = new TH1F(Form("frame_%s",cname), "", 24, xlo, xhi);
    frame->SetMaximum(yhi);
    frame->SetMinimum(ylo);
    frame->SetLineColor(kBlack);
    frame->GetYaxis()->SetTitleOffset(0.80);
    frame->GetYaxis()->SetTitleSize(0.075);
    frame->GetXaxis()->SetTitleSize(0.08);
    frame->GetXaxis()->SetTitleOffset(0.90);
    frame->GetXaxis()->SetTitle(xtitle);
    frame->GetYaxis()->SetTitle(ytitle);
    frame->GetXaxis()->SetNdivisions(6);
    frame->GetYaxis()->SetNdivisions(605);
    double lsize = frame->GetLabelSize();
    frame->GetYaxis()->SetLabelSize(lsize*1.2);
    frame->GetXaxis()->SetLabelSize(lsize*1.2);
    frame->Draw();
    return can;
}

void SetGraphStyle(TGraphErrors* g, int marker, int color) {
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.2);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
}

void MakeFigure_Elliptic_KP(int cen=5, int job=0) {

char fname[200];
sprintf(fname,"cen%d.v2pT_kaon_proton_job%d.root",cen,job);
TFile *f = new TFile(fname);

TProfile *p_v2Pt_Kp    = (TProfile*)f->Get("pV2Pt_Kp");
TProfile *p_v2Pt_Km    = (TProfile*)f->Get("pV2Pt_Km");
TProfile *p_v2Pt_K     = (TProfile*)f->Get("pV2Pt_K");
TProfile *p_v2Pt_K_PP  = (TProfile*)f->Get("pV2Pt_K_PP");
TProfile *p_v2Pt_p     = (TProfile*)f->Get("pV2Pt_p");
TProfile *p_v2Pt_pbar  = (TProfile*)f->Get("pV2Pt_pbar");
TProfile *p_v2Pt_pr    = (TProfile*)f->Get("pV2Pt_pr");
TProfile *p_v2Pt_pr_PP = (TProfile*)f->Get("pV2Pt_pr_PP");
TProfile *p_v2Eta_K    = (TProfile*)f->Get("pV2Eta_K");
TProfile *p_v2Eta_pr   = (TProfile*)f->Get("pV2Eta_pr");

TH1D *h_Pt_K    = (TH1D*)f->Get("Hist_Pt_K");
TH1D *h_Pt_pr   = (TH1D*)f->Get("Hist_Pt_pr");
TH1D *h_dPhi_K  = (TH1D*)f->Get("Hist_dPhi_K");
TH1D *h_dPhi_pr = (TH1D*)f->Get("Hist_dPhi_pr");

    gStyle->SetPalette(1);
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
    gStyle->SetPadGridX( false );
    gStyle->SetPadGridY( false );
    gStyle->SetLabelSize(.05,"X");
    gStyle->SetLabelSize(.05,"Y");
    gStyle->SetTitleSize(.06,"X");
    gStyle->SetTitleSize(.06,"Y");

  gStyle->SetPadBorderMode (0);
  gStyle->SetPadColor (kWhite);
  gStyle->SetCanvasColor (kWhite);
  gStyle->SetFrameBorderMode (0);
  gStyle->SetCanvasBorderMode (0);
gStyle->SetTitleBorderSize (0);

TH1F* frame = 0;
TLatex* tex = 0;

//========================= canvas 1a: v2(pT) K+ / K- =========================
TCanvas* can1a = MakeFrame("Flow_v2pT_K", frame, 0, 3.0, -0.02, 0.14, "p_{T} (GeV/c)", "v_{2}");

TGraphErrors* g_Kp = ProfileToGraphKP(p_v2Pt_Kp, -0.008);
  SetGraphStyle(g_Kp, kOpenCircle, 2);
  g_Kp->Draw("pe1");

TGraphErrors* g_Km = ProfileToGraphKP(p_v2Pt_Km, 0.008);
  SetGraphStyle(g_Km, kOpenSquare, 4);
  g_Km->Draw("pe1");

TLegend* leg1a = new TLegend(0.70,0.20,0.92,0.36);
  leg1a->SetBorderSize(0);
  leg1a->SetFillStyle(0);
  leg1a->SetTextSize(0.055);
  leg1a->AddEntry(g_Kp,"K^{+}","pe");
  leg1a->AddEntry(g_Km,"K^{-}","pe");
  leg1a->Draw();

   tex = new TLatex(0.15,0.125,"30 - 40% Au+Au 200 GeV (EPOS)");
   tex->SetTextSize(0.065);
   tex->SetTextColor(1);
   tex->Draw();

   tex = new TLatex(0.15,0.110,"v_{2}{RP}, |#eta| < 1");
   tex->SetTextSize(0.055);
   tex->SetTextColor(1);
   tex->Draw();

can1a->SaveAs("v2_pT_kaon.png");

//========================= canvas 1b: v2(pT) p / pbar =========================
TCanvas* can1b = MakeFrame("Flow_v2pT_p", frame, 0, 3.0, -0.02, 0.18, "p_{T} (GeV/c)", "v_{2}");

TGraphErrors* g_p = ProfileToGraphKP(p_v2Pt_p, -0.008);
  SetGraphStyle(g_p, kFullCircle, kGreen+2);
  g_p->Draw("pe1");

TGraphErrors* g_pbar = ProfileToGraphKP(p_v2Pt_pbar, 0.008);
  SetGraphStyle(g_pbar, kFullSquare, kMagenta+1);
  g_pbar->Draw("pe1");

TLegend* leg1b = new TLegend(0.70,0.20,0.92,0.36);
  leg1b->SetBorderSize(0);
  leg1b->SetFillStyle(0);
  leg1b->SetTextSize(0.055);
  leg1b->AddEntry(g_p,"p","pe");
  leg1b->AddEntry(g_pbar,"#bar{p}","pe");
  leg1b->Draw();

   tex = new TLatex(0.15,0.160,"30 - 40% Au+Au 200 GeV (EPOS)");
   tex->SetTextSize(0.065);
   tex->SetTextColor(1);
   tex->Draw();

   tex = new TLatex(0.15,0.140,"v_{2}{RP}, |#eta| < 1");
   tex->SetTextSize(0.055);
   tex->SetTextColor(1);
   tex->Draw();

can1b->SaveAs("v2_pT_proton.png");

//========================= canvas 2a: v2(pT) RP vs PP, K =========================
TCanvas* can2a = MakeFrame("Flow_v2pT_K_plane", frame, 0, 3.0, -0.02, 0.14, "p_{T} (GeV/c)", "v_{2} (K^{#pm})");

TGraphErrors* g_K_RP = ProfileToGraphKP(p_v2Pt_K, -0.008);
  SetGraphStyle(g_K_RP, kOpenCircle, 4);
  g_K_RP->Draw("pe1");

TGraphErrors* g_K_PP = ProfileToGraphKP(p_v2Pt_K_PP, 0.008);
  SetGraphStyle(g_K_PP, kFullCircle, 2);
  g_K_PP->Draw("pe1");

TLegend* leg2a = new TLegend(0.58,0.20,0.92,0.40);
  leg2a->SetBorderSize(0);
  leg2a->SetFillStyle(0);
  leg2a->SetTextSize(0.05);
  leg2a->AddEntry(g_K_RP,"v_{2}{RP} (#Psi = #phi_{ev})","pe");
  leg2a->AddEntry(g_K_PP,"v_{2}{PP} (#Psi = #phi_{ev}+#phi_{r})","pe");
  leg2a->Draw();

   tex = new TLatex(0.15,0.125,"30 - 40% Au+Au 200 GeV (EPOS)");
   tex->SetTextSize(0.065);
   tex->SetTextColor(1);
   tex->Draw();

   tex = new TLatex(0.15,0.110,"K^{+}+K^{-}, |#eta| < 1");
   tex->SetTextSize(0.055);
   tex->SetTextColor(1);
   tex->Draw();

can2a->SaveAs("v2_pT_kaon_RPvsPP.png");

//========================= canvas 2b: v2(pT) RP vs PP, p =========================
TCanvas* can2b = MakeFrame("Flow_v2pT_p_plane", frame, 0, 3.0, -0.02, 0.18, "p_{T} (GeV/c)", "v_{2} (p+#bar{p})");

TGraphErrors* g_pr_RP = ProfileToGraphKP(p_v2Pt_pr, -0.008);
  SetGraphStyle(g_pr_RP, kOpenSquare, 4);
  g_pr_RP->Draw("pe1");

TGraphErrors* g_pr_PP = ProfileToGraphKP(p_v2Pt_pr_PP, 0.008);
  SetGraphStyle(g_pr_PP, kFullSquare, 2);
  g_pr_PP->Draw("pe1");

TLegend* leg2b = new TLegend(0.58,0.20,0.92,0.40);
  leg2b->SetBorderSize(0);
  leg2b->SetFillStyle(0);
  leg2b->SetTextSize(0.05);
  leg2b->AddEntry(g_pr_RP,"v_{2}{RP} (#Psi = #phi_{ev})","pe");
  leg2b->AddEntry(g_pr_PP,"v_{2}{PP} (#Psi = #phi_{ev}+#phi_{r})","pe");
  leg2b->Draw();

   tex = new TLatex(0.15,0.160,"30 - 40% Au+Au 200 GeV (EPOS)");
   tex->SetTextSize(0.065);
   tex->SetTextColor(1);
   tex->Draw();

   tex = new TLatex(0.15,0.140,"p+#bar{p}, |#eta| < 1");
   tex->SetTextSize(0.055);
   tex->SetTextColor(1);
   tex->Draw();

can2b->SaveAs("v2_pT_proton_RPvsPP.png");

//========================= canvas 3a: v2(eta) K =========================
TCanvas* can3a = MakeFrame("Flow_v2eta_K", frame, -6, 6, -0.02, 0.06, "#eta", "v_{2} (K^{#pm})");

TGraphErrors* g_etaK = ProfileToGraphKP(p_v2Eta_K, 0.);
  SetGraphStyle(g_etaK, kFullCircle, 4);
  g_etaK->Draw("pe1");

   tex = new TLatex(-5.5,0.053,"30 - 40% Au+Au 200 GeV (EPOS)");
   tex->SetTextSize(0.065);
   tex->SetTextColor(1);
   tex->Draw();

   tex = new TLatex(-5.5,0.046,"K^{#pm}, v_{2}{RP}, p_{T} > 0.05 GeV/c");
   tex->SetTextSize(0.055);
   tex->SetTextColor(1);
   tex->Draw();

can3a->SaveAs("v2_eta_kaon.png");

//========================= canvas 3b: v2(eta) p =========================
TCanvas* can3b = MakeFrame("Flow_v2eta_p", frame, -6, 6, -0.02, 0.08, "#eta", "v_{2} (p+#bar{p})");

TGraphErrors* g_etaPr = ProfileToGraphKP(p_v2Eta_pr, 0.);
  SetGraphStyle(g_etaPr, kFullSquare, kGreen+2);
  g_etaPr->Draw("pe1");

   tex = new TLatex(-5.5,0.070,"30 - 40% Au+Au 200 GeV (EPOS)");
   tex->SetTextSize(0.065);
   tex->SetTextColor(1);
   tex->Draw();

   tex = new TLatex(-5.5,0.061,"p+#bar{p}, v_{2}{RP}, p_{T} > 0.05 GeV/c");
   tex->SetTextSize(0.055);
   tex->SetTextColor(1);
   tex->Draw();

can3b->SaveAs("v2_eta_proton.png");

//========================= canvas 4: QA distributions =========================
gStyle->SetOptTitle(1);	//re-enable pad titles for the QA panels
TCanvas *c4 = new TCanvas("c4", "QA (kaons/protons)", 1200, 800);
c4->Divide(2, 2);

c4->cd(1);
gPad->SetLogy();
h_Pt_K->SetTitle("Transverse Momentum K^{#pm} (pT)");
h_Pt_K->Draw();

c4->cd(2);
h_dPhi_K->SetTitle("#phi - #Psi_{RP} (K^{#pm})");
h_dPhi_K->SetMinimum(0);
h_dPhi_K->Draw();

c4->cd(3);
gPad->SetLogy();
h_Pt_pr->SetTitle("Transverse Momentum p+#bar{p} (pT)");
h_Pt_pr->Draw();

c4->cd(4);
h_dPhi_pr->SetTitle("#phi - #Psi_{RP} (p+#bar{p})");
h_dPhi_pr->SetMinimum(0);
h_dPhi_pr->Draw();

c4->SaveAs("v2_QA_kaon_proton.png");

return;
}
