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
TGraphErrors* ProfileToGraph(TProfile* p, float xshift = 0.) {
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

void MakeFigure_Elliptic(int cen=5, int job=0) {

char fname[200];
sprintf(fname,"cen%d.v2pT_pion_job%d.root",cen,job);
TFile *f = new TFile(fname);

TProfile *p_v2Pt_pip   = (TProfile*)f->Get("pV2Pt_pip");
TProfile *p_v2Pt_pim   = (TProfile*)f->Get("pV2Pt_pim");
TProfile *p_v2Pt_pi    = (TProfile*)f->Get("pV2Pt_pi");
TProfile *p_v2Pt_pi0   = (TProfile*)f->Get("pV2Pt_pi0");
TProfile *p_v2Pt_pi_PP = (TProfile*)f->Get("pV2Pt_pi_PP");
TProfile *p_v2Eta_pi   = (TProfile*)f->Get("pV2Eta_pi");

TH1D *h_Pt_pi   = (TH1D*)f->Get("Hist_Pt_pi");
TH1D *h_dPhi_pi = (TH1D*)f->Get("Hist_dPhi_pi");
TH1D *h_b       = (TH1D*)f->Get("Hist_b");
TH1D *h_Psi     = (TH1D*)f->Get("Hist_Psi");

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

//========================= canvas 1: v2(pT) pi+ / pi- / pi0 =========================
TCanvas* can1 = new TCanvas("Flow_v2pT", "Flow_v2pT", 40,40,740,580);
 can1-> SetTopMargin(0.06);
 can1-> SetRightMargin(0.02);
 can1-> SetBottomMargin(0.17);
 can1-> SetLeftMargin(0.13);
 can1->Draw();

  gPad->SetGridx(0);
  gPad->SetGridy(0);
   gPad->SetTickx();
   gPad->SetTicky();
  TH1F* histGraph = new TH1F("Fl", "", 24, 0, 3.0);
  histGraph->SetMaximum(0.14);
  histGraph->SetMinimum(-0.01);
  histGraph->SetLineColor(kBlack);
  histGraph->GetYaxis()->SetTitleOffset(0.80);
  histGraph->GetYaxis()->SetTitleSize(0.075);
  histGraph->GetXaxis()->SetTitleSize(0.08);
  histGraph->GetXaxis()->SetTitleOffset(0.90);
  histGraph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histGraph->GetYaxis()->SetTitle("v_{2}");
  histGraph->GetXaxis()->SetNdivisions(6);
  histGraph->GetYaxis()->SetNdivisions(605);
  double lsize=histGraph->GetLabelSize();
  histGraph->GetYaxis()->SetLabelSize(lsize*1.2);
  histGraph->GetXaxis()->SetLabelSize(lsize*1.2);
  histGraph->Draw();

TGraphErrors* g_pip = ProfileToGraph(p_v2Pt_pip, -0.008);
  g_pip->SetMarkerStyle(kOpenCircle);
  g_pip->SetMarkerSize(1.2);
  g_pip->SetMarkerColor(2);
  g_pip->SetLineColor(2);
  g_pip->SetLineWidth(2);
  g_pip->Draw("pe1");

TGraphErrors* g_pim = ProfileToGraph(p_v2Pt_pim, 0.008);
  g_pim->SetMarkerStyle(kOpenSquare);
  g_pim->SetMarkerSize(1.2);
  g_pim->SetMarkerColor(4);
  g_pim->SetLineColor(4);
  g_pim->SetLineWidth(2);
  g_pim->Draw("pe1");

TGraphErrors* g_pi0 = ProfileToGraph(p_v2Pt_pi0, 0.);
  g_pi0->SetMarkerStyle(kOpenDiamond);
  g_pi0->SetMarkerSize(1.4);
  g_pi0->SetMarkerColor(kGreen+2);
  g_pi0->SetLineColor(kGreen+2);
  g_pi0->SetLineWidth(2);
  g_pi0->Draw("pe1");

TLegend* leg1 = new TLegend(0.62,0.20,0.92,0.42);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetTextSize(0.055);
  leg1->AddEntry(g_pip,"#pi^{+}","pe");
  leg1->AddEntry(g_pim,"#pi^{-}","pe");
  leg1->AddEntry(g_pi0,"#pi^{0}","pe");
  leg1->Draw();

   TLatex *tex = new TLatex(0.15,0.125,"30 - 40% Au+Au 200 GeV (EPOS)");
   tex->SetTextSize(0.065);
   tex->SetTextColor(1);
   tex->Draw();

   tex = new TLatex(0.15,0.110,"v_{2}{RP}, |#eta| < 1");
   tex->SetTextSize(0.055);
   tex->SetTextColor(1);
   tex->Draw();

can1->SaveAs("v2_pT_pion.png");

//========================= canvas 2: v2(pT) RP vs PP =========================
TCanvas* can2 = new TCanvas("Flow_v2pT_plane", "Flow_v2pT_plane", 40,40,740,580);
 can2-> SetTopMargin(0.06);
 can2-> SetRightMargin(0.02);
 can2-> SetBottomMargin(0.17);
 can2-> SetLeftMargin(0.13);
 can2->Draw();

  gPad->SetGridx(0);
  gPad->SetGridy(0);
   gPad->SetTickx();
   gPad->SetTicky();
  TH1F* histGrap = new TH1F("Flo", "", 24, 0, 3.0);
  histGrap->SetMaximum(0.14);
  histGrap->SetMinimum(-0.01);
  histGrap->SetLineColor(kBlack);
  histGrap->GetYaxis()->SetTitleOffset(0.80);
  histGrap->GetYaxis()->SetTitleSize(0.075);
  histGrap->GetXaxis()->SetTitleSize(0.08);
  histGrap->GetXaxis()->SetTitleOffset(0.90);
  histGrap->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histGrap->GetYaxis()->SetTitle("v_{2} (#pi^{#pm})");
  histGrap->GetXaxis()->SetNdivisions(6);
  histGrap->GetYaxis()->SetNdivisions(605);
  lsize=histGrap->GetLabelSize();
  histGrap->GetYaxis()->SetLabelSize(lsize*1.2);
  histGrap->GetXaxis()->SetLabelSize(lsize*1.2);
  histGrap->Draw();

TGraphErrors* g_RP = ProfileToGraph(p_v2Pt_pi, -0.008);
  g_RP->SetMarkerStyle(kOpenCircle);
  g_RP->SetMarkerSize(1.2);
  g_RP->SetMarkerColor(4);
  g_RP->SetLineColor(4);
  g_RP->SetLineWidth(2);
  g_RP->Draw("pe1");

TGraphErrors* g_PP = ProfileToGraph(p_v2Pt_pi_PP, 0.008);
  g_PP->SetMarkerStyle(kFullCircle);
  g_PP->SetMarkerSize(1.2);
  g_PP->SetMarkerColor(2);
  g_PP->SetLineColor(2);
  g_PP->SetLineWidth(2);
  g_PP->Draw("pe1");

TLegend* leg2 = new TLegend(0.58,0.20,0.92,0.40);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.05);
  leg2->AddEntry(g_RP,"v_{2}{RP} (#Psi = #phi_{ev})","pe");
  leg2->AddEntry(g_PP,"v_{2}{PP} (#Psi = #phi_{ev}+#phi_{r})","pe");
  leg2->Draw();

   tex = new TLatex(0.15,0.125,"30 - 40% Au+Au 200 GeV (EPOS)");
   tex->SetTextSize(0.065);
   tex->SetTextColor(1);
   tex->Draw();

   tex = new TLatex(0.15,0.110,"#pi^{+}+#pi^{-}, |#eta| < 1");
   tex->SetTextSize(0.055);
   tex->SetTextColor(1);
   tex->Draw();

can2->SaveAs("v2_pT_pion_RPvsPP.png");

//========================= canvas 3: v2(eta) =========================
TCanvas* can3 = new TCanvas("Flow_v2eta", "Flow_v2eta", 40,40,740,580);
 can3-> SetTopMargin(0.06);
 can3-> SetRightMargin(0.02);
 can3-> SetBottomMargin(0.17);
 can3-> SetLeftMargin(0.13);
 can3->Draw();

  gPad->SetGridx(0);
  gPad->SetGridy(0);
   gPad->SetTickx();
   gPad->SetTicky();
  TH1F* histGra = new TH1F("Flow", "", 24, -6, 6);
  histGra->SetMaximum(0.06);
  histGra->SetMinimum(-0.01);
  histGra->SetLineColor(kBlack);
  histGra->GetYaxis()->SetTitleOffset(0.80);
  histGra->GetYaxis()->SetTitleSize(0.075);
  histGra->GetXaxis()->SetTitleSize(0.08);
  histGra->GetXaxis()->SetTitleOffset(0.90);
  histGra->GetXaxis()->SetTitle("#eta");
  histGra->GetYaxis()->SetTitle("v_{2} (#pi^{#pm})");
  histGra->GetXaxis()->SetNdivisions(6);
  histGra->GetYaxis()->SetNdivisions(605);
  lsize=histGra->GetLabelSize();
  histGra->GetYaxis()->SetLabelSize(lsize*1.2);
  histGra->GetXaxis()->SetLabelSize(lsize*1.2);
  histGra->Draw();

TGraphErrors* g_eta = ProfileToGraph(p_v2Eta_pi, 0.);
  g_eta->SetMarkerStyle(kFullCircle);
  g_eta->SetMarkerSize(1.2);
  g_eta->SetMarkerColor(4);
  g_eta->SetLineColor(4);
  g_eta->SetLineWidth(2);
  g_eta->Draw("pe1");

   tex = new TLatex(-5.5,0.053,"30 - 40% Au+Au 200 GeV (EPOS)");
   tex->SetTextSize(0.065);
   tex->SetTextColor(1);
   tex->Draw();

   tex = new TLatex(-5.5,0.046,"v_{2}{RP}, p_{T} > 0.05 GeV/c");
   tex->SetTextSize(0.055);
   tex->SetTextColor(1);
   tex->Draw();

can3->SaveAs("v2_eta_pion.png");

//========================= canvas 4: QA distributions =========================
gStyle->SetOptTitle(1);	//re-enable pad titles for the QA panels
TCanvas *c4 = new TCanvas("c4", "QA (pions)", 1200, 800);
c4->Divide(2, 2);

c4->cd(1);
gPad->SetLogy();
h_Pt_pi->SetTitle("Transverse Momentum #pi^{#pm} (pT)");
h_Pt_pi->Draw();

c4->cd(2);
h_dPhi_pi->SetTitle("#phi - #Psi_{RP} (#pi^{#pm})");
h_dPhi_pi->SetMinimum(0);
h_dPhi_pi->Draw();

c4->cd(3);
h_b->SetTitle("Impact Parameter (b)");
h_b->Draw();

c4->cd(4);
h_Psi->SetTitle("#Psi_{RP}");
h_Psi->SetMinimum(0);
h_Psi->Draw();

c4->SaveAs("v2_QA_pion.png");

return;
}
