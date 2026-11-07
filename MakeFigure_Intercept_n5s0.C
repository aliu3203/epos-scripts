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
#include "TCanvas.h"
#include "TGraphErrors.h"

void MakeFigure_Intercept_n5s0(int cen=5) {
TF1 *fun = new TF1("fun","[0]",0,2);

char fname2[200];
sprintf(fname2, "cen0.b.root");
TFile *f2 = new TFile(fname2);


char fname[200];
sprintf(fname,"cen5.AVFD_Parity112_EPD_pT02_lab_job0.root");
//sprintf(fname,"cen%d.AVFD_Parity112_EPD_pT02_lab_ns0.parentv2_57_new.root",cen);
//sprintf(fname,"cen%d.AVFD_Isobar_Parity112_EPD_pT02_lab_ns0.parentv2_57.root",cen);
//sprintf(fname,"Test_Yunshan_00.root");


TFile *f = new TFile(fname);

TH1D *h_Mult_All = (TH1D*)f->Get("hMult_All");
TH1D *h_Mult_ptCut = (TH1D*)f->Get("hMult_ptCut");
TH1D *h_Mult_ptEtaCut = (TH1D*)f->Get("hMult_ptEtaCut");

TH1D *h_Pt           = (TH1D*)f->Get("Hist_Pt");
TH1D *h_b            = (TH1D*)f2->Get("Hist_b");     // The one you added
TH1D *h_Phi          = (TH1D*)f->Get("Hist_Phi");   // The one you added
TH2D *h_EtaPt        = (TH2D*)f->Get("EtaPtDist");  // Existing 2D Eta-Pt
TH1D *h_y            = (TH1D*)f->Get("Hist_Y");

// To get a 1D distribution of Eta from the 2D EtaPtDist histogram:
TH1D *h_Eta = h_EtaPt->ProjectionX("h_Eta");

// 2. Plot them on a new Canvas
TCanvas *c1 = new TCanvas("c1", "Distributions for pT > 0.05", 1200, 800);
c1->Divide(3, 2); // Split canvas into 6 pads

c1->cd(1);
h_b->SetTitle("Impact Parameter (b)");
h_b->Draw();
/*
c1->cd(2);
h_Multiplicity->SetTitle("Multiplicity");
h_Multiplicity->Draw();
*/
c1->cd(2);
gPad->SetLogy();
h_Pt->GetXaxis()->SetRangeUser(0.1, h_Pt->GetXaxis()->GetXmax());
h_Pt->SetTitle("Transverse Momentum (pT)");
h_Pt->Draw();

c1->cd(3);
h_Eta->SetTitle("Pseudorapidity (Eta)");
h_Eta->Draw();

c1->cd(4);
h_Phi->SetTitle("Azimuthal Angle (Phi)");
h_Phi->SetMinimum(0);
h_Phi->Draw();

c1->cd(5);
h_y->SetTitle("Rapidity (y)");
h_y->GetXaxis()->SetRangeUser(-6.0, 6.0); // Usually matches Eta range conceptually
h_y->Draw();

c1->SaveAs("Distributions.png");

TCanvas *c2 = new TCanvas("c2", "Multiplicity (All Particles)", 1200, 800);

c2->cd(1);
h_Mult_All->SetTitle("Multiplicity (All Particles)");
//gPad->SetLogy();
h_Mult_All->GetXaxis()->SetRangeUser(0,1500); 

h_Mult_All->Draw();

c2->SaveAs("Distributions2.png");

TCanvas *c3 = new TCanvas("c3", "Multiplicity (pT > 0.05)", 1200, 800);

c3->cd(1);
h_Mult_ptCut->GetXaxis()->SetRangeUser(0, 1200); 
h_Mult_ptCut->SetTitle("Multiplicity (pT > 0.05)");
h_Mult_ptCut->Draw();

c3->SaveAs("Distributions3.png");

TCanvas *c4 = new TCanvas("c4", "Multiplicity (pT > 0.05, Eta < 0.5)", 1200, 500);

c4->cd(1);
gPad->SetLogy();
h_Mult_ptEtaCut->GetXaxis()->SetRangeUser(0, 200); 
h_Mult_ptEtaCut->SetTitle("Multiplicity (pT > 0.05, Eta < 0.5)");
h_Mult_ptEtaCut->Draw();

c4->SaveAs("Distributions4.png");
TH1F *hSpeciesAll = (TH1F*)f->Get("hSpeciesAll");
if (!hSpeciesAll) {
    std::cout << "Error: hSpeciesAll could not be found in the ROOT file!" << std::endl;
    return; // Stop execution before the seg fault happens
}

TCanvas *cPie = new TCanvas("cPie", "Particle Species Distribution", 800, 600);
cPie->cd(1);
// Clean up the histogram before drawing
hSpeciesAll->LabelsDeflate("X"); // Trims any unused empty bins
hSpeciesAll->LabelsOption("d");  // Sorts the slices from largest to smallest

//hSpeciesAll->GetXaxis()->SetRange(15, 30);
hSpeciesAll->Draw("PIE L");
// Draw as a pie chart with labels
cPie->SaveAs("PIE.png");

return;







TProfile *Hist_cos = (TProfile*)f->Get("Hist_cos");
Hist_cos->GetEntries();
float a1_positive = Hist_cos->GetBinContent(3);
float a1_positive_err = Hist_cos->GetBinError(3);
float a1_negative = Hist_cos->GetBinContent(4);
float a1_negative_err = Hist_cos->GetBinError(4);
cout<<"a1(pi+) = "<<a1_positive<<" +/- "<<a1_positive_err<<endl;
cout<<"a1(pi-) = "<<a1_negative<<" +/- "<<a1_negative_err<<endl;
cout<<endl;

TProfile2D *Parity_Q_R   = (TProfile2D*)f->Get("Parity_Q");
TH1* Parity_Q_ss_R = (TH1*)Parity_Q_R->ProjectionY("Parity_Q_ss_R",3,3);
TH1* Parity_Q_os_R = (TH1*)Parity_Q_R->ProjectionY("Parity_Q_os_R",4,4);
TProfile2D *Parity_Q_Rs = (TProfile2D*)f->Get("Parity_Q2");
TH1* Parity_Q_ss_Rs = (TH1*)Parity_Q_Rs->ProjectionY("Parity_Q_ss_Rs",3,3);
TH1* Parity_Q_os_Rs = (TH1*)Parity_Q_Rs->ProjectionY("Parity_Q_os_Rs",4,4);

//single Q
TProfile* p_a1_Q2_p = (TProfile*)f->Get("p_a1_Q2_p");
TProfile* p_a1_Q2_n = (TProfile*)f->Get("p_a1_Q2_n");
TProfile *p_v2p_Qs  = (TProfile*)f->Get("p_v2_p_Q2");
TProfile *p_v2s_Qs  = (TProfile*)f->Get("p_v2_Q2");

//pair Q
TProfile* p_a1_QQ_p = (TProfile*)f->Get("p_a1_QQ_p");
TProfile* p_a1_QQ_n = (TProfile*)f->Get("p_a1_QQ_n");
TProfile *p_v2s_Qp  = (TProfile*)f->Get("p_v2_s_Q");
TProfile *p_v2p_Qp  = (TProfile*)f->Get("p_v2_Q");

const int N_bins = 50;
float v2s_QQ[N_bins], v2s_QQ_err[N_bins], v2s_Q2[N_bins], v2s_Q2_err[N_bins];
float v2p_QQ[N_bins], v2p_QQ_err[N_bins], v2p_Q2[N_bins], v2p_Q2_err[N_bins];
float dg_QQ[N_bins], dg_QQ_err[N_bins], dg_Q2[N_bins], dg_Q2_err[N_bins];
float a1_QQ[N_bins], a1_QQ_err[N_bins], a1_Q2[N_bins], a1_Q2_err[N_bins];
for(int i=0;i<N_bins;i++) {
	v2s_QQ[i] = p_v2s_Qp->GetBinContent(i+1) / 100.;
	v2s_QQ_err[i] = p_v2s_Qp->GetBinError(i+1) / 100.;
	v2s_Q2[i] = p_v2s_Qs->GetBinContent(i+1) / 100.;
        v2s_Q2_err[i] = p_v2s_Qs->GetBinError(i+1) / 100.;

	v2p_QQ[i] = p_v2p_Qp->GetBinContent(i+1) / 100.;
        v2p_QQ_err[i] = p_v2p_Qp->GetBinError(i+1) / 100.;
        v2p_Q2[i] = p_v2p_Qs->GetBinContent(i+1) / 100.;
        v2p_Q2_err[i] = p_v2p_Qs->GetBinError(i+1) / 100.;

	dg_QQ[i] = (Parity_Q_os_R->GetBinContent(i+1) - Parity_Q_ss_R->GetBinContent(i+1)) / 100.;
	dg_QQ_err[i] = sqrt(pow(Parity_Q_ss_R->GetBinError(i+1),2) + pow(Parity_Q_os_R->GetBinError(i+1),2))  / 100.;
	dg_Q2[i] = (Parity_Q_os_Rs->GetBinContent(i+1) - Parity_Q_ss_Rs->GetBinContent(i+1)) / 100.;
	dg_Q2_err[i] = sqrt(pow(Parity_Q_ss_Rs->GetBinError(i+1),2) + pow(Parity_Q_os_Rs->GetBinError(i+1),2))  / 100.;

	float pos = p_a1_Q2_p->GetBinContent(i+1);
	float pos_err = p_a1_Q2_p->GetBinError(i+1);
        float neg = p_a1_Q2_n->GetBinContent(i+1);
        float neg_err = p_a1_Q2_n->GetBinError(i+1);

        a1_Q2[i] = 0.5*(pos*pos + neg*neg) - pos*neg;
        a1_Q2_err[i] = sqrt(pos*pos*pos_err*pos_err + neg*neg*neg_err*neg_err + pos*pos*neg_err*neg_err + neg*neg*pos_err*pos_err);

	pos = p_a1_QQ_p->GetBinContent(i+1);
	pos_err = p_a1_QQ_p->GetBinError(i+1);
	neg = p_a1_QQ_n->GetBinContent(i+1);
	neg_err = p_a1_QQ_n->GetBinError(i+1);

	a1_QQ[i] = 0.5*(pos*pos + neg*neg) - pos*neg;
        a1_QQ_err[i] = sqrt(pos*pos*pos_err*pos_err + neg*neg*neg_err*neg_err + pos*pos*neg_err*neg_err + neg*neg*pos_err*pos_err);
}

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

TCanvas* can3 = new TCanvas("Flow_v", "Flow_v", 40,40,740,580);
 can3-> SetTopMargin(0.06);
 can3-> SetRightMargin(0.01);
 can3-> SetBottomMargin(0.17);
 can3-> SetLeftMargin(0.11);
 can3->Draw();

  gPad->SetGridx(0);
  gPad->SetGridy(0);
   gPad->SetTickx();
   gPad->SetTicky();
  TString* histGraphName = new TString("Fl");
  TH1F* histGraph = new TH1F(histGraphName->Data(), "", 24, 0, 0.13);
  histGraph->SetMaximum(0.00059);
  histGraph->SetMinimum(-0.00003);
  histGraph->SetLineColor(kBlack);
  histGraph->GetYaxis()->SetTitleOffset(0.67);
  histGraph->GetYaxis()->SetTitleSize(0.075);
  histGraph->GetXaxis()->SetTitleSize(0.08);
  histGraph->GetXaxis()->SetTitleOffset(0.90);
//  histGraph->GetYaxis()->CenterTitle();
  histGraph->GetXaxis()->SetTitle("single v_{2}{single q^{2}_{2}}");
//  histGraph->GetXaxis()->CenterTitle();
  histGraph->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  histGraph->GetXaxis()->SetNdivisions(6);
  histGraph->GetYaxis()->SetNdivisions(605);
  double lsize=histGraph->GetLabelSize();
  histGraph->GetYaxis()->SetLabelSize(lsize*1.2);
  histGraph->GetXaxis()->SetLabelSize(lsize*1.2);
  histGraph->Draw();

TGraphErrors* a1_v2s_Qs = new TGraphErrors(N_bins,v2s_Q2,a1_Q2,v2s_Q2_err,a1_Q2_err);
  a1_v2s_Qs->SetMarkerStyle(kOpenDiamond);
  a1_v2s_Qs->SetMarkerSize(1.5);
  a1_v2s_Qs->SetMarkerColor(2);
  a1_v2s_Qs->SetLineColor(2);
  a1_v2s_Qs->SetFillColor(2);
  a1_v2s_Qs->SetLineStyle(1);
  a1_v2s_Qs->SetLineWidth(2);
//  a1_v2s_Qs->Draw("pe1");

//TF1 *fit_a1_v2s_Qs = new TF1("fit_a1_v2s_Qs","[0]+[1]*x+[2]*x*x",0,0.13);
TF1 *fit_a1_v2s_Qs = new TF1("fit_a1_v2s_Qs","[0]+[1]*x",0,0.13);
fit_a1_v2s_Qs->SetParameters(0,0.01,0);
fit_a1_v2s_Qs->SetLineStyle(4);
fit_a1_v2s_Qs->SetLineColor(2);
a1_v2s_Qs->Fit("fit_a1_v2s_Qs","0E","",-0.01,0.13);
//fit_a1_v2s_Qs->Draw("same");

   TLatex *tex = new TLatex(0.04,0.000025,Form("intercept(2#font[12]{a}_{1}^{2}) = (%3.2f#pm%3.2f)#times10^{-5}",100000*fit_a1_v2s_Qs->GetParameter(0),100000*fit_a1_v2s_Qs->GetParError(0)));
   tex->SetTextSize(0.035);
   tex->SetTextColor(2);
//   tex->Draw();

TGraphErrors* g112_v2s_Qs = new TGraphErrors(N_bins,v2s_Q2,dg_Q2,v2s_Q2_err,dg_Q2_err);
  g112_v2s_Qs->SetMarkerStyle(kOpenCircle);
  g112_v2s_Qs->SetMarkerSize(1.5);
  g112_v2s_Qs->SetMarkerColor(4);
  g112_v2s_Qs->SetLineColor(4);
  g112_v2s_Qs->SetFillColor(4);
  g112_v2s_Qs->SetLineStyle(1);
  g112_v2s_Qs->SetLineWidth(2);
  g112_v2s_Qs->Draw("pe1");

//TF1 *fit_dg_v2s_Qs = new TF1("fit_dg_v2s_Qs","[0]+[1]*x+[2]*x*x",0,0.13);
TF1 *fit_dg_v2s_Qs = new TF1("fit_dg_v2s_Qs","[0]+[1]*x",0,0.13);
fit_dg_v2s_Qs->SetParameters(0,0.01,0);
fit_dg_v2s_Qs->SetLineStyle(4);
fit_dg_v2s_Qs->SetLineColor(4);
g112_v2s_Qs->Fit("fit_dg_v2s_Qs","0E","",-0.01,0.13);
fit_dg_v2s_Qs->Draw("same");

   tex = new TLatex(0.003,0.00035,Form("y-int(#Delta#gamma_{112}) = (%3.2f#pm%3.2f)#times10^{-5}",100000*fit_dg_v2s_Qs->GetParameter(0),100000*fit_dg_v2s_Qs->GetParError(0)));
   tex->SetTextSize(0.065);
   tex->SetTextColor(4);
   tex->Draw();

   tex = new TLatex(0.003,0.00052,"30 - 40% Au+Au 200 GeV (AVFD)");
   tex->SetTextSize(0.08);
   tex->SetTextColor(1);
   tex->Draw();

   tex = new TLatex(0.003,0.00045,"n_{5}/s = 0");
   tex->SetTextSize(0.08);
   tex->SetTextColor(1);
   tex->Draw();

   tex = new TLatex(0.118,0.00002,"(a)");
   tex->SetTextSize(0.08);
   tex->SetTextColor(1);
   tex->Draw();

TCanvas* can1 = new TCanvas("Flow_v1", "Flow_v1", 40,40,740,580);
 can1-> SetTopMargin(0.06);
 can1-> SetRightMargin(0.01);
 can1-> SetBottomMargin(0.17);
 can1-> SetLeftMargin(0.11);
 can1->Draw();

  gPad->SetGridx(0);
  gPad->SetGridy(0);
   gPad->SetTickx();
   gPad->SetTicky();
  TString* histGraphNam = new TString("Flo");
  TH1F* histGrap = new TH1F(histGraphNam->Data(), "", 24, 0, 0.13);
  histGrap->SetMaximum(0.00059);
  histGrap->SetMinimum(-0.00003);
  histGrap->SetLineColor(kBlack);
  histGrap->GetYaxis()->SetTitleOffset(0.67);
  histGrap->GetYaxis()->SetTitleSize(0.075);
  histGrap->GetXaxis()->SetTitleSize(0.08);
  histGrap->GetXaxis()->SetTitleOffset(0.90);
//  histGrap->GetYaxis()->CenterTitle();
  histGrap->GetXaxis()->SetTitle("pair v_{2}{single q^{2}_{2}}");
//  histGrap->GetXaxis()->CenterTitle();
  histGrap->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  histGrap->GetXaxis()->SetNdivisions(6);
  histGrap->GetYaxis()->SetNdivisions(605);
  lsize=histGrap->GetLabelSize();
  histGrap->GetYaxis()->SetLabelSize(lsize*1.2);
  histGrap->GetXaxis()->SetLabelSize(lsize*1.2);
  histGrap->Draw();

TGraphErrors* a1_v2p_Qs = new TGraphErrors(N_bins,v2p_Q2,a1_Q2,v2p_Q2_err,a1_Q2_err);
  a1_v2p_Qs->SetMarkerStyle(kOpenDiamond);
  a1_v2p_Qs->SetMarkerSize(1.5);
  a1_v2p_Qs->SetMarkerColor(2);
  a1_v2p_Qs->SetLineColor(2);
  a1_v2p_Qs->SetFillColor(2);
  a1_v2p_Qs->SetLineStyle(1);
  a1_v2p_Qs->SetLineWidth(2);
//  a1_v2p_Qs->Draw("pe1");

//TF1 *fit_a1_v2p_Qs = new TF1("fit_a1_v2p_Qs","[0]+[1]*x+[2]*x*x",0,0.13);
TF1 *fit_a1_v2p_Qs = new TF1("fit_a1_v2p_Qs","[0]+[1]*x",0,0.13);
fit_a1_v2p_Qs->SetParameters(0,0.01,0);
fit_a1_v2p_Qs->SetLineStyle(4);
fit_a1_v2p_Qs->SetLineColor(2);
a1_v2p_Qs->Fit("fit_a1_v2p_Qs","0E","",-0.01,0.13);
//fit_a1_v2p_Qs->Draw("same");

   tex = new TLatex(0.04,0.000025,Form("intercept(2#font[12]{a}_{1}^{2}) = (%3.2f#pm%3.2f)#times10^{-5}",100000*fit_a1_v2p_Qs->GetParameter(0),100000*fit_a1_v2p_Qs->GetParError(0)));
   tex->SetTextSize(0.035);
   tex->SetTextColor(2);
//   tex->Draw();

TGraphErrors* g112_v2p_Qs = new TGraphErrors(N_bins,v2p_Q2,dg_Q2,v2p_Q2_err,dg_Q2_err);
  g112_v2p_Qs->SetMarkerStyle(kOpenCircle);
  g112_v2p_Qs->SetMarkerSize(1.5);
  g112_v2p_Qs->SetMarkerColor(4);
  g112_v2p_Qs->SetLineColor(4);
  g112_v2p_Qs->SetFillColor(4);
  g112_v2p_Qs->SetLineStyle(1);
  g112_v2p_Qs->SetLineWidth(2);
  g112_v2p_Qs->Draw("pe1");

//TF1 *fit_dg_v2p_Qs = new TF1("fit_dg_v2p_Qs","[0]+[1]*x+[2]*x*x",0,0.13);
TF1 *fit_dg_v2p_Qs = new TF1("fit_dg_v2p_Qs","[0]+[1]*x",0,0.13);
fit_dg_v2p_Qs->SetParameters(0,0.01,0);
fit_dg_v2p_Qs->SetLineStyle(4);
fit_dg_v2p_Qs->SetLineColor(4);
g112_v2p_Qs->Fit("fit_dg_v2p_Qs","0E","",-0.01,0.13);
fit_dg_v2p_Qs->Draw("same");

   tex = new TLatex(0.003,0.00045,Form("y-int(#Delta#gamma_{112}) = (%3.2f#pm%3.2f)#times10^{-5}",100000*fit_dg_v2p_Qs->GetParameter(0),100000*fit_dg_v2p_Qs->GetParError(0)));
   tex->SetTextSize(0.065);
   tex->SetTextColor(4);
   tex->Draw();

   tex = new TLatex(0.003,0.00052,"30 - 40% Au+Au 200 GeV (AVFD)");
   tex->SetTextSize(0.08);
   tex->SetTextColor(1);
//   tex->Draw();

   tex = new TLatex(0.003,0.00042,"n_{5}/s = 0");
   tex->SetTextSize(0.08);
   tex->SetTextColor(1);
//   tex->Draw();

   tex = new TLatex(0.118,0.00002,"(d)");
   tex->SetTextSize(0.08);
   tex->SetTextColor(1);
   tex->Draw();

TCanvas* can2 = new TCanvas("Flow_v2", "Flow_v2", 40,40,740,580);
 can2-> SetTopMargin(0.06);
 can2-> SetRightMargin(0.01);
 can2-> SetBottomMargin(0.17);
 can2-> SetLeftMargin(0.11);
 can2->Draw();

  gPad->SetGridx(0);
  gPad->SetGridy(0);
   gPad->SetTickx();
   gPad->SetTicky();
  TString* histGraphNa = new TString("Flow");
  TH1F* histGra = new TH1F(histGraphNa->Data(), "", 24, 0, 0.13);
  histGra->SetMaximum(0.00059);
  histGra->SetMinimum(-0.00003);
  histGra->SetLineColor(kBlack);
  histGra->GetYaxis()->SetTitleOffset(0.67);
  histGra->GetYaxis()->SetTitleSize(0.075);
  histGra->GetXaxis()->SetTitleSize(0.08);
  histGra->GetXaxis()->SetTitleOffset(0.90);
//  histGra->GetYaxis()->CenterTitle();
  histGra->GetXaxis()->SetTitle("pair v_{2}{pair q^{2}_{2}}");
//  histGra->GetXaxis()->CenterTitle();
  histGra->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  histGra->GetXaxis()->SetNdivisions(6);
  histGra->GetYaxis()->SetNdivisions(605);
  lsize=histGra->GetLabelSize();
  histGra->GetYaxis()->SetLabelSize(lsize*1.2);
  histGra->GetXaxis()->SetLabelSize(lsize*1.2);
  histGra->Draw();

TGraphErrors* a1_v2p_Qp = new TGraphErrors(N_bins,v2p_QQ,a1_QQ,v2p_QQ_err,a1_QQ_err);
  a1_v2p_Qp->SetMarkerStyle(kOpenDiamond);
  a1_v2p_Qp->SetMarkerSize(1.5);
  a1_v2p_Qp->SetMarkerColor(2);
  a1_v2p_Qp->SetLineColor(2);
  a1_v2p_Qp->SetFillColor(2);
  a1_v2p_Qp->SetLineStyle(1);
  a1_v2p_Qp->SetLineWidth(2);
//  a1_v2p_Qp->Draw("pe1");

//TF1 *fit_a1_v2p_Qp = new TF1("fit_a1_v2p_Qp","[0]+[1]*x+[2]*x*x",0,0.112);
TF1 *fit_a1_v2p_Qp = new TF1("fit_a1_v2p_Qp","[0]+[1]*x",0,0.112);
fit_a1_v2p_Qp->SetParameters(0,0.01,0);
fit_a1_v2p_Qp->SetLineStyle(4);
fit_a1_v2p_Qp->SetLineColor(2);
a1_v2p_Qp->Fit("fit_a1_v2p_Qp","0E","",-0.01,0.13);
//fit_a1_v2p_Qp->Draw("same");

   tex = new TLatex(0.04,0.000025,Form("intercept(2#font[12]{a}_{1}^{2}) = (%3.2f#pm%3.2f)#times10^{-5}",100000*fit_a1_v2p_Qp->GetParameter(0),100000*fit_a1_v2p_Qp->GetParError(0)));
   tex->SetTextSize(0.035);
   tex->SetTextColor(2);
//   tex->Draw();

TGraphErrors* g112_v2p_Qp = new TGraphErrors(N_bins,v2p_QQ,dg_QQ,v2p_QQ_err,dg_QQ_err);
  g112_v2p_Qp->SetMarkerStyle(kOpenCircle);
  g112_v2p_Qp->SetMarkerSize(1.5);
  g112_v2p_Qp->SetMarkerColor(4);
  g112_v2p_Qp->SetLineColor(4);
  g112_v2p_Qp->SetFillColor(4);
  g112_v2p_Qp->SetLineStyle(1);
  g112_v2p_Qp->SetLineWidth(2);
  g112_v2p_Qp->Draw("pe1");

//TF1 *fit_dg_v2p_Qp = new TF1("fit_dg_v2p_Qp","[0]+[1]*x+[2]*x*x",0,0.112);
TF1 *fit_dg_v2p_Qp = new TF1("fit_dg_v2p_Qp","[0]+[1]*x",0,0.112);
fit_dg_v2p_Qp->SetParameters(0,0.01,0);
fit_dg_v2p_Qp->SetLineStyle(4);
fit_dg_v2p_Qp->SetLineColor(4);
g112_v2p_Qp->Fit("fit_dg_v2p_Qp","0E","",-0.01,0.13);
fit_dg_v2p_Qp->Draw("same");

   tex = new TLatex(0.003,0.00045,Form("y-int(#Delta#gamma_{112}) = (%3.2f#pm%3.2f)#times10^{-5}",100000*fit_dg_v2p_Qp->GetParameter(0),100000*fit_dg_v2p_Qp->GetParError(0)));
   tex->SetTextSize(0.065);
   tex->SetTextColor(4);
   tex->Draw();

   tex = new TLatex(0.003,0.00052,"30 - 40% Au+Au 200 GeV (AVFD)");
   tex->SetTextSize(0.08);
   tex->SetTextColor(1);
//   tex->Draw();

   tex = new TLatex(0.003,0.00042,"n_{5}/s = 0");
   tex->SetTextSize(0.08);
   tex->SetTextColor(1);
//   tex->Draw();

   tex = new TLatex(0.118,0.00002,"(b)");
   tex->SetTextSize(0.08);
   tex->SetTextColor(1);
   tex->Draw();

TCanvas* can4 = new TCanvas("Flow_v4", "Flow_v4", 40,40,740,580);
 can4-> SetTopMargin(0.06);
 can4-> SetRightMargin(0.01);
 can4-> SetBottomMargin(0.17);
 can4-> SetLeftMargin(0.11);
 can4->Draw();

  gPad->SetGridx(0);
  gPad->SetGridy(0);
   gPad->SetTickx();
   gPad->SetTicky();
  TString* histGraphN = new TString("Flow4");
  TH1F* histGr = new TH1F(histGraphN->Data(), "", 24, 0, 0.13);
  histGr->SetMaximum(0.00059);
  histGr->SetMinimum(-0.00003);
  histGr->SetLineColor(kBlack);
  histGr->GetYaxis()->SetTitleOffset(0.67);
  histGr->GetYaxis()->SetTitleSize(0.075);
  histGr->GetXaxis()->SetTitleSize(0.08);
  histGr->GetXaxis()->SetTitleOffset(0.90);
//  histGr->GetYaxis()->CenterTitle();
  histGr->GetXaxis()->SetTitle("single v_{2}{pair q^{2}_{2}}");
//  histGr->GetXaxis()->CenterTitle();
  histGr->GetYaxis()->SetTitle("#Delta#gamma_{112}");
  histGr->GetXaxis()->SetNdivisions(6);
  histGr->GetYaxis()->SetNdivisions(605);
  lsize=histGr->GetLabelSize();
  histGr->GetYaxis()->SetLabelSize(lsize*1.2);
  histGr->GetXaxis()->SetLabelSize(lsize*1.2);
  histGr->Draw();

TGraphErrors* a1_v2s_Qp = new TGraphErrors(N_bins,v2s_QQ,a1_QQ,v2s_QQ_err,a1_QQ_err);
  a1_v2s_Qp->SetMarkerStyle(kOpenDiamond);
  a1_v2s_Qp->SetMarkerSize(1.5);
  a1_v2s_Qp->SetMarkerColor(2);
  a1_v2s_Qp->SetLineColor(2);
  a1_v2s_Qp->SetFillColor(2);
  a1_v2s_Qp->SetLineStyle(1);
  a1_v2s_Qp->SetLineWidth(2);
//  a1_v2s_Qp->Draw("pe1");

//TF1 *fit_a1_v2s_Qp = new TF1("fit_a1_v2s_Qp","[0]+[1]*x+[2]*x*x",0,0.12);
TF1 *fit_a1_v2s_Qp = new TF1("fit_a1_v2s_Qp","[0]+[1]*x",0,0.12);
fit_a1_v2s_Qp->SetParameters(0,0.01,0);
fit_a1_v2s_Qp->SetLineStyle(4);
fit_a1_v2s_Qp->SetLineColor(2);
a1_v2s_Qp->Fit("fit_a1_v2s_Qp","0E","",-0.01,0.13);
//fit_a1_v2s_Qp->Draw("same");

   tex = new TLatex(0.04,0.000025,Form("intercept(2#font[12]{a}_{1}^{2}) = (%3.2f#pm%3.2f)#times10^{-5}",100000*fit_a1_v2s_Qp->GetParameter(0),100000*fit_a1_v2s_Qp->GetParError(0)));
   tex->SetTextSize(0.035);
   tex->SetTextColor(2);
//   tex->Draw();

TGraphErrors* g112_v2s_Qp = new TGraphErrors(N_bins,v2s_QQ,dg_QQ,v2s_QQ_err,dg_QQ_err);
  g112_v2s_Qp->SetMarkerStyle(kOpenCircle);
  g112_v2s_Qp->SetMarkerSize(1.5);
  g112_v2s_Qp->SetMarkerColor(4);
  g112_v2s_Qp->SetLineColor(4);
  g112_v2s_Qp->SetFillColor(4);
  g112_v2s_Qp->SetLineStyle(1);
  g112_v2s_Qp->SetLineWidth(2);
  g112_v2s_Qp->Draw("pe1");

//TF1 *fit_dg_v2s_Qp = new TF1("fit_dg_v2s_Qp","[0]+[1]*x+[2]*x*x",0,0.12);
TF1 *fit_dg_v2s_Qp = new TF1("fit_dg_v2s_Qp","[0]+[1]*x",0,0.12);
fit_dg_v2s_Qp->SetParameters(0,0.01,0);
fit_dg_v2s_Qp->SetLineStyle(4);
fit_dg_v2s_Qp->SetLineColor(4);
g112_v2s_Qp->Fit("fit_dg_v2s_Qp","0E","",-0.01,0.13);
fit_dg_v2s_Qp->Draw("same");

   tex = new TLatex(0.003,0.00045,Form("y-int(#Delta#gamma_{112}) = (%3.2f#pm%3.2f)#times10^{-5}",100000*fit_dg_v2s_Qp->GetParameter(0),100000*fit_dg_v2s_Qp->GetParError(0)));
   tex->SetTextSize(0.065);
   tex->SetTextColor(4);
   tex->Draw();

   tex = new TLatex(0.003,0.00052,"30 - 40% Au+Au 200 GeV (AVFD)");
   tex->SetTextSize(0.08);
   tex->SetTextColor(1);
//   tex->Draw();

   tex = new TLatex(0.003,0.00042,"n_{5}/s = 0");
   tex->SetTextSize(0.08);
   tex->SetTextColor(1);
//   tex->Draw();

   tex = new TLatex(0.118,0.00002,"(c)");
   tex->SetTextSize(0.08);
   tex->SetTextColor(1);
   tex->Draw();
}
