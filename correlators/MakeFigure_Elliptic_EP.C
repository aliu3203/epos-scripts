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

// Event-plane resolution R(chi) = sqrt(pi)/2 chi exp(-chi^2/4) [I0(chi^2/4) + I1(chi^2/4)]
double ResEP(double chi) {
    double x = chi*chi/4.;
    return sqrt(TMath::Pi())/2. * chi * exp(-x) * (TMath::BesselI0(x) + TMath::BesselI1(x));
}

// invert R(chi) by bisection on chi in [0,4]
double ChiFromRes(double R) {
    double lo = 0., hi = 4.;
    for(int it = 0; it < 60; it++) {
        double mid = 0.5*(lo+hi);
        if(ResEP(mid) < R) lo = mid; else hi = mid;
    }
    return 0.5*(lo+hi);
}

// Turn a TProfile into a TGraphErrors scaled by 1/R, skipping empty bins
TGraphErrors* ProfileToGraphEP(TProfile* p, double R, float xshift = 0.) {
    TGraphErrors* g = new TGraphErrors();
    int ip = 0;
    for(int i = 1; i <= p->GetNbinsX(); i++) {
        if(p->GetBinEntries(i) < 1) continue;
        g->SetPoint(ip, p->GetBinCenter(i) + xshift, p->GetBinContent(i)/R);
        g->SetPointError(ip, 0, p->GetBinError(i)/R);
        ip++;
    }
    return g;
}

// Styled canvas + axis frame
TCanvas* MakeFrameEP(const char* cname, float xlo, float xhi,
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
    TH1F* frame = new TH1F(Form("frame_%s",cname), "", 24, xlo, xhi);
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

void SetGraphStyleEP(TGraphErrors* g, int marker, int color) {
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.2);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
}

void MakeFigure_Elliptic_EP(int cen=5, int job=0) {

char fname[200];
sprintf(fname,"cen%d.v2pT_EP_job%d.root",cen,job);
TFile *f = new TFile(fname);

TProfile* pRes = (TProfile*)f->Get("pRes");

//---- event-plane resolutions ----
double cosAB     = pRes->GetBinContent(1);
double cosAB_err = pRes->GetBinError(1);
double R_sub = sqrt(cosAB);
double R_sub_err = cosAB_err/(2.*R_sub);

double chi_sub = ChiFromRes(R_sub);
double R_full = ResEP(sqrt(2.)*chi_sub);

cout<<endl;
cout<<"<cos2(PsiA-PsiB)>  = "<<cosAB<<" +/- "<<cosAB_err<<endl;
cout<<"R_sub  = sqrt(<cos2(PsiA-PsiB)>) = "<<R_sub<<" +/- "<<R_sub_err<<endl;
cout<<"chi_sub = "<<chi_sub<<",  R_full = R(sqrt(2)*chi_sub) = "<<R_full<<endl;
cout<<"cross-check with the true RP:  <cos2(PsiFull-PsiRP)> = "<<pRes->GetBinContent(2)
    <<",  <cos2(PsiA-PsiRP)> = "<<pRes->GetBinContent(3)
    <<",  <cos2(PsiB-PsiRP)> = "<<pRes->GetBinContent(4)<<endl;
cout<<endl;

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

//---- one canvas per species: v2{EP}, v2{sub-EP}, v2{true RP} ----
const int nSpec = 3;
const char* specKey[nSpec]   = {"pi","K","pr"};
const char* specLabel[nSpec] = {"#pi^{+}+#pi^{-}","K^{+}+K^{-}","p+#bar{p}"};
const char* specFile[nSpec]  = {"v2_pT_pion_EP.png","v2_pT_kaon_EP.png","v2_pT_proton_EP.png"};
const float specYmax[nSpec]  = {0.14, 0.14, 0.18};

for(int is = 0; is < nSpec; is++) {

    TProfile* p_EP  = (TProfile*)f->Get(Form("pV2Pt_%s_EP",specKey[is]));
    TProfile* p_sub = (TProfile*)f->Get(Form("pV2Pt_%s_sub",specKey[is]));
    TProfile* p_RP  = (TProfile*)f->Get(Form("pV2Pt_%s_RP",specKey[is]));

    TCanvas* can = MakeFrameEP(Form("Flow_v2pT_EP_%s",specKey[is]),
                               0, 3.0, -0.02, specYmax[is], "p_{T} (GeV/c)", "v_{2}");

    TGraphErrors* g_EP = ProfileToGraphEP(p_EP, R_full, -0.010);
      SetGraphStyleEP(g_EP, kFullCircle, 2);
      g_EP->Draw("pe1");

    TGraphErrors* g_sub = ProfileToGraphEP(p_sub, R_sub, 0.010);
      SetGraphStyleEP(g_sub, kOpenSquare, 4);
      g_sub->Draw("pe1");

    TGraphErrors* g_RP = ProfileToGraphEP(p_RP, 1., 0.);	//truth, no correction
      SetGraphStyleEP(g_RP, kOpenCross, kGray+2);
      g_RP->Draw("pe1");

    TLegend* leg = new TLegend(0.55,0.20,0.92,0.42);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.045);
      leg->AddEntry(g_EP, "v_{2}{EP} (full, /R_{full})","pe");
      leg->AddEntry(g_sub,"v_{2}{#eta-sub EP} (/R_{sub})","pe");
      leg->AddEntry(g_RP, "v_{2}{RP} (true plane)","pe");
      leg->Draw();

    TLatex* tex = new TLatex(0.15, specYmax[is]*0.89, "30 - 40% Au+Au 200 GeV (EPOS)");
      tex->SetTextSize(0.065);
      tex->SetTextColor(1);
      tex->Draw();

    tex = new TLatex(0.15, specYmax[is]*0.78, Form("%s, |#eta| < 1",specLabel[is]));
      tex->SetTextSize(0.055);
      tex->SetTextColor(1);
      tex->Draw();

    tex = new TLatex(0.15, specYmax[is]*0.68, Form("R_{full} = %4.3f, R_{sub} = %4.3f",R_full,R_sub));
      tex->SetTextSize(0.045);
      tex->SetTextColor(kGray+2);
      tex->Draw();

    can->SaveAs(specFile[is]);
}

//---- one canvas per method: pi / K / p compared ----
const int nMeth = 3;
const char* methKey[nMeth]   = {"EP","sub","RP"};
const char* methLabel[nMeth] = {"v_{2}{EP} (full, /R_{full})","v_{2}{#eta-sub EP} (/R_{sub})","v_{2}{RP} (true plane)"};
const char* methFile[nMeth]  = {"v2_pT_species_EPfull.png","v2_pT_species_EPsub.png","v2_pT_species_RP.png"};
const double methR[nMeth]    = {R_full, R_sub, 1.};
const float  methYmax[nMeth] = {0.22, 0.22, 0.18};

const int specMarker[nSpec] = {kFullCircle, kOpenSquare, kOpenDiamond};
const int specColor[nSpec]  = {2, 4, kGreen+2};

for(int im = 0; im < nMeth; im++) {

    TCanvas* can = MakeFrameEP(Form("Flow_v2pT_species_%s",methKey[im]),
                               0, 3.0, -0.02, methYmax[im], "p_{T} (GeV/c)", "v_{2}");

    TLegend* leg = new TLegend(0.72,0.20,0.92,0.42);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->SetTextSize(0.055);

    for(int is = 0; is < nSpec; is++) {
        TProfile* prof = (TProfile*)f->Get(Form("pV2Pt_%s_%s",specKey[is],methKey[im]));
        TGraphErrors* g = ProfileToGraphEP(prof, methR[im], 0.010*(is-1));
        SetGraphStyleEP(g, specMarker[is], specColor[is]);
        g->Draw("pe1");
        leg->AddEntry(g, specLabel[is], "pe");
    }
    leg->Draw();

    TLatex* tex = new TLatex(0.15, methYmax[im]*0.92, "30 - 40% Au+Au 200 GeV (EPOS)");
      tex->SetTextSize(0.065);
      tex->SetTextColor(1);
      tex->Draw();

    tex = new TLatex(0.15, methYmax[im]*0.84, Form("%s, |#eta| < 1",methLabel[im]));
      tex->SetTextSize(0.05);
      tex->SetTextColor(1);
      tex->Draw();

    can->SaveAs(methFile[im]);
}

//---- QA canvas: EP angle distributions ----
gStyle->SetOptTitle(1);	//re-enable pad titles for the QA panels
TH1D* h_PsiFull = (TH1D*)f->Get("Hist_PsiFull");
TH1D* h_PsiAB   = (TH1D*)f->Get("Hist_PsiAB");
TH1D* h_nEP     = (TH1D*)f->Get("Hist_nEP");
TH1D* h_b       = (TH1D*)f->Get("Hist_b");

TCanvas *c4 = new TCanvas("c4", "QA (event planes)", 1200, 800);
c4->Divide(2, 2);

c4->cd(1);
h_PsiFull->SetTitle("#Psi_{2}^{full} - #Psi_{RP}");
h_PsiFull->SetMinimum(0);
h_PsiFull->Draw();

c4->cd(2);
h_PsiAB->SetTitle("#Psi_{2}^{A} - #Psi_{2}^{B}");
h_PsiAB->SetMinimum(0);
h_PsiAB->Draw();

c4->cd(3);
h_nEP->SetTitle("N tracks in full Q-vector");
h_nEP->Draw();

c4->cd(4);
h_b->SetTitle("Impact Parameter (b)");
h_b->Draw();

c4->SaveAs("v2_QA_EP.png");

return;
}
