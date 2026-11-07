#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <functional>
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TSystem.h"

// Compare the outputs of Elliptic.C / Elliptic_KP.C / Elliptic_EP.C / Elliptic_ESS.C between two
// datasets (e.g. auau200-7 vs an independently generated EPOS sample with the same config).
//
// Each dataset dir holds the merged cen5.*_job0.root and, optionally, chunks/c*/cen5.*_job0.root
// (independent subsamples written by scripts/flow_compare.sh). With chunks, every error is the
// subsample spread: EPOS oversamples each hydro event 20x (nfreeze 20 -- the 20 events share b,
// Psi_RP and the initial state), so the per-particle TProfile errors are too small. Without
// chunks the TProfile/histogram errors are used.
//
// EP v2 is resolution-corrected and ESS is fitted per dataset (and per chunk) with the same
// recipe as MakeFigure_Elliptic_EP.C and MakeFigure_Elliptic_ESS.C.
//
// Writes compare_*.png and compare_summary.txt to outDir.
// Run: root -l -b -q 'Compare_Datasets.C+("refDir","newDir","outDir","ref label","new label")'

const int colRef = kBlack;
const int colNew = kRed+1;
const int nKey = 4;
const char* fileKey[nKey] = {"v2pT_pion","v2pT_kaon_proton","v2pT_EP","v2q2_ESS"};
FILE* gSum = 0;

struct Dataset {
    TString dir, label;
    TFile* tot[nKey];
    std::vector<TFile*> ch[nKey];	// chunk files, same chunk order for every key
};

struct BinStat { double x, y, e, eNaive, n; };

// --- EP resolution, identical to MakeFigure_Elliptic_EP.C ---
double ResEP_c(double chi) {
    double x = chi*chi/4.;
    return sqrt(TMath::Pi())/2. * chi * exp(-x) * (TMath::BesselI0(x) + TMath::BesselI1(x));
}
double ChiFromRes_c(double R) {
    double lo = 0., hi = 4.;
    for(int it = 0; it < 60; it++) {
        double mid = 0.5*(lo+hi);
        if(ResEP_c(mid) < R) lo = mid; else hi = mid;
    }
    return 0.5*(lo+hi);
}
double RsubOf(TFile* f)  { return sqrt(((TProfile*)f->Get("pRes"))->GetBinContent(1)); }
double RfullOf(TFile* f) { return ResEP_c(sqrt(2.)*ChiFromRes_c(RsubOf(f))); }
double One(TFile*)       { return 1.; }

// events in a file: Hist_b (Elliptic, KP, EP) or hq2s (ESS)
double EventsIn(TFile* f) {
    TH1* h = (TH1*)f->Get("Hist_b");
    if(!h) h = (TH1*)f->Get("hq2s");
    return h ? h->GetEntries() : 0;
}

bool OpenDataset(Dataset& d) {
    for(int k = 0; k < nKey; k++) {
        d.tot[k] = TFile::Open(Form("%s/cen5.%s_job0.root", d.dir.Data(), fileKey[k]));
        if(!d.tot[k]) { cout << "missing " << fileKey[k] << " in " << d.dir << endl; return false; }
    }
    std::vector<std::string> names;
    void* dp = gSystem->OpenDirectory(Form("%s/chunks", d.dir.Data()));
    if(dp) {
        while(const char* e = gSystem->GetDirEntry(dp)) if(e[0] == 'c') names.push_back(e);
        gSystem->FreeDirectory(dp);
    }
    std::sort(names.begin(), names.end());
    for(auto& n : names) {
        TFile* f[nKey]; bool ok = true;
        for(int k = 0; k < nKey && ok; k++) {
            TString path = Form("%s/chunks/%s/cen5.%s_job0.root", d.dir.Data(), n.c_str(), fileKey[k]);
            f[k] = gSystem->AccessPathName(path) ? 0 : TFile::Open(path);
            ok = f[k] != 0;
        }
        if(!ok) { cout << "skipping incomplete chunk " << n << endl; continue; }
        for(int k = 0; k < nKey; k++) d.ch[k].push_back(f[k]);
    }
    cout << d.label << ": " << d.ch[0].size() << " chunks" << endl;
    return true;
}

// subsample variance of a weighted mean: K/(K-1) * Sum_k (w_k/W)^2 (x_k - xbar)^2
double SubsampleErr(const std::vector<double>& x, const std::vector<double>& w, double xbar) {
    int K = x.size();
    if(K < 2) return 0;
    double W = 0; for(double v : w) W += v;
    if(W <= 0) return 0;
    double var = 0;
    for(int c = 0; c < K; c++) if(w[c] > 0) var += pow(w[c]/W, 2) * pow(x[c]-xbar, 2);
    return sqrt(var*K/(K-1.));
}

// per-bin value and error of profile `name` (file key k), scaled by 1/scale(file)
std::vector<BinStat> ProfStats(Dataset& d, int k, const char* name, double (*scale)(TFile*) = One) {
    TProfile* p = (TProfile*)d.tot[k]->Get(name);
    double s = scale(d.tot[k]);
    int K = d.ch[k].size();
    std::vector<TProfile*> pc; std::vector<double> sc;
    for(int c = 0; c < K; c++) { pc.push_back((TProfile*)d.ch[k][c]->Get(name)); sc.push_back(scale(d.ch[k][c])); }
    std::vector<BinStat> out;
    for(int i = 1; i <= p->GetNbinsX(); i++) {
        BinStat b;
        b.x = p->GetBinCenter(i); b.n = p->GetBinEntries(i);
        b.y = p->GetBinContent(i)/s; b.eNaive = p->GetBinError(i)/s; b.e = b.eNaive;
        if(K >= 2 && b.n > 0) {
            std::vector<double> xs, ws;
            for(int c = 0; c < K; c++) { ws.push_back(pc[c]->GetBinEntries(i)); xs.push_back(ws.back() > 0 ? pc[c]->GetBinContent(i)/sc[c] : 0); }
            b.e = SubsampleErr(xs, ws, b.y);
        }
        out.push_back(b);
    }
    return out;
}

// per-bin value and error of histogram `name`, normalised per event and bin width (perEvent)
// or to unit integral
std::vector<BinStat> HistStats(Dataset& d, int k, const char* name, bool perEvent, int rebin = 1) {
    auto norm = [&](TFile* f, TH1* h) { return perEvent ? EventsIn(f)*h->GetBinWidth(1) : h->Integral(); };
    TH1* h = (TH1*)d.tot[k]->Get(name)->Clone(Form("%s_%s_tot", name, d.label.Data()));
    h->Rebin(rebin);
    double s = norm(d.tot[k], h);
    int K = d.ch[k].size();
    std::vector<TH1*> hc; std::vector<double> sc, wev;
    for(int c = 0; c < K; c++) {
        TH1* hh = (TH1*)d.ch[k][c]->Get(name)->Clone(Form("%s_%s_%d", name, d.label.Data(), c));
        hh->Rebin(rebin);
        hc.push_back(hh); sc.push_back(norm(d.ch[k][c], hh)); wev.push_back(EventsIn(d.ch[k][c]));
    }
    std::vector<BinStat> out;
    for(int i = 1; i <= h->GetNbinsX(); i++) {
        BinStat b;
        b.x = h->GetBinCenter(i); b.n = h->GetBinContent(i);
        b.y = h->GetBinContent(i)/s; b.eNaive = h->GetBinError(i)/s; b.e = b.eNaive;
        if(K >= 2) {
            std::vector<double> xs;
            for(int c = 0; c < K; c++) xs.push_back(hc[c]->GetBinContent(i)/sc[c]);
            b.e = SubsampleErr(xs, wev, b.y);
        }
        out.push_back(b);
    }
    return out;
}

// a scalar observable evaluated on the merged file and on each chunk; fn returns {value, naive error}
typedef std::function<std::pair<double,double>(TFile*)> ScalarFn;
std::pair<double,double> Scalar(Dataset& d, int k, ScalarFn fn) {
    std::pair<double,double> t = fn(d.tot[k]);
    int K = d.ch[k].size();
    if(K < 2) return t;
    std::vector<double> xs, ws;
    for(int c = 0; c < K; c++) { xs.push_back(fn(d.ch[k][c]).first); ws.push_back(EventsIn(d.ch[k][c])); }
    return std::make_pair(t.first, SubsampleErr(xs, ws, t.first));
}

void Row(const char* name, std::pair<double,double> a, std::pair<double,double> b) {
    double den = sqrt(a.second*a.second + b.second*b.second);
    double sig = den > 0 ? (b.first-a.first)/den : 0;
    fprintf(gSum, "VAL  %-30s ref = %10.5f +/- %-9.5f new = %10.5f +/- %-9.5f  ratio = %.4f  diff = %+5.2f sigma\n",
            name, a.first, a.second, b.first, b.second, a.first != 0 ? b.first/a.first : 0, sig);
}

// bins usable in both: >= 20 entries per subsample on average
bool Usable(const BinStat& b, int K) { return b.n >= 20.*std::max(K,1) && b.e > 0; }

void ChiCompare(const char* name, const std::vector<BinStat>& a, const std::vector<BinStat>& b,
                int Ka, int Kb, double xlo, double xhi, bool countsBins = false) {
    double chi2 = 0; int ndf = 0; double maxPull = 0, maxPullX = 0;
    for(size_t i = 0; i < a.size(); i++) {
        if(a[i].x < xlo || a[i].x > xhi) continue;
        if(countsBins ? (a[i].e <= 0 || b[i].e <= 0 || a[i].n < 20 || b[i].n < 20) : (!Usable(a[i],Ka) || !Usable(b[i],Kb))) continue;
        double pull = (b[i].y - a[i].y)/sqrt(a[i].e*a[i].e + b[i].e*b[i].e);
        chi2 += pull*pull; ndf++;
        if(fabs(pull) > fabs(maxPull)) { maxPull = pull; maxPullX = a[i].x; }
    }
    fprintf(gSum, "CHI2 %-30s chi2/ndf = %7.1f / %-3d  p = %.3f  max pull = %+5.2f at x = %.3f\n",
            name, chi2, ndf, ndf > 0 ? TMath::Prob(chi2, ndf) : 0, maxPull, maxPullX);
}

// median ratio of subsample to naive error over usable bins
double DesignFactor(const std::vector<BinStat>& a, int K, double xlo, double xhi) {
    std::vector<double> r;
    for(auto& b : a) if(b.x >= xlo && b.x <= xhi && Usable(b,K) && b.eNaive > 0) r.push_back(b.e/b.eNaive);
    if(r.empty()) return 0;
    std::sort(r.begin(), r.end());
    return r[r.size()/2];
}

TGraphErrors* StatGraph(const std::vector<BinStat>& s, int K, double xshift, bool countsBins = false) {
    TGraphErrors* g = new TGraphErrors();
    for(auto& b : s) {
        if(countsBins ? b.n < 20 : !Usable(b,K)) continue;
        int ip = g->GetN();
        g->SetPoint(ip, b.x + xshift, b.y);
        g->SetPointError(ip, 0, b.e);
    }
    return g;
}

// overlay pad (top) + new/ref ratio pad (bottom) in column `col` of `ncol`
void DrawOverlay(TCanvas* c, int col, int ncol, const std::vector<BinStat>& a, const std::vector<BinStat>& b,
                 int Ka, int Kb, double xlo, double xhi, double ylo, double yhi,
                 const char* xtitle, const char* ytitle, const char* id, const char* label, int marker,
                 bool legend, const char* refLabel, const char* newLabel, bool logy = false, bool countsBins = false) {
    double x0 = (double)col/ncol, x1 = (double)(col+1)/ncol;
    c->cd();
    TPad* top = new TPad(Form("top_%s",id), "", x0, 0.32, x1, 1.0);
    TPad* bot = new TPad(Form("bot_%s",id), "", x0, 0.0,  x1, 0.32);
    top->SetBottomMargin(0.0); top->SetTopMargin(0.06); top->SetLeftMargin(0.17); top->SetRightMargin(0.03);
    bot->SetTopMargin(0.0);    bot->SetBottomMargin(0.32); bot->SetLeftMargin(0.17); bot->SetRightMargin(0.03);
    top->SetTickx(); top->SetTicky(); bot->SetTickx(); bot->SetTicky();
    top->SetLogy(logy);
    top->Draw(); bot->Draw();

    top->cd();
    TH1F* fr = new TH1F(Form("fr_%s",id), "", 24, xlo, xhi);
    fr->SetMinimum(ylo); fr->SetMaximum(yhi);
    fr->GetYaxis()->SetTitle(ytitle);
    fr->GetYaxis()->SetTitleSize(0.065); fr->GetYaxis()->SetTitleOffset(1.20);
    fr->GetYaxis()->SetLabelSize(0.055); fr->GetYaxis()->SetNdivisions(505);
    fr->Draw();

    double dx = 0.006*(xhi-xlo);
    TGraphErrors* gA = StatGraph(a, Ka, -dx, countsBins);
    TGraphErrors* gB = StatGraph(b, Kb,  dx, countsBins);
    gA->SetMarkerStyle(marker);   gA->SetMarkerColor(colRef); gA->SetLineColor(colRef); gA->SetMarkerSize(0.9);
    gB->SetMarkerStyle(marker+4); gB->SetMarkerColor(colNew); gB->SetLineColor(colNew); gB->SetMarkerSize(0.9);
    gA->Draw("pe same");
    gB->Draw("pe same");

    TLatex* tex = new TLatex();
    tex->SetNDC(); tex->SetTextSize(0.065);
    tex->DrawLatex(0.22, 0.85, label);
    if(legend) {
        TLegend* leg = new TLegend(0.40, 0.05, 0.96, 0.22);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.050);
        leg->AddEntry(gA, refLabel, "pe");
        leg->AddEntry(gB, newLabel, "pe");
        leg->Draw();
    }

    bot->cd();
    TH1F* fr2 = new TH1F(Form("fr2_%s",id), "", 24, xlo, xhi);
    fr2->SetMinimum(0.85); fr2->SetMaximum(1.15);
    fr2->GetXaxis()->SetTitle(xtitle); fr2->GetYaxis()->SetTitle("new / ref");
    fr2->GetXaxis()->SetTitleSize(0.13); fr2->GetXaxis()->SetTitleOffset(1.0); fr2->GetXaxis()->SetLabelSize(0.11);
    fr2->GetYaxis()->SetTitleSize(0.11); fr2->GetYaxis()->SetTitleOffset(0.70); fr2->GetYaxis()->SetLabelSize(0.10);
    fr2->GetYaxis()->SetNdivisions(504);
    fr2->Draw();
    TLine* one = new TLine(xlo, 1., xhi, 1.); one->SetLineStyle(2); one->Draw();
    TGraphErrors* gR = new TGraphErrors();
    for(size_t i = 0; i < a.size(); i++) {
        if(a[i].x < xlo || a[i].x > xhi || fabs(a[i].y) < 1e-12) continue;
        if(countsBins ? (a[i].n < 20 || b[i].n < 20) : (!Usable(a[i],Ka) || !Usable(b[i],Kb))) continue;
        double r = b[i].y/a[i].y;
        int ip = gR->GetN();
        gR->SetPoint(ip, a[i].x, r);
        gR->SetPointError(ip, 0, fabs(r)*sqrt(pow(a[i].e/a[i].y,2) + pow(b[i].e/b[i].y,2)));
    }
    gR->SetMarkerStyle(marker); gR->SetMarkerColor(colNew); gR->SetLineColor(colNew); gR->SetMarkerSize(0.8);
    gR->Draw("pe same");
}

// entry-weighted profile average over [xlo,xhi] (pT-integrated v2) with naive error
std::pair<double,double> ProfAverage(TProfile* p, double xlo, double xhi) {
    double sw = 0, swy = 0, swe2 = 0;
    for(int i = 1; i <= p->GetNbinsX(); i++) {
        double x = p->GetBinCenter(i);
        if(x < xlo || x > xhi) continue;
        double w = p->GetBinEntries(i);
        sw += w; swy += w*p->GetBinContent(i); swe2 += pow(w*p->GetBinError(i),2);
    }
    return std::make_pair(sw > 0 ? swy/sw : 0, sw > 0 ? sqrt(swe2)/sw : 0);
}

// ESS linear fit, same recipe as MakeFigure_Elliptic_ESS.C: all bins with >= 20 events
TF1* FitESS(TProfile* p, const char* name) {
    TGraphErrors* g = new TGraphErrors();
    for(int i = 1; i <= p->GetNbinsX(); i++) {
        if(p->GetBinEntries(i) < 20) continue;
        int ip = g->GetN();
        g->SetPoint(ip, p->GetBinCenter(i), p->GetBinContent(i));
        g->SetPointError(ip, 0, p->GetBinError(i));
    }
    TF1* f = new TF1(name, "[0]+[1]*x", 0, 4);
    g->Fit(f, "Q0");
    delete g;
    return f;
}

void Compare_Datasets(const char* refDir = ".", const char* newDir = ".", const char* outDir = ".",
                      const char* refLabel = "reference", const char* newLabel = "new") {

    gStyle->SetOptStat(0); gStyle->SetOptTitle(0); gStyle->SetOptFit(0);
    gStyle->SetFrameBorderMode(0); gStyle->SetCanvasColor(0); gStyle->SetPadColor(0);
    gStyle->SetPadBorderMode(0); gStyle->SetCanvasBorderMode(0); gStyle->SetEndErrorSize(0);
    TH1::AddDirectory(kFALSE);

    Dataset R, N;
    R.dir = refDir; R.label = "ref";
    N.dir = newDir; N.label = "new";
    if(!OpenDataset(R) || !OpenDataset(N)) return;
    int KR = R.ch[0].size(), KN = N.ch[0].size();

    gSum = fopen(Form("%s/compare_summary.txt", outDir), "w");
    fprintf(gSum, "# ref = %s  (%s), %d subsamples\n# new = %s  (%s), %d subsamples\n",
            refLabel, refDir, KR, newLabel, newDir, KN);
    fprintf(gSum, "# errors: %s\n", (KR >= 2 && KN >= 2) ? "subsample spread (independent chunks)" : "TProfile/histogram (no chunks for at least one dataset)");
    fprintf(gSum, "NEV  ref = %.0f  new = %.0f   (events passing 8.4<b<9.2)\n", EventsIn(R.tot[0]), EventsIn(N.tot[0]));

    //================= scalars =================
    auto histMean = [](const char* h) { return ScalarFn([h](TFile* f) {
        TH1* x = (TH1*)f->Get(h); return std::make_pair(x->GetMean(), x->GetMeanError()); }); };
    auto histRMS  = [](const char* h) { return ScalarFn([h](TFile* f) {
        TH1* x = (TH1*)f->Get(h); return std::make_pair(x->GetRMS(), x->GetRMSError()); }); };
    auto yield    = [](const char* h) { return ScalarFn([h](TFile* f) {
        TH1* x = (TH1*)f->Get(h); double n = EventsIn(f), s = x->Integral(); return std::make_pair(s/n, sqrt(s)/n); }); };
    auto v2int    = [](const char* p) { return ScalarFn([p](TFile* f) {
        return ProfAverage((TProfile*)f->Get(p), 0.2, 2.0); }); };
    auto resBin   = [](int i) { return ScalarFn([i](TFile* f) {
        TProfile* p = (TProfile*)f->Get("pRes"); return std::make_pair(p->GetBinContent(i), p->GetBinError(i)); }); };
    auto essConst = [](int i) { return ScalarFn([i](TFile* f) {
        return std::make_pair(((TH1*)f->Get("hConst"))->GetBinContent(i), 0.); }); };
    auto essFit   = [](const char* p, int par) { return ScalarFn([p,par](TFile* f) {
        TF1* fit = FitESS((TProfile*)f->Get(p), "fit_tmp");
        std::pair<double,double> r(fit->GetParameter(par), fit->GetParError(par)); delete fit; return r; }); };

    Row("<b> [fm]",                 Scalar(R,0,histMean("Hist_b")),  Scalar(N,0,histMean("Hist_b")));
    Row("rms(b) [fm]",              Scalar(R,0,histRMS("Hist_b")),   Scalar(N,0,histRMS("Hist_b")));
    Row("dN/ev pi+- |eta|<1",       Scalar(R,0,yield("Hist_Pt_pi")), Scalar(N,0,yield("Hist_Pt_pi")));
    Row("dN/ev K+- |eta|<1",        Scalar(R,1,yield("Hist_Pt_K")),  Scalar(N,1,yield("Hist_Pt_K")));
    Row("dN/ev p+pbar |eta|<1",     Scalar(R,1,yield("Hist_Pt_pr")), Scalar(N,1,yield("Hist_Pt_pr")));
    Row("<pT> pi+- [GeV]",          Scalar(R,0,histMean("Hist_Pt_pi")), Scalar(N,0,histMean("Hist_Pt_pi")));
    Row("<pT> K+- [GeV]",           Scalar(R,1,histMean("Hist_Pt_K")),  Scalar(N,1,histMean("Hist_Pt_K")));
    Row("<pT> p+pbar [GeV]",        Scalar(R,1,histMean("Hist_Pt_pr")), Scalar(N,1,histMean("Hist_Pt_pr")));
    Row("<N> EP Q-vector tracks",   Scalar(R,2,histMean("Hist_nEP")), Scalar(N,2,histMean("Hist_nEP")));
    Row("rms(N) EP Q-vector",       Scalar(R,2,histRMS("Hist_nEP")),  Scalar(N,2,histRMS("Hist_nEP")));
    Row("v2{RP} pi+- 0.2<pT<2",     Scalar(R,0,v2int("pV2Pt_pi")),  Scalar(N,0,v2int("pV2Pt_pi")));
    Row("v2{RP} K+- 0.2<pT<2",      Scalar(R,1,v2int("pV2Pt_K")),   Scalar(N,1,v2int("pV2Pt_K")));
    Row("v2{RP} p+pbar 0.2<pT<2",   Scalar(R,1,v2int("pV2Pt_pr")),  Scalar(N,1,v2int("pV2Pt_pr")));
    Row("v2{RP} pi0 0.2<pT<2",      Scalar(R,0,v2int("pV2Pt_pi0")), Scalar(N,0,v2int("pV2Pt_pi0")));
    Row("v2{PP} pi+- 0.2<pT<2",     Scalar(R,0,v2int("pV2Pt_pi_PP")), Scalar(N,0,v2int("pV2Pt_pi_PP")));
    const char* resLab[4] = {"<cos2(PsiA-PsiB)>","<cos2(PsiFull-PsiRP)>","<cos2(PsiA-PsiRP)>","<cos2(PsiB-PsiRP)>"};
    for(int i = 1; i <= 4; i++) Row(resLab[i-1], Scalar(R,2,resBin(i)), Scalar(N,2,resBin(i)));
    ScalarFn rsub  = [](TFile* f) { TProfile* p = (TProfile*)f->Get("pRes"); double r = sqrt(p->GetBinContent(1));
                                    return std::make_pair(r, p->GetBinError(1)/(2*r)); };
    ScalarFn rfull = [](TFile* f) { TProfile* p = (TProfile*)f->Get("pRes"); double r = sqrt(p->GetBinContent(1));
                                    double e = p->GetBinError(1)/(2*r); double R0 = ResEP_c(sqrt(2.)*ChiFromRes_c(r));
                                    return std::make_pair(R0, fabs(ResEP_c(sqrt(2.)*ChiFromRes_c(r+e)) - R0)); };
    Row("R_sub", Scalar(R,2,rsub), Scalar(N,2,rsub));
    Row("R_full", Scalar(R,2,rfull), Scalar(N,2,rfull));
    Row("v2{2} (ESS, pi+-)",        Scalar(R,3,essConst(1)), Scalar(N,3,essConst(1)));
    Row("v2pair{2} (ESS, pi+-)",    Scalar(R,3,essConst(2)), Scalar(N,3,essConst(2)));
    Row("<single q2^2>",            Scalar(R,3,histMean("hq2s")), Scalar(N,3,histMean("hq2s")));
    Row("<pair q2^2>",              Scalar(R,3,histMean("hq2p")), Scalar(N,3,histMean("hq2p")));

    const char* essKey[4] = {"pV2s_q2s","pV2p_q2s","pV2s_q2p","pV2p_q2p"};
    const char* essTag[4] = {"(a) single q_{2}^{2}, single v_{2}","(b) single q_{2}^{2}, pair v_{2}",
                             "(c) pair q_{2}^{2}, single v_{2}","(d) pair q_{2}^{2}, pair v_{2}"};
    const char* essLab[4] = {"ESS (a) s-q2 -> s-v2","ESS (b) s-q2 -> p-v2","ESS (c) p-q2 -> s-v2","ESS (d) p-q2 -> p-v2"};
    std::pair<double,double> essR[4][2], essN[4][2];
    for(int k = 0; k < 4; k++) {
        for(int par = 0; par < 2; par++) {
            essR[k][par] = Scalar(R,3,essFit(essKey[k],par));
            essN[k][par] = Scalar(N,3,essFit(essKey[k],par));
        }
        Row(Form("%s slope", essLab[k]),     essR[k][1], essN[k][1]);
        Row(Form("%s intercept", essLab[k]), essR[k][0], essN[k][0]);
    }

    //================= figures + bin-by-bin chi2 =================
    const char* sLab[3] = {"#pi^{+}+#pi^{-}","K^{+}+K^{-}","p+#bar{p}"};
    const char* sKey[3] = {"pi","K","pr"};
    const char* sShort[3] = {"pi","K","p"};
    int sFile[3] = {0,1,1};
    double sYmax[3] = {0.14, 0.14, 0.20};

    // 1) true-RP and participant-plane v2(pT)
    const char* plane[2] = {"RP","PP"};
    for(int ip = 0; ip < 2; ip++) {
        TCanvas* c = new TCanvas(Form("c_%s",plane[ip]), "", 1500, 620);
        for(int s = 0; s < 3; s++) {
            TString name = ip == 0 ? Form("pV2Pt_%s",sKey[s]) : Form("pV2Pt_%s_PP",sKey[s]);
            auto a = ProfStats(R, sFile[s], name);
            auto b = ProfStats(N, sFile[s], name);
            DrawOverlay(c, s, 3, a, b, KR, KN, 0, 3.0, -0.005, sYmax[s], "p_{T} (GeV/c)",
                        ip == 0 ? "v_{2}{true RP}" : "v_{2}{PP}", Form("%s%d",plane[ip],s), sLab[s], 20,
                        s == 0, refLabel, newLabel);
            ChiCompare(Form("v2(pT) %s %s", plane[ip], sShort[s]), a, b, KR, KN, 0.0, 3.0);
            if(ip == 0 && s == 0)
                fprintf(gSum, "DEFF v2(pT) RP pi: median subsample/naive error = ref %.2f, new %.2f\n",
                        DesignFactor(a, KR, 0, 3), DesignFactor(b, KN, 0, 3));
        }
        c->SaveAs(Form("%s/compare_v2pT_%s.png", outDir, plane[ip]));
    }

    // 2) v2(eta), true RP
    {
        TCanvas* c = new TCanvas("c_eta", "", 1500, 620);
        for(int s = 0; s < 3; s++) {
            TString name = Form("pV2Eta_%s", sKey[s]);
            auto a = ProfStats(R, sFile[s], name);
            auto b = ProfStats(N, sFile[s], name);
            DrawOverlay(c, s, 3, a, b, KR, KN, -6, 6, -0.005, 0.08, "#eta", "v_{2}{true RP}",
                        Form("eta%d",s), sLab[s], 20, s == 0, refLabel, newLabel);
            ChiCompare(Form("v2(eta) RP %s", sShort[s]), a, b, KR, KN, -6, 6);
        }
        c->SaveAs(Form("%s/compare_v2eta_RP.png", outDir));
    }

    // 3) resolution-corrected v2{EP} and v2{eta-sub EP}; each dataset/chunk uses its own resolution
    {
        const char* meth[2] = {"EP","sub"};
        for(int im = 0; im < 2; im++) {
            TCanvas* c = new TCanvas(Form("c_%s",meth[im]), "", 1500, 620);
            for(int s = 0; s < 3; s++) {
                TString name = Form("pV2Pt_%s_%s", sKey[s], meth[im]);
                auto a = ProfStats(R, 2, name, im == 0 ? RfullOf : RsubOf);
                auto b = ProfStats(N, 2, name, im == 0 ? RfullOf : RsubOf);
                DrawOverlay(c, s, 3, a, b, KR, KN, 0, 3.0, -0.005, sYmax[s]*1.4, "p_{T} (GeV/c)",
                            im == 0 ? "v_{2}{EP}" : "v_{2}{#eta-sub EP}", Form("%s%d",meth[im],s), sLab[s], 21,
                            s == 0, refLabel, newLabel);
                ChiCompare(Form("v2(pT) %s %s", meth[im], sShort[s]), a, b, KR, KN, 0.0, 3.0);
            }
            c->SaveAs(Form("%s/compare_v2pT_%s.png", outDir, im == 0 ? "EPfull" : "EPsub"));
        }
    }

    // 4) pi0 v2(pT) (decayed pi0, ist==1)
    {
        TCanvas* c = new TCanvas("c_pi0", "", 600, 620);
        auto a = ProfStats(R, 0, "pV2Pt_pi0");
        auto b = ProfStats(N, 0, "pV2Pt_pi0");
        DrawOverlay(c, 0, 1, a, b, KR, KN, 0, 3.0, -0.005, 0.14, "p_{T} (GeV/c)", "v_{2}{true RP}",
                    "pi0", "#pi^{0}", 20, true, refLabel, newLabel);
        ChiCompare("v2(pT) RP pi0", a, b, KR, KN, 0.0, 3.0);
        c->SaveAs(Form("%s/compare_v2pT_pi0.png", outDir));
    }

    // 5) ESS 2x2 panels
    {
        TCanvas* c = new TCanvas("c_ess", "", 1200, 1000);
        c->Divide(2, 2, 0.001, 0.001);
        for(int k = 0; k < 4; k++) {
            auto a = ProfStats(R, 3, essKey[k]);
            auto b = ProfStats(N, 3, essKey[k]);
            c->cd(k+1);
            gPad->SetTopMargin(0.04); gPad->SetRightMargin(0.03); gPad->SetBottomMargin(0.14); gPad->SetLeftMargin(0.15);
            gPad->SetTickx(); gPad->SetTicky();
            TH1F* fr = new TH1F(Form("fr_ess%d",k), "", 24, 0, 3.5);
            fr->SetMinimum(0); fr->SetMaximum(0.18);
            fr->GetXaxis()->SetTitle(k < 2 ? "single q_{2}^{2}" : "pair q_{2}^{2}");
            fr->GetYaxis()->SetTitle(k % 2 == 0 ? "single v_{2}" : "pair v_{2}");
            fr->GetXaxis()->SetTitleSize(0.06); fr->GetYaxis()->SetTitleSize(0.06);
            fr->GetYaxis()->SetTitleOffset(1.15); fr->Draw();
            TGraphErrors* gA = StatGraph(a, KR, -0.015);
            TGraphErrors* gB = StatGraph(b, KN,  0.015);
            gA->SetMarkerStyle(20); gA->SetMarkerColor(colRef); gA->SetLineColor(colRef);
            gB->SetMarkerStyle(24); gB->SetMarkerColor(colNew); gB->SetLineColor(colNew);
            gA->Draw("pe same"); gB->Draw("pe same");
            TF1* fA = FitESS((TProfile*)R.tot[3]->Get(essKey[k]), Form("fitR%d",k));
            TF1* fB = FitESS((TProfile*)N.tot[3]->Get(essKey[k]), Form("fitN%d",k));
            fA->SetLineColor(colRef); fA->SetLineStyle(2); fA->SetRange(0,3.5); fA->Draw("same");
            fB->SetLineColor(colNew); fB->SetLineStyle(2); fB->SetRange(0,3.5); fB->Draw("same");
            TLatex* t = new TLatex(); t->SetNDC(); t->SetTextSize(0.050);
            t->DrawLatex(0.19, 0.89, essTag[k]);
            t->SetTextSize(0.040); t->SetTextColor(colRef);
            t->DrawLatex(0.19, 0.81, Form("ref: slope %.4f #pm %.4f, int. %.4f #pm %.4f", essR[k][1].first, essR[k][1].second, essR[k][0].first, essR[k][0].second));
            t->SetTextColor(colNew);
            t->DrawLatex(0.19, 0.75, Form("new: slope %.4f #pm %.4f, int. %.4f #pm %.4f", essN[k][1].first, essN[k][1].second, essN[k][0].first, essN[k][0].second));
            if(k == 0) {
                TLegend* leg = new TLegend(0.45, 0.17, 0.95, 0.31);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.042);
                leg->AddEntry(gA, refLabel, "pe"); leg->AddEntry(gB, newLabel, "pe"); leg->Draw();
            }
            ChiCompare(Form("ESS %s", essLab[k]+4), a, b, KR, KN, 0, 3.5);
        }
        c->SaveAs(Form("%s/compare_ESS.png", outDir));
    }

    // 6) QA: b, pT spectra, EP multiplicity, q2 distributions (counts histograms)
    {
        TCanvas* c = new TCanvas("c_qa", "", 1500, 1300);
        const int nQA = 6;
        int qaFile[nQA]          = {0, 0, 1, 1, 2, 3};
        const char* qaName[nQA]  = {"Hist_b","Hist_Pt_pi","Hist_Pt_K","Hist_Pt_pr","Hist_nEP","hq2s"};
        bool qaPerEv[nQA]        = {false, true, true, true, false, false};
        int qaRebin[nQA]         = {1, 1, 1, 1, 2, 1};
        const char* qaX[nQA]     = {"b (fm)","p_{T} (GeV/c)","p_{T} (GeV/c)","p_{T} (GeV/c)","N_{trk} in EP Q-vector","single q_{2}^{2}"};
        const char* qaY[nQA]     = {"fraction of events","1/N_{ev} dN_{#pi}/dp_{T}","1/N_{ev} dN_{K}/dp_{T}","1/N_{ev} dN_{p}/dp_{T}","fraction of events","fraction of events"};
        const char* qaChi[nQA]   = {"b distribution","pT spectrum pi","pT spectrum K","pT spectrum p","N in EP Q-vector","single q2^2 distribution"};
        const char* qaLab[nQA]   = {"b","#pi^{+}+#pi^{-}","K^{+}+K^{-}","p+#bar{p}","N_{EP}","q_{2}^{2}"};
        double qaLo[nQA] = {8.2, 0, 0, 0, 40, 0}, qaHi[nQA] = {9.4, 3, 3, 3, 240, 6};
        double qaYlo[nQA] = {0, 0.1, 0.05, 0.03, 0, 1e-5}, qaYhi[nQA] = {0.8, 600, 40, 15, 0.45, 0.2};
        bool qaLog[nQA] = {false, true, true, true, false, true};
        // three columns x two rows of (overlay + ratio)
        for(int k = 0; k < nQA; k++) {
            auto a = HistStats(R, qaFile[k], qaName[k], qaPerEv[k], qaRebin[k]);
            auto b = HistStats(N, qaFile[k], qaName[k], qaPerEv[k], qaRebin[k]);
            ChiCompare(qaChi[k], a, b, KR, KN, qaLo[k], qaHi[k], true);
            c->cd();
            TPad* row = new TPad(Form("row%d",k), "", 0, k < 3 ? 0.5 : 0.0, 1, k < 3 ? 1.0 : 0.5);
            row->SetFillStyle(4000); row->Draw();
            row->cd();
            double x0 = (double)(k%3)/3, x1 = (double)(k%3+1)/3;
            TPad* top = new TPad(Form("qtop%d",k), "", x0, 0.32, x1, 1.0);
            TPad* bot = new TPad(Form("qbot%d",k), "", x0, 0.0,  x1, 0.32);
            top->SetBottomMargin(0.0); top->SetTopMargin(0.06); top->SetLeftMargin(0.17); top->SetRightMargin(0.03);
            bot->SetTopMargin(0.0);    bot->SetBottomMargin(0.32); bot->SetLeftMargin(0.17); bot->SetRightMargin(0.03);
            top->SetTickx(); top->SetTicky(); bot->SetTickx(); bot->SetTicky(); top->SetLogy(qaLog[k]);
            top->Draw(); bot->Draw();
            top->cd();
            TH1F* fr = new TH1F(Form("qfr%d",k), "", 24, qaLo[k], qaHi[k]);
            fr->SetMinimum(qaYlo[k]); fr->SetMaximum(qaYhi[k]);
            fr->GetYaxis()->SetTitle(qaY[k]); fr->GetYaxis()->SetTitleSize(0.065);
            fr->GetYaxis()->SetTitleOffset(1.20); fr->GetYaxis()->SetLabelSize(0.055); fr->Draw();
            TGraphErrors* gA = StatGraph(a, KR, 0, true);
            TGraphErrors* gB = StatGraph(b, KN, 0, true);
            gA->SetMarkerStyle(20); gA->SetMarkerColor(colRef); gA->SetLineColor(colRef); gA->SetMarkerSize(0.7);
            gB->SetMarkerStyle(24); gB->SetMarkerColor(colNew); gB->SetLineColor(colNew); gB->SetMarkerSize(0.7);
            gA->Draw("pe same"); gB->Draw("pe same");
            TLatex* t = new TLatex(); t->SetNDC(); t->SetTextSize(0.065); t->DrawLatex(0.70, 0.85, qaLab[k]);
            if(k == 0) {
                TLegend* leg = new TLegend(0.20, 0.72, 0.70, 0.90);
                leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.050);
                leg->AddEntry(gA, refLabel, "pe"); leg->AddEntry(gB, newLabel, "pe"); leg->Draw();
            }
            bot->cd();
            TH1F* fr2 = new TH1F(Form("qfr2%d",k), "", 24, qaLo[k], qaHi[k]);
            fr2->SetMinimum(0.8); fr2->SetMaximum(1.2);
            fr2->GetXaxis()->SetTitle(qaX[k]); fr2->GetYaxis()->SetTitle("new / ref");
            fr2->GetXaxis()->SetTitleSize(0.13); fr2->GetXaxis()->SetTitleOffset(1.0); fr2->GetXaxis()->SetLabelSize(0.11);
            fr2->GetYaxis()->SetTitleSize(0.11); fr2->GetYaxis()->SetTitleOffset(0.70); fr2->GetYaxis()->SetLabelSize(0.10);
            fr2->GetYaxis()->SetNdivisions(504); fr2->Draw();
            TLine* one = new TLine(qaLo[k], 1., qaHi[k], 1.); one->SetLineStyle(2); one->Draw();
            TGraphErrors* gR = new TGraphErrors();
            for(size_t i = 0; i < a.size(); i++) {
                if(a[i].x < qaLo[k] || a[i].x > qaHi[k] || a[i].n < 20 || b[i].n < 20 || a[i].y <= 0) continue;
                double r = b[i].y/a[i].y; int ip = gR->GetN();
                gR->SetPoint(ip, a[i].x, r);
                gR->SetPointError(ip, 0, r*sqrt(pow(a[i].e/a[i].y,2) + pow(b[i].e/b[i].y,2)));
            }
            gR->SetMarkerStyle(20); gR->SetMarkerColor(colNew); gR->SetLineColor(colNew); gR->SetMarkerSize(0.7);
            gR->Draw("pe same");
        }
        c->SaveAs(Form("%s/compare_QA.png", outDir));
    }

    fclose(gSum);
    gSystem->Exec(Form("cat %s/compare_summary.txt", outDir));
    return;
}
