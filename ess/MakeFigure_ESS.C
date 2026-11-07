// MakeFigure_ESS.C -- combine the per-chunk ESS.C outputs and do the ESS projection of arXiv:2307.14997
// (as in the authors' MakeFigure_Intercept_n5s0.C): in each event-shape bin compute Delta gamma112 and the
// elliptic-flow variable, fit Delta gamma112 vs v2 with a straight line in [-0.01, 0.13] and quote the y-intercept.
//
// Per shape bin (from the 8 profiles of ESS.C):
//   single v2, pair v2, resonance v2 = (N_OS v2_OS - N_SS v2_SS)/(N_OS - N_SS)            (Eq. 12)
//   Delta gamma112 = gamma_OS - gamma_SS                                                    (Eq. 3)
//   2a1^2 = (a1+^2 + a1-^2)/2 - a1+ a1-   (true-signal reference; EPOS has no CME -> ~0)    (Eq. 5)
// Central values come from the summed chunks. Errors come from the chunk-to-chunk spread (chunks are disjoint file sets,
// hence independent; per-term errors miss the nfreeze=20 correlations), as in correlators/Combine_Correlators.C:
//   var(x) = K/(K-1) Sum_k (w_k/W)^2 (x_k - x)^2 ,  w_k = events of chunk k in the bin.
// Intercepts are also given x (1 - <v2>)^2, the paper's first-order restoration to the inclusive signal.
//
// Figures (written to <base>): ESS_fig6.png  = Fig. 6/7 (four combinations of single/pair q2^2 and v2),
//                              ESS_fig5_9.png = Fig. 5 (resonance v2) and Fig. 9 (q2,B, linear and 2nd-order fits)
// Table: <base>/ESS_summary.txt
//
// Usage: root -l -b -q 'MakeFigure_ESS.C+("/media/Students/aliu/ess_all","day1 day2 day3 day4")'
using namespace std;
#include <cstdio>
#include <string>
#include <vector>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TMultiGraph.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TProfile2D.h"
#include "TStyle.h"
#include "TSystem.h"

static const int NQ = 5;      // derived quantities
static const char* qTitle[NQ] = {"single v_{2}", "pair v_{2}", "resonance v_{2}", "#Delta#gamma_{112}", "2a_{1}^{2}"};
static const char* kTitle[3]  = {"single q_{2}^{2}", "pair q_{2}^{2}", "q_{2,B}"};
static const char* kTag[3]    = {"q2s", "q2p", "q2B"};
static const double MINEV = 1000;    // minimum events in a shape bin to use it

struct Sample { TProfile2D* p[3]; TH1D* h[3]; };
static Sample tot;
static vector<Sample> chunks;
static double v2incl = 0;
static string summary;

// derived quantities from the 8 observables, summed over shape bins ylo..yhi
static bool Derive(TProfile2D* p, int ylo, int yhi, double* q){
    double m[9], n[9];
    for(int o = 1; o <= 8; o++){
        double sw = 0, swz = 0;
        for(int y = ylo; y <= yhi; y++){ int b = p->GetBin(o, y); sw += p->GetBinEntries(b); swz += p->GetBinEntries(b)*p->GetBinContent(b); }
        if(sw <= 0) return false;
        n[o] = sw; m[o] = swz/sw;
    }
    q[0] = m[1];
    q[1] = m[2];
    q[2] = (n[3]*m[3] - n[4]*m[4])/(n[3] - n[4]);
    q[3] = m[5] - m[6];
    q[4] = 0.5*(m[7]*m[7] + m[8]*m[8]) - m[7]*m[8];
    return true;
}

// value (summed chunks) and chunk-spread error of all derived quantities in shape bins ylo..yhi of variable k
static void ValErr(int k, int ylo, int yhi, double* v, double* e){
    Derive(tot.p[k], ylo, yhi, v);
    double W = tot.h[k]->Integral(ylo, yhi), s[NQ] = {0}, qc[NQ];
    for(auto& c : chunks){
        double w = c.h[k]->Integral(ylo, yhi);
        if(w <= 0 || !Derive(c.p[k], ylo, yhi, qc)) continue;
        for(int i = 0; i < NQ; i++) s[i] += (w/W)*(w/W)*(qc[i]-v[i])*(qc[i]-v[i]);
    }
    int K = chunks.size();
    for(int i = 0; i < NQ; i++) e[i] = (K > 1) ? sqrt(s[i]*K/(K-1.)) : 0;
}

// Delta gamma112 (blue) and 2a1^2 (red) vs elliptic-flow variable xq in bins of shape variable k, with fits
static void Panel(TVirtualPad* pad, int k, int xq, const char* label, bool pol2 = false){
    static int id = 0;
    pad->cd(); pad->SetLeftMargin(0.16); pad->SetBottomMargin(0.13); pad->SetTickx(); pad->SetTicky();
    TGraphErrors* g[2] = {new TGraphErrors(), new TGraphErrors()};
    double xmax = 0, v[NQ], e[NQ];
    for(int b = 1; b <= 50; b++){
        if(tot.h[k]->GetBinContent(b) < MINEV) continue;
        ValErr(k, b, b, v, e);
        for(int j = 0; j < 2; j++){ int n = g[j]->GetN(); g[j]->SetPoint(n, v[xq], v[3+j]); g[j]->SetPointError(n, e[xq], e[3+j]); }
        xmax = max(xmax, v[xq]);
    }
    double fmax = (xq == 2) ? xmax + 0.01 : 0.13;    // reference fit range; resonance v2 extends further
    const int col[2] = {kBlue, kRed}, mk[2] = {kOpenCircle, kOpenDiamond};
    TMultiGraph* mg = new TMultiGraph();
    for(int j = 0; j < 2; j++){ g[j]->SetMarkerStyle(mk[j]); g[j]->SetMarkerColor(col[j]); g[j]->SetLineColor(col[j]); mg->Add(g[j], "P"); }
    mg->SetTitle(Form(";%s{%s};#Delta#gamma_{112}, 2a_{1}^{2}", qTitle[xq], kTitle[k]));
    mg->Draw("A");
    mg->GetXaxis()->SetLimits(0, fmax);
    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.045);
    tex.DrawLatex(0.20, 0.88, label);
    TString plain(label);                                                // label for the text summary
    plain.ReplaceAll("v_{2}", "v2"); plain.ReplaceAll("q_{2}^{2}", "q2^2"); plain.ReplaceAll("q_{2,B}", "q2,B");

    for(int j = 0; j < 2; j++){
        for(int order = 1; order <= (pol2 ? 2 : 1); order++){
            TF1* f = new TF1(Form("fit%d", id++), order == 1 ? "pol1" : "pol2", -0.01, fmax);
            f->SetLineColor(col[j]); f->SetLineStyle(order == 1 ? 1 : 2);
            g[j]->Fit(f, "Q0R");
            f->Draw("same");
            double y0 = f->GetParameter(0), ey0 = f->GetParError(0), sc = (1-v2incl)*(1-v2incl);
            if(order == 1){
                tex.SetTextColor(col[j]);
                tex.DrawLatex(0.20, 0.81 - 0.06*j, Form("y-int(%s) = (%.2f#pm%.2f)#times10^{-5}", qTitle[3+j], 1e5*y0, 1e5*ey0));
            }
            summary += Form("  %-30s %-8s %-5s %+9.2f %8.2f  %+9.2f  %6.2f\n", plain.Data(), j ? "2a1^2" : "Dgamma", order == 1 ? "pol1" : "pol2",
                            1e5*y0, 1e5*ey0, 1e5*y0*sc, f->GetNDF() > 0 ? f->GetChisquare()/f->GetNDF() : 0.);
        }
    }
}

void MakeFigure_ESS(const char* base = "/media/Students/aliu/ess_all", const char* names = "day1 day2 day3 day4"){
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    TH1::AddDirectory(false);

    // ---- read the finished chunks ----
    TObjArray* tok = TString(names).Tokenize(" ");
    double vconst[3] = {0};
    for(int i = 0; i < tok->GetEntries(); i++){
        TString name = ((TObjString*)tok->At(i))->GetString();
        for(int c = 0; c < 1000; c++){
            TString d = Form("%s/%s/chunks/c%03d", base, name.Data(), c);
            if(gSystem->AccessPathName(d)) break;                               // no more chunks
            if(gSystem->AccessPathName(d + "/.done.ESS")){ printf("  skipping unfinished %s\n", d.Data()); continue; }
            TFile* f = TFile::Open(d + "/cen5.ESS_job0.root");
            if(!f || f->IsZombie()){ printf("  cannot open %s\n", d.Data()); continue; }
            Sample s;
            for(int k = 0; k < 3; k++){
                s.p[k] = (TProfile2D*)f->Get(Form("pESS_%s", kTag[k]))->Clone(); s.p[k]->SetDirectory(0);
                s.h[k] = (TH1D*)f->Get(Form("h%s", kTag[k]))->Clone();
            }
            TH1* hc = (TH1*)f->Get("hConst");
            for(int k = 0; k < 3; k++) vconst[k] = hc->GetBinContent(k+1);
            f->Close();
            chunks.push_back(s);
        }
    }
    if(chunks.empty()){ printf("no finished chunks under %s\n", base); return; }
    for(int k = 0; k < 3; k++){
        tot.p[k] = (TProfile2D*)chunks[0].p[k]->Clone(); tot.p[k]->SetDirectory(0);
        tot.h[k] = (TH1D*)chunks[0].h[k]->Clone();
        for(size_t c = 1; c < chunks.size(); c++){ tot.p[k]->Add(chunks[c].p[k]); tot.h[k]->Add(chunks[c].h[k]); }
    }

    // ---- inclusive values (all shape bins incl. overflow), Fig. 4(a) analogue ----
    double v[NQ], e[NQ];
    ValErr(0, 0, 51, v, e);
    v2incl = v[0];
    summary += Form("ESS (arXiv:2307.14997) on EPOS 30-40%% Au+Au 200 GeV: pi+-, 0.2<pT<2, |y|<1, true RP\n");
    summary += Form("%zu chunks, %.0f events\n", chunks.size(), tot.h[0]->Integral(0, 51));
    summary += Form("normalization (pass 1): v2{2} = %.5f, v2pair{2} = %.5f, v2B{2} = %.5f\n", vconst[0], vconst[1], vconst[2]);
    summary += Form("<single q2^2> = %.3f, <pair q2^2> = %.3f, <q2,B> = %.3f (within 0-5)\n\n",
                    tot.h[0]->GetMean(), tot.h[1]->GetMean(), tot.h[2]->GetMean());
    summary += "inclusive (errors: chunk spread)\n";
    const char* qn[NQ] = {"single v2", "pair v2", "resonance v2", "Dgamma112", "2a1^2"};
    for(int i = 0; i < NQ; i++) summary += Form("  %-14s %+.4e +- %.1e\n", qn[i], v[i], e[i]);
    summary += Form("\ny-intercepts in 1e-5 (fit range [-0.01,0.13] in v2; scaled = x (1-<v2>)^2 = x %.4f)\n", (1-v2incl)*(1-v2incl));
    summary += Form("  %-30s %-8s %-5s %9s %8s  %9s  %6s\n", "recipe", "quantity", "fit", "y-int", "error", "scaled", "chi2/ndf");

    // ---- Fig. 6/7: the four combinations ----
    TCanvas* c6 = new TCanvas("c6", "ESS four recipes", 1100, 900);
    c6->Divide(2, 2);
    Panel(c6->cd(1), 0, 0, "(a) single v_{2}{single q_{2}^{2}}");
    Panel(c6->cd(2), 1, 1, "(b) pair v_{2}{pair q_{2}^{2}}");
    Panel(c6->cd(3), 1, 0, "(c) single v_{2}{pair q_{2}^{2}}");
    Panel(c6->cd(4), 0, 1, "(d) pair v_{2}{single q_{2}^{2}}");
    c6->SaveAs(Form("%s/ESS_fig6.png", base));

    // ---- Fig. 5 (resonance v2) and Fig. 9 (q2,B) ----
    TCanvas* c9 = new TCanvas("c9", "ESS resonance v2 and q2,B", 1500, 450);
    c9->Divide(3, 1);
    Panel(c9->cd(1), 0, 2, "Fig.5 resonance v_{2}{single q_{2}^{2}}");
    c9->cd(2); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.13); gPad->SetTickx(); gPad->SetTicky();
    TGraphErrors* gB = new TGraphErrors();
    for(int b = 1; b <= 50; b++){
        if(tot.h[2]->GetBinContent(b) < MINEV) continue;
        ValErr(2, b, b, v, e);
        gB->SetPoint(gB->GetN(), tot.h[2]->GetBinCenter(b), v[0]); gB->SetPointError(gB->GetN()-1, 0, e[0]);
    }
    gB->SetTitle(";q_{2,B};single v_{2}"); gB->SetMarkerStyle(kFullCircle); gB->Draw("AP");
    TLatex tex; tex.SetNDC(); tex.SetTextSize(0.045); tex.DrawLatex(0.20, 0.88, "Fig.9(a)");
    Panel(c9->cd(3), 2, 0, "Fig.9 single v_{2}{q_{2,B}}", true);
    c9->SaveAs(Form("%s/ESS_fig5_9.png", base));

    printf("\n%s", summary.c_str());
    FILE* fs = fopen(Form("%s/ESS_summary.txt", base), "w");
    if(fs){ fputs(summary.c_str(), fs); fclose(fs); }
}
