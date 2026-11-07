using namespace std;

#include "stdio.h"
#include <cmath>
#include <string>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"

// Combine per-chunk Correlators.C outputs (v2, gamma, delta, kappa112 for charged pions vs the true RP)
// into one result per dataset and one over all datasets.
//
// Layout (made by scripts/correlators_all.sh): <base>/<name>/chunks/cNNN/cen5.correlators_job0.root,
// a chunk counting only when cNNN/.done exists. Chunks are disjoint sets of whole files, so they are
// independent subsamples; every error here comes from the chunk-to-chunk spread,
//   var(x) = K/(K-1) Sum_k (w_k/W)^2 (x_k - x)^2 ,  w_k = events used in chunk k,
// which holds the nfreeze-20 correlations (20 events per hydro event) that per-event profile errors
// and Correlators.C's round-robin subgroups miss. kappa112 = Dgamma/(v2*Ddelta) is formed per chunk.
//
// For each set it writes <base>/<set>/cen5.correlators_job0.root: the summed profiles, hMeta, and
// hKappa holding the combined kappa112 with its chunk-spread error, so MakeFigure_Correlators.C and
// MakeFigure_Kappa_STAR.C run on it unmodified. hErrChunk holds the chunk-spread errors of all
// eight observables. The table goes to <base>/summary.txt.
//
// Usage: root -l -b -q 'Combine_Correlators.C+("/media/Students/aliu/correlators_all","auau200-7 alpha_day1")'

static const int NOBS = 8;
static const char* obsName[NOBS]  = {"v2","gamma_OS","gamma_SS","Dgamma","delta_OS","delta_SS","Ddelta","kappa112"};
static const char* profName[NOBS-1] = {"pV2int","pGammaOS","pGammaSS","pDGamma","pDeltaOS","pDeltaSS","pDDelta"};

struct Chunk { double w; double x[NOBS]; double meta[7]; };

struct Result {
    TString name; int K = 0; double W = 0;
    double val[NOBS], err[NOBS], errProf[NOBS];
};

static double KappaOf(const double* x){ return (x[0]*x[6] != 0) ? x[3]/(x[0]*x[6]) : 0.; }

static vector<Chunk> ReadChunks(const TString& dir, vector<TString>& files){
    vector<Chunk> out;
    for(int c = 0; c < 1000; c++){
        TString d = Form("%s/chunks/c%03d", dir.Data(), c);
        if(gSystem->AccessPathName(d)) break;                       // kTRUE = does not exist
        if(gSystem->AccessPathName(d + "/.done")) { printf("  skipping unfinished %s\n", d.Data()); continue; }
        TString f = d + "/cen5.correlators_job0.root";
        TFile* tf = TFile::Open(f);
        if(!tf || tf->IsZombie()){ printf("  cannot open %s\n", f.Data()); continue; }
        Chunk ch;
        TProfile* p0 = (TProfile*)tf->Get(profName[0]);
        ch.w = p0->GetBinEntries(1);
        for(int i = 0; i < NOBS-1; i++) ch.x[i] = ((TProfile*)tf->Get(profName[i]))->GetBinContent(1);
        ch.x[7] = KappaOf(ch.x);
        TH1* hm = (TH1*)tf->Get("hMeta");
        for(int i = 0; i < 7; i++) ch.meta[i] = hm->GetBinContent(i+1);
        tf->Close();
        if(ch.w > 0){ out.push_back(ch); files.push_back(f); }
    }
    return out;
}

// sum the chunk files, compute values + chunk-spread errors, write <outDir>/cen5.correlators_job0.root
static Result Combine(const TString& name, const vector<Chunk>& ch, const vector<TString>& files, const TString& outDir){
    Result r; r.name = name; r.K = ch.size();
    gSystem->mkdir(outDir, true);
    TFile fout(outDir + "/cen5.correlators_job0.root", "RECREATE");

    // summed profiles (all TProfiles in the chunk files; hMeta/hKappa are rebuilt below)
    const char* allProf[] = {"pV2pT","pV2eta","pV2int","pGammaOS","pGammaSS","pDGamma","pDeltaOS","pDeltaSS","pDDelta"};
    vector<TProfile*> sum(9, nullptr);
    for(auto& f : files){
        TFile* tf = TFile::Open(f);
        for(int i = 0; i < 9; i++){
            TProfile* p = (TProfile*)tf->Get(allProf[i]);
            if(!sum[i]){ fout.cd(); sum[i] = (TProfile*)p->Clone(allProf[i]); sum[i]->SetDirectory(&fout); }
            else sum[i]->Add(p);
        }
        tf->Close();
    }
    for(int i = 0; i < NOBS-1; i++){
        TProfile* p = sum[2+i];
        r.val[i] = p->GetBinContent(1);
        r.errProf[i] = p->GetBinError(1);
    }
    r.W = sum[2]->GetBinEntries(1);
    r.val[7] = KappaOf(r.val);
    r.errProf[7] = 0;

    // chunk-spread errors
    for(int i = 0; i < NOBS; i++){
        double v = 0;
        for(auto& c : ch) v += (c.w/r.W)*(c.w/r.W)*(c.x[i]-r.val[i])*(c.x[i]-r.val[i]);
        r.err[i] = (r.K > 1) ? sqrt(v * r.K/(r.K-1.)) : 0.;
    }

    fout.cd();
    // hMeta as in Correlators.C; bins 2-7 are event-weighted chunk averages (bins 4-7, the
    // pair-weighted correlators, are approximate: pairs per event are nearly the same in every chunk)
    TH1D* hMeta = new TH1D("hMeta","meta: 1 Nev, 2 <N+>, 3 <N->, 4 gOS_pw, 5 gSS_pw, 6 dOS_pw, 7 dSS_pw",7,0.5,7.5);
    hMeta->SetBinContent(1, r.W);
    for(int b = 1; b < 7; b++){
        double s = 0; for(auto& c : ch) s += c.w * c.meta[b];
        hMeta->SetBinContent(b+1, s/r.W);
    }
    TH1D* hKappa = new TH1D("hKappa","#kappa_{112}=#Delta#gamma/(v_{2}#Delta#delta); ; value",1,0.5,1.5);
    hKappa->SetBinContent(1, r.val[7]);
    hKappa->SetBinError(1, r.err[7]);
    TH1D* hErr = new TH1D("hErrChunk","chunk-spread values/errors: v2 gOS gSS Dg dOS dSS Dd kappa",NOBS,0.5,NOBS+0.5);
    for(int i = 0; i < NOBS; i++){ hErr->SetBinContent(i+1, r.val[i]); hErr->SetBinError(i+1, r.err[i]); }
    TH1D* hN = new TH1D("hNchunk","number of independent chunks",1,0.5,1.5);
    hN->SetBinContent(1, r.K);
    fout.Write();
    fout.Close();
    return r;
}

void Combine_Correlators(const char* base = "/media/Students/aliu/correlators_all",
                         const char* names = "auau200-7 alpha_day1"){

    TObjArray* tok = TString(names).Tokenize(" ");
    vector<Result> res;
    vector<Chunk> allCh; vector<TString> allFiles;
    for(int i = 0; i < tok->GetEntries(); i++){
        TString n = ((TObjString*)tok->At(i))->GetString();
        vector<TString> files;
        vector<Chunk> ch = ReadChunks(Form("%s/%s", base, n.Data()), files);
        printf("%s: %zu chunks\n", n.Data(), ch.size());
        if(ch.empty()) continue;
        res.push_back(Combine(n, ch, files, Form("%s/%s", base, n.Data())));
        allCh.insert(allCh.end(), ch.begin(), ch.end());
        allFiles.insert(allFiles.end(), files.begin(), files.end());
    }
    if(res.size() > 1) res.push_back(Combine("combined", allCh, allFiles, Form("%s/combined", base)));

    // ---------------- summary ----------------
    string s;
    s += "Correlators.C over all events: charged pions, 0.2<pT<2, |y|<1, true RP, 30-40% (8.4<b<9.2 fm)\n";
    s += "errors: chunk-to-chunk spread (independent file chunks); DEFF = (chunk err / per-event profile err)^2\n\n";
    for(auto& r : res){
        s += Form("== %s: %.0f events, %d chunks ==\n", r.name.Data(), r.W, r.K);
        s += Form("  %-10s %13s %12s %8s %6s\n", "observable", "value", "stat.err", "signif.", "DEFF");
        for(int i = 0; i < NOBS; i++){
            double deff = (r.errProf[i] > 0) ? pow(r.err[i]/r.errProf[i], 2) : 0.;
            s += Form("  %-10s %+13.4e %12.3e %8.1f", obsName[i], r.val[i], r.err[i], r.err[i] > 0 ? fabs(r.val[i])/r.err[i] : 0.);
            s += (deff > 0) ? Form(" %6.2f\n", deff) : "\n";
        }
        s += "\n";
    }
    // consistency of the individual datasets
    int nd = (res.size() > 1) ? res.size()-1 : 0;
    if(nd == 2){
        s += Form("== %s vs %s: (2nd - 1st)/sigma ==\n", res[0].name.Data(), res[1].name.Data());
        for(int i = 0; i < NOBS; i++){
            double d = res[1].val[i]-res[0].val[i], e = hypot(res[0].err[i], res[1].err[i]);
            s += Form("  %-10s %+.2f\n", obsName[i], e > 0 ? d/e : 0.);
        }
    }
    printf("\n%s", s.c_str());
    FILE* fs = fopen(Form("%s/summary.txt", base), "w");
    if(fs){ fputs(s.c_str(), fs); fclose(fs); }
}
