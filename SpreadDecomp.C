// Where does the residual error difference between Dgamma112 and Dgamma132(1ord)
// come from?  Decompose each Delta into its OS and SS pieces and read off the
// per-event RMS of every piece plus the OS-SS covariance, straight from the
// TProfiles written by Correlators_Cent.C.
#include "TFile.h"
#include "TProfile.h"
#include "TMath.h"
#include <cstdio>

static void get(TFile* f, const char* nm, double& mean, double& rms, double& n){
    TProfile* p = (TProfile*)f->Get(nm);
    if(!p){ printf("MISSING %s\n",nm); mean=rms=n=0; return; }
    n    = p->GetBinEntries(1);
    mean = p->GetBinContent(1);
    TProfile* q = (TProfile*)p->Clone(Form("%s_s",nm));
    q->SetErrorOption("s");
    rms  = q->GetBinError(1);          // spread of the filled values, not the error on mean
}

void Spread(const char* fn = "cen.correlators_cent_poi2_job2.root", int ic = 4){
    TFile* f = TFile::Open(fn);
    if(!f || f->IsZombie()){ printf("cannot open %s\n",fn); return; }

    const char* base[] = {"pG112OS","pG112SS","pDg112",
                          "pG132OS1","pG132SS1","pDg132_1",
                          "pG132OS","pG132SS","pDg132"};
    double m[9],r[9],n[9];
    printf("file %s, centrality index %d\n\n", fn, ic);
    printf("%-12s %10s %14s %12s %12s\n","profile","entries","mean","RMS(event)","err(mean)");
    for(int k=0;k<9;k++){
        get(f, Form("%s_c%d",base[k],ic), m[k], r[k], n[k]);
        printf("%-12s %10.0f %+14.6e %12.5e %12.4e\n",base[k],n[k],m[k],r[k],r[k]/sqrt(n[k]));
    }

    // Var(OS-SS) = Var(OS) + Var(SS) - 2Cov  ->  solve for Cov and the correlation
    struct { const char* tag; int os,ss,d; } grp[3] = {
        {"gamma112       ",0,1,2},
        {"gamma132 1ord  ",3,4,5},
        {"gamma132 sym   ",6,7,8} };
    printf("\n%-15s %11s %11s %11s %9s %11s\n",
           "correlator","sig(OS)","sig(SS)","sig(Delta)","corr","cancel");
    for(int g=0;g<3;g++){
        double vo=r[grp[g].os]*r[grp[g].os], vs=r[grp[g].ss]*r[grp[g].ss], vd=r[grp[g].d]*r[grp[g].d];
        double cov = 0.5*(vo+vs-vd);
        double rho = cov/sqrt(vo*vs);
        // "cancel" = how much smaller Delta is than the uncorrelated expectation
        printf("%-15s %11.5e %11.5e %11.5e %9.4f %11.4f\n",
               grp[g].tag, sqrt(vo), sqrt(vs), sqrt(vd), rho, sqrt(vd/(vo+vs)));
    }

    printf("\nerr ratio  Dg132(1ord)/Dg112  = %.4f\n", (r[5]/r[2]));
    printf("err ratio  Dg132(sym) /Dg112  = %.4f\n", (r[8]/r[2]));
    printf("OS-only    132(1ord)/112      = %.4f\n", (r[3]/r[0]));
    printf("SS-only    132(1ord)/112      = %.4f\n", (r[4]/r[1]));
    f->Close();
}
