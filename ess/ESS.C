// ESS.C -- pass 2 of the event-shape-selection analysis (Xu et al., arXiv:2307.14997), ported from the authors'
// macro Parity_EPD_pbp.C to EPOS. Each event is classified by three event-shape variables:
//
//   single q2^2 = |Sum_i e^{2i phi_i}|^2   / [N     (1 + N     v2{2}^2)]       Eq. (9)
//   pair   q2^2 = |Sum_p e^{2i phi_p}|^2   / [Npair (1 + Npair v2pair{2}^2)]   Eq. (11), phi_p = azimuth of p_a+p_b
//   q2,B        = sqrt(Eq. 9 for the 1<|y|<2 particles, with their own v2B{2})   App. A
//
// with the two-particle-correlation constants from pass 1 (V2Const.C, hadd-ed into constFile), as in the paper.
// The reaction plane enters only the observables below. In bins of each variable we record
//
//   1 single v2  <cos(2phi - 2Psi)>            per POI          5 gamma112 OS  <cos(phi_a+phi_b-2Psi)>  per OS pair
//   2 pair v2    <cos(2phi_p - 2Psi)>          per pair (all)   6 gamma112 SS                         per SS pair
//   3 pair v2 OS (for resonance v2, Eq. 12)    per OS pair      7 a1(pi+)      <sin(phi - Psi)>        per pi+
//   4 pair v2 SS                               per SS pair      8 a1(pi-)                              per pi-
//   9 gamma132 OS  <cos(phi_a-3phi_b+2Psi)>    per OS pair     10 gamma132 SS                         per SS pair
//
// Each event fills its per-event mean with weight = number of terms, so the bin content is the same per-pair /
// per-particle average as the reference's per-term fills and GetBinEntries() is the number of terms.
// gamma132 is the reference's correlator1 (Parity_EPD_pbp.C), which takes phi_a from the lower track index; here both
// orders of each pair are averaged, so the result does not depend on the particle ordering in the EPOS output.
//
// Usage (inside a chunk directory): root -l -b -q 'ESS.C+(0,"z-auau_run_*.root","/path/v2const.root",8.4,9.2)'
// (bmin, bmax = impact-parameter window of the sample in fm; must match pass 1)
#include <cstdio>
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "ESS_common.h"

void ESS(int job = 0, const char* filepat = "z-auau_run_*.root", const char* constFile = "v2const.root",
         double bmin = 8.4, double bmax = 9.2){
    ess::bMin = bmin; ess::bMax = bmax;

    // ---- normalization constants from pass 1 ----
    double v2sq[3];
    {
        TFile fc(constFile);
        TProfile* p = (TProfile*)fc.Get("pV2sq");
        if(!p){ printf("cannot read pV2sq from %s\n", constFile); return; }
        for(int k = 0; k < 3; k++) v2sq[k] = p->GetBinContent(k+1);
    }
    printf("constants: v2{2} = %.5f, v2pair{2} = %.5f, v2B{2} = %.5f\n", sqrt(v2sq[0]), sqrt(v2sq[1]), sqrt(v2sq[2]));

    TChain* ch = ess::MakeChain(filepat);
    TFile fout(Form("cen5.ESS_job%d.root", job), "RECREATE");

    // index 0 = single q2^2, 1 = pair q2^2, 2 = q2,B ; all average ~1 with the v2{2} normalization: 0.1-wide bins as in
    // the reference, first 50
    const double qmax[3] = {5, 5, 5};
    const char* tag[3]   = {"q2s", "q2p", "q2B"};
    const char* title[3] = {"single q_{2}^{2}", "pair q_{2}^{2}", "q_{2,B}"};
    TProfile2D* pESS[3]; TH1D* hq[3];
    for(int k = 0; k < 3; k++){
        pESS[k] = new TProfile2D(Form("pESS_%s", tag[k]), Form("observables vs %s;observable;%s", title[k], title[k]),
                                 10, 0.5, 10.5, 50, 0, qmax[k]);
        pESS[k]->Sumw2();                                                // weighted fills
        hq[k]  = new TH1D(Form("h%s", tag[k]), Form("events vs %s;%s", title[k], title[k]), 50, 0, qmax[k]);
    }
    TH1D* hConst = new TH1D("hConst", "constants used: v2{2}, v2pair{2}, v2B{2}", 3, 0.5, 3.5);
    for(int k = 0; k < 3; k++) hConst->SetBinContent(k+1, sqrt(v2sq[k]));

    ess::Event ev;
    std::vector<double> c3, s3;                                          // cos/sin(3phi) of the POI
    Long64_t n = ch->GetEntries(), nUsed = 0;
    printf("%lld entries in chain\n", n);
    for(Long64_t i = 0; i < n; i++){
        if((i+1) % 10000 == 0) printf("Processing entry == %lld == out of %lld.\n", i+1, n);
        if(!ess::ReadEvent(ch, i, ev)) continue;
        const int N = ev.px.size();
        const double C2 = ev.cpsi*ev.cpsi - ev.spsi*ev.spsi, S2 = 2*ev.cpsi*ev.spsi;   // cos/sin(2 Psi)

        // ---- single particles: Q-vector, v2, a1 ----
        double C = 0, S = 0, v2s = 0, a1p = 0, a1m = 0;
        c3.resize(N); s3.resize(N);
        for(int a = 0; a < N; a++){
            double c2 = ev.c[a]*ev.c[a] - ev.s[a]*ev.s[a], s2 = 2*ev.c[a]*ev.s[a];
            c3[a] = ev.c[a]*(4*ev.c[a]*ev.c[a] - 3); s3[a] = ev.s[a]*(3 - 4*ev.s[a]*ev.s[a]);    // cos/sin(3phi)
            C += c2; S += s2;
            v2s += c2*C2 + s2*S2;                                        // cos(2phi - 2Psi)
            double a1 = ev.s[a]*ev.cpsi - ev.c[a]*ev.spsi;                // sin(phi - Psi)
            if(ev.q[a] > 0) a1p += a1; else a1m += a1;
        }

        // ---- pairs: Q-vector, pair v2, gamma112 and gamma132, split OS/SS ----
        double Cp = 0, Sp = 0, vOS = 0, vSS = 0, gOS = 0, gSS = 0, hOS = 0, hSS = 0, nOS = 0, nSS = 0, c2, s2;
        for(int a = 0; a < N; a++) for(int b = a+1; b < N; b++){
            ess::PairDir(ev, a, b, c2, s2);
            Cp += c2; Sp += s2;
            const double ca = ev.c[a], sa = ev.s[a], cb = ev.c[b], sb = ev.s[b];
            double v = c2*C2 + s2*S2;                                                   // cos(2phi_p - 2Psi)
            double g = (ca*cb - sa*sb)*C2 + (sa*cb + ca*sb)*S2;                         // cos(phi_a+phi_b-2Psi)
            double h = 0.5*((ca*c3[b] + sa*s3[b])*C2 - (sa*c3[b] - ca*s3[b])*S2       // cos(phi_a-3phi_b+2Psi)
                          + (cb*c3[a] + sb*s3[a])*C2 - (sb*c3[a] - cb*s3[a])*S2);     // cos(phi_b-3phi_a+2Psi)
            if(ev.q[a]*ev.q[b] < 0){ vOS += v; gOS += g; hOS += h; nOS++; }
            else                   { vSS += v; gSS += g; hSS += h; nSS++; }
        }
        const double Npair = nOS + nSS;

        // ---- event-shape variables ----
        double shape[3];
        shape[0] = (C*C + S*S)/(N*(1 + N*v2sq[0]));
        shape[1] = (Cp*Cp + Sp*Sp)/(Npair*(1 + Npair*v2sq[1]));
        shape[2] = (ev.NB > 0) ? sqrt((ev.QBx*ev.QBx + ev.QBy*ev.QBy)/(ev.NB*(1 + ev.NB*v2sq[2]))) : -1;

        const double mean[10] = {v2s/N, (vOS+vSS)/Npair, vOS/nOS, vSS/nSS, gOS/nOS, gSS/nSS, a1p/ev.Np, a1m/ev.Nm, hOS/nOS, hSS/nSS};
        const double nTerm[10] = {(double)N, Npair, nOS, nSS, nOS, nSS, (double)ev.Np, (double)ev.Nm, nOS, nSS};
        for(int k = 0; k < 3; k++){
            if(shape[k] < 0) continue;                                   // no q2,B particles
            hq[k]->Fill(shape[k]);
            for(int o = 0; o < 10; o++) pESS[k]->Fill(o+1, shape[k], mean[o], nTerm[o]);
        }
        nUsed++;
    }
    printf("<single q2^2> = %.3f, <pair q2^2> = %.3f, <q2,B> = %.3f\n", hq[0]->GetMean(), hq[1]->GetMean(), hq[2]->GetMean());
    printf("events used = %lld\n", nUsed);
    fout.Write();
}
