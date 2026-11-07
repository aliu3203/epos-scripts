// V2Const.C -- pass 1 of the ESS analysis: the ensemble two-particle cumulants that normalize the event-shape
// variables of arXiv:2307.14997, Eqs. (8)-(11), measured on the same events as ESS.C:
//
//   v2{2}^2 = Sum_ev (|Q|^2 - M) / Sum_ev M(M-1),   Q = Sum_j e^{2i phi_j}
//
// for three samples: bin 1 = single POI (M = N), bin 2 = POI pairs (phi_p of p_a+p_b, M = Npair),
// bin 3 = q2,B particles (M = N_B). Filled as a weighted TProfile, so hadd over chunks gives the global value.
//
// Usage (inside a chunk directory): root -l -b -q 'V2Const.C+(0)'  ->  cen5.V2Const_job0.root
#include <cstdio>
#include "TFile.h"
#include "TProfile.h"
#include "ESS_common.h"

void V2Const(int job = 0, const char* filepat = "z-auau_run_*.root"){
    TChain* ch = ess::MakeChain(filepat);
    TFile fout(Form("cen5.V2Const_job%d.root", job), "RECREATE");
    TProfile* pV2sq = new TProfile("pV2sq", "v_{2}{2}^{2}: 1 single, 2 pair, 3 q_{2,B} particles", 3, 0.5, 3.5);

    ess::Event ev;
    Long64_t n = ch->GetEntries(), nUsed = 0;
    printf("%lld entries in chain\n", n);
    for(Long64_t i = 0; i < n; i++){
        if((i+1) % 10000 == 0) printf("Processing entry == %lld == out of %lld.\n", i+1, n);
        if(!ess::ReadEvent(ch, i, ev)) continue;
        const int N = ev.px.size();

        double C = 0, S = 0;
        for(int a = 0; a < N; a++){ C += ev.c[a]*ev.c[a] - ev.s[a]*ev.s[a]; S += 2*ev.c[a]*ev.s[a]; }

        double Cp = 0, Sp = 0, c2, s2;
        for(int a = 0; a < N; a++) for(int b = a+1; b < N; b++){ ess::PairDir(ev, a, b, c2, s2); Cp += c2; Sp += s2; }
        double Npair = 0.5*N*(N-1);

        pV2sq->Fill(1, (C*C + S*S - N)/(N*(N-1.)), N*(N-1.));
        pV2sq->Fill(2, (Cp*Cp + Sp*Sp - Npair)/(Npair*(Npair-1)), Npair*(Npair-1));
        if(ev.NB >= 2) pV2sq->Fill(3, (ev.QBx*ev.QBx + ev.QBy*ev.QBy - ev.NB)/(ev.NB*(ev.NB-1.)), ev.NB*(ev.NB-1.));
        nUsed++;
    }

    const char* nm[3] = {"v2{2}", "v2pair{2}", "v2B{2}"};
    for(int k = 1; k <= 3; k++){
        double x = pV2sq->GetBinContent(k);
        printf("%-10s = %.5f\n", nm[k-1], x > 0 ? sqrt(x) : 0.);
    }
    printf("events used = %lld\n", nUsed);
    fout.Write();
}
