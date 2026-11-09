// V2Const.C -- pass 1 of the ESS analysis: the ensemble elliptic flow that normalizes the event-shape variables of
// arXiv:2307.14997, Eqs. (8)-(11), measured on the same events as ESS.C. EPOS provides the true reaction plane
// (branch "phi"), so v2 is taken directly from its definition instead of a two-particle cumulant:
//
//   v2 = Sum_ev Sum_j cos(2 phi_j - 2 Psi_RP) / Sum_ev M = Sum_ev (Qx cos2Psi + Qy sin2Psi) / Sum_ev M,  Q = Sum_j e^{2i phi_j}
//
// for three samples: bin 1 = single POI (M = N), bin 2 = POI pairs (phi_p of p_a+p_b, M = Npair),
// bin 3 = q2,B particles (M = N_B). Filled as a weighted TProfile, so hadd over chunks gives the global value.
// (A pair cumulant v2pair{2} would be dominated by pairs of pairs sharing a pion -- not flow.)
//
// Usage (inside a chunk directory): root -l -b -q 'V2Const.C+(0)'  ->  cen5.V2Const_job0.root
#include <cstdio>
#include "TFile.h"
#include "TProfile.h"
#include "ESS_common.h"

void V2Const(int job = 0, const char* filepat = "z-auau_run_*.root"){
    TChain* ch = ess::MakeChain(filepat);
    TFile fout(Form("cen5.V2Const_job%d.root", job), "RECREATE");
    TProfile* pV2 = new TProfile("pV2", "v_{2}{RP}: 1 single, 2 pair, 3 q_{2,B} particles", 3, 0.5, 3.5);

    ess::Event ev;
    Long64_t n = ch->GetEntries(), nUsed = 0;
    printf("%lld entries in chain\n", n);
    for(Long64_t i = 0; i < n; i++){
        if((i+1) % 10000 == 0) printf("Processing entry == %lld == out of %lld.\n", i+1, n);
        if(!ess::ReadEvent(ch, i, ev)) continue;
        const int N = ev.px.size();
        const double C2 = ev.cpsi*ev.cpsi - ev.spsi*ev.spsi, S2 = 2*ev.cpsi*ev.spsi;   // cos/sin(2 Psi_RP)

        double C = 0, S = 0;
        for(int a = 0; a < N; a++){ C += ev.c[a]*ev.c[a] - ev.s[a]*ev.s[a]; S += 2*ev.c[a]*ev.s[a]; }

        double Cp = 0, Sp = 0, c2, s2;
        for(int a = 0; a < N; a++) for(int b = a+1; b < N; b++){ ess::PairDir(ev, a, b, c2, s2); Cp += c2; Sp += s2; }
        double Npair = 0.5*N*(N-1);

        pV2->Fill(1, (C*C2 + S*S2)/N, N);
        pV2->Fill(2, (Cp*C2 + Sp*S2)/Npair, Npair);
        if(ev.NB > 0) pV2->Fill(3, (ev.QBx*C2 + ev.QBy*S2)/ev.NB, ev.NB);
        nUsed++;
    }

    const char* nm[3] = {"v2{RP}", "v2pair{RP}", "v2B{RP}"};
    for(int k = 1; k <= 3; k++) printf("%-10s = %.5f\n", nm[k-1], pV2->GetBinContent(k));
    printf("events used = %lld\n", nUsed);
    fout.Write();
}
