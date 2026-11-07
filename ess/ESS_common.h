// ESS_common.h -- event and track selection shared by V2Const.C (pass 1) and ESS.C (pass 2), so that both
// passes of the event-shape-selection analysis (Xu et al., arXiv:2307.14997) see exactly the same events.
//
//   POI          : pi+/- (EPOS id +-120), final state, 0.2 < pT < 2 GeV/c, |y| < 1           (paper Sec. III)
//   q2,B particles: charged pi/K/p (EPOS ids +-120, +-130, +-1120), 0.2 < pT < 2, 1 < |y| < 2  (paper App. A)
//   events       : 8.4 < b < 9.2 fm (30-40%), at least 2 pi+ and 2 pi- among the POI
//   Psi_RP       : the "phi" branch (true reaction plane)
#ifndef ESS_COMMON_H
#define ESS_COMMON_H

#include <cmath>
#include <cstdlib>
#include <vector>
#include "TChain.h"
#include "TLeaf.h"

namespace ess {

const double ptMin = 0.2, ptMax = 2.0;
const double yPOI = 1.0;                 // POI:   |y| < 1
const double yBlo = 1.0, yBhi = 2.0;     // q2,B:  1 < |y| < 2
const double bMin = 8.4, bMax = 9.2;     // 30-40% centrality
const double M_PI_CH = 0.13957, M_K_CH = 0.493677, M_PROTON = 0.938272;

struct Event {
    double cpsi, spsi;                   // cos/sin(Psi_RP)
    std::vector<double> px, py, c, s;    // POI momentum and cos/sin(phi)
    std::vector<int> q;                  // POI charge
    int Np, Nm;                          // number of pi+ / pi-
    double QBx, QBy; int NB;             // q2,B sub-event: Sum cos2phi, Sum sin2phi, multiplicity
};

inline TChain* MakeChain(const char* filepat){
    TChain* ch = new TChain("teposevent");
    ch->Add(filepat);
    ch->SetBranchStatus("*", 0);
    const char* act[] = {"np","bim","phi","px","py","pz","id","ist"};
    for(const char* b : act) ch->SetBranchStatus(b, 1);
    return ch;
}

inline double Rapidity(double m, double px, double py, double pz){
    double E = std::sqrt(m*m + px*px + py*py + pz*pz);
    return 0.5*std::log((E+pz)/(E-pz));
}

// Read entry i; returns false if the event fails the selection.
inline bool ReadEvent(TChain* ch, Long64_t i, Event& ev){
    ch->GetEntry(i);
    double bim = ch->GetLeaf("bim")->GetValue(0);
    if(bim < bMin || bim > bMax) return false;
    double psi = ch->GetLeaf("phi")->GetValue(0);
    ev.cpsi = std::cos(psi); ev.spsi = std::sin(psi);
    ev.px.clear(); ev.py.clear(); ev.c.clear(); ev.s.clear(); ev.q.clear();
    ev.Np = ev.Nm = 0; ev.QBx = ev.QBy = 0; ev.NB = 0;

    int np = (int)ch->GetLeaf("np")->GetValue(0);
    TLeaf *lpx = ch->GetLeaf("px"), *lpy = ch->GetLeaf("py"), *lpz = ch->GetLeaf("pz");
    TLeaf *lid = ch->GetLeaf("id"), *lst = ch->GetLeaf("ist");
    for(int t = 0; t < np; t++){
        if((int)lst->GetValue(t) != 0) continue;                     // final state only
        int id = (int)lid->GetValue(t), aid = std::abs(id);
        double m = (aid == 120) ? M_PI_CH : (aid == 130) ? M_K_CH : (aid == 1120) ? M_PROTON : -1;
        if(m < 0) continue;                                           // not a charged pi/K/p
        double px = lpx->GetValue(t), py = lpy->GetValue(t), pz = lpz->GetValue(t);
        double pt = std::sqrt(px*px + py*py);
        if(pt < ptMin || pt > ptMax) continue;
        double ay = std::fabs(Rapidity(m, px, py, pz));
        double c = px/pt, s = py/pt;                                  // cos(phi), sin(phi)
        if(ay > yBlo && ay < yBhi){                                   // q2,B sub-event
            ev.QBx += c*c - s*s; ev.QBy += 2*c*s; ev.NB++;
        }
        if(aid == 120 && ay < yPOI){                                  // POI
            ev.px.push_back(px); ev.py.push_back(py); ev.c.push_back(c); ev.s.push_back(s);
            ev.q.push_back(id > 0 ? 1 : -1);
            if(id > 0) ev.Np++; else ev.Nm++;
        }
    }
    return ev.Np >= 2 && ev.Nm >= 2;
}

// cos(2 phi_p), sin(2 phi_p) of the pair momentum p_a + p_b (same angle as (lA+lB).Phi() in the reference macro)
inline void PairDir(const Event& ev, int a, int b, double& c2, double& s2){
    double X = ev.px[a] + ev.px[b], Y = ev.py[a] + ev.py[b], r2 = X*X + Y*Y;
    c2 = (X*X - Y*Y)/r2; s2 = 2*X*Y/r2;
}

} // namespace ess
#endif
