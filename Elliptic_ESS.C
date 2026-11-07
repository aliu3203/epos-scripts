using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TChain.h>
#include "TLeaf.h"
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

const float PI = TMath::Pi();

// EPOS particle ids (see analyze.C / idt.dt)
const int ID_pip = 120;    // pi+
const int ID_pim = -120;   // pi-

// POI cuts -- matching arXiv:2307.14997 (Xu et al.): charged pions, |y|<1, 0.2<pT<2, true RP
const float yCut   = 1.0;
const float ptMin  = 0.2;
const float ptMax  = 2.0;

// Event-shape-selection variables of arXiv:2307.14997, Eqs. (9),(10),(11):
//   single q2^2 = [ (Sum sin2phi)^2 + (Sum cos2phi)^2 ] / [ N (1 + N v2^2{2}) ]           (9)
//   pair v2     = <cos(2phi_p - 2Psi_RP)>                                                 (10)
//   pair q2^2   = [ (Sum sin2phi_p)^2 + (Sum cos2phi_p)^2 ] / [ Npair (1 + Npair v2pair^2{2}) ] (11)
//   single v2   = <cos(2phi - 2Psi_RP)>   (Eq. 1)
// phi_p is the azimuthal angle of the momentum sum of each POI pair.
// v2{2}, v2pair{2} are ensemble-averaged two-particle cumulants (static constants).

void Elliptic_ESS(int cen = 5, int job = 0){	//main_function

    TChain* chain = new TChain("teposevent");
    chain->Add("z-*.root");
    //chain->Add("../auau200-6/z-*.root");

    //only read the branches we use
    chain->SetBranchStatus("*",0);
    const char* activeBranches[] = {"np","bim","phi","px","py","pz","id","e","ist"};
    for(int ib = 0; ib < 9; ib++) chain->SetBranchStatus(activeBranches[ib],1);

    char fname_out[200];
    sprintf(fname_out,"cen%d.v2q2_ESS_job%d.root",cen,job);
    TFile fout(fname_out,"RECREATE");

    // per-event scalars stored between the two passes
    std::vector<int>    ev_N;
    std::vector<double> ev_Q2s;     // |Q_single|^2
    std::vector<double> ev_vsum_s;  // Sum cos(2phi   - 2Psi_RP)
    std::vector<double> ev_Npair;
    std::vector<double> ev_Q2p;     // |Q_pair|^2
    std::vector<double> ev_vsum_p;  // Sum cos(2phi_p - 2Psi_RP)

    // globals for the two-particle cumulants v2{2}^2 and v2pair{2}^2
    double num_s2 = 0, den_s2 = 0;  // Sum(|Q|^2 - N)      / Sum N(N-1)
    double num_p2 = 0, den_p2 = 0;  // Sum(|Qp|^2 - Npair) / Sum Npair(Npair-1)

    //================= PASS 1: build Q-vectors, collect per-event scalars =================
    Int_t nentries = chain->GetEntries();
    cout << nentries << "\n";
    for(int i = 0; i < nentries; i++){

        if((i+1)%1000==0) cout << "Pass1 entry == "<< i+1 <<" == out of "<<nentries<<".\n";
        chain->GetEntry(i);

        float bim = chain->GetLeaf("bim")->GetValue(0);
        if(bim < 8.4 || bim > 9.2) continue;	//30-40% centrality

        int NPTracks = (int)chain->GetLeaf("np")->GetValue(0);
        float PsiRP  = chain->GetLeaf("phi")->GetValue(0);	//true reaction plane
        float cos2Psi = cos(2.*PsiRP), sin2Psi = sin(2.*PsiRP);

        TLeaf* leaf_px  = chain->GetLeaf("px");
        TLeaf* leaf_py  = chain->GetLeaf("py");
        TLeaf* leaf_pz  = chain->GetLeaf("pz");
        TLeaf* leaf_id  = chain->GetLeaf("id");
        TLeaf* leaf_m   = chain->GetLeaf("e");	// "e" branch holds the particle mass
        TLeaf* leaf_ist = chain->GetLeaf("ist");

        // collect POI transverse-momentum components for this event
        std::vector<float> vpx, vpy;
        vpx.reserve(256); vpy.reserve(256);

        for(int trk = 0; trk < NPTracks; trk++) {
            if((int)leaf_ist->GetValue(trk) != 0) continue;	//final state
            int pid = (int)leaf_id->GetValue(trk);
            if(!(pid == ID_pip || pid == ID_pim)) continue;

            float px = leaf_px->GetValue(trk);
            float py = leaf_py->GetValue(trk);
            float pz = leaf_pz->GetValue(trk);
            float m  = leaf_m->GetValue(trk);

            float pt = sqrt(px*px+py*py);
            if(pt < ptMin || pt > ptMax) continue;
            float E = sqrt(m*m+px*px+py*py+pz*pz);
            float y = 0.5*log((E+pz)/(E-pz));
            if(fabs(y) > yCut) continue;

            vpx.push_back(px);
            vpy.push_back(py);
        }

        int N = (int)vpx.size();
        if(N < 2) continue;	//need at least one pair

        // single-particle Q-vector and v2 sum
        double Cs = 0, Ss = 0, vsum_s = 0;
        for(int a = 0; a < N; a++){
            float phi = atan2(vpy[a], vpx[a]);
            double c = cos(2.*phi), s = sin(2.*phi);
            Cs += c; Ss += s;
            vsum_s += c*cos2Psi + s*sin2Psi;	// cos(2phi - 2Psi_RP)
        }
        double Q2s = Cs*Cs + Ss*Ss;

        // pair Q-vector and pair-v2 sum: phi_p = angle of (p_a + p_b)
        double Cp = 0, Sp = 0, vsum_p = 0;
        double Npair = 0.5*(double)N*(N-1);
        for(int a = 0; a < N; a++){
            for(int b = a+1; b < N; b++){
                float phip = atan2(vpy[a]+vpy[b], vpx[a]+vpx[b]);
                double c = cos(2.*phip), s = sin(2.*phip);
                Cp += c; Sp += s;
                vsum_p += c*cos2Psi + s*sin2Psi;	// cos(2phi_p - 2Psi_RP)
            }
        }
        double Q2p = Cp*Cp + Sp*Sp;

        // accumulate two-particle cumulant numerators/denominators
        num_s2 += (Q2s - N);          den_s2 += (double)N*(N-1);
        num_p2 += (Q2p - Npair);      den_p2 += Npair*(Npair-1);

        // store per-event scalars for pass 2
        ev_N.push_back(N);
        ev_Q2s.push_back(Q2s);
        ev_vsum_s.push_back(vsum_s);
        ev_Npair.push_back(Npair);
        ev_Q2p.push_back(Q2p);
        ev_vsum_p.push_back(vsum_p);
    }

    double v2_2   = sqrt(num_s2/den_s2 > 0 ? num_s2/den_s2 : 0);	// v2{2}
    double v2p_2  = sqrt(num_p2/den_p2 > 0 ? num_p2/den_p2 : 0);	// v2pair{2}
    double v2_2sq  = v2_2*v2_2;
    double v2p_2sq = v2p_2*v2p_2;

    cout << "\n==== ensemble two-particle cumulants ====\n";
    cout << "v2{2}      = " << v2_2  << "\n";
    cout << "v2pair{2}  = " << v2p_2 << "\n";
    cout << "N events used = " << ev_N.size() << "\n\n";

    //================= profiles: v2 (single/pair) vs q2^2 (single/pair) =================
    // naming pXY: X = q2^2 type on x-axis (S=single,P=pair), Y = v2 type on y-axis
    TProfile* pSS = new TProfile("pV2s_q2s","single v_{2} vs single q_{2}^{2}",40,0,4);
    TProfile* pSP = new TProfile("pV2p_q2s","pair v_{2} vs single q_{2}^{2}",  40,0,4);
    TProfile* pPS = new TProfile("pV2s_q2p","single v_{2} vs pair q_{2}^{2}",  40,0,4);
    TProfile* pPP = new TProfile("pV2p_q2p","pair v_{2} vs pair q_{2}^{2}",    40,0,4);

    // QA: distributions of the event-shape variables
    TH1D* hq2s = new TH1D("hq2s","single q_{2}^{2}",100,0,6);
    TH1D* hq2p = new TH1D("hq2p","pair q_{2}^{2}",  100,0,6);

    // store the constants so the figure macro can label them
    TH1D* hConst = new TH1D("hConst","ESS constants (bin1: v2{2}, bin2: v2pair{2})",2,0.5,2.5);
    hConst->SetBinContent(1, v2_2);
    hConst->SetBinContent(2, v2p_2);

    for(size_t e = 0; e < ev_N.size(); e++){
        int    N     = ev_N[e];
        double Npair = ev_Npair[e];

        double q2s = ev_Q2s[e] / ( N     * (1. + N     * v2_2sq ) );	// Eq. (9)
        double q2p = ev_Q2p[e] / ( Npair * (1. + Npair * v2p_2sq) );	// Eq. (11)

        double v2s = ev_vsum_s[e] / N;		// single v2 = <cos(2phi   - 2Psi_RP)>
        double v2p = ev_vsum_p[e] / Npair;	// pair   v2 = <cos(2phi_p - 2Psi_RP)>  Eq. (10)

        pSS->Fill(q2s, v2s);
        pSP->Fill(q2s, v2p);
        pPS->Fill(q2p, v2s);
        pPP->Fill(q2p, v2p);

        hq2s->Fill(q2s);
        hq2p->Fill(q2p);
    }

    fout.Write();
    return;
}
