// Standard deviation of the PER-PAIR trig term distributions for gamma112, gamma132
// and delta.
//
// Every correlator is an average of one cosine per unordered pair.  This macro puts
// each individual term into a TProfile as it is computed -- one Fill per pair, not one
// Fill per event -- and then reads back the spread of that distribution with
// SetErrorOption("s").  No event averaging, no covariance modelling: just the standard
// deviation of the numbers that go into the sum, and the resulting std/sqrt(Nterms).
//
// Pairs come from the same explicit i<j loop over all N accepted tracks that
// Correlators_Cent.C uses (PAIRLOOP), so the term set is identical.
//
// run: root -l -b -q 'TermSpread.C+(20000,2,20)'
using namespace std;
#include "stdio.h"
#include "TFile.h"
#include <iostream>
#include <vector>
#include <TChain.h>
#include "TLeaf.h"
#include "TMath.h"
#include "TProfile.h"
#include "EposPID.h"

const float MPI_CHG = 0.13957;
const float yCut = 1.0, ptMin = 0.2, ptMax = 2.0;

// 0/1 = OS/SS for each of the four term types
// gamma132's pair term is not symmetric in (a,b), so a pair has two candidate values.
//   _sym    = contribute their average          (shrinks the spread by 1/sqrt(2))
//   _lowidx = contribute the a=i one, always    (plain "one term per pair")
//   _chgord = contribute the a=positive one for OS, a=i for SS  (matches the Q-vector form)
enum { T112=0, T132S=1, T132O=2, T132L=3, TDEL=4, NTYPE=5 };
static const char* typeName[NTYPE] = {"gamma112", "gamma132_sym", "gamma132_chgord",
                                      "gamma132_lowidx", "delta"};

void TermSpread(int nev = 20000, int poiMode = 2, int nfile = 20,
                const char* outname = "cen.termspread.root"){

    if(!EposPID::Load("idt.dt")) return;

    TChain* chain = new TChain("teposevent");
    for(int k=0;k<nfile;k++){
        int idx = 1 + (int)((665.0*k)/nfile);
        chain->Add(Form("z-auau_run_%d.root", idx));
    }
    chain->SetBranchStatus("*",0);
    const char* act[] = {"np","phi","px","py","pz","id","ist"};
    for(int i=0;i<7;i++) chain->SetBranchStatus(act[i],1);

    TFile fout(outname,"RECREATE");

    // one profile per (term type, charge combination); each is filled ONCE PER PAIR
    TProfile* pT[NTYPE][2];
    for(int t=0;t<NTYPE;t++) for(int k=0;k<2;k++)
        pT[t][k] = new TProfile(Form("pTerm_%s_%s", typeName[t], k==0?"OS":"SS"),
                                Form("per-pair term, %s %s", typeName[t], k==0?"OS":"SS"),
                                1,0,1);

    Long64_t nentries = chain->GetEntries();
    if(nev > nentries) nev = nentries;
    printf("%lld entries in chain (%d files), using %d\n", nentries, nfile, nev);

    vector<double> tc1, ts1, tc3, ts3; vector<int> tq;
    long nUsed = 0;
    double nTermTot[NTYPE][2]; memset(nTermTot,0,sizeof(nTermTot));

    for(Long64_t i=0;i<nev;i++){
        if((i+1)%2000==0) cout<<"entry "<<i+1<<" / "<<nev<<"\n";
        chain->GetEntry(i);
        int   NPTracks = (int)chain->GetLeaf("np")->GetValue(0);
        float PsiRP    = chain->GetLeaf("phi")->GetValue(0);
        double C2 = cos(2.*PsiRP), S2 = sin(2.*PsiRP);

        TLeaf* lpx=chain->GetLeaf("px"); TLeaf* lpy=chain->GetLeaf("py");
        TLeaf* lpz=chain->GetLeaf("pz"); TLeaf* lid=chain->GetLeaf("id");
        TLeaf* lst=chain->GetLeaf("ist");

        tc1.clear(); ts1.clear(); tc3.clear(); ts3.clear(); tq.clear();
        for(int trk=0; trk<NPTracks; trk++){
            if((int)lst->GetValue(trk) != 0) continue;
            int pid = (int)lid->GetValue(trk);
            int chg = EposPID::Charge(pid);
            if(chg == 0) continue;
            float px=lpx->GetValue(trk), py=lpy->GetValue(trk), pz=lpz->GetValue(trk);
            float pt2 = px*px+py*py;
            int apid = abs(pid);
            bool isPOI = (poiMode==0)? (apid==120) : (poiMode==2)? (apid==120||apid==130) : true;
            if(!isPOI) continue;
            float pt = sqrt(pt2);
            if(pt < ptMin || pt > ptMax) continue;
            float m  = (poiMode==0) ? MPI_CHG : EposPID::Mass(pid);
            float E  = sqrt(m*m + pt2 + pz*pz);
            float y  = 0.5*log((E+pz)/(E-pz));
            if(fabs(y) > yCut) continue;
            float phi = atan2(py,px);
            tc1.push_back(cos(phi));    ts1.push_back(sin(phi));
            tc3.push_back(cos(3.*phi)); ts3.push_back(sin(3.*phi));
            tq.push_back(chg>0 ? +1 : -1);
        }

        int Np=0,Nm=0;
        for(size_t a=0;a<tq.size();a++){ if(tq[a]>0) Np++; else Nm++; }
        if(Np < 2 || Nm < 2) continue;

        const int N = (int)tc1.size();
        for(int a=0;a<N;a++){
            const double ci=tc1[a], si=ts1[a], Ci=tc3[a], Si=ts3[a];
            const int    qi=tq[a];
            for(int b=a+1;b<N;b++){
                const double cj=tc1[b], sj=ts1[b], Cj=tc3[b], Sj=ts3[b];
                const int k = (qi*tq[b] < 0) ? 0 : 1;

                double Rp = ci*cj - si*sj, Ip = si*cj + ci*sj;
                double t112 = Rp*C2 + Ip*S2;              // cos(pa + pb - 2Psi)
                double tdel = ci*cj + si*sj;              // cos(pa - pb)

                double Xij = ci*Cj + si*Sj, Yij = si*Cj - ci*Sj;
                double tij = Xij*C2 - Yij*S2;             // cos(pa - 3pb + 2Psi)
                double Xji = cj*Ci + sj*Si, Yji = sj*Ci - cj*Si;
                double tji = Xji*C2 - Yji*S2;             // reversed ordering

                double t132s = 0.5*(tij + tji);
                double t132o = (k==0) ? ((qi>0)? tij : tji) : tij;

                // ONE Fill per pair per term type -- this is the distribution whose
                // standard deviation is being asked for
                pT[T112 ][k]->Fill(0.5, t112 );
                pT[T132S][k]->Fill(0.5, t132s);
                pT[T132O][k]->Fill(0.5, t132o);
                pT[T132L][k]->Fill(0.5, tij  );   // a = i, no charge rule at all
                pT[TDEL ][k]->Fill(0.5, tdel );
            }
        }
        nUsed++;
    }

    //------------------- report -------------------
    printf("\n================ per-pair term distributions ================\n");
    printf("events used = %ld\n\n", nUsed);
    printf("%-16s %-4s %14s %14s %13s %13s\n",
           "term", "chg", "N terms", "mean", "std dev", "std/sqrt(N)");

    double mean[NTYPE][2], sdev[NTYPE][2], nter[NTYPE][2];
    for(int t=0;t<NTYPE;t++) for(int k=0;k<2;k++){
        TProfile* p = pT[t][k];
        nter[t][k] = p->GetBinEntries(1);
        mean[t][k] = p->GetBinContent(1);
        TProfile* q = (TProfile*)p->Clone(Form("%s_s", p->GetName()));
        q->SetErrorOption("s");                       // -> GetBinError is the std dev
        sdev[t][k] = q->GetBinError(1);
        nTermTot[t][k] = nter[t][k];
        printf("%-16s %-4s %14.0f %+14.6e %13.6f %13.4e\n",
               typeName[t], k==0?"OS":"SS", nter[t][k], mean[t][k], sdev[t][k],
               nter[t][k]>0 ? sdev[t][k]/sqrt(nter[t][k]) : 0.0);
    }

    printf("\n---- Delta = OS - SS, errors from the per-term std devs ----\n");
    printf("%-16s %14s %13s\n", "term", "Delta", "err (indep terms)");
    for(int t=0;t<NTYPE;t++){
        double d  = mean[t][0]-mean[t][1];
        double e  = sqrt( sdev[t][0]*sdev[t][0]/nter[t][0]
                        + sdev[t][1]*sdev[t][1]/nter[t][1] );
        printf("%-16s %+14.6e %13.4e\n", typeName[t], d, e);
    }

    printf("\n---- std dev ratios, each gamma132 variant over gamma112 ----\n");
    for(int t=1;t<TDEL;t++)
        printf("%-16s / gamma112 :  OS %.5f   SS %.5f\n", typeName[t],
               sdev[t][0]/sdev[T112][0], sdev[t][1]/sdev[T112][1]);

    // scale the per-term errors up to a production-sized sample
    const double NEV_PROD = 676000.0;
    printf("\n---- scaled to %.0f events (terms/event held fixed) ----\n", NEV_PROD);
    printf("%-16s %13s\n", "term", "err(Delta)");
    for(int t=0;t<NTYPE;t++){
        double scale = NEV_PROD/(double)nUsed;
        double e = sqrt( sdev[t][0]*sdev[t][0]/(nter[t][0]*scale)
                       + sdev[t][1]*sdev[t][1]/(nter[t][1]*scale) );
        printf("%-16s %13.4e\n", typeName[t], e);
    }

    fout.Write();
    printf("\nwrote %s\n", outname);
}
