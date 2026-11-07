// Pair-count audit for Correlators_Cent.C: how many pairs does each correlator
// actually average over, on the PRODUCTION selection (poiMode=2, all files)?
// Prints the mean and RMS of Np, Nm and every denominator the macro divides by.
using namespace std;
#include "stdio.h"
#include <iostream>
#include <TChain.h>
#include "TLeaf.h"
#include "TMath.h"
#include "EposPID.h"

const float MPI_CHG = 0.13957;
const float yCut = 1.0, ptMin = 0.2, ptMax = 2.0;

void PairCount(int nev = 20000, int poiMode = 2, int nfile = 20){
    if(!EposPID::Load("idt.dt")) return;
    // Contiguous reads from files spread over the whole production range -- striding
    // one entry at a time across 665 files thrashes the disk for no extra information,
    // since every file covers the same b range.
    TChain* chain = new TChain("teposevent");
    for(int k=0;k<nfile;k++){
        int idx = 1 + (int)((665.0*k)/nfile);
        chain->Add(Form("z-auau_run_%d.root", idx));
    }
    chain->SetBranchStatus("*",0);
    const char* act[] = {"np","phi","px","py","pz","id","ist"};
    for(int i=0;i<7;i++) chain->SetBranchStatus(act[i],1);

    Long64_t nentries = chain->GetEntries();
    printf("%lld entries in chain (%d files), sampling %d\n", nentries, nfile, nev);
    if(nev > nentries) nev = nentries;
    Long64_t step = 1;

    const int NQ = 8;
    const char* nm[NQ] = {"Np (+)","Nm (-)","N tot","NpOS unord","NpSS unord",
                          "NoOS ord","NoSS ord","NpOS+NpSS"};
    double s[NQ]={0}, s2[NQ]={0};
    long n=0;
    for(Long64_t i=0;i<nentries && n<nev;i+=step){
        chain->GetEntry(i);
        int   NPTracks = (int)chain->GetLeaf("np")->GetValue(0);
        TLeaf* lpx=chain->GetLeaf("px"); TLeaf* lpy=chain->GetLeaf("py");
        TLeaf* lpz=chain->GetLeaf("pz"); TLeaf* lid=chain->GetLeaf("id");
        TLeaf* lst=chain->GetLeaf("ist");
        int Np=0, Nm=0;
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
            if(chg > 0) Np++; else Nm++;
        }
        if(Np < 2 || Nm < 2) continue;
        double NpOS = (double)Np*Nm;
        double NpSS = 0.5*Np*(Np-1) + 0.5*Nm*(Nm-1);
        double v[NQ] = {(double)Np,(double)Nm,(double)(Np+Nm),NpOS,NpSS,
                        2.0*Np*Nm,(double)Np*(Np-1)+(double)Nm*(Nm-1),NpOS+NpSS};
        for(int k=0;k<NQ;k++){ s[k]+=v[k]; s2[k]+=v[k]*v[k]; }
        n++;
    }
    printf("\nevents sampled = %ld\n\n", n);
    printf("%-14s %14s %14s\n","quantity","mean","RMS");
    for(int k=0;k<NQ;k++)
        printf("%-14s %14.1f %14.1f\n", nm[k], s[k]/n, sqrt(s2[k]/n-(s[k]/n)*(s[k]/n)));
    double mp=s[0]/n, mm=s[1]/n, mos=s[3]/n, mss=s[4]/n;
    printf("\nSS - OS = %+.1f   (equals Var(N) - <N> per charge if Np,Nm independent)\n", mss-mos);
    printf("ratio SS/OS = %.4f\n", mss/mos);
    printf("check: <Np><Nm> = %.1f  vs  <Np*Nm> = %.1f  (differ if Np,Nm correlated)\n", mp*mm, mos);
}
