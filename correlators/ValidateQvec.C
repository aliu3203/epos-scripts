// Brute-force validation of the Q-vector algebra used in Correlators_Cent.C for
// gamma112 / gamma132 / delta, plus a direct test of the gamma132 pair-ordering
// convention.
//
// For a subset of events every correlator is computed twice: (a) exactly as the
// production macro does it, from Q1/Q2/Q3, and (b) from an explicit double loop over
// pairs.  Any disagreement in the per-event value, or in the number of pairs used, is
// a bug in the production macro.
//
// The macro also builds gamma132 THREE ways over the IDENTICAL unordered pair set
// that gamma112 uses (so the pair counts are equal by construction):
//   132_sym  : (1/2)[cos(pa-3pb+2Psi) + cos(pb-3pa+2Psi)]   <- what production uses
//   132_1ord : cos(pa-3pb+2Psi) only, alpha = (+) for OS, lower index for SS
//   132_idty : cos(pa+pb-2Psi) * cos(2(pa-pb))              <- algebraic identity
// 132_sym and 132_idty must agree term by term; 132_1ord must agree in the MEAN but
// carries the extra antisymmetric zero-mean term.  Comparing the per-pair RMS of these
// against gamma112's is what explains the error ratio between Dgamma112 and Dgamma132.
using namespace std;
#include "stdio.h"
#include "TFile.h"
#include <iostream>
#include <vector>
#include <TChain.h>
#include "TLeaf.h"
#include "TMath.h"
#include "EposPID.h"

const float MPI_CHG = 0.13957;
const float yCut = 1.0, ptMin = 0.2, ptMax = 2.0;

void ValidateQvec(int nev = 2000, int poiMode = 2, const char* filepat = "z-auau_run_1.root"){

    if(!EposPID::Load("idt.dt")) return;

    TChain* chain = new TChain("teposevent");
    chain->Add(filepat);
    chain->SetBranchStatus("*",0);
    const char* act[] = {"np","bim","phi","px","py","pz","id","ist"};
    for(int i=0;i<8;i++) chain->SetBranchStatus(act[i],1);

    Long64_t nentries = chain->GetEntries();
    if(nev > nentries) nev = nentries;
    printf("validating on %d of %lld events, poiMode=%d\n", nev, nentries, poiMode);

    // running means, Q-vector version vs brute-force version
    double sQ[7]={0}, sB[7]={0};   // 0:g112OS 1:g112SS 2:g132OS 3:g132SS 4:dOS 5:dSS 6:v2
    double sQ2[2]={0}, sB2[2]={0}; // sum of squares of Dg112, Dg132 (event spread)
    long   nUsed=0;
    double worst[7]={0};
    long   nPairMismatch=0;
    // pair-count bookkeeping: PAIRS (not cosine terms) each correlator averages over
    double totPairs[6]={0};

    // ---- per-PAIR statistics, accumulated over every pair of every event ----
    // these decide the error ratio between the correlators
    double pSum[5]={0}, pSum2[5]={0}; long pN[5]={0};   // 0:112 1:132sym 2:132_1ord 3:delta 4:132idty
    double maxIdentityDiff = 0;
    // single-ordering event averages, over the SAME pair set as gamma112
    double s1ordOS=0, s1ordSS=0, s1ord2=0;

    for(Long64_t i=0;i<nev;i++){
        chain->GetEntry(i);
        int   NPTracks = (int)chain->GetLeaf("np")->GetValue(0);
        float PsiRP    = chain->GetLeaf("phi")->GetValue(0);
        double C2 = cos(2.*PsiRP), S2 = sin(2.*PsiRP);

        TLeaf* lpx=chain->GetLeaf("px"); TLeaf* lpy=chain->GetLeaf("py");
        TLeaf* lpz=chain->GetLeaf("pz"); TLeaf* lid=chain->GetLeaf("id");
        TLeaf* lst=chain->GetLeaf("ist");

        double pP=0,qP=0,p2P=0,q2P=0,p3P=0,q3P=0; int Np=0;
        double pM=0,qM=0,p2M=0,q2M=0,p3M=0,q3M=0; int Nm=0;
        vector<double> phiP, phiM;

        for(int trk=0; trk<NPTracks; trk++){
            if((int)lst->GetValue(trk) != 0) continue;
            int pid = (int)lid->GetValue(trk);
            int chg = EposPID::Charge(pid);
            if(chg == 0) continue;
            float px=lpx->GetValue(trk), py=lpy->GetValue(trk), pz=lpz->GetValue(trk);
            float pt2 = px*px+py*py;
            int apid = abs(pid);
            bool isPOI;
            if      (poiMode==0) isPOI = (apid==120);
            else if (poiMode==2) isPOI = (apid==120 || apid==130);
            else                 isPOI = true;
            if(!isPOI) continue;
            float pt = sqrt(pt2);
            if(pt < ptMin || pt > ptMax) continue;
            float m  = (poiMode==0) ? MPI_CHG : EposPID::Mass(pid);
            float E  = sqrt(m*m + pt2 + pz*pz);
            float y  = 0.5*log((E+pz)/(E-pz));
            if(fabs(y) > yCut) continue;
            float phi = atan2(py,px);
            double c1=cos(phi),s1=sin(phi),c2=cos(2.*phi),s2=sin(2.*phi),c3=cos(3.*phi),s3=sin(3.*phi);
            if(chg > 0){ pP+=c1; qP+=s1; p2P+=c2; q2P+=s2; p3P+=c3; q3P+=s3; Np++; phiP.push_back(phi); }
            else       { pM+=c1; qM+=s1; p2M+=c2; q2M+=s2; p3M+=c3; q3M+=s3; Nm++; phiM.push_back(phi); }
        }
        if(Np < 2 || Nm < 2) continue;

        //================= Q-vector version (verbatim from Correlators_Cent.C) =========
        double NpOS = (double)Np*Nm;
        double NpSS = 0.5*Np*(Np-1) + 0.5*Nm*(Nm-1);
        double NoOS = 2.0*Np*Nm;
        double NoSS = (double)Np*(Np-1) + (double)Nm*(Nm-1);

        double g112OS = ((pP*pM - qP*qM)*C2 + (pP*qM + qP*pM)*S2) / NpOS;
        double g112SS = ( 0.5*((pP*pP-qP*qP-p2P)*C2 + (2.*pP*qP-q2P)*S2)
                        + 0.5*((pM*pM-qM*qM-p2M)*C2 + (2.*pM*qM-q2M)*S2) ) / NpSS;
        double Xpm = pP*p3M + qP*q3M,  Ypm = qP*p3M - pP*q3M;
        double Xmp = pM*p3P + qM*q3P,  Ymp = qM*p3P - pM*q3P;
        double g132OS = ((Xpm*C2 - Ypm*S2) + (Xmp*C2 - Ymp*S2)) / NoOS;
        double XppSelf = p2P, YppSelf = -q2P;
        double XmmSelf = p2M, YmmSelf = -q2M;
        double Xpp = pP*p3P + qP*q3P,  Ypp = qP*p3P - pP*q3P;
        double Xmm = pM*p3M + qM*q3M,  Ymm = qM*p3M - pM*q3M;
        double g132SS = ( ((Xpp-XppSelf)*C2 - (Ypp-YppSelf)*S2)
                        + ((Xmm-XmmSelf)*C2 - (Ymm-YmmSelf)*S2) ) / NoSS;
        double dOS = (pP*pM + qP*qM) / NpOS;
        double dSS = ( 0.5*(pP*pP+qP*qP-Np) + 0.5*(pM*pM+qM*qM-Nm) ) / NpSS;
        double v2ev = ((p2P+p2M)*C2 + (q2P+q2M)*S2) / (double)(Np+Nm);

        //================= brute-force pair loops =====================================
        // Every loop below runs over UNORDERED pairs -- the identical pair set for all
        // correlators.  gamma132 simply gets a different function of that same pair.
        double bg112OS=0,bg112SS=0,bg132OS=0,bg132SS=0,bdOS=0,bdSS=0,bv2=0;
        double b1ordOS=0,b1ordSS=0;
        long   nOSpair=0,nSSpair=0;

        for(size_t a=0;a<phiP.size();a++) for(size_t b=0;b<phiM.size();b++){
            double t112 = cos(phiP[a]+phiM[b]-2.*PsiRP);
            double tf   = cos(phiP[a]-3.*phiM[b]+2.*PsiRP);   // alpha = (+)
            double tr   = cos(phiM[b]-3.*phiP[a]+2.*PsiRP);   // alpha = (-)
            double tsym = 0.5*(tf+tr);
            double tid  = t112*cos(2.*(phiP[a]-phiM[b]));     // algebraic identity
            double td   = cos(phiP[a]-phiM[b]);
            if(fabs(tsym-tid)>maxIdentityDiff) maxIdentityDiff=fabs(tsym-tid);
            bg112OS += t112; bg132OS += tsym; b1ordOS += tf; bdOS += td;
            pSum[0]+=t112; pSum2[0]+=t112*t112;
            pSum[1]+=tsym; pSum2[1]+=tsym*tsym;
            pSum[2]+=tf;   pSum2[2]+=tf*tf;
            pSum[3]+=td;   pSum2[3]+=td*td;
            pSum[4]+=tid;  pSum2[4]+=tid*tid;
            pN[0]++; pN[1]++; pN[2]++; pN[3]++; pN[4]++;
            nOSpair++;
        }
        const vector<double>* same[2] = {&phiP,&phiM};
        for(int c=0;c<2;c++){
            const vector<double>& v = *same[c];
            for(size_t a=0;a<v.size();a++) for(size_t b=a+1;b<v.size();b++){
                double t112 = cos(v[a]+v[b]-2.*PsiRP);
                double tf   = cos(v[a]-3.*v[b]+2.*PsiRP);     // alpha = lower index
                double tr   = cos(v[b]-3.*v[a]+2.*PsiRP);
                double tsym = 0.5*(tf+tr);
                double tid  = t112*cos(2.*(v[a]-v[b]));
                if(fabs(tsym-tid)>maxIdentityDiff) maxIdentityDiff=fabs(tsym-tid);
                bg112SS += t112; bg132SS += tsym; b1ordSS += tf; bdSS += cos(v[a]-v[b]);
                nSSpair++;
            }
        }
        for(size_t a=0;a<phiP.size();a++) bv2 += cos(2.*(phiP[a]-PsiRP));
        for(size_t a=0;a<phiM.size();a++) bv2 += cos(2.*(phiM[a]-PsiRP));
        bv2 /= (double)(Np+Nm);

        // do the unordered pair counts match the denominators the macro divides by?
        if(nOSpair!=(long)NpOS || nSSpair!=(long)NpSS) nPairMismatch++;

        bg112OS/=nOSpair; bg112SS/=nSSpair; bg132OS/=nOSpair; bg132SS/=nSSpair;
        b1ordOS/=nOSpair; b1ordSS/=nSSpair; bdOS/=nOSpair; bdSS/=nSSpair;

        double Q[7]={g112OS,g112SS,g132OS,g132SS,dOS,dSS,v2ev};
        double B[7]={bg112OS,bg112SS,bg132OS,bg132SS,bdOS,bdSS,bv2};
        for(int k=0;k<7;k++){
            sQ[k]+=Q[k]; sB[k]+=B[k];
            double d=fabs(Q[k]-B[k]); if(d>worst[k]) worst[k]=d;
        }
        double A1=Q[0]-Q[1], A3=Q[2]-Q[3];
        sQ2[0]+=A1*A1; sQ2[1]+=A3*A3;
        double bA1=B[0]-B[1], bA3=B[2]-B[3];
        sB2[0]+=bA1*bA1; sB2[1]+=bA3*bA3;
        double A31 = b1ordOS-b1ordSS;
        s1ordOS+=b1ordOS; s1ordSS+=b1ordSS; s1ord2+=A31*A31;
        double T[6]={(double)nOSpair,(double)nSSpair,(double)nOSpair,(double)nSSpair,
                     (double)nOSpair,(double)nSSpair};
        for(int k=0;k<6;k++) totPairs[k]+=T[k];
        nUsed++;
    }

    const char* nm[7]={"g112_OS","g112_SS","g132_OS","g132_SS","delta_OS","delta_SS","v2"};
    printf("\nevents used = %ld,  events with a pair-count/denominator mismatch = %ld\n\n",nUsed,nPairMismatch);
    printf("%-10s %14s %14s %12s %12s\n","quantity","Q-vector","brute force","diff","max|ev diff|");
    for(int k=0;k<7;k++)
        printf("%-10s %+14.6e %+14.6e %+12.2e %12.2e\n",nm[k],sQ[k]/nUsed,sB[k]/nUsed,
               (sQ[k]-sB[k])/nUsed,worst[k]);

    double m1=(sQ[0]-sQ[1])/nUsed, m3=(sQ[2]-sQ[3])/nUsed;
    double m31=(s1ordOS-s1ordSS)/nUsed;
    printf("\nDgamma112              = %+.6e\n",m1);
    printf("Dgamma132 (symmetrised)= %+.6e\n",m3);
    printf("Dgamma132 (one order)  = %+.6e\n",m31);
    printf("per-event RMS: Dg112 %.5e   Dg132sym %.5e   Dg132_1ord %.5e\n",
           sqrt(sQ2[0]/nUsed-m1*m1), sqrt(sQ2[1]/nUsed-m3*m3), sqrt(s1ord2/nUsed-m31*m31));
    printf("   RMS ratio  Dg112/Dg132sym  = %.3f\n",
           sqrt(sQ2[0]/nUsed-m1*m1)/sqrt(sQ2[1]/nUsed-m3*m3));
    printf("   RMS ratio  Dg112/Dg132_1ord= %.3f   <-- identical filling\n",
           sqrt(sQ2[0]/nUsed-m1*m1)/sqrt(s1ord2/nUsed-m31*m31));

    const char* tn[6]={"112_OS","112_SS","132_OS","132_SS","del_OS","del_SS"};
    printf("\nmean number of unordered PAIRS averaged per event:\n");
    for(int k=0;k<6;k++) printf("   %-8s %10.1f\n",tn[k],totPairs[k]/nUsed);

    printf("\nmax |sym - identity| over every pair = %.3e   (identity: 132 term = 112 term * cos2(dphi))\n",
           maxIdentityDiff);
    const char* pn[5]={"112 term","132 sym term","132 1-order","delta term","112*cos2dphi"};
    printf("\nper-PAIR statistics over %ld OS pairs (this is what sets the errors):\n",pN[0]);
    printf("   %-14s %12s %12s\n","term","mean","RMS");
    for(int k=0;k<5;k++){
        double mu=pSum[k]/pN[k], rms=sqrt(pSum2[k]/pN[k]-mu*mu);
        printf("   %-14s %+12.5e %12.5f\n",pn[k],mu,rms);
    }
}
