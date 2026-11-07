using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "EposPID.h"

// -----------------------------------------------------------------------------
// Centrality calibration by the reference-multiplicity percentile method.
//
//   1. build the refmult distribution (final-state charged, |eta|<1, pT>0.15)
//      over a MINIMUM-BIAS sample;
//   2. integrate from the high-refmult end: the cut where the integral reaches
//      5% of all events is the 0-5% boundary, 10% gives 5-10%, and so on down
//      to 80%;
//   3. write the 9 thresholds to cen_cuts.txt for Correlators_Cent.C to read.
//
// IMPORTANT: this is only meaningful on a minimum-bias sample (EPOS optns with
// the full impact-parameter range). Run it on a b-restricted sample and the
// percentiles refer to that slice, not to the Au+Au cross section -- the macro
// checks the bim range and refuses to write cuts if it looks b-restricted.
// -----------------------------------------------------------------------------

void CentralityCalib(int job = 0, const char* filepat = "z-auau_run_*.root"){

    if(!EposPID::Load("idt.dt")) return;

    TChain* chain = new TChain("teposevent");
    chain->Add(filepat);
    chain->SetBranchStatus("*",0);
    const char* act[] = {"np","bim","px","py","pz","id","ist","npartproj","nparttarg"};
    for(int i=0;i<9;i++) chain->SetBranchStatus(act[i],1);

    Long64_t nentries = chain->GetEntries();
    cout << nentries << " entries in chain\n";
    if(nentries<=0){ printf("no entries\n"); return; }

    TH1D* hRef   = new TH1D("hRef","refmult (charged, |#eta|<1, p_{T}>0.15);refmult;events",2000,-0.5,1999.5);
    TH1D* hBim   = new TH1D("hBim","impact parameter;b (fm);events",400,0,20);
    TProfile* pNpartVsRef = new TProfile("pNpartVsRef","<N_{part}> vs refmult;refmult;N_{part}",200,0,2000);
    TProfile* pBimVsRef   = new TProfile("pBimVsRef","<b> vs refmult;refmult;b (fm)",200,0,2000);

    double bmin=1e9, bmax=-1e9;

    for(Long64_t i=0;i<nentries;i++){
        if((i+1)%10000==0) cout<<"Processing entry == "<<i+1<<" == out of "<<nentries<<".\n";
        chain->GetEntry(i);

        int   NPTracks = (int)chain->GetLeaf("np")->GetValue(0);
        float bim      = chain->GetLeaf("bim")->GetValue(0);
        int   Npart    = (int)chain->GetLeaf("npartproj")->GetValue(0)
                       + (int)chain->GetLeaf("nparttarg")->GetValue(0);
        if(bim<bmin) bmin=bim;
        if(bim>bmax) bmax=bim;

        TLeaf* lpx=chain->GetLeaf("px"); TLeaf* lpy=chain->GetLeaf("py");
        TLeaf* lpz=chain->GetLeaf("pz"); TLeaf* lid=chain->GetLeaf("id");
        TLeaf* lst=chain->GetLeaf("ist");

        int refmult = 0;
        for(int trk=0; trk<NPTracks; trk++){
            if((int)lst->GetValue(trk) != 0) continue;              // final state only
            if(EposPID::Charge((int)lid->GetValue(trk)) == 0) continue;
            float px=lpx->GetValue(trk), py=lpy->GetValue(trk), pz=lpz->GetValue(trk);
            float pt2 = px*px+py*py;
            if(pt2 < REFMULT_PTMIN*REFMULT_PTMIN) continue;
            if(pz*pz >= SINH_ETACUT_SQ*pt2) continue;               // |eta|<1 without log/atan
            refmult++;
        }

        hRef->Fill(refmult);
        hBim->Fill(bim);
        pNpartVsRef->Fill(refmult, Npart);
        pBimVsRef  ->Fill(refmult, bim);
    }

    // ---------------- percentile cuts, integrating from the high end ----------------
    double total = hRef->Integral(0, hRef->GetNbinsX()+1);
    printf("\n=========== centrality calibration ===========\n");
    printf("events = %.0f   bim range in sample = [%.3f, %.3f] fm\n", total, bmin, bmax);

    bool minbias = (bmin < 1.0 && bmax > 13.0);
    if(!minbias){
        printf("\n  *** WARNING: this sample is NOT minimum bias (b in [%.2f,%.2f] fm). ***\n", bmin, bmax);
        printf("  *** Percentile cuts derived here would be meaningless, so none are written. ***\n");
        printf("  *** Regenerate with the full b range (see notes in the message/CLAUDE.md). ***\n");
    }

    int cut[NCENT];
    double frac[NCENT];
    for(int ic=0; ic<NCENT; ic++){
        double target = CENT_HI[ic]/100.0 * total;   // integral from the top
        double run = 0; int thr = 0;
        for(int b=hRef->GetNbinsX()+1; b>=0; b--){
            run += hRef->GetBinContent(b);
            if(run >= target){ thr = (int)hRef->GetBinLowEdge(b); break; }
        }
        cut[ic]  = thr;
        frac[ic] = run/total*100.0;
    }

    printf("\n  %-9s %14s %12s %12s %10s\n","centrality","refmult range","achieved %","<N_part>","<b> fm");
    for(int ic=0; ic<NCENT; ic++){
        // mean Npart / b for events in this bin
        double sN=0,sB=0,sW=0;
        int lo = cut[ic], hi = (ic==0)? 100000 : cut[ic-1]-1;
        for(int b=1;b<=pNpartVsRef->GetNbinsX();b++){
            double x = pNpartVsRef->GetBinCenter(b);
            if(x<lo || x>hi) continue;
            double w = pNpartVsRef->GetBinEntries(b);
            sN += pNpartVsRef->GetBinContent(b)*w;
            sB += pBimVsRef  ->GetBinContent(b)*w;
            sW += w;
        }
        char rng[32];
        if(ic==0) snprintf(rng,32,">= %d",cut[ic]);
        else      snprintf(rng,32,"%d - %d",cut[ic],cut[ic-1]-1);
        printf("  %-9s %14s %11.2f %12.1f %10.2f\n", CentLabel(ic), rng, frac[ic],
               sW>0? sN/sW:0.0, sW>0? sB/sW:0.0);
    }

    if(minbias){
        ofstream out("cen_cuts.txt");
        out << "# refmult lower thresholds, most central first (0-5%,5-10%,...,70-80%)\n";
        out << "# generated by CentralityCalib.C from " << (long long)total << " min-bias events\n";
        for(int ic=0; ic<NCENT; ic++) out << cut[ic] << "\n";
        out.close();
        printf("\nwrote cen_cuts.txt\n");
    }

    char fo[200]; sprintf(fo,"cen_calib_job%d.root",job);
    TFile fout(fo,"RECREATE");
    hRef->Write(); hBim->Write(); pNpartVsRef->Write(); pBimVsRef->Write();
    fout.Close();
    printf("wrote %s\n", fo);
}
