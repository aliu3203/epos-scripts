using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TTree.h"
#include "TMath.h"

// Count events in the 30-40% centrality window (b = 8.4-9.2 fm) and
// diagnose where events are "lost" when the b-cut is applied.
//
// Checks performed per event:
//   1. float b in [8.4, 9.2]          -- the correct cut
//   2. int   b in [8.4, 9.2]          -- the truncation bug in analyze.C
//      (int bim = leaf_b->GetValue(0); 8.7 -> 8 -> fails the cut)
//   3. events outside the window are printed with entry number, file, and b

const float bMin = 8.4;
const float bMax = 9.2;

void CountCentrality(int cen = 5, int job = 0){

    TChain* chain = new TChain("teposevent");
    chain->Add("z-auau_run_*.root");

    //only read the branch we use -- everything else stays compressed on disk
    chain->SetBranchStatus("*",0);
    chain->SetBranchStatus("bim",1);

    Long64_t nentries = chain->GetEntries();
    cout << nentries << " total entries in chain\n";

    Long64_t nInFloat = 0;      // correct float comparison
    Long64_t nInInt   = 0;      // buggy int-truncated comparison
    Long64_t nBelow   = 0;      // b < 8.4
    Long64_t nAbove   = 0;      // b > 9.2
    Long64_t nEdge    = 0;      // outside but within 1e-3 fm of a boundary
    float bLo =  1e9, bHi = -1e9;

    TH1D* Hist_b = new TH1D("Hist_b","b, all events",200,8.0,9.6);

    for(Long64_t i = 0; i < nentries; i++){

        if((i+1)%1000==0) cout << "Processing entry == "<< i+1 <<" == out of "<<nentries<<".\n";
        chain->GetEntry(i);

        TLeaf* leaf_b = chain->GetLeaf("bim");
        float b = leaf_b->GetValue(0);

        Hist_b->Fill(b);
        if(b < bLo) bLo = b;
        if(b > bHi) bHi = b;

        //1. correct cut
        bool inFloat = !(b < bMin || b > bMax);
        if(inFloat) nInFloat++;

        //2. the analyze.C truncation bug: int bim = leaf_b->GetValue(0);
        int bInt = leaf_b->GetValue(0);
        bool inInt = !(bInt < bMin || bInt > bMax);
        if(inInt) nInInt++;

        //3. report anything genuinely outside the window
        if(!inFloat){
            if(b < bMin) nBelow++;
            else         nAbove++;
            if(fabs(b - bMin) < 1e-3 || fabs(b - bMax) < 1e-3) nEdge++;
            cout << "OUTSIDE: entry " << i
                 << "  b = " << b
                 << "  file = " << chain->GetFile()->GetName() << "\n";
        }
    }

    cout << "\n==== centrality count, b in [" << bMin << ", " << bMax << "] fm ====\n";
    cout << "total events              = " << nentries << "\n";
    cout << "in window (float cut)     = " << nInFloat << "\n";
    cout << "in window (int cut, bug)  = " << nInInt
         << "   <-- what analyze.C's 'int bim' would keep\n";
    cout << "below " << bMin << " fm            = " << nBelow << "\n";
    cout << "above " << bMax << " fm            = " << nAbove << "\n";
    cout << "of those, within 1e-3 of a boundary = " << nEdge
         << "  (float rounding of the generator limits)\n";
    cout << "b range in data           = [" << bLo << ", " << bHi << "]\n";

    char foname[200];
    sprintf(foname,"cen%d.bcount_job%d.root",cen,job);
    TFile* fout = new TFile(foname,"recreate");
    Hist_b->Write();
    fout->Close();
}
