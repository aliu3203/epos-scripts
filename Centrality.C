#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
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

const int Nbin = 50;
void Centrality(int job=5) {

TH1D* Hist_b = new TH1D("Hist_b","Hist_b",Nbin,0,20);


TChain *chain = new TChain("teposevent");
char fname[200],fname2[200];
chain->Add("z-auau_run_*.root");
//chain->Add("../auau200-6/z-*.root");

int Nevent = chain->GetEntries();
cout << Nevent;
        for(int i=0;i<Nevent;i++) {
                if((i+1)%1000==0) cout<<"Processing entry == "<< i+1 <<" == out of "<<Nevent<<".\n";
                chain->GetEntry(i);
                TLeaf* leaf_b   = chain->GetLeaf("bim");
		TLeaf* leaf_Np_p= chain->GetLeaf("npartproj");
                TLeaf* leaf_Np_t= chain->GetLeaf("nparttarg");
		TLeaf* leaf_Mult= chain->GetLeaf("np");
                float b = leaf_b->GetValue(0);
                if(false && (b < 8.4 || b > 9.2)){
                        continue;
                }
		int Np = leaf_Np_p->GetValue(0) + leaf_Np_t->GetValue(0);
		if(Np<3) continue;
		int Mult = leaf_Mult->GetValue(0);
                Hist_b->Fill(b);
/*
		TLeaf* leaf_PID   = chain->GetLeaf("PID");
		int netCharge = 0;
          	for(int trk = 0; trk < Mult; trk++) {
                	int PID   = leaf_PID->GetValue(trk);
			if(fabs(PID)==531||fabs(PID)==533||fabs(PID)==513||fabs(PID)==511||fabs(PID)==4312||fabs(PID)==4332||fabs(PID)==441||fabs(PID)==4132||fabs(PID)==4112||fabs(PID)==421||fabs(PID)==423||fabs(PID)==3322||fabs(PID)==313||fabs(PID)==423||fabs(PID)==2112||fabs(PID)==111||fabs(PID)==22||fabs(PID)==311||fabs(PID)==3122) continue;
			if(PID==211) netCharge += 1;   //pion+
			else if(PID==-211) netCharge += -1;
			else if(PID==321) netCharge += 1;   //kaon+
                        else if(PID==-321) netCharge += -1;   
                        else if(PID==2212) netCharge += 1;   //proton
                        else if(PID==-2212) netCharge += -1;
                        else if(PID==11||PID==13) netCharge += -1;  //e- ,mu- 
                        else if(PID==-11||PID==-13) netCharge += 1;  
                        else if(PID==42) netCharge += 1;  //deuteron
                        else if(PID==-42) netCharge += -1;
			else if(PID==-3222||PID==3112||PID==3114||PID==3312||PID==3314||PID==3334) netCharge += -1;
                        else if(PID==3222||PID==-3112||PID==-3114||PID==-3312||PID==-3314||PID==-3334) netCharge += 1;
			else if(PID==411||PID==413||PID==431 ||PID==433) netCharge += 1;  //D+
                        else if(PID==-411||PID==-413||PID==-431 ||PID==-433) netCharge += -1;
                        else if(PID==4122||PID==4232||PID==4212||PID==4322||PID==5222||PID==-5132||PID==4412) netCharge += 1;  //Lambda_c+
                        else if(PID==-4122||PID==-4232||PID==-4212||PID==-4322||PID==-5222||PID==5132||PID==-4412) netCharge +=-1;  
                        else if(PID==4222||PID==4422) netCharge += 2; //cascade_c++
                        else if(PID==-4222||PID==-4422) netCharge +=- 2;
                        else if(PID==521||PID==523) netCharge += 1; //B0
                        else if(PID==-521||PID==-523) netCharge +=-1;
			else cout<<"PID = "<<PID<<",";
		}
		if(netCharge != 2*79) {cout<<endl;cout<<"net charge = "<<netCharge - 2*79<<endl;}
*/
        }

        float sum = (float)Hist_b->Integral();
        int j=1;
        for(int i=0;i<Nbin;i++) {
                float add = (float)Hist_b->Integral(1,i+1);
                float center = Hist_b->GetBinCenter(i+1);
                center += float(20)/Nbin/2;
                if(add > sum*0.05*j) {cout<<0.05*j<<"  "<<center<<"  "<<(add-sum*0.05*j)/(sum*0.05*j)<<endl;j++;}


        }


char foname[200];
sprintf(foname,"cen0.b.root");
TFile *fout = new TFile(foname,"recreate");
Hist_b->Write();
fout->Close();
}
