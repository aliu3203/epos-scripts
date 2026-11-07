using namespace std;

#include "stdio.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <string>
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
#include "TVector3.h"
#include "TLorentzVector.h"

const float PI = TMath::Pi();
const int nHar = 2;
const int opt_boost = 0;

void analyze(int cen = 5, int job = 0, const Char_t* inFileName = "test.list"){	//main_function

        

    float Eweight = 1;
    TChain* chain = new TChain("teposevent");

    chain->Add("z-*.root");
    //chain->Add("../auau200-6/z-*.root");
    //std::cout << "Number of events in chain: " << chain->GetEntries() << endl;

    char fname_out[200];
    sprintf(fname_out,"cen%d.AVFD_Parity1%d%d_EPD_pT02_lab_job%d.root",cen,nHar-1,nHar,job);
    if(opt_boost) sprintf(fname_out,"cen%d.AVFD_Parity1%d%d_EPD_pT02_Boost_job%d.root",cen,nHar-1,nHar,job);
    TFile fout(fname_out,"RECREATE");
    //defining variables
    Int_t   Centrality, NPTracks, fParticles_, QB_index;	//
    Float_t psi, b, Charge, Eta, Theta, Phi, Pt, Mass, yy;		//track info	
    Float_t Px, Py, Pz, PID, energy;

    // Maps to store charge and names, keyed by Particle ID
    std::map<int, float> pCharge;
    std::map<int, std::string> pName;
    pCharge[-1120] = -1.0;
    pName[-1120] = "pbar";

    // Read the idt.dt file that contains table of particle ids/charges
    std::ifstream file("idt.dt");
    std::string line;
    while(std::getline(file, line)) {
        // Skip empty lines or comments starting with '!'
        if(line.empty() || line[0] == '!') continue;
        
        std::istringstream iss(line);
        int id_epos, id_pdg, id_qgsjet, id_gheisha, id_sibyll;
        std::string name;
        int ifl1, ifl2, ifl3, counter;
        float mass, charge;
        
        // Parse up to the 12th column (Charge)
        if(iss >> id_epos >> id_pdg >> id_qgsjet >> id_gheisha >> id_sibyll 
            >> name >> ifl1 >> ifl2 >> ifl3 >> counter >> mass >> charge) {
            
            // NOTE: If your PIDAsso uses EPOS IDs rather than PDG IDs, 
            // change `id_pdg` to `id_epos` inside the brackets below.
            pCharge[id_epos] = charge; 
            pName[id_epos] = name;
        }
    }
    TH1F *hSpeciesAll = new TH1F("hSpeciesAll", "All Charged Particles", 1, 0, 1);
    hSpeciesAll->SetCanExtend(TH1::kAllAxes);

    //defining histograms
    TH1D *hCentrality = new TH1D("hCentrality","hCentrality",10,0,10);
    //TH1D *hMult = new T1D("hMult","hMult",1000,-0.5,999.5);
    TH1D *hMult = new TH1D("hMult","hMult",100,-0.5,99.5);

    TH1D *hMult_All = new TH1D("hMult_All", "Multiplicity (All Particles)", 6000, -0.5, 3000); // Revert to 999.5 tmrw
    //TH1D *hPt_All = new TH1D("hPt_All", "Pt (All Particles)", 100, -0.5, 99.5);
    //TH1D *hPhi_All = new TH1D("hPt_All", )

    TH1D *hMult_ptCut = new TH1D("hMult_ptCut", "Multiplicity (pT > 0.05)", 6000, -0.5, 3000);
    TH1D *hMult_ptEtaCut = new TH1D("hMult_ptEtaCut", "Multiplicity (pT > 0.05 & |eta| < 0.5)", 2000, -0.5, 1000);


    TH1D* Hist_Y = new TH1D("Hist_Y", "Hist_Y", 26, -6, 6);
    TH2D *hEtaPtDist = new TH2D("EtaPtDist","EtaPtDist",26, -6, 6,300,0,15);
    TH1D* Hist_Pt = new TH1D("Hist_Pt","Hist_Pt",60,0,3);
    TH1D* Hist_Phi = new TH1D("Hist_Phi","Hist_Phi",72,-PI,PI);

    Int_t nentries = chain->GetEntries();
    cout << nentries << "\n";
    for(int i = 0; i < nentries; i++){

        if((i+1)%1000==0) cout << "Processing entry == "<< i+1 <<" == out of "<<nentries<<".\n";
        chain->GetEntry(i);

        //		TLeaf* leaf_psi = chain->GetLeaf("Psi");
        TLeaf* leaf_b   = chain->GetLeaf("bim");
        float bim = leaf_b->GetValue(0);
        if(false && (bim < 8.4 || bim > 9.2)){
            continue;
        }
        TLeaf* leaf_NoTracks = chain->GetLeaf("np");
        TLeaf* leaf_Np_p= chain->GetLeaf("npartproj");
        TLeaf* leaf_Np_t= chain->GetLeaf("nparttarg");
        int Np = leaf_Np_p->GetValue(0) + leaf_Np_t->GetValue(0);
        NPTracks= (int)leaf_NoTracks->GetValue(0);

        TLeaf* leaf_PxV0       = chain->GetLeaf("px");
        TLeaf* leaf_PyV0       = chain->GetLeaf("py");
        TLeaf* leaf_PzV0       = chain->GetLeaf("pz");
        TLeaf* leaf_PIDV0    = chain->GetLeaf("id");
        TLeaf* leaf_e           = chain->GetLeaf("e");

        TLeaf* leaf_ist        = chain->GetLeaf("ist");

        int mult_all = 0;
        int mult_ptCut = 0;
        int mult_ptEtaCut = 0;
        for(int trk = 0; trk < NPTracks; trk++) {

            int istAsso = (int)leaf_ist->GetValue(trk);
            if (istAsso != 0) continue;

            float PxAsso    = leaf_PxV0->GetValue(trk);
            float PyAsso    = leaf_PyV0->GetValue(trk); //if((i+1)%2==0) PyAsso *= -1;
            float PzAsso    = leaf_PzV0->GetValue(trk);
            float PIDAsso   = leaf_PIDV0->GetValue(trk);
            float MassAsso    = leaf_e->GetValue(trk);
            //TLorentzVector part_vec;
            //part_vec.SetPxPyPzE(PxAsso, PyAsso, PzAsso, energy);
            //std::cout << "First Particle ID: " << PIDAsso << std::endl;
            //float MassAsso = part_vec.M();
            //float MassAsso  = leaf_MassV0->GetValue(trk);
            //float ChargeAsso = 0;
            //if(PIDAsso== ID_pip || PIDAsso== ID_Kp || PIDAsso== ID_pp) ChargeAsso = 1;
            //if(PIDAsso== ID_pim || PIDAsso== ID_Km || PIDAsso== ID_pb) ChargeAsso =-1;
            float PtAsso = sqrt(PxAsso*PxAsso+PyAsso*PyAsso);
            float PhiAsso = atan2(PyAsso,PxAsso);
            float ThetaAsso = atan2(PtAsso,PzAsso);
            float EtaAsso = -log(tan(ThetaAsso/2.));
            yy = 0.5*log((sqrt(MassAsso*MassAsso+PxAsso*PxAsso+PyAsso*PyAsso+PzAsso*PzAsso)+PzAsso)/(sqrt(MassAsso*MassAsso+PxAsso*PxAsso+PyAsso*PyAsso+PzAsso*PzAsso)-PzAsso));

            int pid = (int)PIDAsso;
            // Check if the particle exists in our table and if it has a non-zero charge
            if (pCharge.find(pid) == pCharge.end() || pCharge[pid] == 0.0) {
                continue;
            }

            // Dynamically assign the charge
            float ChargeAsso = pCharge[pid]; 

            // Fill the pie chart using the particle's string name!
            hSpeciesAll->Fill(pName[pid].c_str(), 1.0);

            if(!(pid == 130 || pid == -130 || pid == 120 || pid == -120 || pid == 1120 || pid == -1120)){
                continue;
            }

            if(ChargeAsso <= 0.0001 && ChargeAsso >= -0.0001) {
                continue;
            }
            mult_all++;
            if(PtAsso > 0.05){
                mult_ptCut++;
                if(EtaAsso < 0.5 && EtaAsso > -0.5){
                    mult_ptEtaCut++;
                }

                hEtaPtDist->Fill(EtaAsso,PtAsso,Eweight);
                Hist_Pt->Fill(PtAsso,Eweight);
                Hist_Y->Fill(yy);
                Hist_Phi->Fill(PhiAsso);
            }
        }
        hMult_All->Fill(mult_all);
        hMult_ptCut->Fill(mult_ptCut);
        hMult_ptEtaCut->Fill(mult_ptEtaCut);
    }
    fout.Write();
    return;
}
