#ifndef EPOSPID_H
#define EPOSPID_H
// Charge / mass lookup for EPOS particle ids, parsed from idt.dt at runtime.
//   idt.dt column 1  = id_EPOS,  column 11 = mass (GeV),  column 12 = charge (e)
// Most final-state ids are small (pi=120, K=130, p=1120, Xi=2330, ...), so a flat
// array gives O(1) lookup; the handful of huge ids (nuclei/remnants, up to 8e8)
// fall back to a map. The flat path matters: refmult loops every track of every
// event (~9500 x 646k), so a std::map lookup there would dominate the runtime.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <fstream>
#include <sstream>

namespace EposPID {

  const int  FLATMAX = 100000;
  static signed char  gChg[2*FLATMAX+1];   // offset by FLATMAX; charge in units of e
  static float        gMass[2*FLATMAX+1];
  static std::map<int,signed char> gChgBig;
  static std::map<int,float>       gMassBig;
  static bool gLoaded = false;

  inline bool Load(const char* fname = "idt.dt"){
    if(gLoaded) return true;
    memset(gChg, 0, sizeof(gChg));
    memset(gMass,0, sizeof(gMass));
    std::ifstream in(fname);
    if(!in){ printf("EposPID: cannot open %s\n", fname); return false; }
    std::string line; int nread=0;
    while(std::getline(in,line)){
      // skip comments / blanks
      size_t p = line.find_first_not_of(" \t");
      if(p==std::string::npos || line[p]=='!') continue;
      std::istringstream ss(line);
      std::string tok; std::string t[20]; int nt=0;
      while(nt<20 && (ss>>tok)) t[nt++]=tok;
      if(nt < 12) continue;
      int   id = atoi(t[0].c_str());
      float m  = atof(t[10].c_str());
      float q  = atof(t[11].c_str());
      signed char qc = (signed char)((q>0)? (q+0.5) : (q-0.5));   // round toward sign
      if(abs(id) <= FLATMAX){ gChg[id+FLATMAX]=qc; gMass[id+FLATMAX]=m; }
      else { gChgBig[id]=qc; gMassBig[id]=m; }
      nread++;
    }
    gLoaded = true;
    printf("EposPID: loaded %d entries from %s\n", nread, fname);
    return true;
  }

  // charge in units of e (0 for neutral / unknown)
  inline int Charge(int id){
    if(abs(id) <= FLATMAX) return gChg[id+FLATMAX];
    std::map<int,signed char>::const_iterator it = gChgBig.find(id);
    if(it != gChgBig.end()) return it->second;
    // antiparticle not listed explicitly -> negate the particle's charge
    it = gChgBig.find(-id);
    return (it != gChgBig.end()) ? -it->second : 0;
  }

  inline float Mass(int id){
    if(abs(id) <= FLATMAX) return gMass[id+FLATMAX];
    std::map<int,float>::const_iterator it = gMassBig.find(id);
    if(it != gMassBig.end()) return it->second;
    it = gMassBig.find(-id);
    return (it != gMassBig.end()) ? it->second : 0.f;
  }
}

// ---------------------------------------------------------------------------
// Reference multiplicity (the centrality estimator), defined to mimic STAR's
// TPC refmult: final-state (ist==0) CHARGED particles with |eta|<1 and pT>0.15.
// The |eta|<1 test is done as |pz| < sinh(1)*pT so no log/atan is needed --
// this is evaluated on every track of every event, so it must stay cheap.
const float REFMULT_PTMIN  = 0.15;
const float REFMULT_ETACUT = 1.0;
const float SINH_ETACUT    = 1.1752012;          // sinh(1.0)
const float SINH_ETACUT_SQ = 1.3810978;          // sinh(1.0)^2

// ---------------------------------------------------------------------------
// Centrality bins (STAR convention, 9 bins spanning 0-80%)
const int   NCENT = 9;
const double CENT_LO[NCENT] = { 0,  5, 10, 20, 30, 40, 50, 60, 70};
const double CENT_HI[NCENT] = { 5, 10, 20, 30, 40, 50, 60, 70, 80};
const double CENT_MID[NCENT]= {2.5,7.5,15, 25, 35, 45, 55, 65, 75};
inline const char* CentLabel(int i){
  static char buf[16]; snprintf(buf,16,"%g-%g%%",CENT_LO[i],CENT_HI[i]); return buf;
}

// Assign a centrality bin from refmult given the 9 lower thresholds
// (cut[i] = minimum refmult to be in bin i; cut[0] is the most central).
// Returns -1 for events more peripheral than 80%.
inline int CentBinFromRefmult(int refmult, const int* cut){
  for(int i=0;i<NCENT;i++) if(refmult >= cut[i]) return i;
  return -1;
}

#endif
