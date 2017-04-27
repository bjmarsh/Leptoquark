#ifndef QCDLOOPER_h
#define QCDLOOPER_h

// C++ includes
/* #include <string> */
#include <vector>

// ROOT includes
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "Math/LorentzVector.h"

//MT2
#include "/home/users/bemarsh/analysis/mt2/current/MT2Analysis/MT2CORE/mt2tree.h"
#include "/home/users/bemarsh/analysis/mt2/current/MT2Analysis/MT2CORE/SR.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class Looper {

 public:

    Looper();
    ~Looper();

    void SetSignalRegions();
    void loop(TChain* chain, std::string output_name = "test.root");
    void fillHistos(std::map<std::string, TH1*>& h_1d, const std::string& dirname, const std::string& s = "");

 private:

    TFile * outfile_;
    mt2tree t;
    float evtweight_;
    std::vector<SR> SRVec;
      
    float M1,M2,M1p,M2p,DM,DMp,AvgM,AvgMp;
    float mLL, Mbb, deltaPhiLL, deltaPhibb, deltaPhiMETJet;
    bool isMuMu, isEMu;
};

LorentzVector GetLorentzVector(float pt, float eta, float phi);

#endif

