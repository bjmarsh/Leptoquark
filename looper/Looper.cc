// C++
#include <iostream>
#include <vector>
#include <set>
#include <cmath>

// ROOT
#include "TDirectory.h"
#include "TTreeCache.h"
#include "Math/VectorUtil.h"
#include "TVector2.h"
#include "TBenchmark.h"

// Tools
#include "../CORE/Tools/utils.h"
#include "../CORE/Tools/goodrun.h"
#include "../CORE/Tools/dorky/dorky.h"

// header
#include "Looper.h"

//MT2
#include "../MT2Analysis/MT2CORE/Plotting/PlotUtilities.h"

#define PI 3.14159265359

using namespace std;
using namespace duplicate_removal;

class mt2tree;
class SR;

// turn on to apply json file to data
bool applyJSON = true;

const float lumi = 36.26;

int nbins_DMAvgM_x = 2;
int nbins_DMAvgM_y = 5;
float bins_DMAvgM_x[3] = {0, 100, 1000};
float bins_DMAvgM_y[6] = {0, 200, 400, 600, 800, 2000};

//_______________________________________
Looper::Looper(){
}
//_______________________________________
Looper::~Looper(){

};

// samples all have their xsec set to 1.0 pb. Fix here
float SignalXsecCorrection(const string &name){

    if(name.find("T2tt_RPV_700") != string::npos)
        return 0.0670476;
    if(name.find("T2tt_RPV_900") != string::npos)
        return 0.0128895;
    if(name.find("T2tt_RPV_1100") != string::npos)
        return 0.00307413;

    return 1.;
}

//_______________________________________
void Looper::SetSignalRegions(){

    SR sr;

    sr.SetName("mumu_baseline");
    sr.SetVar("passZmass", 0, 2);
    sr.SetVar("isMuMu", 1, 2);
    sr.SetVar("isEMu", 0, 1);
    sr.SetVar("DM",0,-1);
    sr.SetVar("AvgM",0,-1);
    sr.SetVar("bidx2",0,-1);
    SRVec.push_back(sr);
    outfile_->mkdir(sr.GetName().c_str());

    sr.SetName("emu_baseline");
    sr.SetVar("passZmass", 0, 2);
    sr.SetVar("isMuMu", 0, 1);
    sr.SetVar("isEMu", 1, 2);
    sr.SetVar("DM",0,-1);
    sr.SetVar("AvgM",0,-1);
    sr.SetVar("bidx2",0,-1);
    SRVec.push_back(sr);
    outfile_->mkdir(sr.GetName().c_str());

    sr.SetName("mumu_zmass");
    sr.SetVar("passZmass", 1, 2);
    sr.SetVar("isMuMu", 1, 2);
    sr.SetVar("isEMu", 0, 1);
    sr.SetVar("DM",0,-1);
    sr.SetVar("AvgM",0,-1);
    sr.SetVar("bidx2",0,-1);
    SRVec.push_back(sr);
    outfile_->mkdir(sr.GetName().c_str());

    sr.SetName("emu_zmass");
    sr.SetVar("passZmass", 1, 2);
    sr.SetVar("isMuMu", 0, 1);
    sr.SetVar("isEMu", 1, 2);
    sr.SetVar("DM",0,-1);
    sr.SetVar("AvgM",0,-1);
    sr.SetVar("bidx2",0,-1);
    SRVec.push_back(sr);
    outfile_->mkdir(sr.GetName().c_str());

    sr.SetName("mumu_sr");
    sr.SetVar("passZmass", 1, 2);
    sr.SetVar("isMuMu", 1, 2);
    sr.SetVar("isEMu", 0, 1);
    sr.SetVar("DM",0,100);
    sr.SetVar("AvgM",0,-1);
    sr.SetVar("bidx2",0,-1);
    SRVec.push_back(sr);
    outfile_->mkdir(sr.GetName().c_str());

    sr.SetName("mumu_crdm");
    sr.SetVar("passZmass", 1, 2);
    sr.SetVar("isMuMu", 1, 2);
    sr.SetVar("isEMu", 0, 1);
    sr.SetVar("DM",100,-1);
    sr.SetVar("AvgM",0,-1);
    sr.SetVar("bidx2",0,-1);
    SRVec.push_back(sr);
    outfile_->mkdir(sr.GetName().c_str());

    sr.SetName("emu_sr");
    sr.SetVar("passZmass", 1, 2);
    sr.SetVar("isMuMu", 0, 1);
    sr.SetVar("isEMu", 1, 2);
    sr.SetVar("DM",0,100);
    sr.SetVar("AvgM",0,-1);
    sr.SetVar("bidx2",0,-1);
    SRVec.push_back(sr);
    outfile_->mkdir(sr.GetName().c_str());

    sr.SetName("emu_crdm");
    sr.SetVar("passZmass", 1, 2);
    sr.SetVar("isMuMu", 0, 1);
    sr.SetVar("isEMu", 1, 2);
    sr.SetVar("DM",100,-1);
    sr.SetVar("AvgM",0,-1);
    sr.SetVar("bidx2",0,-1);
    SRVec.push_back(sr);
    outfile_->mkdir(sr.GetName().c_str());

    sr.SetName("mumu_sr_AvgM500");
    sr.SetVar("passZmass", 1, 2);
    sr.SetVar("isMuMu", 1, 2);
    sr.SetVar("isEMu", 0, 1);
    sr.SetVar("DM",0,100);
    sr.SetVar("AvgM",500,-1);
    sr.SetVar("bidx2",0,-1);
    SRVec.push_back(sr);
    outfile_->mkdir(sr.GetName().c_str());

    sr.SetName("mumu_sr_AvgM500_Btop2");
    sr.SetVar("passZmass", 1, 2);
    sr.SetVar("isMuMu", 1, 2);
    sr.SetVar("isEMu", 0, 1);
    sr.SetVar("DM",0,100);
    sr.SetVar("AvgM",500,-1);
    sr.SetVar("bidx2",1, 2);
    SRVec.push_back(sr);
    outfile_->mkdir(sr.GetName().c_str());

    sr.SetName("emu_sr_AvgM500");
    sr.SetVar("passZmass", 1, 2);
    sr.SetVar("isMuMu", 0, 1);
    sr.SetVar("isEMu", 1, 2);
    sr.SetVar("DM",0,100);
    sr.SetVar("AvgM",500,-1);
    sr.SetVar("bidx2",0,-1);
    SRVec.push_back(sr);
    outfile_->mkdir(sr.GetName().c_str());

    sr.SetName("emu_sr_AvgM500_Btop2");
    sr.SetVar("passZmass", 1, 2);
    sr.SetVar("isMuMu", 0, 1);
    sr.SetVar("isEMu", 1, 2);
    sr.SetVar("DM",0,100);
    sr.SetVar("AvgM",500,-1);
    sr.SetVar("bidx2",1, 2);
    SRVec.push_back(sr);
    outfile_->mkdir(sr.GetName().c_str());

}


//_______________________________________
void Looper::loop(TChain* chain, std::string output_name){

    gROOT->cd();
    cout << "[Looper::loop] creating output file: " << output_name << endl;

    outfile_ = new TFile(output_name.c_str(),"RECREATE") ; 

    const char* json_file = "/home/users/bemarsh/analysis/mt2/current/MT2Analysis/babymaker/jsons/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_snt.txt";
    if (applyJSON) {
        cout << "Loading json file: " << json_file << endl;
        set_goodrun_file(json_file);
    }

    cout << "[Looper::loop] setting up histos" << endl;

    SetSignalRegions();

    // File Loop
    int nDuplicates = 0;
    int nEvents = chain->GetEntries();
    unsigned int nEventsChain = nEvents;
    cout << "[Looper::loop] running on " << nEventsChain << " events" << endl;
    unsigned int nEventsTotal = 0;
    TObjArray *listOfFiles = chain->GetListOfFiles();
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;
    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        cout << "[Looper::loop] running on file: " << currentFile->GetTitle() << endl;

        // Get File Content
        TFile f( currentFile->GetTitle() );
        TTree *tree = (TTree*)f.Get("mt2");
        TTreeCache::SetLearnEntries(10);
        tree->SetCacheSize(128*1024*1024);
    
        t.Init(tree);

        // Event Loop
        unsigned int nEventsTree = tree->GetEntriesFast();
        for( unsigned int event = 0; event < nEventsTree; ++event) {
      
            t.GetEntry(event);

            //---------------------
            // bookkeeping and progress report
            //---------------------
            ++nEventsTotal;
            if (nEventsTotal%10000==0) {
                ULong64_t i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
                if (isatty(1)) {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                           "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
                    fflush(stdout);
                }else {
                    cout<<i_permille/10.<<" ";
                    if (nEventsTotal%100000==0) cout<<endl;
                }
            }

            //---------------------
            // skip duplicates -- needed when running on mutiple streams in data
            //---------------------
            if( t.isData ) {
                DorkyEventIdentifier id(t.run, t.evt, t.lumi);
                if (is_duplicate(id) ){
                    ++nDuplicates;
                    continue;
                }
            }

            outfile_->cd();
            evtweight_ = 1.;

            // apply relevant weights to MC
            if (!t.isData) {
                evtweight_ = t.evt_scale1fb * lumi;
            } // !isData

            // fix the xsec of signal samples
            // it is set to 1 pb in the ntuples
            evtweight_ *= SignalXsecCorrection(currentFile->GetTitle());

            plot1D("h_weight",  evtweight_, 1, h_1d_global, ";evtweight", 100, 0, 1);

            //---------------------
            // basic event selection and cleaning
            //---------------------

            if( applyJSON && t.isData && !goodrun(t.run, t.lumi) ) continue;
      
            // MET filters (first 2 data only)
            if (t.isData) {
                if (!t.Flag_globalTightHalo2016Filter) continue; 
                if (!t.Flag_badMuonFilter) continue;
            }
            if (!t.Flag_goodVertices) continue;
            if (!t.Flag_eeBadScFilter) continue;
            if (!t.Flag_HBHENoiseFilter) continue;
            if (!t.Flag_HBHENoiseIsoFilter) continue;
            if (!t.Flag_EcalDeadCellTriggerPrimitiveFilter) continue;
            if (!t.Flag_badChargedHadronFilter) continue;

            // trigger requirement in data
            if (t.isData && !(t.HLT_DoubleMu || t.HLT_DoubleMu_NonIso || t.HLT_SingleMu_NonIso))
                continue;
      

            // some simple baseline selections
            if (t.nVert == 0) continue;

            //count the number of loose bjets with pt>40
            int nBJet40 = 0;
            bidx1 = -1;
            bidx2 = -1;
            for(int i=0; i<t.nJet30; i++){
                if(t.jet_pt[i] < 40)
                    break;

                // is loose bjet
                if(t.jet_btagCSV[i] > 0.5426){
                    nBJet40++;
                    if(nBJet40==1)
                        bidx1 = i;
                    if(nBJet40==2)
                        bidx2 = i;
                }
            }

            if(nBJet40<2) continue;

            // check if we have 2 OS leptons with pt>25/25
            // for muons, MT2 babies already apply loose ID and iso<0.2
            // for elecs, MT2 babies apply veto ID and iso<0.1
            if(t.nlep<2) continue;
            // // just muons for now
            // if(abs(t.lep_pdgId[0])!=13 || abs(t.lep_pdgId[1])!=13)
            //     continue;
            // opposite sign
            if(t.lep_charge[0]*t.lep_charge[1] > 0)
                continue;
            // both must have pt>25
            if(t.lep_pt[1]<25)
                continue;

            //require electrons to be loose (already applied at baby-level for mus)
            if(abs(t.lep_pdgId[0])==11 && t.lep_tightId[0]>=1) continue;
            if(abs(t.lep_pdgId[1])==11 && t.lep_tightId[0]>=1) continue;

            isMuMu = false;
            isEMu = false;
            if(abs(t.lep_pdgId[0])==13 && abs(t.lep_pdgId[1])==13)
                isMuMu = true;
            if((abs(t.lep_pdgId[0])==11 && abs(t.lep_pdgId[1])==13) || 
               (abs(t.lep_pdgId[0])==13 && abs(t.lep_pdgId[1])==11))
                isEMu = true;

            // look at invariant masses of both pairings of b/mu's
            LorentzVector b1_p4 = GetLorentzVector(t.jet_pt[bidx1], t.jet_eta[bidx1], t.jet_phi[bidx1]);
            LorentzVector b2_p4 = GetLorentzVector(t.jet_pt[bidx2], t.jet_eta[bidx2], t.jet_phi[bidx2]);
            LorentzVector l1_p4 = GetLorentzVector(t.lep_pt[0], t.lep_eta[0], t.lep_phi[0]);
            LorentzVector l2_p4 = GetLorentzVector(t.lep_pt[1], t.lep_eta[1], t.lep_phi[1]);

            mLL = (l1_p4+l2_p4).M();
            mbb = (b1_p4+b2_p4).M();

            M1 = (b1_p4+l1_p4).M();
            M2 = (b2_p4+l2_p4).M();
            
            // the "primed" variables (trailing p) are the second possible lep-b pairing.
            // check below which one yields a lower DM, and use those
            M1p = (b2_p4+l1_p4).M();
            M2p = (b1_p4+l2_p4).M();

            DM = fabs(M1-M2);
            DMp = fabs(M1p-M2p);

            AvgM = (M1+M2)/2;
            AvgMp = (M1p+M2p)/2;            

            ST = b1_p4.pt() + b2_p4.pt() + l1_p4.pt() + l2_p4.pt();
            
            deltaPhiLL = fabs(t.lep_phi[0] - t.lep_phi[1]);
            if(deltaPhiLL > PI)
                deltaPhiLL = 2*PI - deltaPhiLL;

            deltaPhibb = fabs(t.jet_phi[bidx1] - t.jet_phi[bidx2]);
            if(deltaPhibb > PI)
                deltaPhibb = 2*PI - deltaPhibb;

            float deltaPhiMETJet1 = fabs(t.met_phi - t.jet_phi[bidx1]);
            if(deltaPhiMETJet1 > PI)
                deltaPhiMETJet1 = 2*PI - deltaPhiMETJet1;
            float deltaPhiMETJet2 = fabs(t.met_phi - t.jet_phi[bidx2]);
            if(deltaPhiMETJet2 > PI)
                deltaPhiMETJet2 = 2*PI - deltaPhiMETJet2;
            deltaPhiMETJet = min(deltaPhiMETJet1, deltaPhiMETJet2);
            
            std::map<std::string, float> values;
            values["passZmass"] = (mLL<91.19-15 || mLL>91.19+15);
            values["isMuMu"] = isMuMu;
            values["isEMu"] = isEMu;
            if(DM<DMp){
                values["DM"] = DM;
                values["AvgM"] = AvgM;
            }else{
                values["DM"] = DMp;
                values["AvgM"] = AvgMp;
            }
            values["bidx2"] = bidx2;
            for(unsigned int i=0; i<SRVec.size(); i++){
                if(SRVec.at(i).PassesSelection(values))
                    fillHistos(SRVec.at(i).srHistMap, SRVec.at(i).GetName().c_str());
            }
      
        }//end loop on events in a file
  
        delete tree;
        f.Close();
    }//end loop on files
  
    cout << "[Looper::loop] processed " << nEventsTotal << " events" << endl;
    if ( nEventsChain != nEventsTotal ) {
        std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    }

    cout << nDuplicates << " duplicate events were skipped." << endl;

    //---------------------
    // Save Plots
    //---------------------

    outfile_->cd();
    for(unsigned int i=0; i<SRVec.size(); i++)
        savePlotsDir(SRVec.at(i).srHistMap,outfile_,SRVec.at(i).GetName().c_str());

    //---------------------
    // Write and Close file
    //---------------------
    outfile_->Write();
    outfile_->Close();
    delete outfile_;

    return;
}

void Looper::fillHistos(std::map<std::string, TH1*>& h_1d, const std::string& dirname, const std::string& s) {

    TDirectory * dir = (TDirectory*)outfile_->Get(dirname.c_str());
    if (dir == 0) {
        dir = outfile_->mkdir(dirname.c_str());
           
    } 
    dir->cd();

    // use primed variables if second lep-b pairing is better
    if(DMp<DM){
        M1 = M1p;
        M2 = M2p;
        DM = DMp;
        AvgM = AvgMp;
    }

    // if(AvgM>1000 && evtweight_>1.0)
    //     cout << t.evt << " " << AvgM << " " << evtweight_ << endl;
    
    plot1D("h_Events"+s,  1, 1, h_1d, ";Events, Unweighted", 1, 0, 2);
    plot1D("h_Events_w"+s,  1,   evtweight_, h_1d, ";Events, Weighted", 1, 0, 2);
    plot1D("h_M1"+s,        M1,   evtweight_, h_1d, "b/mu Inv Mass",100,0,2000);
    plot1D("h_M2"+s,        M2,   evtweight_, h_1d, "b/mu Inv Mass",100,0,2000);
    plot1D("h_MLL"+s,      mLL,   evtweight_, h_1d, ";m_{LL}",100,0,1000);
    plot1D("h_Mbb"+s,      mbb,   evtweight_, h_1d, ";m_{bb}",100,0,1000);
    plot1D("h_ST"+s,      ST,   evtweight_, h_1d, ";S_{T}",100,0,3000);
    plot1D("h_deltaPhiLL"+s,      deltaPhiLL,   evtweight_, h_1d, ";#Delta#phi_{LL}",100,0,PI);
    plot1D("h_deltaPhibb"+s,      deltaPhibb,   evtweight_, h_1d, ";#Delta#phi_{bb}",100,0,PI);
    plot1D("h_deltaPhiMETJet"+s,      deltaPhiMETJet,   evtweight_, h_1d, ";#Delta#phi_{MET-Jet}",100,0,PI);
    plot1D("h_bidx1"+s,         bidx1,     evtweight_, h_1d, ";bidx1",5,0,5);
    plot1D("h_bidx2"+s,         bidx2,     evtweight_, h_1d, ";bidx2",5,0,5);
    plot2D("h_M1M2"+s, M1, M2, evtweight_, h_1d, "b/mu Inv Masses;M1;M2",100,0,2000,100,0,2000);

    plot1D("h_DM"+s,        DM,   evtweight_, h_1d, "min(#Delta M)",100,0,1000);
    plot1D("h_AvgM"+s,       AvgM,   evtweight_, h_1d, "Avg M for min(#Delta M)",100,0,2000);
    plot2D("h_DMAvgM"+s, DM, AvgM, evtweight_, h_1d, ";#Delta M;Avg M",100,0,1000,100,0,2000);
    plot2D("h_DMAvgM_binned"+s, DM, AvgM, evtweight_, h_1d, ";#Delta M;Avg M", nbins_DMAvgM_x, bins_DMAvgM_x, nbins_DMAvgM_y, bins_DMAvgM_y);

    plot1D("h_met"+s,        t.met_pt,   evtweight_, h_1d, "min(#Delta M)",100,0,1000);


    outfile_->cd();
    return;
}

LorentzVector GetLorentzVector(float pt, float eta, float phi){

    float theta = 2*atan(exp(-eta));
    float magp = pt/sin(theta);
    float px = pt*cos(phi);
    float py = pt*sin(phi);
    float pz = magp*cos(theta);

    LorentzVector p4(px,py,pz,magp);
        return p4;             
}
