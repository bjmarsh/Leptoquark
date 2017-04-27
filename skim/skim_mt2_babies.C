#include <fstream>
#include <sstream>
#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

using namespace std;

void skim_mt2_babies(string inpath = "/nfs-6/userdata/mt2/V00-00-03", string outpath = "/nfs-6/userdata/mt2/V00-00-03_skim_nj2_ht450_mt2gt50", string tag = "ttall_msdecays") {
  
    //--------------------------------------------------
    // cut for output files
    //--------------------------------------------------
 
    const char* sel = "nJet30>=2 && jet_pt[1]>=40 && nlep>=2 && lep_charge[0]*lep_charge[1]<0";

    cout << "Skimming with selection : "<<sel<<endl;

    //--------------------------------------------------
    // input and output file
    //--------------------------------------------------
  
    const char* infilename = Form("%s/%s*.root",inpath.c_str(),tag.c_str());
    const char* outfilename = Form("%s/%s.root",outpath.c_str(),tag.c_str());
  
    //--------------------------------------------------
    // cout stuff
    //--------------------------------------------------
  
    cout << "Reading in : " << infilename << endl;
    cout << "Writing to : " << outfilename << endl;
    cout << "Selection : " << sel << endl;
  
    //--------------------------------------------------
    // read input file, write to output files
    //--------------------------------------------------
  
    long long max_tree_size = 5000000000LL; // 5GB
    TTree::SetMaxTreeSize(max_tree_size);
  
    TChain *chain = new TChain("mt2");
    chain->Add(infilename);

    unsigned int input_entries = chain->GetEntries();
    cout << "Input tree has entries: " << input_entries << endl;
  
    //-------------------
    // skim
    //-------------------
  
    TFile *outfile = TFile::Open(outfilename, "RECREATE");
    assert( outfile != 0 );
    TTree* outtree = chain->CopyTree( sel );
    unsigned int output_entries = outtree->GetEntries();
    cout << "Output tree has entries: " << output_entries << endl
         << "Reduction factor of: " << double(input_entries)/double(output_entries) << endl;
    outtree->Write();
    outfile->Close();

}
