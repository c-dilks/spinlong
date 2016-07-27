#ifndef RunInfo_
#define RunInfo_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <map>
#include <stdexcept>

#include "TSystem.h"
#include "TObject.h"
#include "TTree.h"
#include "TFile.h"

#include "Environ.h"

const Int_t NRUNS = 2000;


class RunInfo : public TObject
{
  public:
    RunInfo();
    Int_t Index(Int_t runnum0);
    Int_t IndexCounts(Int_t runnum0);
    Int_t IndexPol(Int_t runnum0);
    Int_t GetFill(Int_t runnum0);
    Float_t Rellum(Int_t runnum0, Int_t rellumi, char * detector);
    Float_t RellumErr(Int_t runnum0, Int_t rellumi, char * detector);
    Bool_t RellumConsistent(Int_t runnum0); // true if rellumi consistent (defined for R3)
    Float_t BluePol(Int_t runnum0);
    Float_t YellPol(Int_t runnum0);
    Float_t BluePolErr(Int_t runnum0);
    Float_t YellPolErr(Int_t runnum0);
    Int_t BlueSpin(Int_t runnum0, Int_t bXing);
    Int_t YellSpin(Int_t runnum0, Int_t bXing);
    Bool_t Kicked(Int_t runnum0, Int_t bXing);
    Int_t Pattern(Int_t runnum0);

    Environ * env;

  private:

    // branches -- rtree tree
    Int_t i_rtree; // run index
    Int_t runnum;
    Int_t fill_rtree;
    Double_t t; // run time
    Float_t R[3][10]; // rellum [detector] [which rellum 1-10]
    Float_t R_err[3][10]; // rellum uncertainty [detector] [which rellum 1-10]
    Bool_t isConsistent; // true iff zdc & vpd agree

    // branches -- counts tree
    Int_t i_counts; // run index
    Int_t runnum_counts; // run number
    Int_t bx; // bXing # @ STAR
    Int_t blue; // blue spin
    Int_t yell; // yellow spin
    Bool_t kicked; // true if bXing empty according to scaler data

    // branches -- pol tree
    Int_t fill_pol,runnum_pol;
    Float_t b_pol,y_pol; // polarization
    Float_t b_pol_err,y_pol_err; // polarization uncertainty
    Float_t b_pol_avg,y_pol_avg; // beam current-weighted fill-by-fill polarization
    Int_t pattern_no;


    // data arrays
    Int_t blue_spin_arr[NRUNS][120]; // [max no. runs (assumed)] [no. bXings] 
    Int_t yell_spin_arr[NRUNS][120]; 
    Bool_t kicked_bx_map[NRUNS][120]; 

    // index maps
    std::map<Int_t, Int_t> idx_rellum_map; // runnum => tree entry number in rtree.root : TTree * rellum
    std::map<Int_t, Int_t> idx_counts_map; // runnum => run index - 1 in counts.root
    std::map<Int_t, Int_t> idx_pol_map; // runnum => tree entry number in pol.root : TTree * poltr
    std::map<Int_t, Int_t>::iterator iter;
    Int_t idx;

    std::map<std::string, Int_t> sca_idx; // scaler detector name --> enumerator
    std::map<Int_t,std::string> sca_det; // enumerator --> scaler detector name


    // TFiles
    TFile * counts;
    TFile * rtree;
    TFile * pol;
    TTree * counts_tr;
    TTree * rtree_tr;
    TTree * pol_tr;


    ClassDef(RunInfo,2);


};

#endif
