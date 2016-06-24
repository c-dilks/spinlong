#ifndef KinBounds_
#define KinBounds_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

#include "Environ.h"

class KinBounds : public TObject
{
  public:
    KinBounds(Environ * env0);

    Float_t PtThreshLow(Int_t run_index,Int_t class_index,Int_t trig_index);
    Float_t PtThreshHigh(Int_t run_index,Int_t class_index,Int_t trig_index);
    Bool_t PtInRange(Float_t pt0, Int_t runnum0 ,Int_t class_index0, Int_t trig_index0);

    Float_t EnThreshLow(Int_t run_index,Int_t class_index,Int_t trig_index);
    Float_t EnThreshHigh(Int_t run_index,Int_t class_index,Int_t trig_index);
    Bool_t EnInRange(Float_t en0, Int_t runnum0 ,Int_t class_index0, Int_t trig_index0);

    Int_t GetRunIdx(Int_t runnum0); // get kinvarvsrun runindex from given runnum


    // assume maxima here
    enum { kMaxNTRIGS=30, kMaxNCLASSES=10, kMaxNRUNS=2000 };


  private: 
    Environ * env;
    TFile * setdepfile;
    TTree * thtr;

    Int_t th_runnum,th_index,th_class,th_trig;
    Float_t thresh,thresh_err;
    char which_thresh[32];

    Int_t NTRIGS,NCLASSES;

    Float_t pt_arr[kMaxNTRIGS][kMaxNCLASSES][kMaxNRUNS];
    Float_t ptE_arr[kMaxNTRIGS][kMaxNCLASSES][kMaxNRUNS];
    Float_t en_arr[kMaxNTRIGS][kMaxNCLASSES][kMaxNRUNS];
    Float_t enE_arr[kMaxNTRIGS][kMaxNCLASSES][kMaxNRUNS];

    std::map<Int_t, Int_t> runnum_map; // runnum --> kinvarvsrun runindex

    ClassDef(KinBounds,1);
};

#endif
