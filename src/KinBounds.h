#ifndef KinBounds_
#define KinBounds_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <map>
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

#include "Environ.h"
#include "LevelTwo.h"
#include "EventClass.h"

class KinBounds : public TObject
{
  public:
    KinBounds(Environ * env0, 
              LevelTwo * LT0,
              EventClass * ev0);
    Float_t PtThreshLow(Int_t run_index,Int_t class_index,Int_t trig_index);
    Float_t PtThreshLowErr(Int_t run_index,Int_t class_index,Int_t trig_index);
    Float_t PtThreshHigh(Int_t run_index,Int_t class_index,Int_t trig_index);
    Float_t PtThreshHighErr(Int_t run_index,Int_t class_index,Int_t trig_index);
    Bool_t PtIsGood(Float_t pt0, Int_t run_index0 ,Int_t class_index0, Int_t trig_index0);


    enum { kMaxNTRIGS=30, kMaxNCLASSES=10, kMaxNRUNS=2000 };


  private: 
    Environ * env;
    LevelTwo * LT;
    EventClass * ev;
    TFile * setdepfile;
    TTree * thtr;

    Int_t th_runnum,th_index,th_class,th_trig;
    Float_t Tpt,TptE;

    Int_t NTRIGS,NCLASSES;

    Float_t pt_arr[kMaxNTRIGS][kMaxNCLASSES][kMaxNRUNS];
    Float_t ptE_arr[kMaxNTRIGS][kMaxNCLASSES][kMaxNRUNS];

    ClassDef(KinBounds,1);
};

#endif
