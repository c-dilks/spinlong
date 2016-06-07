#include "KinBounds.h"

ClassImp(KinBounds)

using namespace std;

KinBounds::KinBounds(Environ * env0, 
                     LevelTwo * LT0,
                     EventClass * ev0) {
  env = env0;
  LT = LT0;
  ev = ev0;

  NTRIGS = LT->N;
  NCLASSES = ev->N;


  // open up threshtr and set branch addresses
  char filename[1024];
  sprintf(filename,"%s/diagset/setdep.root",env->SpinDir);
  setdepfile = new TFile(filename,"READ");

  thtr = (TTree*) setdepfile->Get("threshtr");
  thtr->SetBranchAddress("runnum",&th_runnum);
  thtr->SetBranchAddress("index",&th_index);
  thtr->SetBranchAddress("class",&th_class);
  thtr->SetBranchAddress("trig",&th_trig);
  thtr->SetBranchAddress("ptthresh",&Tpt);
  thtr->SetBranchAddress("ptthresh_err",&TptE);


  // initialize arrays
  for(int t=0; t<kMaxNTRIGS; t++) {
    for(int c=0; c<kMaxNCLASSES; c++) {
      for(int r=0; r<kMaxNRUNS; r++) {
        pt_arr[t][c][r] = 0;
        ptE_arr[t][c][r] = 0;
      };
    };
  };


  // fill arrays, indexed by run_index
  for(int x=0; x<thtr->GetEntries(); x++) {
    thtr->GetEntry(x);
    pt_arr[th_trig][th_class][th_index] = Tpt;
    ptE_arr[th_trig][th_class][th_index] = TptE;
  };
};



Float_t KinBounds::PtThreshLow(Int_t run_index,Int_t class_index,Int_t trig_index) {
  if(run_index>=0 && run_index<kMaxNRUNS &&
     class_index>=0 && class_index<kMaxNCLASSES &&
     trig_index>=0 && trig_index<kMaxNTRIGS) {
    return pt_arr[trig_index][class_index][run_index];
  }
  else return -1;
};


Float_t KinBounds::PtThreshLowErr(Int_t run_index,Int_t class_index,Int_t trig_index) {
  if(run_index>=0 && run_index<kMaxNRUNS &&
     class_index>=0 && class_index<kMaxNCLASSES &&
     trig_index>=0 && trig_index<kMaxNTRIGS) {
    return ptE_arr[trig_index][class_index][run_index];
  }
  else return -1;
};


Float_t KinBounds::PtThreshHigh(Int_t run_index,Int_t class_index,Int_t trig_index) {
  return 10;
};


Float_t KinBounds::PtThreshHighErr(Int_t run_index,Int_t class_index,Int_t trig_index) {
  return 0;
};


Bool_t KinBounds::PtIsGood(Float_t pt0, Int_t run_index0 ,Int_t class_index0, Int_t trig_index0) {
  if(pt0 >= PtThreshLow(run_index0,class_index0,trig_index0) &&
     pt0 <= PtThreshHigh(run_index0,class_index0,trig_index0))
      return true;
  else return false;
};
