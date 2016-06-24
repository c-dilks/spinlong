#include "KinBounds.h"

ClassImp(KinBounds)

using namespace std;

KinBounds::KinBounds(Environ * env0) {
  env = env0;


  // open up threshtr and set branch addresses
  char filename[1024];
  sprintf(filename,"%s/diagset/setdep.root",env->SpinDir);
  setdepfile = new TFile(filename,"READ");

  thtr = (TTree*) setdepfile->Get("threshtr");
  thtr->SetBranchAddress("runnum",&th_runnum);
  thtr->SetBranchAddress("index",&th_index);
  thtr->SetBranchAddress("class",&th_class);
  thtr->SetBranchAddress("trig",&th_trig);
  thtr->SetBranchAddress("thresh",&thresh);
  thtr->SetBranchAddress("thresh_err",&thresh_err);
  thtr->SetBranchAddress("which_thresh",which_thresh);


  // initialize arrays
  for(int t=0; t<kMaxNTRIGS; t++) {
    for(int c=0; c<kMaxNCLASSES; c++) {
      for(int r=0; r<kMaxNRUNS; r++) {
        pt_arr[t][c][r] = 0;
        ptE_arr[t][c][r] = 0;
        en_arr[t][c][r] = 0;
        enE_arr[t][c][r] = 0;
      };
    };
  };


  // fill arrays, indexed by kinvarvsrun run_index
  // also fill runnum_map, which maps runnum to kinvarvsrun run_index
  for(int x=0; x<thtr->GetEntries(); x++) {
    thtr->GetEntry(x);
    if(th_class-1>kMaxNCLASSES || th_trig-1>kMaxNTRIGS || th_index-1>kMaxNRUNS) {
      fprintf(stderr,"ERROR: KinBounds::KinBounds -- assumed max is too low, segfault will occur\n");
    };
    if(!strcmp(which_thresh,"pt")) {
      pt_arr[th_trig][th_class][th_index] = thresh;
      ptE_arr[th_trig][th_class][th_index] = thresh_err;
      runnum_map.insert(std::pair<Int_t,Int_t>(th_runnum,th_index)); // only need to fill it once...
    }
    else if(!strcmp(which_thresh,"en")) {
      en_arr[th_trig][th_class][th_index] = thresh;
      enE_arr[th_trig][th_class][th_index] = thresh_err;
    };
  };
};


//------------------ pT ------------------------------------------------------------
Float_t KinBounds::PtThreshLow(Int_t run_index,Int_t class_index,Int_t trig_index) {
  if(run_index>=0 && run_index<kMaxNRUNS &&
     class_index>=0 && class_index<kMaxNCLASSES &&
     trig_index>=0 && trig_index<kMaxNTRIGS) {
    return pt_arr[trig_index][class_index][run_index];
  }
  else return -1;
};


Float_t KinBounds::PtThreshHigh(Int_t run_index,Int_t class_index,Int_t trig_index) {
  return env->PtHigh;
};


Bool_t KinBounds::PtInRange(Float_t pt0, Int_t runnum0, Int_t class_index0, Int_t trig_index0) {
  Int_t idx = GetRunIdx(runnum0);
  if(idx<0) return false;
  if(pt0 >= PtThreshLow(idx,class_index0,trig_index0) &&
     pt0 <= PtThreshHigh(idx,class_index0,trig_index0))
      return true;
  else return false;
};



//------------------ En ------------------------------------------------------------
Float_t KinBounds::EnThreshLow(Int_t run_index,Int_t class_index,Int_t trig_index) {
  if(run_index>=0 && run_index<kMaxNRUNS &&
     class_index>=0 && class_index<kMaxNCLASSES &&
     trig_index>=0 && trig_index<kMaxNTRIGS) {
    return en_arr[trig_index][class_index][run_index];
  }
  else return -1;
};


Float_t KinBounds::EnThreshHigh(Int_t run_index,Int_t class_index,Int_t trig_index) {
  return env->EnHigh;
};


Bool_t KinBounds::EnInRange(Float_t en0, Int_t runnum0, Int_t class_index0, Int_t trig_index0) {
  Int_t idx = GetRunIdx(runnum0);
  if(idx<0) return false;
  if(en0 >= EnThreshLow(idx,class_index0,trig_index0) &&
     en0 <= EnThreshHigh(idx,class_index0,trig_index0))
      return true;
  else return false;
};


//----------------------------------------------------------------------------
Int_t KinBounds::GetRunIdx(Int_t runnum0) {
  Int_t retval;
  try { retval = runnum_map.at(runnum0); }
  catch(const std::out_of_range& e) {
    fprintf(stderr,"ERROR: KinBounds::runnum_map does not contain run %d; cannot proceed, will probably segfault\n",runnum0);
    retval=-1;
  };
  return retval;
};
