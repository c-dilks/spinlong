#include "RunInfo.h"


ClassImp(RunInfo)

namespace {
  enum scaler_enum {kBBC,kZDC,kVPD};
};

RunInfo::RunInfo()
{
  env = new Environ();


  // read trees
  counts = new TFile(env->counts_file,"READ");
  rtree = new TFile(env->rtree_file,"READ");
  pol = new TFile(env->pol_file,"READ");
  counts_tr = (TTree*) counts->Get("sca");
  rtree_tr = (TTree*) rtree->Get("rellum");
  pol_tr = (TTree*) pol->Get("poltr");
  

  // set rtree branch addresses
  rtree_tr->SetBranchAddress("i",&i_rtree);
  rtree_tr->SetBranchAddress("runnum",&runnum);
  rtree_tr->SetBranchAddress("fill",&fill_rtree);
  rtree_tr->SetBranchAddress("t",&t);
  rtree_tr->SetBranchAddress("isConsistent",&isConsistent);
  rtree_tr->SetBranchAddress("pattern",&pattern_no);

  sca_idx.insert(std::pair<std::string,Int_t>(std::string("bbc"),kBBC));
  sca_idx.insert(std::pair<std::string,Int_t>(std::string("zdc"),kZDC));
  sca_idx.insert(std::pair<std::string,Int_t>(std::string("vpd"),kVPD));

  std::map<std::string,Int_t>::iterator sca_iter;
  for(sca_iter = sca_idx.begin(); sca_iter!=sca_idx.end(); sca_iter++) {
    sca_det.insert(std::pair<Int_t,std::string>(sca_iter->second,sca_iter->first));
  };

  char R_name[3][10][16]; // [detector] [1-10]
  char R_err_name[3][10][16]; // [detector] [1-10]

  for(Int_t r=1; r<10; r++)
  {
    for(Int_t dd=0; dd<3; dd++) {
      sprintf(R_name[dd][r],"R%d_%srsc",r,std::string(sca_det.at(dd)).data());
      sprintf(R_err_name[dd][r],"R%d_%s_rsc_err",r,std::string(sca_det.at(dd)).data());

      rtree_tr->SetBranchAddress(R_name[dd][r],&(R[dd][r]));
      rtree_tr->SetBranchAddress(R_err_name[dd][r],&(R_err[dd][r]));
    };
  };


  // set counts tree branch addresses
  counts_tr->SetBranchAddress("i",&i_counts);
  counts_tr->SetBranchAddress("runnum",&runnum_counts);
  counts_tr->SetBranchAddress("bx",&bx);
  counts_tr->SetBranchAddress("blue",&blue);
  counts_tr->SetBranchAddress("yell",&yell);
  counts_tr->SetBranchAddress("kicked",&kicked);


  // set pol tree branch addresses
  pol_tr->SetBranchAddress("fill",&fill_pol);
  pol_tr->SetBranchAddress("runnum",&runnum_pol);
  pol_tr->SetBranchAddress("b_pol",&b_pol);
  pol_tr->SetBranchAddress("y_pol",&y_pol);
  pol_tr->SetBranchAddress("b_pol_avg",&b_pol_avg);
  pol_tr->SetBranchAddress("y_pol_avg",&y_pol_avg);
  //pol_tr->SetBranchAddress("b_pol_err",&b_pol_err);
  //pol_tr->SetBranchAddress("y_pol_err",&y_pol_err);
  b_pol_err=0; // for now..
  y_pol_err=0; // for now..



  // build spin maps --> gives blue/yell spin + kicked status for 
  //                     any given run index and bXing no.
  for(Int_t q=0; q<NRUNS; q++)
  {
    for(Int_t qq=0; qq<120; qq++)
    {
      blue_spin_arr[q][qq]=0;
      yell_spin_arr[q][qq]=0;
      kicked_bx_map[q][qq]=false;
    }
  };
  for(Int_t q=0; q<counts_tr->GetEntries(); q++)
  {
    counts_tr->GetEntry(q);
    if(idx_counts_map.find(runnum_counts)==idx_counts_map.end())
    {
      idx_counts_map.insert(std::pair<Int_t,Int_t>(runnum_counts,i_counts-1));
    };

    blue_spin_arr[i_counts-1][bx] = blue;
    yell_spin_arr[i_counts-1][bx] = yell;
    kicked_bx_map[i_counts-1][bx] = kicked;
  };
  printf("spin maps built\n");


  // build runnum => rellum table entry # hash
  for(Int_t q=0; q<rtree_tr->GetEntries(); q++)
  {
    rtree_tr->GetEntry(q);
    if(idx_rellum_map.find(runnum)==idx_rellum_map.end()) {
      idx_rellum_map.insert(std::pair<Int_t,Int_t>(runnum,q));
      //printf("inserting %d => %d\n",runnum,idx_rellum_map.at(runnum));
    }
    else {
      fprintf(stderr,"ERROR: duplicate run %d found in rtree.root, ignoring!\n",runnum);
    }; 
  };


  // build runnum => pol tree entry 
  TString equality;
  for(Int_t q=0; q<pol_tr->GetEntries(); q++) {
    pol_tr->GetEntry(q);
    if(idx_pol_map.find(runnum_pol)==idx_pol_map.end()) {
      idx_pol_map.insert(std::pair<Int_t,Int_t>(runnum_pol,q));
    }
    else {
      fprintf(stderr,"ERROR: duplicate run %d found in pol.root, ignoring!\n",runnum_pol);
      equality = Form("runnum==%d",runnum_pol);
      fprintf(stderr,"       entries=%lld\n",pol_tr->GetEntries(equality.Data()));
    };
  };

};


// return run index in rtree.root file; otherwise returns -1
Int_t RunInfo::Index(Int_t runnum0) {
  iter = idx_rellum_map.find(runnum0);
  if(iter != idx_rellum_map.end()) {
    return iter->second;
  }
  else return -1;
};

// returns run index in counts.root file; otherwise returns -1
Int_t RunInfo::IndexCounts(Int_t runnum0) {
  iter = idx_counts_map.find(runnum0);
  if(iter != idx_counts_map.end()) {
    return iter->second;
  }
  else return -1;
};

// returns run index in pol.root file; otherwise returns -1
Int_t RunInfo::IndexPol(Int_t runnum0) {
  iter = idx_pol_map.find(runnum0);
  if(iter != idx_pol_map.end()) {
    return iter->second;
  }
  else return -1;
};


// return fill for given runnumber; if not in rellum tree returns 0
Int_t RunInfo::GetFill(Int_t runnum0) {
  idx = Index(runnum0);
  if(idx>=0) {
    rtree_tr->GetEntry(idx);
    return fill_rtree;
  }
  else return 0;
};



// return relative luminosity for given run, rellum number, and scaler detector
Float_t RunInfo::Rellum(Int_t runnum0, Int_t rellumi, char * detector) {
  Int_t scai;
  try { scai = sca_idx.at(std::string(detector)); }
  catch ( const std::out_of_range& e) {
    fprintf(stderr,"ERROR: %s is not a valid scaler detector\n",detector);
    return 0;
  };
  
  idx = Index(runnum0);
  if(idx>=0) {
    rtree_tr->GetEntry(idx);
    return R[scai][rellumi];
  }
  else {
    fprintf(stderr,"WARNING: RunInfo::Rellum returning 0 for missing run %d\n",runnum0);
    return 0;
  };
};


// return relative luminosity error for given run, rellum number, and scaler detector
Float_t RunInfo::RellumErr(Int_t runnum0, Int_t rellumi, char * detector) {
  Int_t scai;
  try { scai = sca_idx.at(std::string(detector)); }
  catch ( const std::out_of_range& e) {
    fprintf(stderr,"ERROR: %s is not a valid scaler detector\n",detector);
    return 0;
  };
  
  idx = Index(runnum0);
  if(idx>=0) {
    rtree_tr->GetEntry(idx);
    return R_err[scai][rellumi];
  }
  else {
    fprintf(stderr,"WARNING: RunInfo::RellumErr returning 0 for missing run %d\n",runnum0);
    return 0;
  };
};


// return relative luminosity consistency boolean for given run
Bool_t RunInfo::RellumConsistent(Int_t runnum0)
{
  idx = Index(runnum0);
  if(idx>=0) {
    rtree_tr->GetEntry(idx);
    return isConsistent;
  }
  else {
    fprintf(stderr,"WARNING RunInfo::RellumConsistent returning false for missing run %d\n",runnum0);
    return false;
  };
};


// return blue beam polarization
Float_t RunInfo::BluePol(Int_t runnum0)
{
  idx = IndexPol(runnum0);
  if(idx>=0) {
    pol_tr->GetEntry(idx);
    return b_pol/100.0;
    //return b_pol_avg/100.0;
  }
  else {
    fprintf(stderr,"WARNING: RunInfo::BluePol returning 0 for missing run %d\n",runnum0);
    return 0;
  };
};

// return yellow beam polarization
Float_t RunInfo::YellPol(Int_t runnum0)
{
  idx = IndexPol(runnum0);
  if(idx>=0) {
    pol_tr->GetEntry(idx);
    return y_pol/100.0;
    //return y_pol_avg/100.0;
  }
  else {
    fprintf(stderr,"WARNING: RunInfo::YellPol returning 0 for missing run %d\n",runnum0);
    return 0;
  };
};

// return blue beam polarization error // NOTE -- CURRENTLY RETURNS ONLY ZERO!!!!!!!!!!!!!!11
//  -- polarization error is not needed downstream of asymmetry analysis...
Float_t RunInfo::BluePolErr(Int_t runnum0)
{
  idx = IndexPol(runnum0);
  if(idx>=0) {
    pol_tr->GetEntry(idx);
    return b_pol_err;
  }
  else {
    fprintf(stderr,"WARNING: RunInfo::BluePolErr returning 0 for missing run %d\n",runnum0);
    return 0;
  };
};

// return yellow beam polarization error // NOTE -- CURRENTLY RETURNS ONLY ZERO!!!!!!!!!!!!!!11
//  -- polarization error is not needed downstream of asymmetry analysis...
Float_t RunInfo::YellPolErr(Int_t runnum0)
{
  idx = IndexPol(runnum0);
  if(idx>=0) {
    pol_tr->GetEntry(idx);
    return y_pol_err;
  }
  else {
    fprintf(stderr,"WARNING: RunInfo::YellPolErr returning 0 for missing run %d\n",runnum0);
    return 0;
  };
};


// return blue spin
Int_t RunInfo::BlueSpin(Int_t runnum0, Int_t bXing)
{
  idx = IndexCounts(runnum0);
  if(idx>=0) {
    if(bXing>=0 && bXing<120) return blue_spin_arr[idx][bXing];
    else
    {
      fprintf(stderr,"bXing out of range\n");
      return 0;
    };
  }
  else return 0;
};

// return yellow spin
Int_t RunInfo::YellSpin(Int_t runnum0, Int_t bXing)
{
  idx = IndexCounts(runnum0);
  if(idx>=0) {
    if(bXing>=0 && bXing<120) return yell_spin_arr[idx][bXing];
    else
    {
      fprintf(stderr,"bXing out of range\n");
      return 0;
    };
  }
  else return 0;
};

// return kicked boolean
Bool_t RunInfo::Kicked(Int_t runnum0, Int_t bXing)
{
  idx = IndexCounts(runnum0);
  if(idx>=0) {
    if(bXing>=0 && bXing<120) return kicked_bx_map[idx][bXing];
    else
    {
      fprintf(stderr,"bXing out of range\n");
      return false;
    };
  }
  else return false;
};

// return spin pattern
Int_t RunInfo::Pattern(Int_t runnum0)
{
  idx = Index(runnum0);
  if(idx>=0) {
    rtree_tr->GetEntry(idx);
    return pattern_no;
  }
  else {
    fprintf(stderr,"WARNING RunInfo::Pattern returning 0 for missing run %d\n",runnum0);
    return 0;
  };
};
