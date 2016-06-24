void DrawThreshes(const char * chosen_thresh="pt") {
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * LT = new LevelTwo(RD->env);
  EventClass * ev = new EventClass(RD->env);


  TFile * infile = new TFile("diagset/setdep.root","READ");
  TTree * tr = (TTree*) infile->Get("threshtr");
  Int_t index,cls,trg;
  Float_t th,the;
  char which_thresh[32];
  tr->SetBranchAddress("index",&index);
  tr->SetBranchAddress("class",&cls);
  tr->SetBranchAddress("trig",&trg);
  tr->SetBranchAddress("thresh",&th);
  tr->SetBranchAddress("thresh_err",&the);
  tr->SetBranchAddress("which_thresh",which_thresh);

  Int_t NCLASSES_tmp = ev->N;
  Int_t NTRIGS_tmp = LT->N;
  
  const Int_t NCLASSES = NCLASSES_tmp;
  const Int_t NTRIGS = NTRIGS_tmp;

  TFile * outfile = new TFile("diagset/thresh.root","RECREATE");
  TGraphErrors * gr[NCLASSES][NTRIGS];
  TString gr_n[NCLASSES][NTRIGS];
  Int_t nn[NCLASSES][NTRIGS];

  for(int c=0; c<NCLASSES; c++) {
    for(int t=0; t<NTRIGS; t++) {
      gr[c][t] = new TGraphErrors();
      gr_n[c][t] = Form("thresh_vs_i__%s__%s",ev->Name(c),(LT->Name(t)).Data());
      gr[c][t]->SetName(gr_n[c][t].Data());
      gr[c][t]->SetTitle(gr_n[c][t].Data());
      gr[c][t]->SetMarkerStyle(kFullCircle);
      gr[c][t]->SetMarkerSize(0.5);
      gr[c][t]->SetMarkerColor(kAzure);
      nn[c][t]=0;
    };
  };

  for(int x=0; x<tr->GetEntries(); x++) {
    tr->GetEntry(x);
    if(!strcmp(chosen_thresh,which_thresh)) {
      gr[cls][trg]->SetPoint(nn[cls][trg],index,th);
      gr[cls][trg]->SetPointError(nn[cls][trg],0,the);
      nn[cls][trg]++;
    };
  };

  for(int c=0; c<NCLASSES; c++) {
    for(int t=0; t<NTRIGS; t++) {
      gr[c][t]->Write();
    };
  };

};

