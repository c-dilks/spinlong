void TestLevelTwo(TString setn="080Ba") {
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * T = new LevelTwo(RD->env);
  EventClass * ev = new EventClass(RD->env);

  TString filename = TString(RD->env->output_dir)+"/Outputset"+setn+".root";
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("TwoTr");

  Int_t Rnum;
  Int_t L2sum[2];
  UInt_t lastdsm[8];

  tr->SetBranchAddress("Rnum",&Rnum);
  tr->SetBranchAddress("L2sum",L2sum);
  tr->SetBranchAddress("lastdsm",lastdsm);


  // print triggers
  tr->GetEntry(0);
  T->PrintTrigIds(Rnum);
};
