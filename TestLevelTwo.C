void TestLevelTwo(TString setn="080Ba") {
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * T = new LevelTwo(RD->env);
  EventClass * ev = new EventClass(RD->env);

  T->debug = true;
  T->tcu->debug = false;

  TString filename = TString(RD->env->output_dir)+"/Outputset"+setn+".root";
  TFile * infile = new TFile(filename,"READ");
  TTree * tr = (TTree*) infile->Get("TwoTr");

  Int_t L2sum[2];
  UInt_t lastdsm[8];
  Int_t evnum,ClIndex;

  tr->SetBranchAddress("Rnum",&(T->runnum));
  tr->SetBranchAddress("L2sum",T->L2sum);
  tr->SetBranchAddress("lastdsm",T->tcu->lastdsm);
  tr->SetBranchAddress("Evnum",&evnum);
  tr->SetBranchAddress("ClIndex",&ClIndex);



  // event loop
  Int_t ENT = tr->GetEntries();
  ENT = 100;
  for(int i=0; i<ENT; i++) {
    tr->GetEntry(i);
    if(ClIndex==0) {
      printf("\n");
      T->PrintTrigIds();
      printf("\n");
      T->PrintVars();
      printf("evnum=%d\n",evnum);
      for(int n=0; n<T->N; n++) T->Fired(n);
    };
  };
};
