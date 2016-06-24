// takes kinematic variable and event class and plots the distribution
// of that variable for each run in a 2d histogram; plots are normalised
// to equalise the scale

void KinVarVsRun(char * var="Pt", char * chosen_trig="All")
{
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * T = new LevelTwo(RD->env);
  EventClass * ev = new EventClass(RD->env);
  KinBounds * KB = new KinBounds(RD->env,T,ev);


  const Int_t N_BINS=100;

  Int_t N_CLASS_tmp = ev->N;
  const Int_t N_CLASS = N_CLASS_tmp;

  TChain * tc = new TChain("str");
  TString files = TString(RD->env->redset_dir)+"/Red*.root";
  tc->Add(files.Data());
  //tc->Add("redset/RedOutputset080s1.root");

  Float_t E12,Pt,Eta,Phi,M12,Z,N12,ClIndex;
  Int_t runnum;
  tc->SetBranchAddress("runnum",&runnum);
  tc->SetBranchAddress("E12",&E12);
  tc->SetBranchAddress("Pt",&Pt);
  tc->SetBranchAddress("Eta",&Eta);
  tc->SetBranchAddress("Phi",&Phi);
  tc->SetBranchAddress("M12",&M12);
  tc->SetBranchAddress("Z",&Z);
  tc->SetBranchAddress("N12",&N12);
  tc->SetBranchAddress("ClIndex",&ClIndex);
  tc->SetBranchAddress("L2sum",T->L2sum);
  tc->SetBranchAddress("lastdsm",T->tcu->lastdsm);


  Float_t bin_low,bin_high;
  Float_t * kinvar;

  Bool_t use_threshold = false;

  if(!strcmp(var,"Eta"))
  {
    bin_low = RD->env->EtaLow;
    bin_high = RD->env->EtaHigh;
    kinvar = &Eta;
  }
  else if(!strcmp(var,"Phi"))
  {
    bin_low = RD->env->PhiLow;
    bin_high = RD->env->PhiHigh;
    kinvar = &Phi;
  }
  else if(!strcmp(var,"Pt"))
  {
    bin_low = RD->env->PtLow;
    bin_high = RD->env->PtHigh;
    kinvar = &Pt;
    use_threshold = true;
  }
  else if(!strcmp(var,"E12"))
  {
    bin_low = RD->env->EnLow;
    bin_high = RD->env->EnHigh;
    kinvar = &E12;
    use_threshold = true;
  }
  else if(!strcmp(var,"M12"))
  {
    bin_low = 0;
    bin_high = 1;
    kinvar = &M12;
  }
  else if(!strcmp(var,"Z"))
  {
    bin_low = 0;
    bin_high = 1;
    kinvar = &Z;
  }

  Int_t count=0;
  const Int_t MAX_RUNS = 2000;

  Int_t run_idx[MAX_RUNS]; // "count" --> run index
  Int_t runnum_arr[MAX_RUNS]; // "count" --> run number
  Int_t max_idx=0; // max run index
  Int_t runnum_tmp=0;

  TH1D * h[N_CLASS][MAX_RUNS]; // assumes max number of runs
  char h_name[N_CLASS][MAX_RUNS][32];

  Float_t threshold[N_CLASS][MAX_RUNS];
  Int_t trig_index = T->Index(TString(chosen_trig));


  for(int c=0; c<N_CLASS; c++) {
    for(int r=0; r<MAX_RUNS; r++) {
      threshold[c][r] = -1;
    };
  };

  Bool_t exclude_run;

  //for(Int_t i=0; i<(tc->GetEntries())/100; i++) {
  for(Int_t i=0; i<(tc->GetEntries())/10; i++) {
  //for(Int_t i=0; i<tc->GetEntries(); i++) {
    tc->GetEntry(i);
    T->runnum = runnum;

    if(RD->RellumConsistent(runnum) && RD->BluePol(runnum)>0 && RD->YellPol(runnum)>0)
    {
      ev->SetKinematics(runnum,E12,Pt,Eta,Phi,M12,Z,N12,ClIndex);
      if(!(ev->ExcludedRun()))
      {
        if(runnum!=runnum_tmp)
        {
          if(runnum_tmp!=0) count++;
          for(Int_t c=0; c<N_CLASS; c++)
          {
            sprintf(h_name[c][count],"h_%d_%s",runnum,ev->Name(c));
            h[c][count] = new TH1D(h_name[c][count],h_name[c][count],N_BINS,bin_low,bin_high);
          };
          printf("%s %d %d\n",var,count,runnum);
          runnum_tmp = runnum;
          run_idx[count] = RD->Index(runnum);
          runnum_arr[count] = runnum;
          max_idx = (run_idx[count]>max_idx) ? run_idx[count]:max_idx;
          if(use_threshold) {
            for(Int_t c=0; c<N_CLASS; c++) {
              if(!strcmp(var,"Pt"))
                threshold[c][run_idx[count]] = KB->PtThreshLow(run_idx[count],c,trig_index);
              else if(!strcmp(var,"E12")) 
                threshold[c][run_idx[count]] = KB->EnThreshLow(run_idx[count],c,trig_index);

              //if(c==1) printf("%d %d %f\n",runnum,run_idx[count],threshold[c][run_idx[count]]);
            };
          };
        };
        for(Int_t c=0; c<N_CLASS; c++)
        {
          if(ev->Valid(c) && T->Fired(TString(chosen_trig)))
            h[c][count]->Fill(*kinvar);
        };
      };
    };
  };

  char h2_t[N_CLASS][64];
  char h2_n[N_CLASS][32];
  TH2D * h2[N_CLASS];
  TGraphErrors * g[N_CLASS];
  char g_n[N_CLASS][64];
  char g_t[N_CLASS][64];
  Int_t g_c[N_CLASS];
  TGraphErrors * gt[N_CLASS];
  char gt_n[N_CLASS][64];
  char gt_t[N_CLASS][64];
  Int_t gt_c[N_CLASS];
  Float_t mean,rms;
  for(Int_t c=0; c<N_CLASS; c++)
  {
    sprintf(h2_t[c],"%s %s vs. run index",ev->Title(c),var);
    sprintf(h2_n[c],"h2_%s",ev->Name(c));
    //h2[c] = new TH2D(h2_n[c],h2_t[c],count+1,0,count+1,N_BINS,bin_low,bin_high);
    h2[c] = new TH2D(h2_n[c],h2_t[c],max_idx+1,0,max_idx+1,N_BINS,bin_low,bin_high);

    sprintf(g_t[c],"%s average %s vs. run index",ev->Title(c),var);
    sprintf(g_n[c],"g_%s",ev->Name(c));
    g[c] = new TGraphErrors();
    g_c[c] = 0;

    g[c]->SetName(g_n[c]);
    g[c]->SetTitle(g_t[c]);
    g[c]->SetMarkerStyle(kFullCircle);
    g[c]->SetMarkerColor(kViolet);

    sprintf(gt_t[c],"%s %s threshold vs. run index",ev->Title(c),var);
    sprintf(gt_n[c],"gt_%s",ev->Name(c));
    gt[c] = new TGraphErrors();
    gt_c[c] = 0;

    gt[c]->SetName(gt_n[c]);
    gt[c]->SetTitle(gt_t[c]);
    gt[c]->SetMarkerStyle(kFullCircle);
    gt[c]->SetMarkerColor(kAzure);

    Float_t binpos;
    Int_t binn;
    for(Int_t k=0; k<count+1; k++)
    {
      if(h[c][k]!=NULL) {
        h[c][k]->Scale(1/h[c][k]->Integral());
        if(h[c][k]->GetEntries()>0) {
          mean = h[c][k]->GetMean();
          rms = h[c][k]->GetRMS();
          if(!(mean>=bin_low && mean<=bin_high)) { mean=10000; rms=10000; };
        }
        else { mean=10000; rms=10000; };
        if(mean<10000) {
          g[c]->SetPoint(g_c[c],k,mean);
          //g[c]->SetPoint(g_c[c],run_idx[k],mean);
          g[c]->SetPointError(g_c[c],0,rms);
          g_c[c]++;
        };

        if(threshold[c][run_idx[k]]>=0) {
          gt[c]->SetPoint(gt_c[c],k,threshold[c][run_idx[k]]);
          //gt[c]->SetPoint(gt_c[c],run_idx[k],threshold[c][run_idx[k]]);
          gt[c]->SetPointError(gt_c[c],0,0);
          gt_c[c]++;
        };

        for(Int_t b=0; b<h[c][k]->GetNbinsX(); b++)
        {
          binpos = h[c][k]->GetBinCenter(b);
          binn = h2[c]->FindBin(k,binpos);
          //binn = h2[c]->FindBin(run_idx[k],binpos);
          h2[c]->SetBinContent(binn,h[c][k]->GetBinContent(b));
          
          printf("runnum=%d h[c=%d][k=%d]->GetEntries()=%d\n",
            runnum_arr[k],c,k,h[c][k]->GetEntries());
        };
      };
    };
  };

  // build tree with kinvar_vs_run's local run index, to be used if you
  // see something strange in the distributions and want to characterize
  // which run it is
  TString outfile_n = Form("kinvarset/%s_%s_vs_run.root",chosen_trig,var);
  TFile * outfile = new TFile(outfile_n.Data(),"RECREATE");
  TTree * runindex_tr = new TTree();
  Int_t i0,runnum0;
  runindex_tr->Branch("idx",&i0,"idx/I");
  runindex_tr->Branch("runnum",&runnum0,"runnum/I");
  for(Int_t k=0; k<count+1; k++) {
    i0 = k;
    runnum0 = runnum_arr[k];
    runindex_tr->Fill();
  };


  TCanvas * canv = new TCanvas("canv","canv",1000,N_CLASS*400);
  gStyle->SetOptStat(0);
  canv->Divide(1,N_CLASS);
  for(Int_t c=0; c<N_CLASS; c++)
  {
    canv->cd(c+1);
    h2[c]->Draw("colz");
    h2[c]->Write();
  };

  canv->Write();

  for(Int_t c=0; c<N_CLASS; c++) g[c]->Write();
  if(use_threshold) { for(Int_t c=0; c<N_CLASS; c++) gt[c]->Write(); };
  runindex_tr->Write("runindex_tr");
};

