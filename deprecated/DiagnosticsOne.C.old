// plots kinematic correlation and trigger overlaps for all FMS triggers
//
// this script only stream in one redset file and is to be run
// within condor

void DiagnosticsOne(const char * infile_name = "RedOutputset122Ca.root")
{
  // if enableOverlap=true, fill overlap matrices; may cause slow-down
  const Bool_t enableOverlap = true;

  const Int_t NBINS=75; // NUMBER OF BINS (default 400)
  const Int_t NBINS_RDIST=100; // number of bins for variable vs. run index plots (default 100)
  const Int_t MAXRUNS=50; // arbitrary max number of runs in redset file 

  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * LT = new LevelTwo(RD->env);
  EventClass * ev = new EventClass(RD->env,true);


  // open tree
  char infile_full_name[256];
  sprintf(infile_full_name,"%s/%s",RD->env->redset_dir,infile_name);
  TFile * infile = new TFile(infile_full_name,"READ");
  TTree * tr = (TTree*) infile->Get("str");

  Float_t E12,Pt,Eta,Phi,M12,Z,b_pol,y_pol;
  Int_t ClIndex;
  Bool_t kicked,isConsistent;
  Int_t runnum,bx;
  Float_t N12;

  Float_t E12_min,Pt_min,Eta_min,Phi_min;
  Float_t E12_max,Pt_max,Eta_max,Phi_max;
  E12_min=Pt_min=Eta_min=Phi_min=1000;
  E12_max=Pt_max=Eta_max=Phi_max=0;

  tr->SetBranchAddress("runnum",&runnum); // LT->runnum = runnum called directly after tr->GetEntry()
  tr->SetBranchAddress("Bunchid7bit",&bx);
  tr->SetBranchAddress("E12",&E12);
  tr->SetBranchAddress("Pt",&Pt);
  tr->SetBranchAddress("Eta",&Eta);
  tr->SetBranchAddress("Phi",&Phi);
  tr->SetBranchAddress("M12",&M12);
  tr->SetBranchAddress("Z",&Z);
  tr->SetBranchAddress("N12",&N12);
  tr->SetBranchAddress("ClIndex",&ClIndex);
  tr->SetBranchAddress("L2sum",LT->L2sum);
  tr->SetBranchAddress("lastdsm",LT->tcu->lastdsm);



  // trigger list
  Int_t N_TRIG_tmp = LT->N;
  const Int_t N_TRIG = N_TRIG_tmp;

  // event classes
  Int_t N_CLASS_tmp = ev->N;
  const Int_t N_CLASS = N_CLASS_tmp;


  // get kinematic ranges using the current binning set by Bin_Splitter.C
  Phi_min = RD->env->PhiLow;
  Phi_max = RD->env->PhiHigh;
  Eta_min = RD->env->EtaLow;
  Eta_max = RD->env->EtaHigh;
  E12_min = RD->env->EnLow;
  E12_max = RD->env->EnHigh;
  Pt_min = RD->env->PtLow;
  Pt_max = RD->env->PtHigh;

  // override binning defaults
  E12_min = 0;
  E12_max = 100;
  Pt_min = 0; 
  Pt_max = 10;

  
  // define histograms
  TH2D * pt_vs_eta[N_CLASS][N_TRIG];
  TH2D * en_vs_eta[N_CLASS][N_TRIG];
  TH2D * pt_vs_phi[N_CLASS][N_TRIG];
  TH2D * en_vs_phi[N_CLASS][N_TRIG];
  TH2D * eta_vs_phi[N_CLASS][N_TRIG];
  TH2D * pt_vs_en[N_CLASS][N_TRIG];
  TH2D * y_vs_x[N_CLASS][N_TRIG];

  TH2D * z_vs_eta[N_CLASS][N_TRIG];
  TH2D * z_vs_phi[N_CLASS][N_TRIG];
  TH2D * mass_vs_en[N_CLASS][N_TRIG];
  TH2D * mass_vs_pt[N_CLASS][N_TRIG];
  TH1D * mass_dist[N_CLASS][N_TRIG];
  TH1D * z_dist[N_CLASS][N_TRIG];

  TH1D * pt_rdist[N_CLASS][N_TRIG][MAXRUNS];
  TH1D * en_rdist[N_CLASS][N_TRIG][MAXRUNS];
  TH1D * eta_rdist[N_CLASS][N_TRIG][MAXRUNS];
  TH1D * phi_rdist[N_CLASS][N_TRIG][MAXRUNS];
  TH1D * mass_rdist[N_CLASS][N_TRIG][MAXRUNS];
  TH1D * z_rdist[N_CLASS][N_TRIG][MAXRUNS];

  char pt_vs_eta_n[N_CLASS][N_TRIG][64];
  char en_vs_eta_n[N_CLASS][N_TRIG][64];
  char pt_vs_phi_n[N_CLASS][N_TRIG][64];
  char en_vs_phi_n[N_CLASS][N_TRIG][64];
  char eta_vs_phi_n[N_CLASS][N_TRIG][64];
  char pt_vs_en_n[N_CLASS][N_TRIG][64];
  char y_vs_x_n[N_CLASS][N_TRIG][64];
  char z_vs_eta_n[N_CLASS][N_TRIG][64];
  char z_vs_phi_n[N_CLASS][N_TRIG][64];
  char mass_vs_en_n[N_CLASS][N_TRIG][64];
  char mass_vs_pt_n[N_CLASS][N_TRIG][64];
  char mass_dist_n[N_CLASS][N_TRIG][64];
  char z_dist_n[N_CLASS][N_TRIG][64];

  char pt_vs_eta_t[N_CLASS][N_TRIG][256];
  char en_vs_eta_t[N_CLASS][N_TRIG][256];
  char pt_vs_phi_t[N_CLASS][N_TRIG][256];
  char en_vs_phi_t[N_CLASS][N_TRIG][256];
  char eta_vs_phi_t[N_CLASS][N_TRIG][256];
  char pt_vs_en_t[N_CLASS][N_TRIG][256];
  char y_vs_x_t[N_CLASS][N_TRIG][256];
  char z_vs_eta_t[N_CLASS][N_TRIG][256];
  char z_vs_phi_t[N_CLASS][N_TRIG][256];
  char mass_vs_en_t[N_CLASS][N_TRIG][256];
  char mass_vs_pt_t[N_CLASS][N_TRIG][256];
  char mass_dist_t[N_CLASS][N_TRIG][256];
  char z_dist_t[N_CLASS][N_TRIG][256];

  char pt_rdist_n[N_CLASS][N_TRIG][MAXRUNS][64];
  char en_rdist_n[N_CLASS][N_TRIG][MAXRUNS][64];
  char eta_rdist_n[N_CLASS][N_TRIG][MAXRUNS][64];
  char phi_rdist_n[N_CLASS][N_TRIG][MAXRUNS][64];
  char mass_rdist_n[N_CLASS][N_TRIG][MAXRUNS][64];
  char z_rdist_n[N_CLASS][N_TRIG][MAXRUNS][64];

  char pt_rdist_n_tmp[N_CLASS][N_TRIG][MAXRUNS][64];
  char en_rdist_n_tmp[N_CLASS][N_TRIG][MAXRUNS][64];
  char eta_rdist_n_tmp[N_CLASS][N_TRIG][MAXRUNS][64];
  char phi_rdist_n_tmp[N_CLASS][N_TRIG][MAXRUNS][64];
  char mass_rdist_n_tmp[N_CLASS][N_TRIG][MAXRUNS][64];
  char z_rdist_n_tmp[N_CLASS][N_TRIG][MAXRUNS][64];

  Float_t low_mass,high_mass;

  for(Int_t t=0; t<N_TRIG; t++)
  {
    for(Int_t c=0; c<N_CLASS; c++)
    {
      low_mass=0; // default
      high_mass=1; // default
      if(c==ev->Idx("thr")) 
      {
        low_mass=0.7; // to match cutoff in EventClass.cxx
        high_mass=6;
      }
      else if(c==ev->Idx("etm")) 
      {
        low_mass=0.3;
      }
      /*
      else if(c==ev->Idx("jps"))
      {
        low_mass=0;
        high_mass=6;
      };
      */

      sprintf(pt_vs_eta_n[c][t],"%s_%s_pt_vs_eta",(LT->Name(t)).Data(),ev->Name(c));
      sprintf(en_vs_eta_n[c][t],"%s_%s_en_vs_eta",(LT->Name(t)).Data(),ev->Name(c));
      sprintf(pt_vs_phi_n[c][t],"%s_%s_pt_vs_phi",(LT->Name(t)).Data(),ev->Name(c));
      sprintf(en_vs_phi_n[c][t],"%s_%s_en_vs_phi",(LT->Name(t)).Data(),ev->Name(c));
      sprintf(eta_vs_phi_n[c][t],"%s_%s_eta_vs_phi",(LT->Name(t)).Data(),ev->Name(c));
      sprintf(pt_vs_en_n[c][t],"%s_%s_pt_vs_en",(LT->Name(t)).Data(),ev->Name(c));
      sprintf(y_vs_x_n[c][t],"%s_%s_y_vs_x",(LT->Name(t)).Data(),ev->Name(c));
      sprintf(z_vs_eta_n[c][t],"%s_%s_z_vs_eta",(LT->Name(t)).Data(),ev->Name(c));
      sprintf(z_vs_phi_n[c][t],"%s_%s_z_vs_phi",(LT->Name(t)).Data(),ev->Name(c));
      sprintf(mass_vs_en_n[c][t],"%s_%s_mass_vs_en",(LT->Name(t)).Data(),ev->Name(c));
      sprintf(mass_vs_pt_n[c][t],"%s_%s_mass_vs_pt",(LT->Name(t)).Data(),ev->Name(c));
      sprintf(mass_dist_n[c][t],"%s_%s_mass_dist",(LT->Name(t)).Data(),ev->Name(c));
      sprintf(z_dist_n[c][t],"%s_%s_z_dist",(LT->Name(t)).Data(),ev->Name(c));

      sprintf(pt_vs_eta_t[c][t],"%s %s --- p_{T} vs. #eta;#eta;p_{T}",
        (LT->Name(t)).Data(),ev->Title(c));
      sprintf(en_vs_eta_t[c][t],"%s %s --- E vs. #eta;#eta;E",
        (LT->Name(t)).Data(),ev->Title(c));
      sprintf(pt_vs_phi_t[c][t],"%s %s --- p_{T} vs. #phi;#phi;p_{T}",
        (LT->Name(t)).Data(),ev->Title(c));
      sprintf(en_vs_phi_t[c][t],"%s %s --- E vs. #phi;#phi;E",
        (LT->Name(t)).Data(),ev->Title(c));
      sprintf(eta_vs_phi_t[c][t],"%s %s --- #eta vs. #phi;#phi;#eta",
        (LT->Name(t)).Data(),ev->Title(c));
      sprintf(pt_vs_en_t[c][t],"%s %s --- p_{T} vs. E;E;p_{T}",
        (LT->Name(t)).Data(),ev->Title(c));
      sprintf(y_vs_x_t[c][t],"%s %s --- y vs. x;x;y",
        (LT->Name(t)).Data(),ev->Title(c));

      //printf("here\n");
      sprintf(z_vs_eta_t[c][t],"%s %s --- Z vs. #eta (%s cuts w/o Z-cut);#eta;Z",
        (LT->Name(t)).Data(),ev->Title(c),ev->Title(c));
      sprintf(z_vs_phi_t[c][t],"%s %s --- Z vs. #phi (%s cuts w/o Z-cut);#phi;Z",
        (LT->Name(t)).Data(),ev->Title(c),ev->Title(c));
      sprintf(mass_vs_en_t[c][t],"%s %s --- M vs. E (%s cuts w/o M-cut);E;M",
        (LT->Name(t)).Data(),ev->Title(c),ev->Title(c));
      sprintf(mass_vs_pt_t[c][t],"%s %s --- M vs. p_{T} (%s cuts w/o M-cut);p_{T};M",
        (LT->Name(t)).Data(),ev->Title(c),ev->Title(c));
      //printf("VVVVV\n");
      sprintf(mass_dist_t[c][t],"%s %s --- M distribution (%s cuts w/o M-cut);M",
        (LT->Name(t)).Data(),ev->Title(c),ev->Title(c));
      //printf("%s\n",mass_dist_t[c][t]);
      sprintf(z_dist_t[c][t],"%s %s --- Z distribution (%s cuts w/o Z-cut);Z",
        (LT->Name(t)).Data(),ev->Title(c),ev->Title(c));
      //printf("^^^^^\n");

      pt_vs_eta[c][t] = new TH2D(pt_vs_eta_n[c][t],pt_vs_eta_t[c][t],
        NBINS,Eta_min,Eta_max,NBINS,Pt_min,Pt_max);
      en_vs_eta[c][t] = new TH2D(en_vs_eta_n[c][t],en_vs_eta_t[c][t],
        NBINS,Eta_min,Eta_max,NBINS,E12_min,E12_max);
      pt_vs_phi[c][t] = new TH2D(pt_vs_phi_n[c][t],pt_vs_phi_t[c][t],
        NBINS,Phi_min,Phi_max,NBINS,Pt_min,Pt_max);
      en_vs_phi[c][t] = new TH2D(en_vs_phi_n[c][t],en_vs_phi_t[c][t],
        NBINS,Phi_min,Phi_max,NBINS,E12_min,E12_max);
      eta_vs_phi[c][t] = new TH2D(eta_vs_phi_n[c][t],eta_vs_phi_t[c][t],
        NBINS,Phi_min,Phi_max,NBINS,Eta_min,Eta_max);
      pt_vs_en[c][t] = new TH2D(pt_vs_en_n[c][t],pt_vs_en_t[c][t],
        NBINS,E12_min,E12_max,NBINS,Pt_min,Pt_max);
      y_vs_x[c][t] = new TH2D(y_vs_x_n[c][t],y_vs_x_t[c][t],
        NBINS,-110,110,NBINS,-110,110);
      z_vs_eta[c][t] = new TH2D(z_vs_eta_n[c][t],z_vs_eta_t[c][t],
        NBINS,Eta_min,Eta_max,NBINS,0,1);
      z_vs_phi[c][t] = new TH2D(z_vs_phi_n[c][t],z_vs_phi_t[c][t],
        NBINS,Phi_min,Phi_max,NBINS,0,1);
      mass_vs_en[c][t] = new TH2D(mass_vs_en_n[c][t],mass_vs_en_t[c][t],
        NBINS,E12_min,E12_max,NBINS,low_mass,high_mass);
      mass_vs_pt[c][t] = new TH2D(mass_vs_pt_n[c][t],mass_vs_pt_t[c][t],
        NBINS,Pt_min,Pt_max,NBINS,low_mass,high_mass);

      z_dist[c][t] = new TH1D(z_dist_n[c][t],z_dist_t[c][t],NBINS,0,1);
      mass_dist[c][t] = new TH1D(mass_dist_n[c][t],mass_dist_t[c][t],NBINS,low_mass,high_mass);

      //printf("----- c=%d t=%d %s %s\n",c,t,ev->Name(c),(LT->Name(t)).Data());
      for(Int_t ru=0; ru<MAXRUNS; ru++)
      {
        sprintf(pt_rdist_n[c][t][ru],"%s_%s_pt_rdist",(LT->Name(t)).Data(),ev->Name(c));
        sprintf(en_rdist_n[c][t][ru],"%s_%s_en_rdist",(LT->Name(t)).Data(),ev->Name(c));
        sprintf(eta_rdist_n[c][t][ru],"%s_%s_eta_rdist",(LT->Name(t)).Data(),ev->Name(c));
        sprintf(phi_rdist_n[c][t][ru],"%s_%s_phi_rdist",(LT->Name(t)).Data(),ev->Name(c));
        sprintf(mass_rdist_n[c][t][ru],"%s_%s_mass_rdist",(LT->Name(t)).Data(),ev->Name(c));
        sprintf(z_rdist_n[c][t][ru],"%s_%s_z_rdist",(LT->Name(t)).Data(),ev->Name(c));

        sprintf(pt_rdist_n_tmp[c][t][ru],"%s_%d",pt_rdist_n[c][t][ru],ru);
        sprintf(en_rdist_n_tmp[c][t][ru],"%s_%d",en_rdist_n[c][t][ru],ru);
        sprintf(eta_rdist_n_tmp[c][t][ru],"%s_%d",eta_rdist_n[c][t][ru],ru);
        sprintf(phi_rdist_n_tmp[c][t][ru],"%s_%d",phi_rdist_n[c][t][ru],ru);
        sprintf(mass_rdist_n_tmp[c][t][ru],"%s_%d",mass_rdist_n[c][t][ru],ru);
        sprintf(z_rdist_n_tmp[c][t][ru],"%s_%d",z_rdist_n[c][t][ru],ru);

        pt_rdist[c][t][ru] = new TH1D(pt_rdist_n_tmp[c][t][ru],pt_rdist_n_tmp[c][t][ru],
          NBINS_RDIST,Pt_min,Pt_max);
        en_rdist[c][t][ru] = new TH1D(en_rdist_n_tmp[c][t][ru],en_rdist_n_tmp[c][t][ru],
          NBINS_RDIST,E12_min,E12_max);
        eta_rdist[c][t][ru] = new TH1D(eta_rdist_n_tmp[c][t][ru],eta_rdist_n_tmp[c][t][ru],
          NBINS_RDIST,Eta_min,Eta_max);
        phi_rdist[c][t][ru] = new TH1D(phi_rdist_n_tmp[c][t][ru],phi_rdist_n_tmp[c][t][ru],
          NBINS_RDIST,Phi_min,Phi_max);
        mass_rdist[c][t][ru] = new TH1D(mass_rdist_n_tmp[c][t][ru],mass_rdist_n_tmp[c][t][ru],
          NBINS_RDIST,low_mass,high_mass);
        z_rdist[c][t][ru] = new TH1D(z_rdist_n_tmp[c][t][ru],z_rdist_n_tmp[c][t][ru],
          NBINS_RDIST,0,1);

      };
    };
  };


  // energy-dependent mass distributions for pions 
  TH1D * mass_dist_for_enbin[N_TRIG][10];
  char mass_dist_for_enbin_n[N_TRIG][10][64];
  char mass_dist_for_enbin_t[N_TRIG][10][256];
  for(Int_t t=0; t<N_TRIG; t++)
  {
    for(Int_t ee=0; ee<10; ee++)
    {
      sprintf(mass_dist_for_enbin_n[t][ee],"%s_mass_dist_for_enbin%d",(LT->Name(t)).Data(),ee);
      sprintf(mass_dist_for_enbin_t[t][ee],
        "%s M_{#gamma#gamma} distribution for E_{#gamma#gamma}#in[%d,%d) GeV",
        (LT->Name(t)).Data(),ee*10,(ee+1)*10);
      mass_dist_for_enbin[t][ee] = 
        new TH1D(mass_dist_for_enbin_n[t][ee],mass_dist_for_enbin_t[t][ee],NBINS,0,1);
    };
  };

  // pt-dependent mass distributions for pions
  TH1D * mass_dist_for_ptbin[N_TRIG][10];
  char mass_dist_for_ptbin_n[N_TRIG][10][64];
  char mass_dist_for_ptbin_t[N_TRIG][10][256];
  for(Int_t t=0; t<N_TRIG; t++)
  {
    for(Int_t pp=0; pp<10; pp++)
    {
      sprintf(mass_dist_for_ptbin_n[t][pp],"%s_mass_dist_for_ptbin%d",(LT->Name(t)).Data(),pp);
      sprintf(mass_dist_for_ptbin_t[t][pp],
        "%s M_{#gamma#gamma} distribution for p_{T}#in[%d,%d) GeV/c",
        (LT->Name(t)).Data(),pp,pp+1);
      mass_dist_for_ptbin[t][pp] = 
        new TH1D(mass_dist_for_ptbin_n[t][pp],mass_dist_for_ptbin_t[t][pp],NBINS,0,1);
    };
  };


  // trigger distribution
  char trig_dist_t[64];
  sprintf(trig_dist_t,"Trigger Counts");
  TH1D * trig_dist = new TH1D("trig_dist",trig_dist_t,N_TRIG,0,N_TRIG);
  for(Int_t t=0; t<N_TRIG; t++) trig_dist->GetXaxis()->SetBinLabel(t+1,(LT->Name(t)).Data());


  TH2D * trig_fms_mix[N_CLASS]; // overlap matrix of FMS triggers
                                 
  char trig_fms_mix_n[N_CLASS][64];
  char trig_fms_mix_t[N_CLASS][256];

  for(Int_t c=0; c<N_CLASS; c++)
  {
    sprintf(trig_fms_mix_n[c],"%s_trig_fms_mix",ev->Name(c));

    sprintf(trig_fms_mix_t[c],"%s FMS trigger overlap matrix",ev->Title(c));

    trig_fms_mix[c] = new TH2D(trig_fms_mix_n[c],trig_fms_mix_t[c],
      N_TRIG,0,N_TRIG,N_TRIG,0,N_TRIG);

    for(Int_t t=0; t<N_TRIG; t++) 
    {
      trig_fms_mix[c]->GetXaxis()->SetBinLabel(t+1,(LT->Name(t)).Data());
      trig_fms_mix[c]->GetYaxis()->SetBinLabel(t+1,(LT->Name(t)).Data());
    };
  };


  // fill histograms
  Int_t runnum_tmp=0;
  Int_t runcount=-1;
  Double_t pt_bc,en_bc,eta_bc,phi_bc,z_bc,mass_bc;
  Int_t pt_bn,en_bn,eta_bn,phi_bn,z_bn,mass_bn;

  Int_t ENT = tr->GetEntries();
  //ENT = 100000; // uncomment to do a short loop for testing
  for(Int_t x=0; x<ENT; x++)
  {
    if((x%100000)==0) printf("filling histograms: %.2f%%\n",100*((Float_t)x)/((Float_t)ENT));
    tr->GetEntry(x);

    LT->runnum = runnum;

    kicked = RD->Kicked(runnum,bx);

    // get new polarisation and check rellum
    // also fill kinematic vs. run plots
    if(runnum!=runnum_tmp)
    {
      b_pol = RD->BluePol(runnum);
      y_pol = RD->YellPol(runnum);
      isConsistent = RD->RellumConsistent(runnum);

      runcount++;
      printf(">>> %d <<<\n",runnum);
      for(Int_t t=0; t<N_TRIG; t++)
      {
        for(Int_t c=0; c<N_CLASS; c++)
        {
          sprintf(pt_rdist_n[c][t][runcount],"%s_%d",pt_rdist_n[c][t][runcount],runnum);
          sprintf(en_rdist_n[c][t][runcount],"%s_%d",en_rdist_n[c][t][runcount],runnum);
          sprintf(eta_rdist_n[c][t][runcount],"%s_%d",eta_rdist_n[c][t][runcount],runnum);
          sprintf(phi_rdist_n[c][t][runcount],"%s_%d",phi_rdist_n[c][t][runcount],runnum);
          sprintf(z_rdist_n[c][t][runcount],"%s_%d",z_rdist_n[c][t][runcount],runnum);
          sprintf(mass_rdist_n[c][t][runcount],"%s_%d",mass_rdist_n[c][t][runcount],runnum);

          pt_rdist[c][t][runcount]->SetName(pt_rdist_n[c][t][runcount]);
          en_rdist[c][t][runcount]->SetName(en_rdist_n[c][t][runcount]);
          eta_rdist[c][t][runcount]->SetName(eta_rdist_n[c][t][runcount]);
          phi_rdist[c][t][runcount]->SetName(phi_rdist_n[c][t][runcount]);
          z_rdist[c][t][runcount]->SetName(z_rdist_n[c][t][runcount]);
          mass_rdist[c][t][runcount]->SetName(mass_rdist_n[c][t][runcount]);

          pt_rdist[c][t][runcount]->SetTitle(pt_rdist_n[c][t][runcount]);
          en_rdist[c][t][runcount]->SetTitle(en_rdist_n[c][t][runcount]);
          eta_rdist[c][t][runcount]->SetTitle(eta_rdist_n[c][t][runcount]);
          phi_rdist[c][t][runcount]->SetTitle(phi_rdist_n[c][t][runcount]);
          z_rdist[c][t][runcount]->SetTitle(z_rdist_n[c][t][runcount]);
          mass_rdist[c][t][runcount]->SetTitle(mass_rdist_n[c][t][runcount]);
        };
      };
      runnum_tmp=runnum;
    }

    // rellum / pol cut
    if( kicked==0 && isConsistent==1 && b_pol>0 && y_pol>0)
    {
      // fill fms trigger plot
      for(Int_t t=0; t<N_TRIG; t++)
      {
        if(LT->Fired(t))
          trig_dist->Fill(t);
      };

      // below here we have if statements for each event class; within each
      // if block, we have trigger loops, since it's more efficient to first pick the event and 
      // then loop through all triggers

      ev->SetKinematics(runnum,E12,Pt,Eta,Phi,M12,Z,N12,ClIndex);
      for(Int_t c=0; c<N_CLASS; c++)
      {
        if(c==ev->Idx("dpi")) continue; // don't bother with this yet
        // fill main kinematic correlation plots
        if(ev->Valid(c))
        {
          ev->FiducialGeom(Eta,Phi,0); // compute x and y coordinates
          for(Int_t t=0; t<N_TRIG; t++)
          {
            if(LT->Fired(t))
            {
              pt_vs_eta[c][t]->Fill(Eta,Pt);
              en_vs_eta[c][t]->Fill(Eta,E12);
              pt_vs_phi[c][t]->Fill(Phi,Pt);
              en_vs_phi[c][t]->Fill(Phi,E12);
              eta_vs_phi[c][t]->Fill(Phi,Eta);
              pt_vs_en[c][t]->Fill(E12,Pt);
              y_vs_x[c][t]->Fill(ev->Xd,ev->Yd);

              pt_rdist[c][t][runcount]->Fill(Pt);
              en_rdist[c][t][runcount]->Fill(E12);
              eta_rdist[c][t][runcount]->Fill(Eta);
              phi_rdist[c][t][runcount]->Fill(Phi);

              // fill trigger overlap matrices
              if(enableOverlap)
              {
                for(Int_t tt=0; tt<N_TRIG; tt++)
                {
                  if(LT->Fired(tt))
                    trig_fms_mix[c]->Fill(t,tt);
                };
              };
            };
          };
        };

        // fill mass plots
        if(c!=ev->Idx("sph") && ev->ValidWithoutMcut(c))
        {
          for(Int_t t=0; t<N_TRIG; t++)
          {
            if(LT->Fired(t))
            {
              mass_dist[c][t]->Fill(M12);
              mass_vs_en[c][t]->Fill(E12,M12);
              mass_vs_pt[c][t]->Fill(Pt,M12);
              mass_rdist[c][t][runcount]->Fill(M12);
              // fill kinematic-dependent mass distributions for pions
              if(c==ev->Idx("pi0"))
              {
                for(Int_t ee=0; ee<10; ee++)
                {
                  if(E12>=(ee*10) && E12<((ee+1)*10)) mass_dist_for_enbin[t][ee]->Fill(M12);
                };
                for(Int_t pp=0; pp<10; pp++)
                {
                  if(Pt>=pp && Pt<(pp+1)) mass_dist_for_ptbin[t][pp]->Fill(M12);
                };
              };
            };
          };
        };

        // fill z plots
        if((c==ev->Idx("pi0") || c==ev->Idx("etm") /*|| c==ev->Idx("jps")*/) && ev->ValidWithoutZcut(c))
        {
          for(Int_t t=0; t<N_TRIG; t++)
          {
            if(LT->Fired(t))
            {
              z_vs_eta[c][t]->Fill(Eta,Z);
              z_vs_phi[c][t]->Fill(Phi,Z);
              z_dist[c][t]->Fill(Z);
              z_rdist[c][t][runcount]->Fill(Z);
            };
          };
        };
      };
    };
  };


  // build TObjArrays
  TObjArray * pt_vs_eta_arr[N_CLASS];
  TObjArray * en_vs_eta_arr[N_CLASS];
  TObjArray * pt_vs_phi_arr[N_CLASS];
  TObjArray * en_vs_phi_arr[N_CLASS];
  TObjArray * eta_vs_phi_arr[N_CLASS];
  TObjArray * pt_vs_en_arr[N_CLASS];
  TObjArray * y_vs_x_arr[N_CLASS];
  TObjArray * z_vs_eta_arr[N_CLASS];
  TObjArray * z_vs_phi_arr[N_CLASS];
  TObjArray * mass_vs_en_arr[N_CLASS];
  TObjArray * mass_vs_pt_arr[N_CLASS];
  TObjArray * mass_dist_arr[N_CLASS];
  TObjArray * z_dist_arr[N_CLASS];
  TObjArray * mass_dist_for_enbin_arr[10];
  TObjArray * mass_dist_for_ptbin_arr[10];

  TObjArray * pt_rdist_arr[N_CLASS][MAXRUNS];
  TObjArray * en_rdist_arr[N_CLASS][MAXRUNS];
  TObjArray * eta_rdist_arr[N_CLASS][MAXRUNS];
  TObjArray * phi_rdist_arr[N_CLASS][MAXRUNS];
  TObjArray * z_rdist_arr[N_CLASS][MAXRUNS];
  TObjArray * mass_rdist_arr[N_CLASS][MAXRUNS];

  char pt_vs_eta_arr_n[N_CLASS][32];
  char en_vs_eta_arr_n[N_CLASS][32];
  char pt_vs_phi_arr_n[N_CLASS][32];
  char en_vs_phi_arr_n[N_CLASS][32];
  char eta_vs_phi_arr_n[N_CLASS][32];
  char pt_vs_en_arr_n[N_CLASS][32];
  char y_vs_x_arr_n[N_CLASS][32];
  char z_vs_eta_arr_n[N_CLASS][32];
  char z_vs_phi_arr_n[N_CLASS][32];
  char mass_vs_en_arr_n[N_CLASS][32];
  char mass_vs_pt_arr_n[N_CLASS][32];
  char mass_dist_arr_n[N_CLASS][32];
  char z_dist_arr_n[N_CLASS][32];
  char mass_dist_for_enbin_arr_n[10][32];
  char mass_dist_for_ptbin_arr_n[10][32];

  char pt_rdist_arr_n[N_CLASS][MAXRUNS][32];
  char en_rdist_arr_n[N_CLASS][MAXRUNS][32];
  char eta_rdist_arr_n[N_CLASS][MAXRUNS][32];
  char phi_rdist_arr_n[N_CLASS][MAXRUNS][32];
  char z_rdist_arr_n[N_CLASS][MAXRUNS][32];
  char mass_rdist_arr_n[N_CLASS][MAXRUNS][32];


  for(Int_t c=0; c<N_CLASS; c++)
  {
    sprintf(pt_vs_eta_arr_n[c],"%s_pt_vs_eta_arr",ev->Name(c));
    sprintf(en_vs_eta_arr_n[c],"%s_en_vs_eta_arr",ev->Name(c));
    sprintf(pt_vs_phi_arr_n[c],"%s_pt_vs_phi_arr",ev->Name(c));
    sprintf(en_vs_phi_arr_n[c],"%s_en_vs_phi_arr",ev->Name(c));
    sprintf(eta_vs_phi_arr_n[c],"%s_eta_vs_phi_arr",ev->Name(c));
    sprintf(pt_vs_en_arr_n[c],"%s_pt_vs_en_arr",ev->Name(c));
    sprintf(y_vs_x_arr_n[c],"%s_y_vs_x_arr",ev->Name(c));
    sprintf(z_vs_eta_arr_n[c],"%s_z_vs_eta_arr",ev->Name(c));
    sprintf(z_vs_phi_arr_n[c],"%s_z_vs_phi_arr",ev->Name(c));
    sprintf(mass_vs_en_arr_n[c],"%s_mass_vs_en_arr",ev->Name(c));
    sprintf(mass_vs_pt_arr_n[c],"%s_mass_vs_pt_arr",ev->Name(c));
    sprintf(mass_dist_arr_n[c],"%s_mass_dist_arr",ev->Name(c));
    sprintf(z_dist_arr_n[c],"%s_z_dist_arr",ev->Name(c));

    pt_vs_eta_arr[c] = new TObjArray();
    en_vs_eta_arr[c] = new TObjArray();
    pt_vs_phi_arr[c] = new TObjArray();
    en_vs_phi_arr[c] = new TObjArray();
    eta_vs_phi_arr[c] = new TObjArray();
    pt_vs_en_arr[c] = new TObjArray();
    y_vs_x_arr[c] = new TObjArray();
    z_vs_eta_arr[c] = new TObjArray();
    z_vs_phi_arr[c] = new TObjArray();
    mass_vs_en_arr[c] = new TObjArray();
    mass_vs_pt_arr[c] = new TObjArray();
    mass_dist_arr[c] = new TObjArray();
    z_dist_arr[c] = new TObjArray();

    for(Int_t ru=0; ru<MAXRUNS; ru++)
    {
      sprintf(pt_rdist_arr_n[c][ru],"%s_pt_rdist_arr_%d",ev->Name(c),ru);
      sprintf(en_rdist_arr_n[c][ru],"%s_en_rdist_arr_%d",ev->Name(c),ru);
      sprintf(eta_rdist_arr_n[c][ru],"%s_eta_rdist_arr_%d",ev->Name(c),ru);
      sprintf(phi_rdist_arr_n[c][ru],"%s_phi_rdist_arr_%d",ev->Name(c),ru);
      sprintf(z_rdist_arr_n[c][ru],"%s_z_rdist_arr_%d",ev->Name(c),ru);
      sprintf(mass_rdist_arr_n[c][ru],"%s_mass_rdist_arr_%d",ev->Name(c),ru);
      pt_rdist_arr[c][ru] = new TObjArray();
      en_rdist_arr[c][ru] = new TObjArray();
      eta_rdist_arr[c][ru] = new TObjArray();
      phi_rdist_arr[c][ru] = new TObjArray();
      z_rdist_arr[c][ru] = new TObjArray();
      mass_rdist_arr[c][ru] = new TObjArray();
    };
  };
  for(Int_t k=0; k<10; k++)
  {
    sprintf(mass_dist_for_enbin_arr_n[k],"mass_dist_for_enbin_%d_arr",k);
    sprintf(mass_dist_for_ptbin_arr_n[k],"mass_dist_for_ptbin_%d_arr",k);
    mass_dist_for_enbin_arr[k] = new TObjArray();
    mass_dist_for_ptbin_arr[k] = new TObjArray();
  };

  for(Int_t t=0; t<N_TRIG; t++)
  {
    for(Int_t c=0; c<N_CLASS; c++)
    {
      pt_vs_eta_arr[c]->AddLast(pt_vs_eta[c][t]);
      en_vs_eta_arr[c]->AddLast(en_vs_eta[c][t]);
      pt_vs_phi_arr[c]->AddLast(pt_vs_phi[c][t]);
      en_vs_phi_arr[c]->AddLast(en_vs_phi[c][t]);
      eta_vs_phi_arr[c]->AddLast(eta_vs_phi[c][t]);
      pt_vs_en_arr[c]->AddLast(pt_vs_en[c][t]);
      y_vs_x_arr[c]->AddLast(y_vs_x[c][t]);
      z_vs_eta_arr[c]->AddLast(z_vs_eta[c][t]);
      z_vs_phi_arr[c]->AddLast(z_vs_phi[c][t]);
      mass_vs_en_arr[c]->AddLast(mass_vs_en[c][t]);
      mass_vs_pt_arr[c]->AddLast(mass_vs_pt[c][t]);
      mass_dist_arr[c]->AddLast(mass_dist[c][t]);
      z_dist_arr[c]->AddLast(z_dist[c][t]);
      for(Int_t ru=0; ru<=runcount; ru++)
      {
        pt_rdist_arr[c][ru]->AddLast(pt_rdist[c][t][ru]);
        en_rdist_arr[c][ru]->AddLast(en_rdist[c][t][ru]);
        eta_rdist_arr[c][ru]->AddLast(eta_rdist[c][t][ru]);
        phi_rdist_arr[c][ru]->AddLast(phi_rdist[c][t][ru]);
        z_rdist_arr[c][ru]->AddLast(z_rdist[c][t][ru]);
        mass_rdist_arr[c][ru]->AddLast(mass_rdist[c][t][ru]);
      };
    };
    for(Int_t k=0; k<10; k++)
    {
      mass_dist_for_enbin_arr[k]->AddLast(mass_dist_for_enbin[t][k]);
      mass_dist_for_ptbin_arr[k]->AddLast(mass_dist_for_ptbin[t][k]);
    };
  };


  // write output
  char outfilename[256];
  sscanf(infile_name,"RedOutput%s",outfilename);
  sprintf(outfilename,"%s/diag%s",RD->env->diagset_dir,outfilename);
  TFile * outfile = new TFile(outfilename,"RECREATE");
  outfile->cd();

  trig_dist->Write();

  if(enableOverlap)
  {
    outfile->mkdir("overlap_matrices");
    outfile->cd("overlap_matrices");
    for(Int_t c=0; c<N_CLASS; c++) trig_fms_mix[c]->Write();
  };
  outfile->cd();
  outfile->mkdir("rdists");
  outfile->cd("rdists");
  for(Int_t c=0; c<N_CLASS; c++)
  {
    for(Int_t ru=0; ru<=runcount; ru++)
    {
      //printf("ru=%d c=%d %p\n",ru,c,(void*)pt_rdist_arr_n[c][ru]);
      //printf("%s\n",pt_rdist_arr_n[c][ru]);
      pt_rdist_arr[c][ru]->Write(pt_rdist_arr_n[c][ru],TObject::kSingleKey);
      en_rdist_arr[c][ru]->Write(en_rdist_arr_n[c][ru],TObject::kSingleKey);
      eta_rdist_arr[c][ru]->Write(eta_rdist_arr_n[c][ru],TObject::kSingleKey);
      phi_rdist_arr[c][ru]->Write(phi_rdist_arr_n[c][ru],TObject::kSingleKey);
      z_rdist_arr[c][ru]->Write(z_rdist_arr_n[c][ru],TObject::kSingleKey);
      mass_rdist_arr[c][ru]->Write(mass_rdist_arr_n[c][ru],TObject::kSingleKey);
    };
  };
  outfile->cd();
  for(Int_t c=0; c<N_CLASS; c++)
  {
    pt_vs_eta_arr[c]->Write(pt_vs_eta_arr_n[c],TObject::kSingleKey);
    en_vs_eta_arr[c]->Write(en_vs_eta_arr_n[c],TObject::kSingleKey);
    pt_vs_phi_arr[c]->Write(pt_vs_phi_arr_n[c],TObject::kSingleKey);
    en_vs_phi_arr[c]->Write(en_vs_phi_arr_n[c],TObject::kSingleKey);
    eta_vs_phi_arr[c]->Write(eta_vs_phi_arr_n[c],TObject::kSingleKey);
    pt_vs_en_arr[c]->Write(pt_vs_en_arr_n[c],TObject::kSingleKey);
    y_vs_x_arr[c]->Write(y_vs_x_arr_n[c],TObject::kSingleKey);
    z_vs_eta_arr[c]->Write(z_vs_eta_arr_n[c],TObject::kSingleKey);
    z_vs_phi_arr[c]->Write(z_vs_phi_arr_n[c],TObject::kSingleKey);
    mass_vs_en_arr[c]->Write(mass_vs_en_arr_n[c],TObject::kSingleKey);
    mass_vs_pt_arr[c]->Write(mass_vs_pt_arr_n[c],TObject::kSingleKey);
    mass_dist_arr[c]->Write(mass_dist_arr_n[c],TObject::kSingleKey);
    z_dist_arr[c]->Write(z_dist_arr_n[c],TObject::kSingleKey);
  };
  for(Int_t k=0; k<10; k++)
  {
    mass_dist_for_enbin_arr[k]->Write(mass_dist_for_enbin_arr_n[k],TObject::kSingleKey);
  };
  for(Int_t k=0; k<10; k++)
  {
    mass_dist_for_ptbin_arr[k]->Write(mass_dist_for_ptbin_arr_n[k],TObject::kSingleKey);
  };
};
