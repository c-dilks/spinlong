// plots mass distributions, based on environment (WHICH SHOULD BE massenv.sh!)
//

void MassInPTE(const char * infile_name = "RedOutputset087Ba.root")
{

  const Int_t NBINS=75; // NUMBER OF BINS (default 400)
  const Int_t MAXRUNS=50; // arbitrary max number of runs in redset file 

  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * LT = new LevelTwo(RD->env);
  EventClass * ev = new EventClass(RD->env,false);


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
  Float_t E12_bins,Pt_bins,Eta_bins,Phi_bins;

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


  // get kinematic ranges using the current binning set by massenv.sh
  Phi_min = RD->env->PhiLow;
  Phi_max = RD->env->PhiHigh;
  Eta_min = RD->env->EtaLow;
  Eta_max = RD->env->EtaHigh;
  E12_min = RD->env->EnLow;
  E12_max = RD->env->EnHigh;
  Pt_min = RD->env->PtLow;
  Pt_max = RD->env->PtHigh;
  Eta_bins = RD->env->EtaBins;
  E12_bins = RD->env->EnBins;
  Pt_bins = RD->env->PtBins;
  const Int_t N_G = Eta_bins;
  const Int_t N_E = E12_bins;
  const Int_t N_P = Pt_bins;

  
  // bin boundaries
  Double_t low_mass = 0;
  Double_t high_mass = 1;

  Int_t c_idx = ev->Idx("pi0");

  // define histograms
  TH1D * mass_dist[N_G][N_E][N_P][N_TRIG];
  TString mass_dist_n, mass_dist_t;
  TString g_str, e_str, p_str;

  // iterators
  Int_t t,g,e,p;
 

  for(g=0; g<N_G; g++) {
    for(e=0; e<N_E; e++) {
      for(p=0; p<N_P; p++) {
        for(t=0; t<N_TRIG; t++) {
          mass_dist_n = Form("mass_dist_g%d_e%d_p%d_t%s",
            g,e,p,(LT->Name(t)).Data());
          mass_dist_t = "M_{#gamma#gamma} distribution :: ";
          g_str = Form("#eta#in[%.2f,%.2f] :: ",
                        RD->env->EtaDiv(g), RD->env->EtaDiv(g+1));
          e_str = Form("#E_{#gamma#gamma}#in[%.2f,%.2f] :: ",
                        RD->env->EnDiv(e), RD->env->EnDiv(e+1));
          p_str = Form("#p_{T}#in[%.2f,%.2f] :: ",
                        RD->env->PtDiv(p), RD->env->PtDiv(p+1));
          mass_dist_t = mass_dist_t + g_str + e_str + p_str + LT->Name(t);
          mass_dist[g][e][p][t] = new TH1D(mass_dist_n,
                                           mass_dist_t,
                                           NBINS,
                                           low_mass,
                                           high_mass
                                          );
        };
      };
    }; 
  };





  // fill histograms
  Int_t runnum_tmp=0;
  Int_t runcount=-1;
  Int_t g_b,e_b,p_b;

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
      runnum_tmp=runnum;
    }

    // rellum / pol cut
    if( kicked==0 && isConsistent==1 && b_pol>0 && y_pol>0)
    {
      // below here we have if statements for each event class; within each
      // if block, we have trigger loops, since it's more efficient to first pick the event and 
      // then loop through all triggers

      ev->SetKinematics(runnum,E12,Pt,Eta,Phi,M12,Z,N12,ClIndex);

      // fill mass plots
      if(ev->ValidWithoutMcut(c_idx,-1))
      {
        for(Int_t t=0; t<N_TRIG; t++)
        {
          if(LT->Fired(t))
          {
            if(ev->ValidWithoutMcut(c_idx,t)) {
              // determine which kinematic bin
              g_b = e_b = p_b = 0;
              while(Eta > RD->env->EtaDiv(g_b)) g_b++;
              while(E12 > RD->env->EnDiv(g_b))  e_b++;
              while(Pt  > RD->env->PtDiv(g_b))  p_b++;
              mass_dist[g][e][p][t]->Fill(M12);
            };
          };
        };
      }; // eo if ValidWithoutMcut(c_idx,-1)
    }; // eo rellum/pol cut
  }; // eo tree loop


  // acabe de escribir aqui----------------------------------------- 


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
  sprintf(outfilename,"%s_tight/diag%s",RD->env->diagset_dir,outfilename);
  TFile * outfile = new TFile(outfilename,"RECREATE");
  outfile->cd();

  trig_dist->Write();

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
