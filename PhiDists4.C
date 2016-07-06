// -- create phi distributions with various kinematic cuts (eta,pt,E) for each spinbit
//    (kinematic cuts are set as Double_t's below)
// -- phi distributions for "+ -" and "- +" are weighted by relative luminosity
//    A_LL = 1/(P_b*P_y)*(N++ + N-- - RN+- -RN-+)/(N++ + N-- + RN+- + RN-+)
// -- phi distributions are written to phiset/ directory with similar name; 
//    they are named: phi_s[spinbit]_g[eta bin]_p[pt bin]_e[en bin]

void PhiDists4(const char * filename="RedOutputset080Ba.root",Bool_t debug=false)
{
  if(debug) system("rm debug/debug.dat");

  // load libraries
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * LT = new LevelTwo(RD->env);
  EventClass * ev = new EventClass(RD->env);



  // determine index of TRIGGER_TYPE
  Int_t chosen_trig = LT->Index(TString(RD->env->TriggerType));

  // get bins from environment
  Int_t phi_bins0 = RD->env->PhiBins; const Int_t phi_bins = phi_bins0;
  Int_t eta_bins0 = RD->env->EtaBins; const Int_t eta_bins = eta_bins0;
  Int_t en_bins0 = RD->env->EnBins; const Int_t en_bins = en_bins0;
  Int_t pt_bins0 = RD->env->PtBins; const Int_t pt_bins = pt_bins0;

  // get number of event classes
  Int_t N_CLASS_tmp = ev->N;
  const Int_t N_CLASS = N_CLASS_tmp;


  // read redset tree and set output file name
  Int_t runnum,bx,blue,yell,pattern,ClIndex,Nclust,Evnum;
  Float_t M12,N12,E12,Z,Phi,Eta,Pt,b_pol,y_pol;
  Bool_t kicked,isConsistent;
  char setname[32];
  sscanf(filename,"RedOutputset%s",setname);
  char filename2[256];
  sprintf(filename2,"%s/%s",RD->env->redset_dir,filename);
  TFile * infile = new TFile(filename2,"READ");
  char outname[256];
  sprintf(outname,"%s/phi%s",RD->env->phiset_dir,setname);
  TTree * tree = (TTree*) infile->Get("str");
  tree->SetBranchAddress("runnum",&runnum); // LT->runnum = runnum called directly after tree->GetEntry()
  tree->SetBranchAddress("Bunchid7bit",&bx);
  tree->SetBranchAddress("M12",&M12);
  tree->SetBranchAddress("N12",&N12);
  tree->SetBranchAddress("E12",&E12);
  tree->SetBranchAddress("Z",&Z);
  tree->SetBranchAddress("Phi",&Phi);
  tree->SetBranchAddress("Eta",&Eta);
  tree->SetBranchAddress("Pt",&Pt);
  tree->SetBranchAddress("ClIndex",&ClIndex);
  tree->SetBranchAddress("Nclust",&Nclust);
  tree->SetBranchAddress("Evnum",&Evnum);
  tree->SetBranchAddress("L2sum",LT->L2sum);
  tree->SetBranchAddress("lastdsm",LT->tcu->lastdsm);


  // define spinbit strings
  char spinbit_t[4][4];
  char spinbit_n[4][4];
  strcpy(spinbit_t[0],"--"); strcpy(spinbit_n[0],"nn");
  strcpy(spinbit_t[1],"-+"); strcpy(spinbit_n[1],"np");
  strcpy(spinbit_t[2],"+-"); strcpy(spinbit_n[2],"pn");
  strcpy(spinbit_t[3],"++"); strcpy(spinbit_n[3],"pp");


  // get number of runs in the redset file and build array of run numbers
  Int_t runnum_arr[30]; // (n.b. maximum size arbitrarily defined)
  Int_t runnum_tmp=0;
  Int_t NRUNS_tmp = 0;
  if(tree->GetEntries() == 0)
  {
    fprintf(stderr,"ERROR: no str entries for %s\n --> phi file not produced!\n",filename);
    return;
  };
  for(Int_t i=0; i<tree->GetEntries(); i++)
  {
    tree->GetEntry(i);
    if(runnum != runnum_tmp)
    {
      runnum_tmp = runnum;
      ev->runnum=runnum;
      if(!(ev->ExcludedRun()) || !(RD->RellumConsistent(runnum)))
      {
        if(NRUNS_tmp<30) runnum_arr[NRUNS_tmp] = runnum;
        else 
        {
          fprintf(stderr,"ERROR: more than 30 runs in root file; increase arbitrarily defined max\n");
          return;
        };
        NRUNS_tmp++;
      };
    };
  };
  const Int_t NRUNS = NRUNS_tmp;
  for(Int_t r=0; r<NRUNS; r++) printf("%d\n",runnum_arr[r]);
  printf("NRUNS=%d\n",NRUNS);
  
  // catch files with no good runs (not excluded via exclusion_list and isConsistent==true)
  if(NRUNS<1) {
    fprintf(stderr,"[++++++++++++++] No good runs in this file\n");
    return;
  };



  // define kinematic distributions ("wdist" = weighting distribution)
  // -- these are distributions in energy and phi, but with finer binning; they are used 
  //    to weight the horizontal position of points on the asymmetry plots
  // -- one pt wdist for each en bin (and eta bin) 
  // -- one en wdist for each pt bin (and eta bin)
  // -- one invariant mass (mm) wdist for each eta,pt,en bin
  const Int_t NWBINS = 100;
  TH1D * pt_wdist[N_CLASS][eta_bins][en_bins][NRUNS];
  TH1D * en_wdist[N_CLASS][eta_bins][pt_bins][NRUNS];
  TH1D * mm_wdist[N_CLASS][eta_bins][pt_bins][en_bins][NRUNS];
  char pt_wdist_n[N_CLASS][eta_bins][en_bins][NRUNS][64];
  char en_wdist_n[N_CLASS][eta_bins][pt_bins][NRUNS][64];
  char mm_wdist_n[N_CLASS][eta_bins][pt_bins][en_bins][NRUNS][64];
  for(Int_t c=0; c<N_CLASS; c++)
  {
    for(Int_t r=0; r<NRUNS; r++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          sprintf(pt_wdist_n[c][g][e][r],"pt_wdist_%s_g%d_e%d_r%d",ev->Name(c),g,e,runnum_arr[r]);
          pt_wdist[c][g][e][r] = new TH1D(pt_wdist_n[c][g][e][r],pt_wdist_n[c][g][e][r],NWBINS,
            RD->env->PtLow,RD->env->PtHigh);
        };
        for(Int_t p=0; p<pt_bins; p++)
        {
          sprintf(en_wdist_n[c][g][p][r],"en_wdist_%s_g%d_p%d_r%d",ev->Name(c),g,p,runnum_arr[r]);
          en_wdist[c][g][p][r] = new TH1D(en_wdist_n[c][g][p][r],en_wdist_n[c][g][p][r],NWBINS,
            RD->env->EnLow,RD->env->EnHigh);
        };
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            sprintf(mm_wdist_n[c][g][p][e][r],"mm_wdist_%s_g%d_p%d_e%d_r%d",ev->Name(c),g,p,e,runnum_arr[r]);
            mm_wdist[c][g][p][e][r] = new TH1D(mm_wdist_n[c][g][p][e][r],mm_wdist_n[c][g][p][e][r],NWBINS,0,1);
          };
        };
      };
    };
  };


  // define phi distributions  for each spinbit and kinematic bin
  // ---> ROOT is limited in the sense that I cannot make a pointer to
  //      a dim>=6 array... I have thus folded the spinbit and event class
  //      index into one array dimension, 10*(event class) + spinbit
  //      ... if you have a better idea, good. implment it. 
  TH1D * phi_dist[10*N_CLASS+4][eta_bins][pt_bins][en_bins][NRUNS]; // [10*event_class+spinbit] [..] ..
  TString phi_dist_n[10*N_CLASS+4][eta_bins][pt_bins][en_bins][NRUNS];
  TString phi_dist_t[10*N_CLASS+4][eta_bins][pt_bins][en_bins][NRUNS];
  char tmp_str[256];

  for(Int_t r=0; r<NRUNS; r++)
  {
    for(Int_t c=0; c<N_CLASS; c++)
    {
      for(Int_t s=0; s<4; s++)
      {
        for(Int_t g=0; g<eta_bins; g++)
        {
          for(Int_t p=0; p<pt_bins; p++)
          {
            for(Int_t e=0; e<en_bins; e++)
            {
              sprintf(tmp_str,"phi_%s_s%d_g%d_p%d_e%d_r%d",ev->Name(c),
                s,g,p,e,runnum_arr[r]);
              phi_dist_n[10*c+s][g][p][e][r] = TString(tmp_str);

              sprintf(tmp_str,
  "%s #phi distribution :: spin=(%s) #eta#in[%.2f,%.2f) p_{T}#in[%.2f,%.2f) E#in[%.2f,%.2f) :: r%d",
               ev->Title(c),spinbit_t[s],RD->env->EtaDiv(g),RD->env->EtaDiv(g+1),
               RD->env->PtDiv(p),RD->env->PtDiv(p+1),RD->env->EnDiv(e),RD->env->EnDiv(e+1),runnum_arr[r]);
              phi_dist_t[10*c+s][g][p][e][r] = TString(tmp_str);

              phi_dist[10*c+s][g][p][e][r] = new TH1D(phi_dist_n[10*c+s][g][p][e][r].Data(),
                phi_dist_t[10*c+s][g][p][e][r].Data(),phi_bins,RD->env->PhiLow,RD->env->PhiHigh);

              phi_dist[10*c+s][g][p][e][r]->Sumw2();
            };
          }; 
        };
      };
    };
  };



  // fill phi distributions and wdists 
  Int_t ss,gg,pp,ee,rr;
  Int_t Evnum_tmp=0;
  Int_t xx;
  Bool_t EVPinRange;
  rr=-1; runnum_tmp=0;
  printf("fill phi dists...\n");

  //for(Int_t x=0; x<(tree->GetEntries())/10; x++) {
  for(Int_t x=0; x<tree->GetEntries(); x++) {
    if((x%10000)==0) printf("%.2f%%\n",100*((Float_t)x)/((Float_t)tree->GetEntries()));
    ss=gg=pp=ee=-1; // reset 

    tree->GetEntry(x);
    LT->runnum = runnum;

    // check if run is on exclusion_list
    ev->runnum = runnum;
    if(ev->ExcludedRun()) continue;

    // run number --> array index
    if(runnum != runnum_tmp)
    {
      rr=-1;
      for(Int_t r=0; r<NRUNS; r++) { if(runnum_arr[r] == runnum) rr=r; }; 
      runnum_tmp = runnum;
      isConsistent = RD->RellumConsistent(runnum);
      pattern = RD->Pattern(runnum);
      b_pol = RD->BluePol(runnum);
      y_pol = RD->YellPol(runnum);
    };


    // check for rellum consistency and valid polarization and FMS trigger
    if(isConsistent==1 && b_pol>0 && y_pol>0 && LT->Fired(chosen_trig))
    {
      blue = RD->BlueSpin(runnum,bx);
      yell = RD->YellSpin(runnum,bx);
      kicked = RD->Kicked(runnum,bx);

      // check spin and kicked status
      if(abs(blue)==1 && abs(yell)==1 && kicked==0)
      {
        // spin --> array index
        ss = 2 * (blue>0) + (yell>0);

        // kinematic bins --> array indices
        for(Int_t g=0; g<eta_bins; g++) { if(Eta>=RD->env->EtaDiv(g) && Eta<RD->env->EtaDiv(g+1)) { gg=g; break; }; };
        for(Int_t p=0; p<pt_bins;  p++) { if(Pt>=RD->env->PtDiv(p)   && Pt<RD->env->PtDiv(p+1)  ) { pp=p; break; }; };
        for(Int_t e=0; e<en_bins;  e++) { if(E12>=RD->env->EnDiv(e)  && E12<RD->env->EnDiv(e+1) ) { ee=e; break; }; };


        // check for valid array indices (filters out events outside kinematic boundaries)
        //printf("%d %d %d %d %d\n",ss,gg,pp,ee,rr);
        /*if(ss>=0 && gg>=0 && pp>=0 && ee>=0 && rr>=0 && fabs(Phi)>3.1415/2.0)*/ /* north cells only */
        /*if(ss>=0 && gg>=0 && pp>=0 && ee>=0 && rr>=0 && fabs(Phi)<3.1415/2.0)*/ /* south cells only */
        if(ss>=0 && gg>=0 && pp>=0 && ee>=0 && rr>=0)
        {
          // set kinematics variables for event, tcu bits, rp bits
          ev->SetKinematics(runnum,E12,Pt,Eta,Phi,M12,Z,N12,ClIndex);


          // determine event class(es)
          for(Int_t c=0; c<N_CLASS; c++)
          {
            // for di-pions //
            if(c==ev->Idx("dpi") && x+1!=tree->GetEntries()) {
              if(ClIndex!=0) continue; // only use ClIndex==0 pions
              //printf("----------\n");
              Evnum_tmp = Evnum;
              xx = 0;
              while(Evnum==Evnum_tmp && x+xx<tree->GetEntries()) {
                ev->AppendKinematics(E12,Pt,Eta,Phi,M12,Z,N12,xx);
                //printf("x=%d xx=%d Evnum=%d ClIndex=%d Nclust=%d E12=%f Pt=%f\n",
                  //x,xx,Evnum,ClIndex,Nclust,E12,Pt);
                xx++;
                tree->GetEntry(x+xx);
              };
              ev->Nclust = xx;
              tree->GetEntry(x);
            };
            ///////

            if(ev->Valid(c,chosen_trig))
            {
              if(c==ev->Idx("dpi") && x+1!=tree->GetEntries()) {
                phi_dist[10*c+ss][gg][pp][ee][rr]->Fill(ev->delta_phi);
              }
              else {
                phi_dist[10*c+ss][gg][pp][ee][rr]->Fill(Phi);
              }
              pt_wdist[c][gg][ee][rr]->Fill(Pt);
              en_wdist[c][gg][pp][rr]->Fill(E12);
              mm_wdist[c][gg][pp][ee][rr]->Fill(M12);
            };
          };
        };
      };
    };
  };

        

  // make object arrays
  TFile * outfile = new TFile(outname,"RECREATE");
  TObjArray * phi_dist_arr[4][eta_bins][pt_bins][en_bins][N_CLASS];
  TObjArray * pt_wdist_arr[N_CLASS][eta_bins][en_bins];
  TObjArray * en_wdist_arr[N_CLASS][eta_bins][pt_bins];
  TObjArray * mm_wdist_arr[N_CLASS][eta_bins][pt_bins][en_bins];
  char phi_dist_arr_name[4][eta_bins][pt_bins][en_bins][N_CLASS][64];
  char pt_wdist_arr_name[N_CLASS][eta_bins][en_bins][64];
  char en_wdist_arr_name[N_CLASS][eta_bins][pt_bins][64];
  char mm_wdist_arr_name[N_CLASS][eta_bins][pt_bins][en_bins][64];
  for(Int_t c=0; c<N_CLASS; c++)
  {
    for(Int_t e=0; e<en_bins; e++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t s=0; s<4; s++)
          {
            phi_dist_arr[s][g][p][e][c] = new TObjArray();

            sprintf(phi_dist_arr_name[s][g][p][e][c],"phi_dist_%s_s%d_g%d_p%d_e%d",ev->Name(c),s,g,p,e);
            
            for(Int_t r=0; r<NRUNS; r++)
            {
              phi_dist_arr[s][g][p][e][c]->AddLast(phi_dist[10*c+s][g][p][e][r]);
            };
          };
        };
      };
    };
  };
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t e=0; e<en_bins; e++)
    {
      for(Int_t c=0; c<N_CLASS; c++) 
      {
        pt_wdist_arr[c][g][e] = new TObjArray();
        sprintf(pt_wdist_arr_name[c][g][e],"pt_wdist_%s_g%d_e%d",ev->Name(c),g,e);
        for(Int_t r=0; r<NRUNS; r++) pt_wdist_arr[c][g][e]->AddLast(pt_wdist[c][g][e][r]);
      };
    };
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t c=0; c<N_CLASS; c++) 
      {
        en_wdist_arr[c][g][p] = new TObjArray();
        sprintf(en_wdist_arr_name[c][g][p],"en_wdist_%s_g%d_p%d",ev->Name(c),g,p);
        for(Int_t r=0; r<NRUNS; r++) en_wdist_arr[c][g][p]->AddLast(en_wdist[c][g][p][r]);
      };
    };
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        for(Int_t c=0; c<N_CLASS; c++) 
        {
          mm_wdist_arr[c][g][p][e] = new TObjArray();
          sprintf(mm_wdist_arr_name[c][g][p][e],"mm_wdist_%s_g%d_p%d_e%d",ev->Name(c),g,p,e);
          for(Int_t r=0; r<NRUNS; r++) mm_wdist_arr[c][g][p][e]->AddLast(mm_wdist[c][g][p][e][r]);
        };
      };
    };
  };

  
  // write object arrays
  for(Int_t c=0; c<N_CLASS; c++)
  {
    outfile->mkdir(ev->Name(c));
    outfile->cd(ev->Name(c));
    for(Int_t e=0; e<en_bins; e++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t s=0; s<4; s++)
          {
            phi_dist_arr[s][g][p][e][c]->Write(phi_dist_arr_name[s][g][p][e][c],TObject::kSingleKey);
          };
        };
      };
    };
    outfile->cd();
  };

  // write wdists
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t c=0; c<N_CLASS; c++)
    {
      for(Int_t e=0; e<en_bins; e++) pt_wdist_arr[c][g][e]->Write(pt_wdist_arr_name[c][g][e],TObject::kSingleKey);
      for(Int_t p=0; p<pt_bins; p++) en_wdist_arr[c][g][p]->Write(en_wdist_arr_name[c][g][p],TObject::kSingleKey);
      for(Int_t p=0; p<pt_bins; p++) for(Int_t e=0; e<en_bins; e++)
        mm_wdist_arr[c][g][p][e]->Write(mm_wdist_arr_name[c][g][p][e],TObject::kSingleKey);
    };
  };
};
