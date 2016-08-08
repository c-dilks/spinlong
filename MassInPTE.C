// plots mass distributions, based on environment (WHICH SHOULD BE massenv.sh!)
//

void MassInPTE(TString infile_name = "RedOutputset087Ba.root")
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
  sprintf(infile_full_name,"%s/%s",RD->env->redset_dir,infile_name.Data());
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

  // define 2d hists
  TH2D * mass_vs_en[N_G][N_TRIG];
  TH2D * mass_vs_pt[N_G][N_TRIG];
  TH2D * pt_vs_en[N_G][N_TRIG];
  TH2D * mass_vs_eta[N_TRIG];

  TString  mass_vs_en_n[N_G][N_TRIG];
  TString  mass_vs_pt_n[N_G][N_TRIG];
  TString  pt_vs_en_n[N_G][N_TRIG];
  TString  mass_vs_eta_n[N_TRIG];

  TString  mass_vs_en_t[N_G][N_TRIG];
  TString  mass_vs_pt_t[N_G][N_TRIG];
  TString  pt_vs_en_t[N_G][N_TRIG];
  TString  mass_vs_eta_t[N_TRIG];

  TCanvas * mass_vs_en_canv[N_G][N_TRIG];
  TCanvas * mass_vs_pt_canv[N_G][N_TRIG];
  TCanvas * pt_vs_en_canv[N_G][N_TRIG];
  TCanvas * mass_vs_eta_canv[N_TRIG];

  // iterators
  Int_t t,g,e,p;

 
  // initialise distributions
  for(g=0; g<N_G; g++) {
    for(e=0; e<N_E; e++) {
      for(p=0; p<N_P; p++) {
        for(t=0; t<N_TRIG; t++) {
          mass_dist_n = Form("mass_dist_g%d_e%d_p%d_t%s",
            g,e,p,(LT->Name(t)).Data());
          mass_dist_t = "M_{#gamma#gamma} distribution :: ";
          g_str = Form("#eta#in[%.2f,%.2f] :: ",
                        RD->env->EtaDiv(g), RD->env->EtaDiv(g+1));
          e_str = Form("E_{#gamma#gamma}#in[%.2f,%.2f] :: ",
                        RD->env->EnDiv(e), RD->env->EnDiv(e+1));
          p_str = Form("p_{T}#in[%.2f,%.2f] :: ",
                        RD->env->PtDiv(p), RD->env->PtDiv(p+1));
          mass_dist_t = mass_dist_t + g_str + e_str + p_str + LT->Name(t) + " triggers";
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


  // initialise 2d histograms
  for(t=0; t<N_TRIG; t++) {
    for(g=0; g<N_G; g++) {
      mass_vs_en_n[g][t] = Form("mass_vs_en_g%d_%s",g,(LT->Name(t)).Data());
      mass_vs_pt_n[g][t] = Form("mass_vs_pt_g%d_%s",g,(LT->Name(t)).Data());
      pt_vs_en_n[g][t] = Form("pt_vs_en_g%d_%s",g,(LT->Name(t)).Data());

      mass_vs_en_t[g][t] = Form("M_{#gamma#gamma} vs E_{#gamma#gamma} :: #eta#in[%.2f,%.2f] :: %s triggers;E_{#gamma#gamma};M_{#gamma#gamma}",
                                RD->env->EtaDiv(g),RD->env->EtaDiv(g+1),(LT->Name(t)).Data());
      mass_vs_pt_t[g][t] = Form("M_{#gamma#gamma} vs p_{T} :: #eta#in[%.2f,%.2f] :: %s triggers;p_{T};M_{#gamma#gamma}",
                                RD->env->EtaDiv(g),RD->env->EtaDiv(g+1),(LT->Name(t)).Data());
      pt_vs_en_t[g][t] = Form("p_{T} vs E_{#gamma#gamma} :: #eta#in[%.2f,%.2f] :: %s triggers;E_{#gamma#gamma};p_{T}",
                                RD->env->EtaDiv(g),RD->env->EtaDiv(g+1),(LT->Name(t)).Data());

      mass_vs_en[g][t] = new TH2D(mass_vs_en_n[g][t].Data(),mass_vs_en_t[g][t].Data(),NBINS,E12_min,E12_max,NBINS,low_mass,high_mass);
      mass_vs_pt[g][t] = new TH2D(mass_vs_pt_n[g][t].Data(),mass_vs_pt_t[g][t].Data(),NBINS,Pt_min,Pt_max,NBINS,low_mass,high_mass);
      pt_vs_en[g][t] = new TH2D(pt_vs_en_n[g][t].Data(),pt_vs_en_t[g][t].Data(),NBINS,E12_min,E12_max,NBINS,Pt_min,Pt_max);

      mass_vs_en_n[g][t] = mass_vs_en_n[g][t] + "_canv";
      mass_vs_pt_n[g][t] = mass_vs_pt_n[g][t] + "_canv";
      pt_vs_en_n[g][t] = pt_vs_en_n[g][t] + "_canv";

      mass_vs_en_canv[g][t] = new TCanvas(mass_vs_en_n[g][t].Data(),mass_vs_en_n[g][t].Data(),800,800);
      mass_vs_pt_canv[g][t] = new TCanvas(mass_vs_pt_n[g][t].Data(),mass_vs_pt_n[g][t].Data(),800,800);
      pt_vs_en_canv[g][t] = new TCanvas(pt_vs_en_n[g][t].Data(),pt_vs_en_n[g][t].Data(),800,800);
    };

    mass_vs_eta_n[t] = "mass_vs_eta_"+LT->Name(t);
    mass_vs_eta_t[t] = "M_{#gamma#gamma} vs #eta :: "+LT->Name(t)+" triggers;#eta;M_{#gamma#gamma}";
    mass_vs_eta[t] = new TH2D(mass_vs_eta_n[t].Data(),mass_vs_eta_t[t].Data(),NBINS,Eta_min,Eta_max,NBINS,low_mass,high_mass);
    mass_vs_eta_n[t] = mass_vs_eta_n[t] + "_canv";
    mass_vs_eta_canv[t] = new TCanvas(mass_vs_eta_n[t].Data(),mass_vs_eta_n[t].Data(),800,800);
  };


  // bin lines for 2d hists
  TLine * en_line[N_E-1];
  TLine * en2_line[N_E-1]; // (version of en_line, for pt_vs_en)
  TLine * pt_line[N_P-1];
  TLine * pt2_line[N_P-1]; // (version of pt_line, for pt_vs_en)
  TLine * eta_line[N_G-1];

  const Float_t LWIDTH = 2;
  const Int_t LCOLOR = kRed;

  for(e=0; e<N_E-1; e++) {
    en_line[e] = new TLine(RD->env->EnDiv(e+1),low_mass,RD->env->EnDiv(e+1),high_mass);
    en2_line[e] = new TLine(RD->env->EnDiv(e+1),Pt_min,RD->env->EnDiv(e+1),Pt_max);
    en_line[e]->SetLineWidth(LWIDTH);
    en_line[e]->SetLineColor(LCOLOR);
    en2_line[e]->SetLineWidth(LWIDTH);
    en2_line[e]->SetLineColor(LCOLOR);
  };
  for(p=0; p<N_P-1; p++) {
    pt_line[p] = new TLine(RD->env->PtDiv(p+1),low_mass,RD->env->PtDiv(p+1),high_mass);
    pt2_line[p] = new TLine(E12_min,RD->env->PtDiv(p+1),E12_max,RD->env->PtDiv(p+1));
    pt_line[p]->SetLineWidth(LWIDTH);
    pt_line[p]->SetLineColor(LCOLOR);
    pt2_line[p]->SetLineWidth(LWIDTH);
    pt2_line[p]->SetLineColor(LCOLOR);
  };
  for(g=0; g<N_G-1; g++) {
    eta_line[g] = new TLine(RD->env->EtaDiv(g+1),low_mass,RD->env->EtaDiv(g+1),high_mass);
    eta_line[g]->SetLineWidth(LWIDTH);
    eta_line[g]->SetLineColor(LCOLOR);
  };

    


  // fill distributions
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
              g_b = e_b = p_b = -1;
              for(g=0; g<N_G; g++) { if(Eta>=RD->env->EtaDiv(g) && Eta<=RD->env->EtaDiv(g+1)) { g_b=g; break; }; };
              for(e=0; e<N_E; e++) { if(E12>=RD->env->EnDiv(e)  && E12<=RD->env->EnDiv(e+1) ) { e_b=e; break; }; };
              for(p=0; p<N_P; p++) { if(Pt>=RD->env->PtDiv(p)   && Pt<=RD->env->PtDiv(p+1)  ) { p_b=p; break; }; };
              //printf("%d %d %d\n",g_b,e_b,p_b);
              if(g_b>=0 && e_b>=0 && p_b>=0) {
                mass_dist[g_b][e_b][p_b][t]->Fill(M12);
                mass_vs_en[g_b][t]->Fill(E12,M12);
                mass_vs_pt[g_b][t]->Fill(Pt,M12);
                pt_vs_en[g_b][t]->Fill(E12,Pt);
                mass_vs_eta[t]->Fill(Eta,M12);
              };
            };
          };
        };
      }; // eo if ValidWithoutMcut(c_idx,-1)
    }; // eo rellum/pol cut
  }; // eo tree loop


  // draw 2d histograms with bin lines
  for(t=0; t<N_TRIG; t++) {
    for(g=0; g<N_G; g++) {
      mass_vs_en_canv[g][t]->cd();
      mass_vs_en[g][t]->Draw("colz");
      for(int k=0; k<N_E-1; k++) en_line[k]->Draw();

      mass_vs_pt_canv[g][t]->cd();
      mass_vs_pt[g][t]->Draw("colz");
      for(int k=0; k<N_P-1; k++) pt_line[k]->Draw();

      pt_vs_en_canv[g][t]->cd();
      pt_vs_en[g][t]->Draw("colz");
      for(int ke=0; ke<N_E-1; ke++) en2_line[ke]->Draw();
      for(int kp=0; kp<N_P-1; kp++) pt2_line[kp]->Draw();
    };

    mass_vs_eta_canv[t]->cd();
    mass_vs_eta[t]->Draw("colz");
    for(int k=0; k<N_G-1; k++) eta_line[k]->Draw();
  };


  // build obj arrays for output
  TObjArray * mass_dist_arr[N_G][N_E][N_P]; 
  TObjArray * mass_vs_en_arr[N_G];
  TObjArray * mass_vs_pt_arr[N_G];
  TObjArray * pt_vs_en_arr[N_G];
  TObjArray * mass_vs_eta_arr;

  TString mass_dist_arr_n[N_G][N_E][N_P];
  TString mass_vs_en_arr_n[N_G];
  TString mass_vs_pt_arr_n[N_G];
  TString pt_vs_en_arr_n[N_G];
  TString mass_vs_eta_arr_n;

  for(g=0; g<N_G; g++) {
    for(e=0; e<N_E; e++) {
      for(p=0; p<N_P; p++) {
        mass_dist_arr_n[g][e][p] = Form("mass_dist_arr_g%d_e%d_p%d",g,e,p);
        mass_dist_arr[g][e][p] = new TObjArray();
        for(t=0; t<N_TRIG; t++) {
          mass_dist_arr[g][e][p]->AddLast(mass_dist[g][e][p][t]);
        };
      };
    };
  };

  mass_vs_eta_arr_n = "mass_vs_eta_arr";
  mass_vs_eta_arr = new TObjArray();
  for(g=0; g<N_G; g++) {
    mass_vs_en_arr_n[g] = Form("mass_vs_en_arr_g%d",g);
    mass_vs_pt_arr_n[g] = Form("mass_vs_pt_arr_g%d",g);
    pt_vs_en_arr_n[g] = Form("pt_vs_en_arr_g%d",g);

    mass_vs_en_arr[g] = new TObjArray();
    mass_vs_pt_arr[g] = new TObjArray();
    pt_vs_en_arr[g] = new TObjArray();

    for(t=0; t<N_TRIG; t++) {
      mass_vs_en_arr[g]->AddLast(mass_vs_en[g][t]);
      mass_vs_pt_arr[g]->AddLast(mass_vs_pt[g][t]);
      pt_vs_en_arr[g]->AddLast(pt_vs_en[g][t]);
    };
  };
  for(t=0; t<N_TRIG; t++) mass_vs_eta_arr->AddLast(mass_vs_eta[t]);



  // write output
  char mkdircmd[2048];
  sprintf(mkdircmd,".! mkdir -p %s",RD->env->massset_dir);
  gROOT->ProcessLine(mkdircmd);
  TString outfilename = infile_name;
  TRegexp outfilename_re("RedOutputset");
  TString outfilename_rep = Form("%s/mass",RD->env->massset_dir);
  outfilename(outfilename_re) = outfilename_rep;
  TFile * outfile = new TFile(outfilename.Data(),"RECREATE");
  outfile->cd();
  for(g=0; g<N_G; g++) {
    for(e=0; e<N_E; e++) {
      for(p=0; p<N_P; p++) {
        mass_dist_arr[g][e][p]->Write(mass_dist_arr_n[g][e][p].Data(),TObject::kSingleKey);
      };
    };
  };

  for(g=0; g<N_G; g++) mass_vs_en_arr[g]->Write(mass_vs_en_arr_n[g],TObject::kSingleKey);
  for(g=0; g<N_G; g++) mass_vs_pt_arr[g]->Write(mass_vs_pt_arr_n[g],TObject::kSingleKey);
  for(g=0; g<N_G; g++) pt_vs_en_arr[g]->Write(pt_vs_en_arr_n[g],TObject::kSingleKey);
  mass_vs_eta_arr->Write(mass_vs_eta_arr_n,TObject::kSingleKey);
};
