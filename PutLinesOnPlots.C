// subroutine of hadd_MassInPTE; this just puts the bin lines on the kinematic
// plots in massset*/all.root (or any other file in that directory that you want it to)

void PutLinesOnPlots(TString infile_n = "all.root") {
  // load libs
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * LT = new LevelTwo(RD->env);
  EventClass * ev = new EventClass(RD->env,false);

  // open file
  TString dir(RD->env->massset_dir);
  infile_n = dir+"/"+infile_n;
  TFile * infile = new TFile(infile_n.Data(),"UPDATE");


  // trigger list
  Int_t N_TRIG_tmp = LT->N;
  const Int_t N_TRIG = N_TRIG_tmp;

  // iterators
  Int_t t,g,e,p;


  // get kinematic ranges using the current binning set by massenv.sh
  Float_t E12_min,Pt_min,Eta_min,Phi_min;
  Float_t E12_max,Pt_max,Eta_max,Phi_max;
  Float_t E12_bins,Pt_bins,Eta_bins,Phi_bins;
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
  Double_t low_mass = 0;
  Double_t high_mass = 1;


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


  // define 2d hists
  TH2D * mass_vs_en[N_G][N_TRIG];
  TH2D * mass_vs_pt[N_G][N_TRIG];
  TH2D * pt_vs_en[N_G][N_TRIG];
  TH2D * mass_vs_eta[N_TRIG];

  // define canvases
  TCanvas * mass_vs_en_canv[N_G][N_TRIG];
  TCanvas * mass_vs_pt_canv[N_G][N_TRIG];
  TCanvas * pt_vs_en_canv[N_G][N_TRIG];
  TCanvas * mass_vs_eta_canv[N_TRIG];

  TString  mass_vs_en_n[N_G][N_TRIG];
  TString  mass_vs_pt_n[N_G][N_TRIG];
  TString  pt_vs_en_n[N_G][N_TRIG];
  TString  mass_vs_eta_n[N_TRIG];

  // obtain 2d histograms
  for(t=0; t<N_TRIG; t++) {
    for(g=0; g<N_G; g++) {
      mass_vs_en_n[g][t] = Form("mass_vs_en_g%d_%s_canv",g,(LT->Name(t)).Data());
      mass_vs_pt_n[g][t] = Form("mass_vs_pt_g%d_%s_canv",g,(LT->Name(t)).Data());
      pt_vs_en_n[g][t] = Form("pt_vs_en_g%d_%s_canv",g,(LT->Name(t)).Data());

      mass_vs_en_canv[g][t] = new TCanvas(mass_vs_en_n[g][t].Data(),mass_vs_en_n[g][t].Data(),800,800);
      mass_vs_pt_canv[g][t] = new TCanvas(mass_vs_pt_n[g][t].Data(),mass_vs_pt_n[g][t].Data(),800,800);
      pt_vs_en_canv[g][t] = new TCanvas(pt_vs_en_n[g][t].Data(),pt_vs_en_n[g][t].Data(),800,800);
    };

    mass_vs_eta_n[t] = "mass_vs_eta_"+LT->Name(t)+"_canv";
    mass_vs_eta_canv[t] = new TCanvas(mass_vs_eta_n[t].Data(),mass_vs_eta_n[t].Data(),800,800);
  };


  // read input arrays
  TObjArray * mass_vs_en_arr[N_G];
  TObjArray * mass_vs_pt_arr[N_G];
  TObjArray * pt_vs_en_arr[N_G];
  TObjArray * mass_vs_eta_arr;

  TString mass_vs_en_arr_n[N_G];
  TString mass_vs_pt_arr_n[N_G];
  TString pt_vs_en_arr_n[N_G];
  TString mass_vs_eta_arr_n;

  mass_vs_eta_arr_n = "mass_vs_eta_arr";
  mass_vs_eta_arr = (TObjArray*) infile->Get(mass_vs_eta_arr_n.Data());
  for(g=0; g<N_G; g++) {
    mass_vs_en_arr_n[g] = Form("mass_vs_en_arr_g%d",g);
    mass_vs_pt_arr_n[g] = Form("mass_vs_pt_arr_g%d",g);
    pt_vs_en_arr_n[g] = Form("pt_vs_en_arr_g%d",g);

    mass_vs_en_arr[g] = (TObjArray*) infile->Get(mass_vs_en_arr_n[g].Data());
    mass_vs_pt_arr[g] = (TObjArray*) infile->Get(mass_vs_pt_arr_n[g].Data());
    pt_vs_en_arr[g] = (TObjArray*) infile->Get(pt_vs_en_arr_n[g].Data());


    for(t=0; t<N_TRIG; t++) {
      mass_vs_en[g][t] = (TH2D*) (mass_vs_en_arr[g]->At(t));
      mass_vs_pt[g][t] = (TH2D*) mass_vs_pt_arr[g]->At(t);
      pt_vs_en[g][t] = (TH2D*) pt_vs_en_arr[g]->At(t);
    };
  };
  for(t=0; t<N_TRIG; t++) mass_vs_eta[t] = (TH2D*) mass_vs_eta_arr->At(t);


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


  // define output arrays
  TObjArray * mass_vs_en_canvarr[N_G];
  TObjArray * mass_vs_pt_canvarr[N_G];
  TObjArray * pt_vs_en_canvarr[N_G];
  TObjArray * mass_vs_eta_canvarr;

  TString mass_vs_en_canvarr_n[N_G];
  TString mass_vs_pt_canvarr_n[N_G];
  TString pt_vs_en_canvarr_n[N_G];
  TString mass_vs_eta_canvarr_n;

  mass_vs_eta_canvarr_n = "mass_vs_eta_canvarr";
  mass_vs_eta_canvarr = new TObjArray();
  for(g=0; g<N_G; g++) {
    mass_vs_en_canvarr_n[g] = Form("mass_vs_en_canvarr_g%d",g);
    mass_vs_pt_canvarr_n[g] = Form("mass_vs_pt_canvarr_g%d",g);
    pt_vs_en_canvarr_n[g] = Form("pt_vs_en_canvarr_g%d",g);

    mass_vs_en_canvarr[g] = new TObjArray();
    mass_vs_pt_canvarr[g] = new TObjArray();
    pt_vs_en_canvarr[g] = new TObjArray();

    for(t=0; t<N_TRIG; t++) {
      mass_vs_en_canvarr[g]->AddLast(mass_vs_en_canv[g][t]);
      mass_vs_pt_canvarr[g]->AddLast(mass_vs_pt_canv[g][t]);
      pt_vs_en_canvarr[g]->AddLast(pt_vs_en_canv[g][t]);
    };
  };
  for(t=0; t<N_TRIG; t++) mass_vs_eta_canvarr->AddLast(mass_vs_eta_canv[t]);



  // write output
  for(g=0; g<N_G; g++) mass_vs_en_canvarr[g]->Write(mass_vs_en_canvarr_n[g],TObject::kSingleKey);
  for(g=0; g<N_G; g++) mass_vs_pt_canvarr[g]->Write(mass_vs_pt_canvarr_n[g],TObject::kSingleKey);
  for(g=0; g<N_G; g++) pt_vs_en_canvarr[g]->Write(pt_vs_en_canvarr_n[g],TObject::kSingleKey);
  mass_vs_eta_canvarr->Write(mass_vs_eta_canvarr_n,TObject::kSingleKey);
};
