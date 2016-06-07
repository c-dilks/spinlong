// draws plots from diag.root

void DrawDiagnostics(const char * filename="diagset/all.root")
{
  TFile * infile = new TFile(filename,"READ");
  Int_t RESX = 400;
  Int_t RESY = 2000;

  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * T = new LevelTwo(RD->env);
  EventClass * ev = new EventClass(RD->env);

  gStyle->SetOptStat(0);

  // event classes
  Int_t N_CLASS_tmp = ev->N;
  const Int_t N_CLASS = N_CLASS_tmp;

  // trigger names
  Int_t N_TRIG_tmp = T->N;
  const Int_t N_TRIG = N_TRIG_tmp;

  // initialise tcanvases
  TCanvas * trig_canv = new TCanvas("trig_canv","trig_canv",500,500);
  trig_canv->SetGrid(1,1);
  TCanvas * corr_canv[N_CLASS];
  TCanvas * corr_MZ_canv[N_CLASS];
  char corr_canv_n[N_CLASS][16];
  char corr_MZ_canv_n[N_CLASS][16];
  for(Int_t c=0; c<N_CLASS; c++)
  {
    sprintf(corr_canv_n[c],"corr_canv_%s",ev->Name(c));
    sprintf(corr_MZ_canv_n[c],"corr_MZ_canv_%s",ev->Name(c));
    corr_canv[c] = new TCanvas(corr_canv_n[c],corr_canv_n[c],RESX,RESY);
    corr_MZ_canv[c] = new TCanvas(corr_MZ_canv_n[c],corr_MZ_canv_n[c],RESX,RESY);
    corr_canv[c]->Divide(1,6);
    corr_MZ_canv[c]->Divide(1,6);
    for(Int_t cc=1; cc<=6; cc++) 
    {
      corr_canv[c]->GetPad(cc)->SetLogz();
      if(c<=4) corr_MZ_canv[c]->GetPad(cc)->SetLogz();
      corr_canv[c]->GetPad(cc)->SetGrid(1,1);
      corr_MZ_canv[c]->GetPad(cc)->SetGrid(1,1);
    };
  };


  // read object arrays
  TObjArray * pt_vs_eta_arr[N_CLASS];
  TObjArray * en_vs_eta_arr[N_CLASS];
  TObjArray * pt_vs_phi_arr[N_CLASS];
  TObjArray * en_vs_phi_arr[N_CLASS];
  TObjArray * eta_vs_phi_arr[N_CLASS];
  TObjArray * pt_vs_en_arr[N_CLASS];
  TObjArray * z_vs_eta_arr[N_CLASS];
  TObjArray * z_vs_phi_arr[N_CLASS];
  TObjArray * mass_vs_en_arr[N_CLASS];
  TObjArray * mass_vs_pt_arr[N_CLASS];
  TObjArray * mass_dist_arr[N_CLASS];
  TObjArray * z_dist_arr[N_CLASS];
  char pt_vs_eta_arr_n[N_CLASS][32];
  char en_vs_eta_arr_n[N_CLASS][32];
  char pt_vs_phi_arr_n[N_CLASS][32];
  char en_vs_phi_arr_n[N_CLASS][32];
  char eta_vs_phi_arr_n[N_CLASS][32];
  char pt_vs_en_arr_n[N_CLASS][32];
  char z_vs_eta_arr_n[N_CLASS][32];
  char z_vs_phi_arr_n[N_CLASS][32];
  char mass_vs_en_arr_n[N_CLASS][32];
  char mass_vs_pt_arr_n[N_CLASS][32];
  char mass_dist_arr_n[N_CLASS][32];
  char z_dist_arr_n[N_CLASS][32];
  for(Int_t c=0; c<N_CLASS; c++)
  {
    sprintf(pt_vs_eta_arr_n[c],"%s_pt_vs_eta_arr",ev->Name(c));
    sprintf(en_vs_eta_arr_n[c],"%s_en_vs_eta_arr",ev->Name(c));
    sprintf(pt_vs_phi_arr_n[c],"%s_pt_vs_phi_arr",ev->Name(c));
    sprintf(en_vs_phi_arr_n[c],"%s_en_vs_phi_arr",ev->Name(c));
    sprintf(eta_vs_phi_arr_n[c],"%s_eta_vs_phi_arr",ev->Name(c));
    sprintf(pt_vs_en_arr_n[c],"%s_pt_vs_en_arr",ev->Name(c));
    sprintf(z_vs_eta_arr_n[c],"%s_z_vs_eta_arr",ev->Name(c));
    sprintf(z_vs_phi_arr_n[c],"%s_z_vs_phi_arr",ev->Name(c));
    sprintf(mass_vs_en_arr_n[c],"%s_mass_vs_en_arr",ev->Name(c));
    sprintf(mass_vs_pt_arr_n[c],"%s_mass_vs_pt_arr",ev->Name(c));
    sprintf(mass_dist_arr_n[c],"%s_mass_dist_arr",ev->Name(c));
    sprintf(z_dist_arr_n[c],"%s_z_dist_arr",ev->Name(c));
    pt_vs_eta_arr[c] = (TObjArray*) infile->Get(pt_vs_eta_arr_n[c]);
    en_vs_eta_arr[c] = (TObjArray*) infile->Get(en_vs_eta_arr_n[c]);
    pt_vs_phi_arr[c] = (TObjArray*) infile->Get(pt_vs_phi_arr_n[c]);
    en_vs_phi_arr[c] = (TObjArray*) infile->Get(en_vs_phi_arr_n[c]);
    eta_vs_phi_arr[c] = (TObjArray*) infile->Get(eta_vs_phi_arr_n[c]);
    pt_vs_en_arr[c] = (TObjArray*) infile->Get(pt_vs_en_arr_n[c]);
    z_vs_eta_arr[c] = (TObjArray*) infile->Get(z_vs_eta_arr_n[c]);
    z_vs_phi_arr[c] = (TObjArray*) infile->Get(z_vs_phi_arr_n[c]);
    mass_vs_en_arr[c] = (TObjArray*) infile->Get(mass_vs_en_arr_n[c]);
    mass_vs_pt_arr[c] = (TObjArray*) infile->Get(mass_vs_pt_arr_n[c]);
    mass_dist_arr[c] = (TObjArray*) infile->Get(mass_dist_arr_n[c]);
    z_dist_arr[c] = (TObjArray*) infile->Get(z_dist_arr_n[c]);
  };

  TH1D * trig_dist = (TH1D*) infile->Get("trig_dist");
  trig_dist->SetBarWidth(0.8);
  trig_dist->SetBarOffset(0.3);
  trig_dist->SetFillColor(kBlue);

  char diag_dir[16]; strcpy(diag_dir,"diag_plots");
  char trig_canv_png[64]; sprintf(trig_canv_png,"trig.png");
  char trig_canv_png_print[64]; sprintf(trig_canv_png_print,"%s/trig.png",diag_dir);

  char mkdirstr[512];
  sprintf(mkdirstr,".! mkdir -p diag_plots");
  gROOT->ProcessLine(mkdirstr);
  for(Int_t c=0; c<N_CLASS; c++) {
    sprintf(mkdirstr,".! mkdir -p diag_plots/%s_corr",ev->Name(c));
    gROOT->ProcessLine(mkdirstr);
    sprintf(mkdirstr,".! mkdir -p diag_plots/%s_corr_MZ",ev->Name(c));
    gROOT->ProcessLine(mkdirstr);
  };


  // png names
  char corr_canv_png[N_CLASS][N_TRIG][64];
  char corr_canv_png_print[N_CLASS][N_TRIG][64];
  char corr_MZ_canv_png[N_CLASS][N_TRIG][64];
  char corr_MZ_canv_png_print[N_CLASS][N_TRIG][64];
  for(Int_t t=0; t<N_TRIG; t++)
  {
    for(Int_t c=0; c<N_CLASS; c++) 
    {
      sprintf(corr_MZ_canv_png[c][t],"%s_corr_MZ/%s_corr_MZ_%s.png",
        ev->Name(c),ev->Name(c),(T->Name(t)).Data());
      sprintf(corr_MZ_canv_png_print[c][t],"%s/%s_corr_MZ/%s_corr_MZ_%s.png",
        diag_dir,ev->Name(c),ev->Name(c),(T->Name(t)).Data());
      sprintf(corr_canv_png[c][t],"%s_corr/%s_corr_%s.png",
       ev->Name(c),ev->Name(c),(T->Name(t)).Data());
      sprintf(corr_canv_png_print[c][t],"%s/%s_corr/%s_corr_%s.png",
       diag_dir,ev->Name(c),ev->Name(c),(T->Name(t)).Data());
    };
  };


  // print pdfs & pngs
  trig_canv->cd();
  trig_dist->Draw("bar2");
  trig_canv->Print(trig_canv_png_print,"png");
  for(Int_t t=0; t<N_TRIG; t++)
  {
    for(Int_t c=0; c<N_CLASS; c++)
    {
      corr_canv[c]->cd(1); pt_vs_eta_arr[c]->At(t)->Draw("colz");
      corr_canv[c]->cd(2); pt_vs_phi_arr[c]->At(t)->Draw("colz");
      corr_canv[c]->cd(3); en_vs_eta_arr[c]->At(t)->Draw("colz");
      corr_canv[c]->cd(4); en_vs_phi_arr[c]->At(t)->Draw("colz");
      corr_canv[c]->cd(5); eta_vs_phi_arr[c]->At(t)->Draw("colz");
      corr_canv[c]->cd(6); pt_vs_en_arr[c]->At(t)->Draw("colz");
      corr_canv[c]->Print(corr_canv_png_print[c][t],"png");
      if(c!=ev->Idx("sph"))
      {
        corr_MZ_canv[c]->cd(1); z_vs_eta_arr[c]->At(t)->Draw("colz");
        corr_MZ_canv[c]->cd(2); z_vs_phi_arr[c]->At(t)->Draw("colz");
        corr_MZ_canv[c]->cd(3); mass_vs_en_arr[c]->At(t)->Draw("colz");
        corr_MZ_canv[c]->cd(4); mass_vs_pt_arr[c]->At(t)->Draw("colz");
        corr_MZ_canv[c]->cd(5); mass_dist_arr[c]->At(t)->Draw();
        corr_MZ_canv[c]->cd(6); z_dist_arr[c]->At(t)->Draw();
        corr_MZ_canv[c]->Print(corr_MZ_canv_png_print[c][t],"png");
      };
    };
  };


  // build html page
  gSystem->RedirectOutput("diag_plots/diag_web.html","w");
  printf("<html>\n<head><title>Kinematic Correlations</title></head>\n<body>\n");
  printf("<img src=%s>\n",trig_canv_png);
  printf("<div style=\"width:%dpx\">\n",RESX*(N_TRIG+1));
  for(Int_t c=0; c<N_CLASS; c++)
  {
    printf("<hr />\n");
    for(Int_t t=0; t<N_TRIG; t++) printf("<img src=%s>\n",corr_canv_png[c][t]);
    if(c!=ev->Idx("sph")) for(Int_t t=0; t<N_TRIG; t++) printf("<img src=%s>\n",corr_MZ_canv_png[c][t]);
  };
  printf("</div>\n</body>\n</html>\n");
  gSystem->RedirectOutput(0);
};
