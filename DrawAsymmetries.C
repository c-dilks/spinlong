// draws asymmetries all on one plot

void DrawAsymmetries(const char * evclass="pi0", const char * filetype="png", const char * asym_file="spin.root")
{
  // ANALYSIS TYPE -- use this swtich to change between a longitudinal and a transverse analysis
  enum atypes {kLong,kTrans};
  const Int_t ANALYSIS_TYPE = kLong;  // ---- SWITCH ---- //

  Int_t NPARAM_tmp;
  switch(ANALYSIS_TYPE) {
    case kLong:
      NPARAM_tmp = 1;
      break;
    case kTrans:
      NPRAM_tmp = 2;
      break;
  };
  const Int_t NPARAM = NPARAM_tmp;


  // open root file
  TFile * asym_tfile = new TFile(asym_file,"READ");

  // environment
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  EventClass * ev = new EventClass(RD->env);


  // asymmetry titles
  const Int_t asym_bins=4;
  char asymmetry[asym_bins][8];
  strcpy(asymmetry[1],"Y-SSA");
  strcpy(asymmetry[2],"B-SSA");
  strcpy(asymmetry[3],"DSA");
  // written out asymmetry titles
  char asymmetry_w[asym_bins][128];
  sprintf(asymmetry_w[1],"yellow single spin asymmetry (%s triggers)",RD->env->TriggerType);
  sprintf(asymmetry_w[2],"blue single spin asymmetry (%s triggers)",RD->env->TriggerType);
  sprintf(asymmetry_w[3],"double spin asymmetry (%s triggers)",RD->env->TriggerType);

  // define asymmetry kinematic dependence plots
  char kindep_name[NPARAM][asym_bins][128];
  char xaxistitle[64];
  Int_t num_bins;
  char dir_title[NPARAM][asym_bins][32];
  if(ANALYSIS_TYPE==kTrans) {
    strcpy(dir_title[0][1],"R_yellow");
    strcpy(dir_title[0][2],"R_blue");
    strcpy(dir_title[0][3],"A_Sigma");
    strcpy(dir_title[1][1],"A_N_yellow");
    strcpy(dir_title[1][2],"A_N_blue");
    strcpy(dir_title[1][3],"A_TT");
  }
  else if(ANALYSIS_TYPE==kLong) {
    strcpy(dir_title[0][1],"A_L_yellow");
    strcpy(dir_title[0][2],"A_L_blue");
    strcpy(dir_title[0][3],"A_LL");
  };
  char kindep_main_title[64];
  if(RD->env->PtBins==1 && RD->env->EnBins!=1 && RD->env->EtaBins==1) 
  {
    for(Int_t z=0; z<NPARAM; z++)
    {
      for(Int_t a=1; a<asym_bins; a++) sprintf(kindep_name[z][a],"/%s/en_dep_z%d_a%d_g0_p0",dir_title[z][a],z,a);
    };
    strcpy(xaxistitle,"E (GeV)");
    sprintf(kindep_main_title,"%s asymmetries vs. E (%s triggers)",ev->Title(evclass),RD->env->TriggerType);
    num_bins = RD->env->EnBins;
  }
  else if(RD->env->PtBins!=1 && RD->env->EnBins==1 && RD->env->EtaBins==1) 
  {
    for(Int_t z=0; z<NPARAM; z++)
    {
      for(Int_t a=1; a<asym_bins; a++) sprintf(kindep_name[z][a],"/%s/pt_dep_z%d_a%d_g0_e0",dir_title[z][a],z,a);
    };
    strcpy(xaxistitle,"p_{#perp}  (GeV/c)");
    sprintf(kindep_main_title,"%s asymmetries vs. p_{T} (%s triggers)",ev->Title(evclass),RD->env->TriggerType);
    num_bins = RD->env->PtBins;
  }
  else if(RD->env->PtBins==1 && RD->env->EnBins==1 && RD->env->EtaBins==1)
  {
    for(Int_t z=0; z<NPARAM; z++)
    {
      for(Int_t a=1; a<asym_bins; a++) sprintf(kindep_name[z][a],"/%s/pt_dep_z%d_a%d_g0_e0",dir_title[z][a],z,a);
    };
    sprintf(kindep_main_title,"%s asymmetries",ev->Title(evclass));
    strcpy(xaxistitle,"single bin");
    num_bins = 1;
  }
  else
  {
    printf("\n<><><><><><><><><><>\n\n");
    printf("pt_bins>1 && en_bins>1 ----------> three.png NOT DRAWN\n");
    printf("spin.root produced\n");
    return;
  };

  // get asym kin dep plots
  TGraphErrors * kindep_gr[NPARAM][asym_bins];
  for(Int_t z=0; z<NPARAM; z++)
  {
    for(Int_t a=1; a<asym_bins; a++) 
    {
      kindep_gr[z][a] = (TGraphErrors*) asym_tfile->Get(kindep_name[z][a]);
      kindep_gr[z][a]->GetXaxis()->SetLabelSize(0.05);
      kindep_gr[z][a]->GetYaxis()->SetLabelSize(0.05);
      kindep_gr[z][a]->GetXaxis()->SetTitleSize(0.06);
      kindep_gr[z][a]->GetYaxis()->SetTitleSize(0.06);
      kindep_gr[z][a]->GetYaxis()->SetTitleOffset(0.8);
      kindep_gr[z][a]->SetMarkerSize(1.3);
      kindep_gr[z][a]->SetMarkerStyle(kFullCircle);
      kindep_gr[z][a]->SetLineWidth(2);
    };
    kindep_gr[z][1]->SetMarkerColor(kBlack);
    kindep_gr[z][2]->SetMarkerColor(kBlack);
    kindep_gr[z][3]->SetMarkerColor(kBlack);
    kindep_gr[z][1]->SetLineColor(kOrange-3);
    kindep_gr[z][2]->SetLineColor(kBlue);
    kindep_gr[z][3]->SetLineColor(kRed);
  };


  // get analysing power vs. phi plots
  // note that for these hists, z=0 == z=1
  const Int_t num_bins0 = num_bins;
  char asym_name[NPARAM][asym_bins][num_bins0][128];
  char asym_title[NPARAM][asym_bins][num_bins0][256];
  TH1D * asym_hist[NPARAM][asym_bins][num_bins0];
  char asym_main_title[asym_bins][100];
  if(RD->env->PtBins==1 && RD->env->EnBins!=1 && RD->env->EtaBins==1) 
  {
    for(Int_t z=0; z<NPARAM; z++)
    {
      for(Int_t a=1; a<asym_bins; a++) 
      {
        for(Int_t e=0; e<RD->env->EnBins; e++)
        {
          sprintf(asym_name[z][a][e],"/%s/asym_a%d_g0_p0_e%d",dir_title[z][a],a,e);
          sprintf(asym_title[z][a][e],"E #in [%.2f,%.2f)",RD->env->EnDiv(e),RD->env->EnDiv(e+1));
          asym_hist[z][a][e] = (TH1D*) asym_tfile->Get(asym_name[z][a][e]);
          asym_hist[z][a][e]->SetTitle(asym_title[z][a][e]);
          asym_hist[z][a][e]->GetXaxis()->SetTitle("#phi");
          asym_hist[z][a][e]->GetYaxis()->SetTitle(asymmetry[a]);
          asym_hist[z][a][e]->GetXaxis()->SetLabelSize(0.05);
          asym_hist[z][a][e]->GetYaxis()->SetLabelSize(0.05);
          asym_hist[z][a][e]->GetXaxis()->SetTitleSize(0.06);
          asym_hist[z][a][e]->GetYaxis()->SetTitleSize(0.06);
          asym_hist[z][a][e]->GetYaxis()->SetTitleOffset(0.8);
        };
        sprintf(asym_main_title[a],"%s %s vs #phi",ev->Title(evclass),asymmetry_w[a]);
      };
    };
  }
  else if(RD->env->PtBins!=1 && RD->env->EnBins==1 && RD->env->EtaBins==1) 
  {
    for(Int_t z=0; z<NPARAM; z++)
    {
      for(Int_t a=1; a<asym_bins; a++) 
      {
        for(Int_t p=0; p<RD->env->PtBins; p++)
        {
          sprintf(asym_name[z][a][p],"/%s/asym_a%d_g0_p%d_e0",dir_title[z][a],a,p);
          sprintf(asym_title[z][a][p],"p_{#perp} #in [%.2f,%.2f)",RD->env->PtDiv(p),RD->env->PtDiv(p+1));
          asym_hist[z][a][p] = (TH1D*) asym_tfile->Get(asym_name[z][a][p]);
          asym_hist[z][a][p]->SetTitle(asym_title[z][a][p]);
          asym_hist[z][a][p]->GetXaxis()->SetTitle("#phi");
          asym_hist[z][a][p]->GetYaxis()->SetTitle(asymmetry[a]);
          asym_hist[z][a][p]->GetXaxis()->SetLabelSize(0.05);
          asym_hist[z][a][p]->GetYaxis()->SetLabelSize(0.05);
          asym_hist[z][a][p]->GetXaxis()->SetTitleSize(0.06);
          asym_hist[z][a][p]->GetYaxis()->SetTitleSize(0.06);
          asym_hist[z][a][p]->GetYaxis()->SetTitleOffset(0.8);
        };
        sprintf(asym_main_title[a],"%s %s vs #phi",ev->Title(evclass),asymmetry_w[a]);
      };
    };
  }
  else if(RD->env->PtBins==1 && RD->env->EnBins==1 && RD->env->EtaBins==1)
  {
    for(Int_t z=0; z<NPARAM; z++)
    {
      for(Int_t a=1; a<asym_bins; a++) 
      {
        sprintf(asym_name[z][a][0],"/%s/asym_a%d_g0_p0_e0",dir_title[z][a],a);
        sprintf(asym_title[z][a][0],"single bin");
        asym_hist[z][a][0] = (TH1D*) asym_tfile->Get(asym_name[z][a][0]);
        asym_hist[z][a][0]->SetTitle(asym_title[z][a][0]);
        asym_hist[z][a][0]->GetXaxis()->SetTitle("#phi");
        asym_hist[z][a][0]->GetYaxis()->SetTitle(asymmetry[a]);
        asym_hist[z][a][0]->GetXaxis()->SetLabelSize(0.05);
        asym_hist[z][a][0]->GetYaxis()->SetLabelSize(0.05);
        asym_hist[z][a][0]->GetXaxis()->SetTitleSize(0.06);
        asym_hist[z][a][0]->GetYaxis()->SetTitleSize(0.06);
        asym_hist[z][a][0]->GetYaxis()->SetTitleOffset(0.8);
        sprintf(asym_main_title[a],"%s %s vs #phi",ev->Title(evclass),asymmetry_w[a]);
      };
    };
  };


  // canvas sizes
  Float_t hsize = 800;
  Float_t vsize = 1000;


  // draw asym kin dep canvas
  TCanvas * kindep_canv[NPARAM];
  char kindep_canv_n[NPARAM][32];
  for(Int_t z=0; z<NPARAM; z++)
  {
    sprintf(kindep_canv_n[z],"canv_kindep_%d",z);
    kindep_canv[z] = new TCanvas(kindep_canv_n[z],kindep_canv_n[z],hsize,vsize);
  };
  TPad * kindep_pad[NPARAM][asym_bins];
  char kindep_pad_n[NPARAM][asym_bins][16];
  TLine * zero_line[NPARAM][asym_bins];
  Float_t padding = 0.05;
  Float_t extra_bottom = 0.04;
  Float_t extra_left = 0.1;
  Float_t interval = (1-2*padding-extra_bottom)/(asym_bins-1);
  TPaveText * kindep_pave = new TPaveText(0.25,0.96,0.75,0.99,"br");
  kindep_pave->AddText(kindep_main_title);
  for(Int_t z=0; z<NPARAM; z++)
  {
    for(Int_t a=1; a<asym_bins; a++)
    {
      zero_line[z][a] = new TLine(kindep_gr[z][a]->GetXaxis()->GetXmin(),0,kindep_gr[z][a]->GetXaxis()->GetXmax(),0);
      zero_line[z][a]->SetLineColor(kCyan+3);
      zero_line[z][a]->SetLineWidth(3);
      zero_line[z][a]->SetLineStyle(2);
      sprintf(kindep_pad_n[z][a],"kpad_%d_%d",z,a);
      if(a==asym_bins-1)
        kindep_pad[z][a] = new TPad(kindep_pad_n[z][a],kindep_pad_n[z][a],
          padding,padding+(asym_bins-a-1)*interval,
          1-padding,padding+extra_bottom+(asym_bins-a)*interval,
          0,0);
      else
        kindep_pad[z][a] = new TPad(kindep_pad_n[z][a],kindep_pad_n[z][a],
          padding,padding+extra_bottom+(asym_bins-a-1)*interval,
          1-padding,padding+extra_bottom+(asym_bins-a)*interval,
          0,0);
      if(a==asym_bins-1) 
      {
        kindep_pad[z][a]->SetTopMargin(0);
        kindep_pad[z][a]->SetBottomMargin(extra_bottom/(interval+extra_bottom));
      }
      else 
      {
        kindep_pad[z][a]->SetBottomMargin(0);
        kindep_pad[z][a]->SetTopMargin(0);
      };
      kindep_canv[z]->cd();
      kindep_pad[z][a]->SetLeftMargin(extra_left);
      kindep_pad[z][a]->SetGrid(1,1);
      kindep_pad[z][a]->Draw();
      kindep_pad[z][a]->cd();
      kindep_gr[z][a]->Draw("APE");
      zero_line[z][a]->Draw();
      kindep_canv[z]->cd();
    };
    kindep_canv[z]->cd();
    kindep_pave->Draw();
  };


  // draw analysing power canvases
  char canv_filename[64];
  sprintf(canv_filename,"asymcanv_%s.root",evclass);
  TFile * outfile = new TFile(canv_filename,"RECREATE");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  TCanvas * asym_canv[NPARAM][asym_bins];
  char asym_canv_n[NPARAM][asym_bins][32];
  TPad * asym_pad[NPARAM][asym_bins][num_bins0];
  char asym_pad_n[NPARAM][asym_bins][num_bins0][8];
  interval = (1-2*padding-extra_bottom)/(num_bins0);
  TPaveText * asym_pave[NPARAM][asym_bins];
  for(Int_t z=0; z<NPARAM; z++)
  {
    for(Int_t a=1; a<asym_bins; a++) 
    {
      sprintf(asym_canv_n[z][a],"canv_z%d_%s",z,asymmetry[a]);
      asym_canv[z][a] = new TCanvas(asym_canv_n[z][a],asym_canv_n[z][a],hsize,num_bins0*vsize/(asym_bins-1));
      asym_pave[z][a] = new TPaveText(0.25,0.96,0.75,0.99,"br");
      asym_pave[z][a]->AddText(asym_main_title[a]);
      for(Int_t n=0; n<num_bins0; n++)
      {
        sprintf(asym_pad_n[z][a][n],"p%d_%d_%d",z,a,n);
        if(n==num_bins0-1)
          asym_pad[z][a][n] = new TPad(asym_pad_n[z][a][n],asym_pad_n[z][a][n],
          padding,padding+(num_bins0-n-1)*interval,
          1-padding,padding+extra_bottom+(num_bins0-n)*interval,
          0,0);
        else
          asym_pad[z][a][n] = new TPad(asym_pad_n[z][a][n],asym_pad_n[z][a][n],
          padding,padding+extra_bottom+(num_bins0-n-1)*interval,
          1-padding,padding+extra_bottom+(num_bins0-n)*interval,
          0,0);
        if(n==num_bins0-1) 
        {
          asym_pad[z][a][n]->SetTopMargin(0);
          asym_pad[z][a][n]->SetBottomMargin(extra_bottom/(interval+extra_bottom));
        }
        else 
        {
          asym_pad[z][a][n]->SetBottomMargin(0);
          asym_pad[z][a][n]->SetTopMargin(0);
        };
        asym_pad[z][a][n]->SetLeftMargin(extra_left);
        asym_pad[z][a][n]->SetGrid(1,1);
        asym_pad[z][a][n]->Draw();
        asym_pad[z][a][n]->cd();
        asym_hist[z][a][n]->Draw();
        asym_canv[z][a]->cd();
      };
      asym_canv[z][a]->cd();
      asym_pave[z][a]->Draw();
    };
  };

  char kindep_canv_png[NPARAM][128];
  for(Int_t z=0; z<NPARAM; z++)
  {
    sprintf(kindep_canv_png[z],"canv_kindep_%d.%s",z,filetype);
    kindep_canv[z]->Write();
    kindep_canv[z]->Print(kindep_canv_png[z],filetype);
  };

  // note that for asym hists, z=0 == z=1
  char asym_canv_png[asym_bins][128];
  for(Int_t a=1; a<asym_bins; a++) 
  {
    sprintf(asym_canv_png[a],"%s.%s",asym_canv_n[0][a],filetype);
    asym_canv[0][a]->Write();
    asym_canv[0][a]->Print(asym_canv_png[a],filetype);
  };
}
