void DrawA_LL(TString dir="output.large", Bool_t onlyALL=false)
{
  TString filename=dir+"/spin_pi0.root";
  TFile * infile = new TFile(filename.Data(),"READ");

  const Int_t N_DIRS = 4;
  TString directory[N_DIRS];
  Bool_t draw_plots[N_DIRS];
  int c,d;

  ///////           DRAW OPTIONS             /////////////
  /******************************************************/
  enum dirs {kALL,kALb,kALy,kBin};
  directory[kALL]="A_LL";            draw_plots[0]=true;
  directory[kALb]="A_L_blue";        draw_plots[1]=false;
  directory[kALy]="A_L_yellow";      draw_plots[2]=false;
  directory[kBin]="bin_weighting";   draw_plots[3]=true;
  /******************************************************/

  // override draw options if onlyALL==true
  if(onlyALL) { 
    for(d=0; d<N_DIRS; d++) draw_plots[d]=false;
    draw_plots[kALL]=true;
  };
  

  TCanvas * canv[N_DIRS][50];
  TString canv_n[N_DIRS][50];
  TGraphErrors * graph[N_DIRS][50];
  TKey * key;
  Float_t fit_lim;

  for(d=0; d<N_DIRS; d++) {
    if(draw_plots[d]) {
      c=0;
      directory[d]="/"+directory[d];
      infile->cd(directory[d].Data());
      TIter next(gDirectory->GetListOfKeys());
      while((key=(TKey*)next())) {
        if(gROOT->GetClass(key->GetClassName())->InheritsFrom("TCanvas"))
        {
          canv[d][c]=(TCanvas*)key->ReadObj();
          canv[d][c]->Draw();
          c++;
        }
        else if(gROOT->GetClass(key->GetClassName())->InheritsFrom("TGraphErrors")) {
          graph[d][c]=(TGraphErrors*)key->ReadObj();
          canv_n[d][c]=Form("%s_canv",key->GetName());
          if(onlyALL && canv_n[d][c].Data()[0]!='p') continue;

          graph[d][c]->SetLineColor(kRed);
          graph[d][c]->SetLineWidth(3);
          graph[d][c]->SetMarkerStyle(kFullCircle);
          graph[d][c]->SetMarkerColor(kBlack);
          graph[d][c]->SetMarkerSize(1.5);

          canv[d][c] = new TCanvas(canv_n[d][c].Data(),
            canv_n[d][c].Data(),1000,800);
          canv[d][c]->SetGrid(1,1);

          gStyle->SetOptFit(1);
          fit_lim = (canv_n[d][c].Data()[0]=='p') ? 1:10;
          graph[d][c]->Fit("pol0","","",fit_lim,10*fit_lim);

          graph[d][c]->Draw("APE");

          c++;
        };
      };
    };
  };
};
