// void draws kinvar vs. runindex for specified file and event class

void DrawKinvarPlot(TString filename="kinvarset/All_Pt_vs_run.root",
                    TString classname="pi0") {
  TFile * infile = new TFile(filename.Data(),"READ");

  char trig[32];
  char typ[32];
  char tmp[128];
  TString tmpstr = filename.ReplaceAll("_"," ");
  sscanf(filename.Data(),"kinvarset/%s %s",trig,typ);
  printf("%s %s\n",trig,typ);

  TString hname = "h2_"+classname;
  TH2D * h = (TH2D*)infile->Get(hname.Data());
  TString gname = "g_"+classname;
  TGraphErrors * g = (TGraphErrors*)infile->Get(gname.Data());
  TString gtname = "gt_"+classname;
  TGraphErrors * gt;
  if(!strcmp(typ,"Pt")) gt = (TGraphErrors*)infile->Get(gtname.Data());
  
  char htitle[256];
  char hnewtitle[512];
  strcpy(htitle,h->GetTitle());
  sprintf(hnewtitle,"%s -- %s triggers",htitle,trig);
  h->SetTitle(hnewtitle);

  h->SetMinimum(0.001);

  if(!strcmp(typ,"Pt")) {
    gt->SetLineWidth(2);
    gt->SetLineColor(kBlack);
  };
  
  TCanvas * c = new TCanvas("c","c",1500,700);
  c->Divide(1,2);
  gStyle->SetOptStat(0);
  for(int x=1; x<=2; x++) c->GetPad(x)->SetGrid(1,1);
  c->cd(1);
  c->GetPad(1)->SetLogz();
  h->Draw("colz");
  if(!strcmp(typ,"Pt")) gt->Draw("LX");
  c->cd(2);
  g->SetFillColor(kGray);
  Float_t xmin = h->GetXaxis()->GetXmin();
  Float_t xmax = h->GetXaxis()->GetXmax();
  g->GetXaxis()->SetRangeUser(xmin,xmax);
  g->Draw("A3");
  g->Draw("PLX");
  
  TString outname = filename.ReplaceAll(".root",".png");
  c->Print(outname.Data(),"png");
};
