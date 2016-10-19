// void draws kinvar vs. runindex for specified file and event class

void DrawKinvarPlot(TString filename="kinvarset/FMSOR_Pt_vs_run.root",
                    TString classname="pi0") {
  TFile * infile = new TFile(filename.Data(),"READ");

  char trig[32];
  char typ[32];
  char tmp[128];
  TString tmpstr = filename;
  tmpstr.ReplaceAll("_"," ");
  TRegexp re("^.*\/");
  tmpstr(re) = "";
  cout << tmpstr << endl;
  sscanf(tmpstr.Data(),"%s %s",trig,typ);
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

  // temporary hack to cut out fluctuations of pt_thresh vs. run to 0
  Double_t xx,yy;
  Int_t gtn = gt->GetN();
  for(int ii=0; ii<gtn; ii++) {
    gt->GetPoint(ii,xx,yy);
    printf("xx=%f yy=%f\n",xx,yy);
    if(yy<0.001) gt->RemovePoint(ii);
  };
   
  // drawing range
  Float_t xmin = h->GetXaxis()->GetXmin();
  Float_t xmax = h->GetXaxis()->GetXmax();
  xmax = 515;
  if(xmax<gtn) { 
    fprintf(stderr,"\nWARNING: xmax < gt->GetN(); graphs will be cut off!\n\n");
  };
  
  TCanvas * c = new TCanvas("c","c",1500,700);
  c->Divide(1,2);
  gStyle->SetOptStat(0);
  for(int x=1; x<=2; x++) c->GetPad(x)->SetGrid(1,1);
  c->cd(1);
  c->GetPad(1)->SetLogz();
  h->GetXaxis()->SetRangeUser(xmin,xmax);
  h->Draw("colz");
  if(!strcmp(typ,"Pt")) gt->Draw("LX");
  c->cd(2);
  g->SetFillColor(kGray);
  printf("%f %f\n",xmin,xmax);
  g->GetXaxis()->SetLimits(xmin,xmax);
  g->GetXaxis()->SetRangeUser(xmin,xmax);
  /*
  g->Draw("A3");
  g->Draw("PLX");
  */
  g->Draw("APLX");
  
  TString outname = filename.ReplaceAll(".root"," "+classname);
  outname = outname+".png";
  c->Print(outname.Data(),"png");
};
