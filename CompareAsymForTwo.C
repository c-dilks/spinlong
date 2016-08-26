// compares A_LL given two files

void CompareAsymForTwo(TString dir0="13_using_lw_pol",
                       TString dir1="13_using_tdep_pol",
                       TString classname0="pi0",
                       TString classname1="pi0"
) {
  TString filename[2];
  filename[0] = "asym_plots/"+dir0+"/spin_"+classname0+".root";
  filename[1] = "asym_plots/"+dir1+"/spin_"+classname1+".root";

  TFile * infile[2];
  TGraphErrors * graph[2];
  int i;
  for(i=0; i<2; i++) {
    infile[i] = new TFile(filename[i].Data(),"READ");
    infile[i]->cd("A_LL");
    graph[i] = (TGraphErrors*) infile[i]->Get("/A_LL/pt_dep_z0_a3_g0_e0");
  };

  graph[0]->SetMarkerColor(kRed);
  graph[1]->SetMarkerColor(kBlue);
  graph[0]->SetLineColor(kRed);
  graph[1]->SetLineColor(kBlue);
  for(i=0;i<2;i++) 
  {
    graph[i]->SetMarkerStyle(kFullCircle);
    graph[i]->GetXaxis()->SetLabelSize(0.08);
    graph[i]->GetYaxis()->SetLabelSize(0.08);
  };

  TMultiGraph * mg = new TMultiGraph();
  for(i=0;i<2;i++) mg->Add(graph[i]);
  TString mg_t = Form("A_{LL} vs. p_{T} -- red for %s %s -- blue for %s %s",
    dir0.Data(),classname0.Data(),dir1.Data(),classname1.Data());
  mg->SetTitle(mg_t.Data());

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(1,2);
  for(k=1;k<3;k++)c1->GetPad(k)->SetGrid(1,1);
  c1->cd(1);
  mg->Draw("APE");
  mg->GetXaxis()->SetLabelSize(0.08);
  mg->GetYaxis()->SetLabelSize(0.08);

  TGraph * diff = new TGraphErrors();
  Double_t x[2];
  Double_t y[2];
  for(Int_t k=0;k<graph[0]->GetN();k++)
  {
    for(i=0;i<2;i++) graph[i]->GetPoint(k,x[i],y[i]);
    diff->SetPoint(k,x[0],y[0]-y[1]);
  };

  diff->SetMarkerStyle(kFullCircle);
  diff->SetMarkerColor(kBlack);
  diff->SetTitle("red A_{LL} minus blue A_{LL}");

  c1->cd(2);
  diff->GetXaxis()->SetLabelSize(0.08);
  diff->GetYaxis()->SetLabelSize(0.08);
  diff->Fit("pol0","","",1,10);
  gStyle->SetOptFit(1);
  diff->Draw("AP");

};
  
