
void BGcorr(TString dir="test_12_large",
            TString classname0="pi0" /*measured signal*/,
            TString classname1="bkg" /*background*/
) {
  const Float_t PURITY = 0.85;
  const Float_t PURITY_ERROR = 0.03;
  const Float_t ASYM_CORR = 0;

  enum sigbkg {kM,kB};

  TString filename[2];
  filename[0] = "asym_plots/"+dir+"/spin_"+classname0+".root";
  filename[1] = "asym_plots/"+dir+"/spin_"+classname1+".root";

  TFile * infile[2];
  TGraphErrors * graph[2];
  int i;
  for(i=0; i<2; i++) {
    infile[i] = new TFile(filename[i].Data(),"READ");
    infile[i]->cd("A_LL");
    graph[i] = (TGraphErrors*) infile[i]->Get("/A_LL/pt_dep_z0_a3_g0_e0");
  };
  TGraphErrors * bgcorr = new TGraphErrors();

  graph[0]->SetMarkerColor(kRed);
  graph[1]->SetMarkerColor(kBlue);
  bgcorr->SetMarkerColor(kGreen+2);
  graph[0]->SetLineColor(kRed);
  graph[1]->SetLineColor(kBlue);
  bgcorr->SetLineColor(kGreen+2);
  for(i=0;i<2;i++) 
  {
    graph[i]->SetMarkerStyle(kFullCircle);
    graph[i]->GetXaxis()->SetLabelSize(0.08);
    graph[i]->GetYaxis()->SetLabelSize(0.08);
  };
  bgcorr->SetMarkerStyle(kFullCircle);
  bgcorr->GetXaxis()->SetLabelSize(0.08);
  bgcorr->GetYaxis()->SetLabelSize(0.08);

  TMultiGraph * mg = new TMultiGraph();
  for(i=0;i<2;i++) mg->Add(graph[i]);
  TString mg_t = Form("A_{LL}^{M} (red) and A_{LL}^{B} (blue) vs. p_{T}^{S(B)}");
  mg->SetTitle(mg_t.Data());

  TMultiGraph * mg2 = new TMultiGraph();
  mg2->Add(graph[kM]);
  mg2->Add(bgcorr);
  TString mg2_t = Form("A_{LL}^{S} (green) and A_{LL}^{M} (red) vs. p_{T}^{S(B)}");
  mg2->SetTitle(mg2_t.Data());


  TGraphErrors * diff = new TGraphErrors();
  TGraphErrors * diff2 = new TGraphErrors();
  Double_t x[2];
  Double_t y[2];
  Double_t ye[2];
  Double_t err;
  for(Int_t k=0;k<graph[0]->GetN();k++)
  {
    for(i=0;i<2;i++) {
      graph[i]->GetPoint(k,x[i],y[i]);
      ye[i] = graph[i]->GetErrorY(k);
    };
    err = TMath::Sqrt(TMath::Power(ye[0],2)+TMath::Power(ye[1],2));
    diff->SetPoint(k,x[0],y[0]-y[1]);
    diff->SetPointError(k,0,err);
  };

  diff->SetMarkerStyle(kFullCircle);
  diff->SetMarkerColor(kBlack);
  diff->SetTitle("A_{LL}^{M}-A_{LL}^{B} vs. p_{T}^{S}");

  diff2->SetMarkerStyle(kFullCircle);
  diff2->SetMarkerColor(kBlack);
  diff2->SetTitle("A_{LL}^{S}-A_{LL}^{M} vs. p_{T}^{S}-p_{T}^{M} :: errors are #sigma^{S}-#sigma^{M}");

  Double_t ab,am,as,dam,dab,f,df,sf,sam,sab,vas,cas,sas;
  Double_t pb,pm,ps,dpm,dpb,spm,spb,vps,cps,sps;

  f = PURITY; // F
  dam = 1/f; // derivative A_S wrt A_M
  dab = (f-1)/f; // derivative A_S wrt A_B
  dpm = dam; // derivative pT_S wrt pT_M
  dpb = dab; // derivative pT_S wrt pT_B
  sf = PURITY_ERROR; // sigma F

  for(Int_t k=0; k<graph[0]->GetN(); k++) {
    graph[kM]->GetPoint(k,pm,am); // measured A and pT
    graph[kB]->GetPoint(k,pb,ab); // bkg A and pT

    as = (1/f) * am - ((1-f)/f) * ab; // corrected A
    ps = (1/f) * pm - ((1-f)/f) * pb; // corrected pT

    sam = graph[kM]->GetErrorY(k); // sigma A_M
    sab = graph[kB]->GetErrorY(k); // sigma A_B
    spm = graph[kM]->GetErrorX(k); // sigma pT_M
    spb = graph[kB]->GetErrorX(k); // sigma pT_B

    df = (ab-am)/(TMath::Power(f,2)); // derivative A_S wrt F

    // diagonal terms
    vas = TMath::Power(dam,2) * TMath::Power(sam,2) +
          TMath::Power(dab,2) * TMath::Power(sab,2) +
          TMath::Power(df,2) * TMath::Power(sf,2);
    vps = TMath::Power(dpm,2) * TMath::Power(spm,2) +
          TMath::Power(dpb,2) * TMath::Power(spb,2) +
          TMath::Power(df,2) * TMath::Power(sf,2);

    // off-diagonal terms
    cas = 2 * dam * dab * ASYM_CORR * sam * sab;
    cps = 2 * dpm * dpb * ASYM_CORR * spm * spb;

    // full error
    sas = TMath::Sqrt(vas + cas);
    sps = TMath::Sqrt(vps + cps);


    // set new point
    bgcorr->SetPoint(k,ps,as);
    bgcorr->SetPointError(k,sps,sas);
    diff2->SetPoint(k,ps-pm,as-am);
    diff2->SetPointError(k,sps-spm,sas-sam);
  };
  bgcorr->SetTitle("background-corrected A_{LL} vs. p_{T}");

  TCanvas * c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);
  for(k=1;k<5;k++)c1->GetPad(k)->SetGrid(1,1);
  c1->cd(1);
  mg->Draw("APE");
  mg->GetXaxis()->SetLabelSize(0.08);
  mg->GetYaxis()->SetLabelSize(0.08);

  c1->cd(2);
  mg2->Draw("APE");
  mg2->GetXaxis()->SetLabelSize(0.08);
  mg2->GetYaxis()->SetLabelSize(0.08);

  c1->cd(3);
  diff->GetXaxis()->SetLabelSize(0.08);
  diff->GetYaxis()->SetLabelSize(0.08);
  diff->Fit("pol0","","",1,10);
  gStyle->SetOptFit(1);
  diff->Draw("AP");

  c1->cd(4);
  diff2->GetXaxis()->SetLabelSize(0.08);
  diff2->GetYaxis()->SetLabelSize(0.08);
  //diff2->Fit("pol0","","",1,10);
  //gStyle->SetOptFit(1);
  diff2->Draw("APE");

};
  
