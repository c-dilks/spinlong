#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooArgusBG.h"
#include "RooGenericPdf.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"

using namespace RooFit;

Double_t mu,sigma,alpha;
Double_t mu_eta,sigma_eta,alpha_eta;
Double_t ch[10];
Double_t pion_purity,background_frac;

void Purity(TString infile_n = "massset_12/all.root",
            Int_t polyOrder = 3,
            Bool_t fitEta = false) {
  // open histogram
  TFile * infile = new TFile(infile_n.Data(),"READ");
  TObjArray * arr = (TObjArray*) infile->Get("mass_dist_arr_g0_e1_p1");
  TH1D * h = (TH1D*) arr->At(1);
  ExeFit(h,polyOrder,fitEta);
  printf("#mu = %.4f GeV/c^{2}\n",mu);
  printf("#sigma = %.4f GeV/c^{2}\n",sigma);
  printf("#alpha = %.4f\n",alpha);
  if(fitEta) {
    printf("\n");
    printf("#mu(#eta) = %.4f GeV/c^{2}\n",mu_eta);
    printf("#sigma(#eta) = %.4f GeV/c^{2}\n",sigma_eta);
    printf("#alpha(#eta) = %.4f\n",alpha_eta);
  };
  printf("\n");
  for(int i=0; i<polyOrder; i++) {
    printf("a_{%d} = %.4f\n",i,ch[i]);
  };
  printf("\n");
  printf("pion_purity = %.2f%%\n",100*pion_purity);
  printf("background_frac = %.2f%%\n",100*background_frac);
  printf("sum = %.2f%%\n",100*(pion_purity+background_frac));
  printf("\n");

};

void ExeFit(TH1D * hh, Int_t polyOrder_, Bool_t fitEta_) {
  if(hh==NULL) { 
    fprintf("ERROR: hh = %p\n",(void*)hh);
    return;
  };

  Int_t chebOrder = polyOrder_;
  Bool_t pionOnly = !fitEta_;

  // pion mass, bounds, and ranges for fit and draw
  const Float_t PION_MASS = 0.135;
  const Float_t ETA_MASS = 0.548;
  const Float_t MASS_LB = 0;
  const Float_t MASS_UB = 1;
  Float_t fit_lb; // fit range
  Float_t fit_ub; // fit range
  Float_t draw_lb; // drawing range
  Float_t draw_ub; // drawing range
  if(pionOnly) {
    fit_lb = 0.1;
    fit_ub = 0.4;
    draw_lb = 0.05;
    draw_ub = 0.5;
  } else {
    fit_lb = 0.1;
    fit_ub = 0.7;
    draw_lb = 0.05;
    draw_ub = 1;
  };
  // mass distribution
  TString hh_t = TString(hh->GetTitle());
  RooRealVar mass("mass","M_{#gamma#gamma}",PION_MASS,MASS_LB,MASS_UB,"GeV/c^{2}");
  RooDataHist massdist("massdist",hh_t.Data(),mass,Import(*hh));


  // SIGNAL FIT MODELS :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  // pi0
  RooRealVar s0("s0","pion_mu",PION_MASS,fit_lb,fit_ub);
  RooRealVar s1("s1","pion_sigma",0.07,0.001,fit_ub-fit_lb);
  RooRealVar s2("s2","pion_alpha",3.2,0.1,20);
  TString formu_gaus = "TMath::Gaus((mass-s0)/s1,0,1)";
  TString formu_skew = "(1+TMath::Erf((1/TMath::Sqrt(2))*s2*((mass-s0)/s1)))";
  TString formu = formu_gaus+"*"+formu_skew;
  RooGenericPdf sig_skew("signal","skewed guassian",formu.Data(),
    RooArgList(mass,s0,s1,s2));


  // eta
  RooRealVar e0("e0","eta_mu",0.55,0.4,0.7);
  RooRealVar e1("e1","eta_sigma",0.5,0.001,2);
  RooRealVar e2("e2","eta_alpha",3.2,0.1,10);
  formu.ReplaceAll("s0","e0");
  formu.ReplaceAll("s1","e1");
  formu.ReplaceAll("s2","e2");
  RooGenericPdf sig_skew_eta("signal_eta","skewed guassian",formu.Data(),
    RooArgList(mass,e0,e1,e2));


  // BACKGROUND FIT MODELS ::::::::::::::::::::::::::::::::::::::::::::::::::

  // chebychev polynomials
  Float_t b_lb = -50.0;
  Float_t b_ub = 50.0;
  RooRealVar b0("b0","b0", -3.8, b_lb, b_ub);
  RooRealVar b1("b1","b1", -16.2, b_lb, b_ub);
  RooRealVar b2("b2","b2", 6.4, b_lb, b_ub);
  RooRealVar b3("b3","b3", b_lb, b_ub);
  RooRealVar b4("b4","b4", b_lb, b_ub);
  switch(chebOrder) {
    case 2:
     RooChebychev bg_cheb("background","bg_cheb",mass,RooArgSet(b0,b1));
     break;
    case 3:
     RooChebychev bg_cheb("background","bg_cheb",mass,RooArgSet(b0,b1,b2));
     break;
    case 4:
     RooChebychev bg_cheb("background","bg_cheb",mass,RooArgSet(b0,b1,b2,b3));
     break;
    case 5:
     RooChebychev bg_cheb("background","bg_cheb",mass,RooArgSet(b0,b1,b2,b3,b4));
     break;
    default:
     fprintf(stderr,"undefined chebychev order\n");
     return;
     break;
  };


  // argus function
  RooRealVar arg_m("arg_m","arg_m", 1.2, 0.8, 2);
  RooRealVar arg_c("arg_c","arg_c", -50, -100, 0);
  RooRealVar arg_p("arg_p","arg_p", 40, 0, 100);
  RooArgusBG bg_argus("background","bg_argus",mass,arg_m,arg_c,arg_p);


  // DEFINE FULL PDF ::::::::::::::::::::::::::::::::::::::::::::::::::::
  //Double_t data_int = hh->Integral(fit_lb,fit_ub);
  //Float_t purity_guess = 0.7;
  //RooRealVar sig_nn("sig_nn","sig_nn",purity_guess*data_int,0,data_int);
  //RooRealVar bg_nn("sig_nn","sig_nn",(1.0-purity_guess)*data_int,0,data_int);
  //RooRealVar eta_nn("eta_nn","eta_nn",0.5*data_int,0,data_int);
  
  RooRealVar frac("frac","background fraction",0.3,0.0,1.0);
  RooRealVar frac2("frac2","background fraction 2",0.3,0.0,1.0);
  if(pionOnly) {
    RooAddPdf model(
      "model",
      "model",
      RooArgList(sig_skew,bg_cheb),
      frac
    );
  } else {
    RooAddPdf model(
      "model",
      "model",
      RooArgList(sig_skew,sig_skew_eta,bg_cheb),
      RooArgSet(frac,frac2)
    );
  };


   
  // perform fit
  mass.setRange("fit_range",fit_lb,fit_ub);
  model.fitTo(massdist,Range("fit_range"));


  // build plot
  TString plotframe_t = hh_t;
  RooPlot * plotframe = mass.frame(Title(plotframe_t.Data()));
  massdist.plotOn(
    plotframe,
    DataError(RooAbsData::SumW2),
    MarkerColor(kRed)
  );
  model.plotOn(
    plotframe,
    Range("Full"),
    LineColor(kBlack),
    LineWidth(3)
  );
  model.plotOn(
    plotframe,
    Range("Full"),
    Components("background"),
    LineStyle(kDashed),
    LineColor(kCyan+1)
  );
  model.plotOn(
    plotframe,
    Range("Full"),
    Components("signal"),
    LineColor(kGreen+1)
  );
  if(!pionOnly) {
    model.plotOn(
      plotframe,
      Range("Full"),
      Components("signal_eta"),
      LineStyle(kDashed),
      LineColor(kRed)
    );
  };
  //model.paramOn(plotframe);
  printf("\n\n");
  model.Print("t");
  printf("\n\n");



  // obtain parameters
  mu = s0.getVal();
  sigma = s1.getVal();
  alpha = s2.getVal();
  if(!pionOnly) {
    mu_eta = e0.getVal();
    sigma_eta = e1.getVal();
    alpha_eta = e2.getVal();
  };
  ch[0] = b0.getVal();
  ch[1] = b1.getVal();
  ch[2] = b2.getVal();
  ch[3] = b3.getVal();
  ch[4] = b4.getVal();





  // determine signal window
  Float_t sig_lb; 
  Float_t sig_ub; 
  Double_t delta = alpha / TMath::Sqrt(1+(alpha*alpha));
  Double_t mean = mu + (sigma*delta);
  Double_t stdev = sigma * TMath::Sqrt(1-(delta*delta));
  Double_t variance = stdev*stdev;
  sig_lb = mean - 3*stdev;
  sig_ub = mean + 3*stdev;
  Double_t sig_lb_bin = hh->FindBin(sig_lb);
  Double_t sig_ub_bin = hh->FindBin(sig_ub);


  // obtain pion purity within signal window
  mass.setRange("int_range",sig_lb,sig_ub);

  RooAbsReal * model_int_roo = model.createIntegral(mass,NormSet(mass),Range("int_range"));
  RooAbsReal * sig_int_roo = sig_skew.createIntegral(mass,NormSet(mass),Range("int_range"));
  RooAbsReal * bg_int_roo =  bg_cheb.createIntegral(mass,NormSet(mass),Range("int_range"));
  
  Double_t model_int = model_int_roo->getVal();
  Double_t sig_int = sig_int_roo->getVal();
  Double_t bg_int = bg_int_roo->getVal();

  printf("model_int -- %.4f\n",model_int);
  printf("sig_int -- %.4f\n",sig_int);
  printf("bg_int -- %.4f\n",bg_int);
  printf("\n");

  pion_purity = frac.getVal() * sig_int / model_int;
  background_frac = (1-frac.getVal()) * bg_int / model_int;
  


  // tlines
  Double_t minim = 1;
  Double_t maxim = 1.2*(hh->GetMaximum());
  TLine * fline_lb = new TLine(fit_lb,minim,fit_lb,maxim);
  TLine * fline_ub = new TLine(fit_ub,minim,fit_ub,maxim);
  fline_lb->SetLineWidth(2); fline_ub->SetLineWidth(2);
  fline_lb->SetLineColor(kMagenta+3); fline_ub->SetLineColor(kMagenta+3);
  TLine * sline_lb = new TLine(sig_lb,0.6*maxim,sig_lb,0.9*maxim);
  TLine * sline_ub = new TLine(sig_ub,0.6*maxim,sig_ub,0.9*maxim);
  TLine * sline_top = new TLine(sig_lb,0.9*maxim,sig_ub,0.9*maxim);
  sline_lb->SetLineWidth(2); sline_ub->SetLineWidth(2);
  sline_top->SetLineWidth(2);
  sline_lb->SetLineColor(kGreen+1); sline_ub->SetLineColor(kGreen+1);
  sline_top->SetLineColor(kGreen+1);


  // tlatex's
  const Int_t MAX_LINES = 20;
  char txt[MAX_LINES][256];
  TLatex * ttxt[MAX_LINES];
  Int_t nt=0;
  sprintf(txt[nt++],"#mu = %.4f GeV/c^{2}",mu);
  sprintf(txt[nt++],"#sigma = %.4f GeV/c^{2}",sigma);
  sprintf(txt[nt++],"#alpha = %.4f",alpha);
  sprintf(txt[nt++],"");
  if(!pionOnly) {
    sprintf(txt[nt++],"#mu(#eta) = %.4f GeV/c^{2}",mu_eta);
    sprintf(txt[nt++],"#sigma(#eta) = %.4f GeV/c^{2}",sigma_eta);
    sprintf(txt[nt++],"#alpha(#eta) = %.4f",alpha_eta);
    sprintf(txt[nt++],"");
  };
  for(int i=0; i<chebOrder; i++) {
    sprintf(txt[nt++],"a_{%d} = %.4f",i,ch[i]);
  };
  sprintf(txt[nt++],"");
  sprintf(txt[nt++],"f = %.4f",frac.getVal());
  sprintf(txt[nt++],"");
  sprintf(txt[nt++],"signal: M_{#gamma#gamma}#in[%.3f, %.3f]",sig_lb,sig_ub);
  sprintf(txt[nt++],"purity: F = %.2f%%  (bg @ %.2f%%)",pion_purity*100,background_frac*100);
  sprintf(txt[nt++],"");
  //sprintf(txt[nt++],"#chi^{2}/NDF=%.4f",plotframe->chiSquare(7));




  // draw canvas
  TCanvas * canv = new TCanvas("canv","canv",1400,800);
  canv->Divide(2,2);
  plotframe->GetXaxis()->SetRangeUser(draw_lb,draw_ub);
  plotframe->GetYaxis()->SetRangeUser(minim,maxim);
  canv->cd(1); 
    canv->GetPad(1)->SetGrid(1,1);
    plotframe->Draw(); fline_lb->Draw(); fline_ub->Draw();
    sline_lb->Draw(); sline_ub->Draw(); sline_top->Draw();
  canv->cd(3); 
    canv->GetPad(3)->SetGrid(1,1);
    canv->GetPad(3)->SetLogy();
    plotframe->Draw(); fline_lb->Draw(); fline_ub->Draw();
  canv->cd(2); 
    canv->GetPad(2)->Range(0,0,1,1);
    Float_t pos=1-0.05;
    for(int k=0; k<nt; k++) {
      ttxt[k] = new TLatex(0.03,pos-=0.06,txt[k]); 
      ttxt[k]->Draw();
    };
};
