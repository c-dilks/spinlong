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
Double_t massL,massM,massH;

// polyOrder is the order of chebychev polynomial for BG
// if polyOrder < 0, we use other trial fit function for BG
// if fitEta = true, attempts to fit skewed gaussian to eta peak as well

void Purity(Int_t eta_bin = 0,
            Int_t en_bin = 0,
            Int_t pt_bin = 0,
            Int_t polyOrder = 5,
            Bool_t fitEta = false,
            TString infile_n = "massset_12/all.root") {
  // load libraries
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;


  // load bin information
  Float_t eta_lb = RD->env->EtaDiv(eta_bin);
  Float_t eta_ub = RD->env->EtaDiv(eta_bin+1);
  Float_t en_lb = RD->env->EnDiv(en_bin);
  Float_t en_ub = RD->env->EnDiv(en_bin+1);
  Float_t pt_lb = RD->env->PtDiv(pt_bin);
  Float_t pt_ub = RD->env->PtDiv(pt_bin+1);


  // open histogram
  TFile * infile = new TFile(infile_n.Data(),"READ");
  TString arrname = Form("mass_dist_arr_g%d_e%d_p%d",eta_bin,en_bin,pt_bin);
  TObjArray * arr = (TObjArray*) infile->Get(arrname.Data());
  TH1D * h = (TH1D*) arr->At(1);

  // execute fit algo
  Int_t hash = 0x100*eta_bin + 0x10*en_bin + pt_bin;
  TCanvas * fitcanv = ExeFit(h,polyOrder,fitEta,hash); 


  // print stats
  printf("\n::::::::::::::::::::::::::\n");
  printf("%.2f <= eta%d < %.2f\n",eta_lb,eta_bin,eta_ub);
  printf("%.2f <= en%d < %.2f\n",en_lb,en_bin,en_ub);
  printf("%.2f <= pt%d < %.2f\n",pt_lb,pt_bin,pt_ub);
  printf("::::::::::::::::::::::::::\n\n");
  printf("#mu = %.4f GeV\n",mu);
  printf("#sigma = %.4f GeV\n",sigma);
  printf("#alpha = %.4f\n",alpha);
  if(fitEta) {
    printf("\n");
    printf("#mu(#eta) = %.4f GeV\n",mu_eta);
    printf("#sigma(#eta) = %.4f GeV\n",sigma_eta);
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
  printf("\n\nbinhash = 0x%03X\n\n",hash);
  printf("\n");

  // output mass cuts
  TString statsfile = Form("%s/stat_g%d_e%d_p%d.dat",RD->env->massset_dir,eta_bin,en_bin,pt_bin);
  gSystem->RedirectOutput(statsfile.Data(),"w");
  printf("%f %f %f %f %f %f %f %f %f %f\n",
    eta_lb, eta_ub,
    en_lb, en_ub,
    pt_lb, pt_ub,
    massL,
    massM,
    massH,
    pion_purity
  );
  gSystem->RedirectOutput(0);


  // print canvas
  TString pdffile = Form("%s/canv_g%d_e%d_p%d.pdf",RD->env->massset_dir,eta_bin,en_bin,pt_bin);
  fitcanv->Print(pdffile.Data(),"pdf");

};



TCanvas * ExeFit(TH1D * hh, 
                 Int_t polyOrder_, 
                 Bool_t fitEta_,
                 Int_t binhash_) {
  // don't fit empty histogram
  if(hh==NULL) { 
    fprintf("ERROR: hh = %p\n",(void*)hh);
    return;
  };


  // control vars
  Int_t chebOrder = polyOrder_;
  Bool_t pionOnly = !fitEta_;
  Int_t binhash = binhash_; // 3 hex-digits 0x[eta_bin][en_bin][pt_bin]
  

  // pion mass, bounds, and ranges for fit and draw
  const Float_t PION_MASS = 0.135;
  const Float_t ETA_MASS = 0.548;
  const Float_t MASS_LB_DEF = 0;
  const Float_t MASS_UB_DEF = 1;
  Float_t fit_lb; // fit range
  Float_t fit_ub; // fit range
  Float_t mass_lb; // mass dist lower bound
  Float_t mass_ub; // mass dist upper bound
  Float_t draw_lb; // drawing range
  Float_t draw_ub; // drawing range
  switch(binhash) {
    case 0x010:
      fit_lb = 0.1;
      fit_ub = 0.6;
      mass_lb = MASS_LB_DEF;
      mass_ub = MASS_UB_DEF;
      draw_lb = 0.05;
      draw_ub = 0.8;
      break;
    default:
      if(pionOnly) {
        fit_lb = 0.1;
        fit_ub = 0.5;
        mass_lb = MASS_LB_DEF;
        mass_ub = MASS_UB_DEF;
        draw_lb = 0.05;
        draw_ub = 0.5;
      } else {
        fit_lb = 0.1;
        fit_ub = 0.7;
        mass_lb = MASS_LB_DEF;
        mass_ub = MASS_UB_DEF;
        draw_lb = 0.05;
        draw_ub = 1;
      };
  };


  // mass distribution
  TString hh_t = TString(hh->GetTitle());
  RooRealVar mass("mass","M_{#gamma#gamma}",PION_MASS,mass_lb,mass_ub,"GeV/c^{2}");
  RooDataHist massdist("massdist",hh_t.Data(),mass,Import(*hh));


  // SIGNAL FIT MODELS :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  // pi0
  switch(binhash) {
    case 0x000:
      RooRealVar s0("s0","pion_mu",PION_MASS,0.1,0.3);
      RooRealVar s1("s1","pion_sigma",0.06,0.01,0.2);
      RooRealVar s2("s2","pion_alpha",3.0,1.0,10.0);
      break;
    case 0x010:
      RooRealVar s0("s0","pion_mu",0.17,0.1,0.3);
      RooRealVar s1("s1","pion_sigma",0.06,0.001,0.2);
      RooRealVar s2("s2","pion_alpha",3.3,1.0,10.0);
      break;
    default:
      RooRealVar s0("s0","pion_mu",PION_MASS,0.1,0.3);
      RooRealVar s1("s1","pion_sigma",0.06,0.01,0.2);
      RooRealVar s2("s2","pion_alpha",3.0,1.0,10.0);
  };
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
  TString formu_bg = "";
  
  // manually entered chebychev polynomials 
  // (for more easily controlling parameter constraints)
  /*
  Float_t b_lb = -50.0; // default upper bound
  Float_t b_ub = 50.0; // default lower bound
  switch(binhash) {
    case 0x000:
      RooRealVar b0("b0","b0", 6.0, b_lb, b_ub); // (not used)
      RooRealVar b1("b1","b1", -20, b_lb, b_ub);
      RooRealVar b2("b2","b2", 20, b_lb, b_ub);
      RooRealVar b3("b3","b3", b_lb, b_ub);
      RooRealVar b4("b4","b4", b_lb, b_ub);
      break;
    default:
      RooRealVar b0("b0","b0", b_lb, b_ub); // (not used)
      RooRealVar b1("b1","b1", b_lb, b_ub);
      RooRealVar b2("b2","b2", b_lb, b_ub);
      RooRealVar b3("b3","b3", b_lb, b_ub);
      RooRealVar b4("b4","b4", b_lb, b_ub);
  };

  TString cheb[5];
  cheb[0]="1"; // 1
  cheb[1]="mass"; // m
  cheb[2]="(2*TMath::Power(mass,2)-1)"; // 2m^2 - 1
  cheb[3]="(4*TMath::Power(mass,3)-3*mass)"; // 4m^3 - 3m
  cheb[4]="(8*TMath::Power(mass,4)-8*TMath::Power(mass,2)+1)"; // 8m^4 - 8m^2 +1
  if(chebOrder>0) {
    // zeroth term has coefficient constrained by other coefficients to
    // ensure that the background is zero at zero mass
    if(chebOrder==3 || chebOrder==4) formu_bg = "b2";
    else if(chebOrder==5) formu_bg = "(b2-b4)";
    else formu_bg = "0";
    formu_bg = "b0"; // OVERRIDE
    for(int co=1; co<chebOrder; co++) {
      formu_bg = Form("%s+b%d*%s",formu_bg.Data(),co,cheb[co].Data());
    };

    switch(chebOrder) {
      case 2:
       RooGenericPdf bg_func("background","bg_func",
                             //formu_bg.Data(),RooArgList(mass,b1));
                             formu_bg.Data(),RooArgList(mass,b0,b1));
       break;
      case 3:
       RooGenericPdf bg_func("background","bg_func",
                             //formu_bg.Data(),RooArgList(mass,b1,b2));
                             formu_bg.Data(),RooArgList(mass,b0,b1,b2));
       break;
      case 4:
       RooGenericPdf bg_func("background","bg_func",
                             //formu_bg.Data(),RooArgList(mass,b1,b2,b3));
                             formu_bg.Data(),RooArgList(mass,b0,b1,b2,b3));
       break;
      case 5:
       RooGenericPdf bg_func("background","bg_func",
                             //formu_bg.Data(),RooArgList(mass,b1,b2,b3,b4));
                             formu_bg.Data(),RooArgList(mass,b0,b1,b2,b3,b4));
       break;
      default:
       return;
    };
  }
  */

  ///*
  // ROOfit chebychev polynomials
  Float_t b_lb = -50.0; // default upper bound
  Float_t b_ub = 50.0; // default lower bound
  switch(binhash) {
    case 0x000:
      RooRealVar b0("b0","b0", 6.0, b_lb, b_ub);
      RooRealVar b1("b1","b1", -20, b_lb, b_ub);
      RooRealVar b2("b2","b2", 20, b_lb, b_ub);
      RooRealVar b3("b3","b3", b_lb, b_ub);
      RooRealVar b4("b4","b4", b_lb, b_ub);
      break;
    default:
      RooRealVar b0("b0","b0", b_lb, b_ub);
      RooRealVar b1("b1","b1", b_lb, b_ub);
      RooRealVar b2("b2","b2", b_lb, b_ub);
      RooRealVar b3("b3","b3", b_lb, b_ub);
      RooRealVar b4("b4","b4", b_lb, b_ub);
  };
  if(chebOrder>0) {
    switch(chebOrder) {
      case 2:
       RooChebychev bg_func("background","bg_func",mass,RooArgSet(b0,b1));
       break;
      case 3:
       RooChebychev bg_func("background","bg_func",mass,RooArgSet(b0,b1,b2));
       break;
      case 4:
       RooChebychev bg_func("background","bg_func",mass,RooArgSet(b0,b1,b2,b3));
       break;
      case 5:
       RooChebychev bg_func("background","bg_func",mass,RooArgSet(b0,b1,b2,b3,b4));
       break;
      default:
       return;
    };
  } 
  //*/
  else {
    // argus function
    RooRealVar arg_m("arg_m","arg_m", 1.2, 0.8, 2);
    RooRealVar arg_c("arg_c","arg_c", -50, -100, 0);
    RooRealVar arg_p("arg_p","arg_p", 40, 0, 100);
    RooArgusBG bg_func("background","bg_func",mass,arg_m,arg_c,arg_p);
  };


  // DEFINE FULL PDF ::::::::::::::::::::::::::::::::::::::::::::::::::::
  switch(binhash) {
    default:
      RooRealVar frac("frac","background fraction",0.7,0.01,1.0);
      RooRealVar frac2("frac2","background fraction 2",0.03,0.0001,1.0); // (for etas)
  };
  if(pionOnly) {
    RooAddPdf model(
      "model",
      "model",
      RooArgList(sig_skew,bg_func),
      frac
    );
  } else {
    RooAddPdf model(
      "model",
      "model",
      RooArgList(sig_skew,sig_skew_eta,bg_func),
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
    Name("p_massdist"),
    DataError(RooAbsData::SumW2),
    MarkerColor(kRed)
  );
  model.plotOn(
    plotframe,
    Name("p_model"),
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
  massL = sig_lb;
  massM = mean;
  massH = sig_ub;


  // obtain pion purity within signal window
  mass.setRange("int_range",sig_lb,sig_ub);

  RooAbsReal * model_int_roo = model.createIntegral(mass,NormSet(mass),Range("int_range"));
  RooAbsReal * sig_int_roo = sig_skew.createIntegral(mass,NormSet(mass),Range("int_range"));
  RooAbsReal * bg_int_roo =  bg_func.createIntegral(mass,NormSet(mass),Range("int_range"));
  if(!pionOnly) 
    RooAbsReal * eta_int_roo = sig_skew_eta.createIntegral(mass,NormSet(mass),Range("int_range"));


  
  Double_t model_int = model_int_roo->getVal();
  Double_t sig_int = sig_int_roo->getVal();
  Double_t bg_int = bg_int_roo->getVal();
  Double_t eta_int=0;
  if(!pionOnly) {
    eta_int = eta_int_roo->getVal();
    bg_int += eta_int;
  };

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
  const Int_t MAX_LINES = 50;
  char txt[MAX_LINES][256];
  TLatex * ttxt[MAX_LINES];
  Int_t nt=0;
  sprintf(txt[nt++],"#mu = %.4f GeV     #mu_{PDF} = %.4f GeV",mu,mean);
  sprintf(txt[nt++],"#sigma = %.4f GeV     #sigma_{PDF} = %.4f GeV",sigma,stdev);
  sprintf(txt[nt++],"#alpha = %.4f",alpha);
  sprintf(txt[nt++],"");
  if(!pionOnly) {
    sprintf(txt[nt++],"#mu(#eta) = %.4f GeV",mu_eta);
    sprintf(txt[nt++],"#sigma(#eta) = %.4f GeV",sigma_eta);
    sprintf(txt[nt++],"#alpha(#eta) = %.4f",alpha_eta);
    sprintf(txt[nt++],"");
  };
  TString units;
  for(int i=0; i<chebOrder; i++) {
    if(i>0) units = Form("GeV^{-%d}",i);
    else units = "";
    sprintf(txt[nt++],"a_{%d} = %.4f %s",i,ch[i],units.Data());
  };
  sprintf(txt[nt++],"");
  sprintf(txt[nt++],"f = %.4f",frac.getVal());
  sprintf(txt[nt++],"");
  sprintf(txt[nt++],"signal: M_{#gamma#gamma}#in[%.3f, %.3f]",sig_lb,sig_ub);
  sprintf(txt[nt++],"");
  sprintf(txt[nt++],"purity: F = %.2f%%  (bg @ %.2f%%)",pion_purity*100,background_frac*100);
  sprintf(txt[nt++],"");
  sprintf(txt[nt++],"#chi^{2}/NDF = %.4f",
    plotframe->chiSquare("p_model","p_massdist"));
  sprintf(txt[nt++],"");
  sprintf(txt[nt++],"binhash = 0x%03X",binhash); 



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
      ttxt[k] = new TLatex(0.03,pos-=0.04,txt[k]); 
      ttxt[k]->Draw();
    };

   printf("\nformu_bg = %s\n",formu_bg.Data());

  return canv;
};
