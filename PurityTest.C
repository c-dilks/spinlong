#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
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

void PurityTest(TString infile_n = "massset_12/all.root") {
  // open histogram
  TFile * infile = new TFile(infile_n.Data(),"READ");
  TObjArray * arr = (TObjArray*) infile->Get("mass_dist_arr_g0_e0_p1");
  TH1D * h = (TH1D*) arr->At(1);
  ExeFit(h);
};

void ExeFit(TH1D * hh) {
  if(hh==NULL) { 
    fprintf("ERROR: hh = %p\n",(void*)hh);
    return;
  };

  // pion mass and bounds
  const Float_t PION_MASS = 0.135;
  const Float_t MASS_LB = 0;
  const Float_t MASS_UB = 1;

  // mass distribution
  TString hh_t = TString(hh->GetTitle());
  RooRealVar mass("mass","M_{#gamma#gamma}",PION_MASS,MASS_LB,MASS_UB,"GeV/c^{2}");
  RooDataHist massdist("massdist",hh_t.Data(),mass,Import(*hh));

  // set ranges
  mass.setRange("signal_range",0.1,0.4); // fit range



  // signal fit function -- gaussian
  RooRealVar sig_mean("sig_mean","sig_mean",PION_MASS,MASS_LB,MASS_UB);
  RooRealVar sig_sigma("sig_sigma","sig_sigma",0.05,MASS_LB,MASS_UB);
  RooGaussian sig_gauss("signal","sig_gauss",mass,sig_mean,sig_sigma);

  // signal fit function -- crystal ball
  RooRealVar sig_alpha("sig_alpha","sig_alpha",-1,-3,0);
  RooRealVar sig_n("sig_n","sig_n",2,0,3);
  RooCBShape sig_cb("signal","sig_cb",mass,sig_mean,sig_sigma,sig_alpha,sig_n);

  // signal fit function -- landau
  RooRealVar lan_mean("lan_mean","lan_mean",PION_MASS,MASS_LB,MASS_UB);
  RooRealVar lan_sigma("lan_sigma","lan_sigma",0.05,MASS_LB,MASS_UB);
  RooLandau sig_landau("signal","sig_landau",mass,lan_mean,lan_sigma);

  // signal fit function -- breit-wigner
  RooRealVar bre_mean("bre_mean","bre_mean",PION_MASS,MASS_LB,MASS_UB);
  RooRealVar bre_sigma("bre_sigma","bre_sigma",0.05,MASS_LB,MASS_UB);
  RooBreitWigner sig_bre("signal","sig_bredau",mass,bre_mean,bre_sigma);

  // signal fit function -- voigtian
  RooRealVar voi_mean("voi_mean","voi_mean",PION_MASS,MASS_LB,MASS_UB);
  RooRealVar voi_sigma("voi_sigma","voi_sigma",0.05,MASS_LB,MASS_UB);
  RooRealVar voi_sigma2("voi_sigma2","voi_sigma2",0.05,MASS_LB,MASS_UB);
  RooVoigtian sig_voi("signal","sig_voigt",mass,voi_mean,voi_sigma,voi_sigma2);

  // signal fit function -- skewed gaussian
  /*
  RooRealVar skew_mean("skew_mean","skew_mean",PION_MASS,MASS_LB,MASS_UB);
  RooRealVar skew_sigma("skew_sigma","skew_sigma",0.05,MASS_LB,MASS_UB);
  RooRealVar skew_alpha("skew_alpha","skew_alpha",0,100);
  RooGenericPdf sig_skew("signal","sig_skew",
    "TMath::Exp(-0.5*TMath::Power((mass-skew_mean)/skew_sigma,2))*(1+TMath::Erf(skew_alpha*mass/TMath::Sqrt(2)))",
    RooArgList(mass,skew_mean,skew_sigma,skew_alpha));
  */
  RooRealVar sp0("sp0","sp0",PION_MASS,MASS_LB,MASS_UB);
  RooRealVar sp1("sp1","sp1",0.05,0.001,1);
  RooRealVar sp2("sp2","sp2",0.0001,100);
  TString formu_gaus = "1/sp1*TMath::Gaus((mass-sp0)/sp1,0,1)";
  TString formu_skew = "(1+TMath::Erf((1/TMath::Sqrt(2))*sp2*((mass-sp0)/sp1)))";
  TString formu = formu_gaus+"*"+formu_skew;
  RooGenericPdf sig_skew("signal","skewed guassian",formu.Data(),
    RooArgList(mass,sp0,sp1,sp2));

  // background fit function -- exponential
  RooRealVar exp_a("exp_a","exp_a",-10,10);
  RooRealVar exp_b("exp_b","exp_b",-10,10);
  RooGenericPdf back_exp("background","back_exp",
    "TMath::Power(mass,exp_a)*TMath::Exp(exp_b*mass)",
    RooArgList(mass,exp_a,exp_b));


  // background fit function -- chebychev
  RooRealVar bg_c0("c0","c0", -20.0, 20.0);
  RooRealVar bg_c1("c1","c1", -20.0, 20.0);
  RooRealVar bg_c2("c2","c2", -20.0, 20.0);
  RooRealVar bg_c3("c3","c3", -20.0, 20.0);
  RooRealVar bg_c4("c4","c4", -20.0, 20.0);
  RooChebychev bg_poly2("background","bg_poly",mass,RooArgSet(bg_c0,bg_c1));
  RooChebychev bg_poly3("background","bg_poly",mass,RooArgSet(bg_c0,bg_c1,bg_c2));
  RooChebychev bg_poly4("background","bg_poly",mass,RooArgSet(bg_c0,bg_c1,bg_c2,bg_c3));
  RooChebychev bg_poly5("background","bg_poly",mass,RooArgSet(bg_c0,bg_c1,bg_c2,bg_c3,bg_c4));

  // background fit function -- argus
  RooRealVar arg_m("arg_m","arg_m", 0.875, 0.8, 1);
  RooRealVar arg_c("arg_c","arg_c", -10, -50, 50);
  RooRealVar arg_p("arg_p","arg_p", 20, 0, 50);
  RooArgusBG bg_argus("background","bg_argus",mass,arg_m,arg_c,arg_p);


  // signal + background //////////////////////////////////////////////////
  RooRealVar frac("frac","background fraction",0.5,0.0,1.0);
  //RooAddPdf model("model","model",RooArgList(bg_poly3,sig_skew),frac); // 0.1-0.4
  RooAddPdf model("model","model",RooArgList(bg_poly3,sig_skew),frac);
  /////////////////////////////////////////////////////////////////////////

  // perform fit
  model.fitTo(massdist,Range("signal_range"));

  // plot
  RooPlot * plotframe = mass.frame(Title("invariant mass distribution with fit"));
  massdist.plotOn(plotframe,DataError(RooAbsData::SumW2));
  model.plotOn(plotframe,Range("Full"),LineColor(kGreen+1));
  model.plotOn(plotframe,Range("Full"),Components("background"),LineStyle(kDashed));
  model.plotOn(plotframe,Range("Full"),Components("signal"),LineStyle(kDashed),LineColor(kRed));
  printf("\n\n");
  model.Print("t");
  printf("\n\n");


  TCanvas * canv = new TCanvas("canv","canv",800,800);
  canv->Divide(1,2);
  canv->GetPad(2)->SetLogy();
  Double_t maxim = hh->GetMaximum();
  plotframe->GetXaxis()->SetRangeUser(0.05,0.5);
  plotframe->GetYaxis()->SetRangeUser(1e2,1.2*maxim);
  canv->cd(1); plotframe->Draw();
  canv->cd(2); plotframe->Draw();
};


  
