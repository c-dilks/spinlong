#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooArgusBG.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"

using namespace RooFit;

void PurityTest(TString infile_n = "massset_12/all.root") {
  // open histogram
  TFile * infile = new TFile(infile_n.Data(),"READ");
  TObjArray * arr = (TObjArray*) infile->Get("mass_dist_arr_g0_e0_p0");
  TH1D * hh = (TH1D*) arr->At(1);

  // pion mass and bounds
  const Float_t PION_MASS = 0.135;
  const Float_t MASS_LB = 0;
  const Float_t MASS_UB = 1;

  // mass distribution
  RooRealVar mass("mass","invariant mass",PION_MASS,MASS_LB,MASS_UB,"GeV/c^2");
  RooDataHist massdist("massdist","invariant mass distribution",mass,Import(*hh));

  // signal fit function -- gaussian
  RooRealVar sig_mean("sig_mean","sig_mean",PION_MASS,MASS_LB,MASS_UB);
  RooRealVar sig_sigma("sig_sigma","sig_sigma",0.05,MASS_LB,MASS_UB);
  RooGaussian sig_gauss("sig_gauss","sig_gauss",mass,sig_mean,sig_sigma);

  // signal fit function -- crystal ball
  RooRealVar sig_alpha("sig_alpha","sig_alpha",-1,-3,0);
  RooRealVar sig_n("sig_n","sig_n",2,0,3);
  RooCBShape sig_cb("sig_cb","sig_cb",mass,sig_mean,sig_sigma,sig_alpha,sig_n);

  // background fit function -- chebychev
  RooRealVar bg_c0("c0","c0", 0.5, -1.0, 1.0);
  RooRealVar bg_c1("c1","c1", 0.5, -1.0, 1.0);
  RooRealVar bg_c2("c2","c2", 0.5, -1.0, 1.0);
  RooRealVar bg_c3("c3","c3", 0.5, -1.0, 1.0);
  RooRealVar bg_c4("c4","c4", 0.5, -1.0, 1.0);
  RooChebychev bg_poly2("bg_poly","bg_poly",mass,RooArgSet(bg_c0,bg_c1));
  RooChebychev bg_poly3("bg_poly","bg_poly",mass,RooArgSet(bg_c0,bg_c1,bg_c2));
  RooChebychev bg_poly4("bg_poly","bg_poly",mass,RooArgSet(bg_c0,bg_c1,bg_c2,bg_c3));
  RooChebychev bg_poly5("bg_poly","bg_poly",mass,RooArgSet(bg_c0,bg_c1,bg_c2,bg_c3,bg_c4));

  // background fit function -- argus
  RooRealVar arg_m("arg_m","arg_m", 0.7, 0.5, 0.9);
  RooRealVar arg_c("arg_c","arg_c", -10, -50, 50);
  RooRealVar arg_p("arg_p","arg_p", 1, 0, 5);
  RooArgusBG bg_argus("bg_argus","bg_argus",mass,arg_m,arg_c,arg_p);


  // signal + background
  RooRealVar frac("frac","background fraction",0.5,0.0,1.0);
  //RooAddPdf model("model","model",RooArgList(bg_poly,sig_gauss),frac);
  RooAddPdf model("model","model",RooArgList(bg_argus,sig_cb),frac);

  // perform fit
  model.fitTo(massdist);

  // plot
  RooPlot * plotframe = mass.frame(Title("invariant mass distribution with fit"));
  massdist.plotOn(plotframe,DataError(RooAbsData::SumW2));
  model.plotOn(plotframe);
  TCanvas * canv = new TCanvas("canv","canv",800,800);
  plotframe->Draw();
};


  
