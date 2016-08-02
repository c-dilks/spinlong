void PeakFinder(Double_t sigma_ = 1,
                Option_t * option_ = "",
                Double_t threshold_ = 0.25,
                TString infile_n="diagset_12_tight/setdep.root",
                TString class_n="pi0",
                TString trig_n="FMSOR",
                TString plot_type="y_vs_x") 
{

  // open array of histos
  TFile * infile = new TFile(infile_n.Data(),"READ");
  TString dirname = class_n+"/"+class_n+"_"+plot_type+"/";
  TString arrname = dirname+class_n+"_"+plot_type+"_"+trig_n;
  TObjArray * arr = (TObjArray*) infile->Get(arrname.Data());
  if(arr==NULL) { fprintf(stderr,"ERROR: %s does not exist in %s\n",arrname.Data(),infile_n.Data()); return; };

  // define outputs
  TString pdfname = "peaks.pdf";
  TString outfile_n = "peaks.root";
  TFile * outfile = new TFile(outfile_n.Data(),"RECREATE");
  gStyle->SetOptStat(0);

  // define looping vars
  Int_t i;
  TH2D * hh;
  TString hh_t,hh_n;
  TRegexp re("^.*set");
  TCanvas * canv = new TCanvas("canv","canv",2000,800);
  Bool_t first = true;
  TSpectrum2 * sp = new TSpectrum2();
  Int_t npeaks;

  // loop over histos
  for(i=0; i<arr->GetEntries(); i++) {
    hh = (TH2D*)(arr->At(i));
    if(hh!=NULL) {
      // execute peak finding algo
      npeaks = sp->Search(hh,sigma_,option_,threshold_);

      // change plot title 
      hh_t = Form("%s -- %d peaks",hh->GetTitle(),npeaks);
      hh->SetTitle(hh_t.Data());

      // draw histos + peaks to canvas
      canv->Clear();
      canv->Divide(2,1);
      for(int k=1;k<3;k++) canv->GetPad(k)->SetGrid(1,1);
      canv->GetPad(2)->SetLogz();
      canv->cd(1); hh->Draw("colz"); sp->Draw("same");
      canv->cd(2); hh->Draw("colz"); sp->Draw("same");

      // add canvas to pdf
      if(first) {
        canv->Print((pdfname+"(").Data(),"pdf");
        first = false;
      }
      else canv->Print(pdfname.Data(),"pdf");

      // add histo to root file
      hh_n = TString(hh->GetName());
      hh_n(re) = "histo";
      hh->Write(hh_n.Data());
    };
  };

  // close pdf
  canv->Print((pdfname+")").Data(),"pdf");
  printf("\n%s printed\n%s built\n",pdfname.Data(),outfile_n.Data());
};


