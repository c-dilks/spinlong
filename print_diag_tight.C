// prints PDFs from diagset_tight/setdep.root

TObjArray * arr;

void print_diag_tight() {
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;

  TString filename = Form("%s_tight/setdep.root",RD->env->diagset_dir);
  TFile * infile = new TFile(filename.Data(),"READ");

  gROOT->ProcessLine(".! rm -rv pdf_kincorr_tight");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  //TIter class_it,plot_it,trig_it;
  TIter class_it(gDirectory->GetListOfKeys());
  TString curr_dir_name,class_name;
  TDirectory * curr_dir;
  TString pdf_name,mkdir_str;
  TCanvas * canv = new TCanvas("canv","canv",200,200);
  while(class_key=(TKey*)class_it()) {
    if(gROOT->GetClass(class_key->GetClassName())->InheritsFrom("TDirectory")) {
      curr_dir = (TDirectory*)class_key->ReadObj();
      class_name = Form("%s",curr_dir->GetName());
      //printf("class_name=%s\n",class_name.Data());
      curr_dir->cd();
      //TIter plot_it(((TDirectory*)class_key)->GetListOfKeys());
      TIter plot_it(curr_dir->GetListOfKeys());
      while(plot_key=(TKey*)plot_it()) {
        if(gROOT->GetClass(plot_key->GetClassName())->InheritsFrom("TDirectory")) {
          curr_dir = (TDirectory*)plot_key->ReadObj();
          curr_dir_name = Form("%s/%s",class_name.Data(),curr_dir->GetName());
          //printf("curr_dir_name=%s\n",curr_dir_name.Data());
          curr_dir->cd();
          TIter trig_it(curr_dir->GetListOfKeys());
          while(trig_key=(TKey*)trig_it()) {
            if(gROOT->GetClass(trig_key->GetClassName())->InheritsFrom("TObjArray")) {
              arr = (TObjArray*)trig_key->ReadObj();
              pdf_name = "pdf_kincorr_tight/"+curr_dir_name+"/"+TString(trig_key->GetName())+".pdf";
              mkdir_str = ".! mkdir -pv pdf_kincorr_tight/"+curr_dir_name;
              //printf("execute: %s\n",mkdir_str.Data());
              gROOT->ProcessLine(mkdir_str.Data());
              for(int x=0; x<arr->GetEntries(); x++) {
                canv->Clear();
                //canv->SetLogz();
                if(!strcmp("TH1D",arr->At(x)->ClassName()))
                  ((TH1D*)(arr->At(x)))->Draw();
                else if(!strcmp("TH2D",arr->At(x)->ClassName()))
                  ((TH2D*)(arr->At(x)))->Draw("colz");
                if(x==0) canv->Print((pdf_name+"(").Data(),"pdf");
                else canv->Print(pdf_name.Data(),"pdf");
              };
              canv->Clear();
              canv->Print((pdf_name+")").Data(),"pdf");
            };
          };
        };
      };
    };
  };
};
