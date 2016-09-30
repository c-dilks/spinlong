// prints PDFs from diagset/setdep.root

TObjArray * arr;

void print_diag(Bool_t useTightCuts=false,
                Int_t whichClass==0) {
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  EventClass * ev = new EventClass(RD->env,useTightCuts?false:true);

  TString tightstr = useTightCuts ? "_tight":"";
  TString filename = Form("%s%s/setdep.root",RD->env->diagset_dir,tightstr.Data());
  TFile * infile = new TFile(filename.Data(),"READ");

  TString yearstr = Form("%d",RD->env->year);
  TString dirstr = "pdf_kincorr_"+yearstr+tightstr;

  TString rmstr = ".! rm -rv "+dirstr;


  if(whichClass==0) gROOT->ProcessLine(rmstr.Data());
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  //TIter class_it,plot_it,trig_it;
  TIter class_it(gDirectory->GetListOfKeys());
  TString curr_dir_name,class_name;
  TDirectory * curr_dir;
  TString pdf_name,mkdir_str;
  TCanvas * canv = new TCanvas("canv","canv",200,200);
  Int_t class_idx;
  while(class_key=(TKey*)class_it()) {
    if(gROOT->GetClass(class_key->GetClassName())->InheritsFrom("TDirectory")) {
      curr_dir = (TDirectory*)class_key->ReadObj();
      class_name = Form("%s",curr_dir->GetName());
      class_idx = ev->Idx(class_name.Data());

      // only execute for class = whichClass; but.. sometimes class_name is not
      // exactly a class (e.g., "overlap_matrices" or "mass"); in that case, 
      // ev->Idx(class_name) will return -1; we do want PDFs for these histograms too,
      // so we print them out along with the class==0 pass
      if( class_idx==whichClass ||
         (class_idx==-1 && whichClass==0)) {


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
                pdf_name = dirstr+"/"+curr_dir_name+"/"+TString(trig_key->GetName())+".pdf";
                mkdir_str = ".! mkdir -pv "+dirstr+"/"+curr_dir_name;
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
        }; // eo key loop
      }; // eo whichClass cut
    };
  }; // eo main key loop
};
