// looks through diagset_*/setdep files for rdists, and plots the kinematics 
// vs. run index, in order to look for overall trends (such as the increase in
// pT threshold due to the radiation damage)
//

void KinVarVsRunFast(Bool_t useTightCuts=false) {
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * LT = new LevelTwo(RD->env);
  EventClass * ev = new EventClass(RD->env,useTightCuts?false:true);

  // open file
  TString filename;
  if(!useTightCuts) filename = Form("%s/setdep.root",       RD->env->diagset_dir);
  else              filename = Form("%s_tight/setdep.root", RD->env->diagset_dir);
  TFile * infile = new TFile(filename.Data(),"READ");


  // list of distributions
  const Int_t N_KIN = 6;
  TString kin_name[N_KIN] = { "pt",
                              "en",
                              "eta",
                              "phi",
                              "z",
                              "mass"
                            };

  // counter maxes
  Int_t N_TRIG_tmp = LT->N;
  const Int_t N_TRIG = N_TRIG_tmp;
  Int_t N_CLASS_tmp = ev->N;
  const Int_t N_CLASS = N_CLASS_tmp;


  // read threshtr
  TTree * tr = (TTree*) infile->Get("threshtr");
  Int_t runnum,index,classnum,trignum;
  Float_t thresh,thresh_err;
  char which_thresh[16];
  tr->SetBranchAddress("runnum",&runnum);
  tr->SetBranchAddress("index",&index);
  tr->SetBranchAddress("class",&classnum);
  tr->SetBranchAddress("trig",&trignum);
  tr->SetBranchAddress("thresh",&thresh);
  tr->SetBranchAddress("thresh_err",&thresh_err);
  tr->SetBranchAddress("which_thresh",which_thresh);
  Float_t thresh_arr[N_CLASS][2][N_TRIG][2000]; // [class] [pt/en] [trig] [assumed max run index]
  TGraph * thresh_gr[N_CLASS][2][N_TRIG]; // [class] [pt/en] [trig]

  Int_t kswitch;
  printf("now filling thresh_arr\n");
  for(Int_t x=0; x<tr->GetEntries(); x++) {
    tr->GetEntry(x); 
    if(!strcmp(which_thresh,"pt")) kswitch=0;
    else if(!strcmp(which_thresh,"en")) kswitch=1;
    else kswitch=-1;

    if(kswitch>0) 
      thresh_arr[classnum][kswitch][trignum][index] = thresh;
  };
  printf("done\n");


  // define vars & pointers
  TH1D * khist;
  TH2D * kvr[N_CLASS][N_KIN][N_TRIG]; // [class]  [kin]  [trigger]
  TString kvr_n[N_CLASS][N_KIN][N_TRIG];
  TString kvr_t[N_CLASS][N_KIN][N_TRIG];
  TObjArray * arr;
  TString tdir_n,arr_n;
  TDirectory * tdir;
  Int_t nbins[N_CLASS][N_KIN][N_TRIG];
  Double_t binlow[N_CLASS][N_KIN][N_TRIG];
  Double_t binhigh[N_CLASS][N_KIN][N_TRIG];
  Double_t integral;
  Double_t bcont,bcent;
  Int_t binn;
  Int_t nruns;
  Int_t count=0;
  Int_t runnum_curr;
  TString title_tmp;


  // main loop
  for(Int_t c=0; c<N_CLASS; c++) {
    for(Int_t k=0; k<N_KIN; k++) {

      tdir_n = Form("/%s/%s_rdist_%s",ev->Name(c),ev->Name(c),kin_name[k].Data());
      tdir = (TDirectory*) infile->Get(tdir_n.Data());
      tdir->cd();

      for(Int_t t=0; t<N_TRIG; t++) {
        printf("c k t = %d %d %d\n",c,k,t);

        arr_n = Form("%s_%s_rdist_%s",(LT->Name(t)).Data(),ev->Name(c),kin_name[k].Data());
        cout << arr_n<<endl;
        arr = (TObjArray*) tdir->Get(arr_n.Data());

        nruns = arr->GetEntries();

        for(Int_t r=0; r<nruns; r++) {
          khist = (TH1D*) arr->At(r);

          title_tmp = Form("%s",khist->GetTitle());
          sscanf(title_tmp(title_tmp.Length()-8,8).Data(),"%d",&runnum_curr);


          // define 2d histogram, if it's the first run
          if(r==0) {
            kvr_n[c][k][t] = Form("%s_%s_%s_vs_run",ev->Name(c),kin_name[k].Data(),(LT->Name(t)).Data());
            kvr_t[c][k][t] = Form("%s_%s_%s_vs_run",ev->Name(c),kin_name[k].Data(),(LT->Name(t)).Data());
            nbins[c][k][t] = khist->GetNbinsX();
            binlow[c][k][t] = khist->GetXaxis()->GetXmin();
            binhigh[c][k][t] = khist->GetXaxis()->GetXmax();
            kvr[c][k][t] = new TH2D(kvr_n[c][k][t].Data(),
                                    kvr_t[c][k][t].Data(),
                                    nruns,0,nruns,
                                    nbins[c][k][t],
                                    binlow[c][k][t],
                                    binhigh[c][k][t]
                                   );
            if(k==0 || k==1) thresh_gr[c][k][t] = new TGraph();
          };

          if(khist->GetEntries() > 0) {


            integral = khist->Integral();
            khist->Scale(1/integral);

            for(Int_t b=1; b<nbins[c][k][t]; b++) {
              bcont = khist->GetBinContent(b);
              bcent = khist->GetBinCenter(b);
              
              binn = kvr[c][k][t]->FindBin(count,bcent);
              kvr[c][k][t]->SetBinContent(binn,bcont);
            };

            if(k==0 || k==1) {
              thresh_gr[c][k][t]->SetPoint(count,count,thresh_arr[c][k][t][RD->Index(runnum_curr)]);
            };

            count++;

          };
        }; // eo for r
      }; // eo for t
    }; // eo for k
  }; // eo for c


  // write output
  TString outfile_n = Form("%s/kinvarvsrun.root",RD->env->diagset_dir);
  TFile * outfile = new TFile(outfile_n.Data(),"RECREATE");
  for(Int_t c=0; c<N_CLASS; c++) {
    for(Int_t k=0; k<N_KIN; k++) {
      for(Int_t t=0; t<N_TRIG; t++) {
        if(kvr[c][k][t]!=NULL) kvr[c][k][t]->Write();
      };
    };
  };

  new TCanvas();
  kvr[1][0][1]->Draw("colz");
  thresh_gr[1][0][1]->Draw("L");
};
