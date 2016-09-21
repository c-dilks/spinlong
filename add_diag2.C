// combines diagsetfiles
//

void add_diag2(Bool_t useTightCuts=false,
               Int_t whichClass=0) {
  gROOT->Reset();

  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * LT = new LevelTwo(RD->env);
  EventClass * ev = new EventClass(RD->env,useTightCuts?false:true);



  // build array of diagset/*.root TFile pointers
  // and get set-names
  const Int_t MAX_NUM_FILES=200;
  TFile * diag_file[MAX_NUM_FILES]; 
  Int_t diag_file_cnt=0;
  TString tight_str = useTightCuts ? "_tight":"";
  TString cmd = Form(".! ls %s%s/diag*.root | sort > toa_files.txt",
    RD->env->diagset_dir,tight_str.Data());
  gROOT->ProcessLine(cmd.Data());
  const Int_t filename_buffer=64;
  char filename[MAX_NUM_FILES][filename_buffer];
  char setname[MAX_NUM_FILES][filename_buffer];
  FILE * toa_files;
  toa_files = fopen("toa_files.txt","r");
  if(toa_files==NULL) {
    fprintf(stderr,"Error opening toa_files.txt\n");
    return;
  }
  else {
    while(!feof(toa_files)) {
      fgets(filename[diag_file_cnt],filename_buffer,toa_files);

      // fgets reads in "returns"; this hack gets rid of them (expects format diagset_${year}/diag*.root
      sscanf(filename[diag_file_cnt],"%s",filename[diag_file_cnt]);
      if(useTightCuts) 
        sscanf(filename[diag_file_cnt],"diagset_%*d_tight/diagset%[^.].root",setname[diag_file_cnt]);
      else
        sscanf(filename[diag_file_cnt],"diagset_%*d/diagset%[^.].root",setname[diag_file_cnt]);
      

      if(strcmp(filename[diag_file_cnt],"")) {
        printf("%d: %s -- %s\n",diag_file_cnt,filename[diag_file_cnt],setname[diag_file_cnt]);
        diag_file[diag_file_cnt] = new TFile(filename[diag_file_cnt],"READ");
        diag_file_cnt++;
      };
    };
  };
  const Int_t NSETS=diag_file_cnt;
  gROOT->ProcessLine(".! rm toa_files.txt");



  /////////////////////////////////////
  // DEFINE CONSTANT NUMBERS
  Int_t NCLASSES_tmp = ev->N; 
  const Int_t NCLASSES = NCLASSES_tmp;

  // count number of kinematic correlation plots
  char keyname[256];
  char temp[256];
  char temp2[256];
  Int_t NPLOTS_tmp=0;
  Bool_t in_one = false;
  diag_file[0]->cd();
  TIter nt(gDirectory->GetListOfKeys());
  while(key=(TKey*)nt()) {
    if(gROOT->GetClass(key->GetClassName())->InheritsFrom("TObjArray")) {
      sprintf(keyname,"%s",key->GetName());
      sscanf(keyname,"%[^_]",temp);
      if(!in_one) {
        strcpy(temp2,temp);
        in_one=true;
      }
      else if(strcmp(temp,temp2)) break;
      NPLOTS_tmp++;
    };
  };
  const Int_t NPLOTS = NPLOTS_tmp;


  // count number of types of rdists
  diag_file[0]->cd();
  TString rdist_types[64]; // assume max number of RTYPES
  Int_t rdist_types_cnt=0;
  TString rdist_type;
  TDirectory * rdir = (TDirectory*)diag_file[0]->Get("rdists");
  rdir->cd();
  TString cls_tmp;
  Int_t int_tmp1,int_tmp2;
  char typ_tmp[32];
  char cls_tmp1[32];
  char cls_tmp2[32];
  in_one=false;
  TIter nt1(rdir->GetListOfKeys());
  Int_t NRTYPES_tmp=0;
  while(rdkey=(TKey*)nt1()) {
    cls_tmp = Form("%s",rdkey->GetName());
    cls_tmp = cls_tmp.ReplaceAll("rdist_arr_","");
    cls_tmp = cls_tmp.ReplaceAll("_"," ");
    sscanf(cls_tmp.Data(),"%s %s %d",cls_tmp1,typ_tmp,&int_tmp1);
    rdist_type = Form("%s",typ_tmp);
    rdist_types[rdist_types_cnt] = rdist_type;
    rdist_types_cnt++;
    //printf("%s -- %s -- %s %s %d\n",rdkey->GetName(),cls_tmp.Data(),cls_tmp1,typ_tmp,int_tmp1);
    if(!in_one) {
      strcpy(cls_tmp2,cls_tmp1);
      int_tmp2 = int_tmp1;
      in_one=true;
    }
    else if(strcmp(cls_tmp1,cls_tmp2) || int_tmp1!=int_tmp2) break;
    NRTYPES_tmp++;
  };
  diag_file[0]->cd();
  const Int_t NRTYPES = NRTYPES_tmp;


  // count number of bins in kinematic-dependent massdistributions
  TIter nt2(gDirectory->GetListOfKeys());
  Int_t mass_cnt=0;
  while(key=(TKey*)nt2()) {
    if(gROOT->GetClass(key->GetClassName())->InheritsFrom("TObjArray")) {
      sprintf(keyname,"%s",key->GetName());
      sscanf(keyname,"%[^_]",temp);
      if(!strcmp(temp,"mass")) mass_cnt++;
    };
  };
  mass_cnt/=2;
  const Int_t NMASSES = mass_cnt;
      

  // other constants
  Int_t NTRIGS_tmp = LT->N;
  const Int_t NTRIGS = NTRIGS_tmp;
  printf("NCLASSES=%d  NPLOTS=%d  NTRIGS=%d  NSETS=%d  NMASSES=%d  NRTYPES=%d\n",
    NCLASSES,NPLOTS,NTRIGS,NSETS,NMASSES,NRTYPES);

  /////////////////////////////////////
  


  // define trigger thresholds tree
  // -- sorry, but this tree has a weird and stupid structure... there is a branch
  //    called which_thresh which tells you which thresholds to believe; e.g., if you
  //    want to plot pt threshold vs. run index, you MUST make the cut which_thresh=="pt";
  //    then plot thresh vs. index
  //    It was done this way because the tree-filling loops below were written well before
  //    this tree was even considered, and it's way too difficult to change the loop structure now
  TString outfile_n = Form("%s%s/setdep.root",RD->env->diagset_dir,tight_str.Data());
  TFile * outfile;
  if(whichClass==0) outfile = new TFile(outfile_n.Data(),"RECREATE");
  else              outfile = new TFile(outfile_n.Data(),"UPDATE");
  TTree * threshtr = new TTree("threshtr","threshtr");
  Int_t th_runnum,th_index,th_class,th_trig;
  char which_thresh[32];
  Float_t thresh;
  Float_t thresh_err;
  threshtr->Branch("runnum",&th_runnum,"runnum/I");
  threshtr->Branch("index",&th_index,"index/I");
  threshtr->Branch("class",&th_class,"class/I");
  threshtr->Branch("trig",&th_trig,"trig/I");
  threshtr->Branch("which_thresh",which_thresh,"which_thresh/C");
  threshtr->Branch("thresh",&thresh,"thresh/F");
  threshtr->Branch("thresh_err",&thresh_err,"thresh_err/F");
  TF1 * curr_fit;



  // define object arrays
  TObjArray * trig_dist_arr = new TObjArray();
    TH1D * trig_dist[NSETS];
  TObjArray * gen_array[NCLASSES][NPLOTS][NSETS];
  TObjArray * new_array[NCLASSES][NPLOTS][NTRIGS];
  TObject * gen_object[NCLASSES][NPLOTS][NTRIGS][NSETS];
    char classname[NCLASSES][NPLOTS][NTRIGS][32];
    char plotname[NCLASSES][NPLOTS][NTRIGS][32];
    char trigname[NCLASSES][NPLOTS][NTRIGS][32];
    char classname_tmp[32];
    char plotname_tmp[32];
    char trigname_tmp[32];
    char char_tmp[32];
    Int_t cc,pp,tt;
    Int_t cc_tmp=-1;
    Int_t pp_tmp=-1;
    Int_t tt_tmp=-1;
    TString plotname_tmp_str;
  char newkeyname[256];
  TObjArray * mass_in[2][NMASSES];
  TObjArray * mass_out[2][NMASSES][NTRIGS];
  TH1D * mass_dist[2][NMASSES][NTRIGS][NSETS];
  Int_t nm,nk;
  char kin_type[8];
  char new_mass_name[256];
  TString mass_out_name[2][NMASSES][NTRIGS];
  TH2D * old_mix[NCLASSES][NSETS];
  TString old_mix_n;
  TString new_mix_n;
  TDirectory * chdir;
  TObjArray * mix_arr[NCLASSES];
  TString mix_arr_name[NCLASSES];
  char rclass[32];
  char rtype[32];
  Int_t rint,rint_tmp;
  Int_t nruns_in_this_set[NSETS];
  for(int ss=0; ss<NSETS; ss++) nruns_in_this_set[ss]=0;
  rint_tmp = -1;
  char rtrig[32];
  const Int_t MAX_NRUNS = 50; // max runs in single runset
  TH1D * rdist[NCLASSES][NRTYPES][MAX_NRUNS][NTRIGS][NSETS];
  TString keyname_tmp;
  TObjArray * rarr[NCLASSES][NRTYPES][MAX_NRUNS][NSETS];
  TObjArray * new_rarr[NCLASSES][NRTYPES][NTRIGS];
  char new_rarr_name[NCLASSES][NRTYPES][NTRIGS][256];
  Int_t rc,rtp,rtg;
  Double_t rmean,rrms;
  char chtmp[4][32];
  TString chtmpstr;
  Int_t runnum;

  Int_t binn;
  Float_t maxx,curr;


  
  ////////////////////////////////
  // DIAGSET FILE LOOP
  ////////////////////////////////
  
  //for(int s=0; s<5; s++) {
  for(int s=0; s<NSETS; s++) {
    //printf("\n");
    printf("processing set %s\n",setname[s]);

    // open diagset file
    diag_file[s]->cd();
    TIter next(gDirectory->GetListOfKeys());

    // loop through keys in the file
    while(key=(TKey*)next()) {
      sprintf(keyname,"%s",key->GetName());
      //printf("[+] %s -- %s\n",setname[s],keyname);

      // trig_dist key sector
      if(!strcmp(keyname,"trig_dist")) {
        sprintf(newkeyname,"%s_%s",keyname,setname[s]);
        trig_dist[s] = (TH1D*) key->ReadObj();
        trig_dist[s]->SetName(newkeyname);
        trig_dist_arr->AddLast(trig_dist[s]);
        //printf("   %s --> %s stored in trig_dist_arr.At(\"%d\")\n",
          //keyname,trig_dist_arr->At(s)->GetName(),s);
      }

      
      // overlap_matrices key sector
      else if(!strcmp(keyname,"overlap_matrices")) {
        chdir = (TDirectory*)key->ReadObj();
        chdir->cd();
        for(int oc=0; oc<NCLASSES; oc++) {
          old_mix_n = Form("/overlap_matrices/%s_trig_fms_mix",ev->Name(oc));
          old_mix[oc][s] = (TH2D*) diag_file[s]->Get(old_mix_n.Data());

          if(s==0) {
            mix_arr[oc] = new TObjArray();
            mix_arr_name[oc] = Form("%s",old_mix[oc][s]->GetName());
          };

          new_mix_n = Form("%s_%s",old_mix[oc][s]->GetName(),setname[s]);
          old_mix[oc][s]->SetName(new_mix_n.Data());

          mix_arr[oc]->AddLast(old_mix[oc][s]);
        };
        diag_file[s]->cd();
      }


      // rdists sector
      else if(!strcmp(keyname,"rdists")) {
        chdir = (TDirectory*)key->ReadObj();
        chdir->cd();
        TIter rnt(chdir->GetListOfKeys());
        while(rkey=(TKey*)rnt()) {
          keyname_tmp = Form("%s",rkey->GetName());
          keyname_tmp = keyname_tmp.ReplaceAll("rdist_arr_","");
          keyname_tmp = keyname_tmp.ReplaceAll("_"," ");
          sscanf(keyname_tmp.Data(),"%s %s %d",rclass,rtype,&rint);

          if(rint > nruns_in_this_set[s]) nruns_in_this_set[s] = rint;

          rc = ev->Idx(rclass);
          for(int qq=0; qq<NRTYPES; qq++) {
            if(!strcmp(rtype,rdist_types[qq].Data())) {
              rtp = qq;
              break;
            };
          };

          rarr[rc][rtp][rint][s] = (TObjArray*) rkey->ReadObj();


          for(int ic=0; ic<rarr[rc][rtp][rint][s]->GetEntries(); ic++) {
            sscanf(((TH1D*)(rarr[rc][rtp][rint][s])->At(ic))->GetName(),"%[^_]",rtrig);
            rtg = LT->Index(TString(rtrig));
            rdist[rc][rtp][rint][rtg][s] = (TH1D*)(rarr[rc][rtp][rint][s]->At(ic));

            // fit pt rdist
            if(!strcmp("pt",rdist_types[rtp].Data())) {
              if(rdist[rc][rtp][rint][rtg][s]->GetEntries()>0) {
                rmean = rdist[rc][rtp][rint][rtg][s]->GetMean();
                rrms = rdist[rc][rtp][rint][rtg][s]->GetRMS();
                rdist[rc][rtp][rint][rtg][s]->Fit("gaus","Q","",rmean-rrms,rmean+rrms);
                chtmpstr = Form("%s",rdist[rc][rtp][rint][rtg][s]->GetName());
                chtmpstr = chtmpstr.ReplaceAll("_"," ");
                sscanf(chtmpstr.Data(),"%s %s %s %s %d",chtmp[0],chtmp[1],chtmp[2],chtmp[3],&runnum);

                th_runnum = runnum;
                th_index = RD->Index(runnum);
                th_class = rc;
                th_trig = rtg;
                curr_fit = rdist[rc][rtp][rint][rtg][s]->GetFunction("gaus");
                if(curr_fit!=NULL) {
                  thresh = curr_fit->GetParameter(1);
                  thresh_err = curr_fit->GetParError(1);
                  sprintf(which_thresh,"%s","pt");

                  // pT thresh at 2/3 of peak height
                  // algorithm written for regular non-tight diagset files
                  // if you run it on diagset_tight files, it *should* produce
                  // the same result... but it doesn't really matter since the 
                  // diagset_tight files are not used for determining thresholds
                  binn = rdist[rc][rtp][rint][rtg][s]->FindBin(thresh);
                  maxx = rdist[rc][rtp][rint][rtg][s]->GetBinContent(binn);
                  curr = maxx;
                  while(curr > (0.66*maxx) && binn>0) {
                    binn--;
                    curr = rdist[rc][rtp][rint][rtg][s]->GetBinContent(binn);
                  };
                  thresh = rdist[rc][rtp][rint][rtg][s]->GetBinCenter(binn);
                  if(binn<=1) {
                    //fprintf(stderr,"WARNING: binn<=1 (event class=%d, trigger=%s)\n",rc,rtrig);
                    thresh=0;
                  };

                  threshtr->Fill();
                };
              };
            }

            // fit en rdist
            else if(!strcmp("en",rdist_types[rtp].Data())) {
              if(rdist[rc][rtp][rint][rtg][s]->GetEntries()>0) {
                rmean = rdist[rc][rtp][rint][rtg][s]->GetMean();
                rrms = rdist[rc][rtp][rint][rtg][s]->GetRMS();
                rdist[rc][rtp][rint][rtg][s]->Fit("gaus","Q","",rmean-rrms,rmean+rrms);
                chtmpstr = Form("%s",rdist[rc][rtp][rint][rtg][s]->GetName());
                chtmpstr = chtmpstr.ReplaceAll("_"," ");
                sscanf(chtmpstr.Data(),"%s %s %s %s %d",chtmp[0],chtmp[1],chtmp[2],chtmp[3],&runnum);

                th_runnum = runnum;
                th_index = RD->Index(runnum);
                th_class = rc;
                th_trig = rtg;
                curr_fit = rdist[rc][rtp][rint][rtg][s]->GetFunction("gaus");
                if(curr_fit!=NULL) {
                  thresh = curr_fit->GetParameter(1);
                  thresh_err = curr_fit->GetParError(1);
                  sprintf(which_thresh,"%s","en");
                  threshtr->Fill();
                };
              };
            };




            if(s==0 && rint==0) { 
              new_rarr[rc][rtp][rtg] = new TObjArray();
              sprintf(new_rarr_name[rc][rtp][rtg],"%s_%s_rdist_%s",rtrig,rclass,rdist_types[rtp].Data());
            };

            new_rarr[rc][rtp][rtg]->AddLast(rdist[rc][rtp][rint][rtg][s]);
          };
        };
      }



      // generic TObjArrays sector (for kinematic correlations plots and mass dists)
      else if(gROOT->GetClass(key->GetClassName())->InheritsFrom("TObjArray")) {

        // get class name and plot name
        sscanf(keyname,"%[^_]_%s_arr",classname_tmp,char_tmp);
        plotname_tmp_str = TString(char_tmp).ReplaceAll("_arr","");
        sprintf(plotname_tmp,"%s",plotname_tmp_str.Data());

        // kinematic correlations sub-sector
        if(strcmp(classname_tmp,"mass")) {
          // get class number cc and plot number pp
          cc = ev->Idx(classname_tmp);


          if(cc!=cc_tmp) {
            cc_tmp = cc;
            pp=0;
          }
          else pp++;
          
          if(cc!=whichClass) continue;  // only do one class at a time for kinematic correlations

          //printf("   %s -- %s->%d -- %s->%d\n",keyname,
            //classname_tmp,cc,plotname_tmp,pp);

          // store pointer to tobjarray in "generic" array of pointers to tobjarrays
          gen_array[cc][pp][s] = (TObjArray*) key->ReadObj();

          // lop through this tobjarray
          for(int e=0; e<gen_array[cc][pp][s]->GetEntries(); e++) {

            // obtain trigger name
            sscanf(gen_array[cc][pp][s]->At(e)->GetName(),"%[^_]",trigname_tmp);
            tt = LT->Index(TString(trigname_tmp));
            //printf("      %s->%d\n",trigname_tmp,tt);

            // store class, plot, and trigger names
            strcpy(classname[cc][pp][tt],classname_tmp);
            strcpy(plotname[cc][pp][tt],plotname_tmp);
            strcpy(trigname[cc][pp][tt],trigname_tmp);

            // if we're reading the first file, instantiate new_array
            if(s==0) new_array[cc][pp][tt] = new TObjArray();

            // store current object, give it a new name, and push it into new_array
            sprintf(newkeyname,"%s_%s_%s_set%s",
              classname[cc][pp][tt],
              trigname[cc][pp][tt],
              plotname[cc][pp][tt],
              setname[s]
            );
            gen_object[cc][pp][tt][s] = gen_array[cc][pp][s]->At(e)->Clone(newkeyname);
            //printf("gen_obj @ %s\n",gen_object[cc][pp][tt][s]->GetName());
            //printf("new_array @ %p\n",new_array[cc][pp][tt]);
            new_array[cc][pp][tt]->AddLast(gen_object[cc][pp][tt][s]);
          };
        }

        // mass distributions sub-sector
        else {
          sscanf(key->GetName(),"mass_dist_for_%2sbin_%d_arr",kin_type,&nm);
          nk = !strcmp(kin_type,"en") ? 0:1;
          mass_in[nk][nm] = (TObjArray*)key->ReadObj();
          for(int xx=0; xx<mass_in[nk][nm]->GetEntries(); xx++) {
            mass_dist[nk][nm][xx][s] = (TH1D*)(mass_in[nk][nm]->At(xx));
            if(s==0) {
              mass_out[nk][nm][xx] = new TObjArray();
              mass_out_name[nk][nm][xx] = Form("%s",mass_dist[nk][nm][xx][s]->GetName());
            };
            sprintf(new_mass_name,"%s_%s",mass_out_name[nk][nm][xx].Data(),setname[s]);
            mass_dist[nk][nm][xx][s]->SetName(new_mass_name);
            mass_out[nk][nm][xx]->AddLast(mass_dist[nk][nm][xx][s]);
          };
        };
      };
    };
  };



  // write trig_dist (only done if whichClass==0) 
  outfile->cd();
  if(whichClass==0) trig_dist_arr->Write("trig_dist_arr",TObject::kSingleKey);

  // write overlap_matrices (only done if whichClass==0)
  TDirectory * om_dir;
  TDirectory * om_classdir[NCLASSES];
  TString om_classname;

  if(whichClass==0) {
    om_dir = outfile->mkdir("overlap_matrices");
    for(int k=0; k<NCLASSES; k++ ) {
      om_dir->cd();
      om_classname = Form("%s_overlap",ev->Name(k));
      om_classdir[k] = om_dir->mkdir(om_classname.Data());
      om_classdir[k]->cd();
      mix_arr[k]->Write(mix_arr_name[k].Data(),TObject::kSingleKey);
    };
    outfile->cd();
  };


  // write en- and pt-dependent mass dists (only if whichClass==0)
  TDirectory * massdir;
  TDirectory * masskindir[2*NMASSES];
  TString kinname;
  TString masskindir_name;


  if(whichClass==0) {
    massdir = outfile->mkdir("mass");
    massdir->cd();

    for(int k=0; k<2; k++ ) {
      kinname = (k==0) ? "en":"pt";
      for(int m=0; m<NMASSES; m++) {
        masskindir_name = Form("%s_%d",kinname.Data(),m);
        masskindir[m+k*NMASSES] = massdir->mkdir(masskindir_name.Data());
        masskindir[m+k*NMASSES]->cd();
        for(int t=0; t<NTRIGS; t++) 
          mass_out[k][m][t]->Write(mass_out_name[k][m][t].Data(),TObject::kSingleKey);
        massdir->cd();
      };
    };
    outfile->cd();
  };

  // write threshtr (only if whichClass==0) 
  if(whichClass==0) threshtr->Write();



  // build directory tree and store Tobjarrays of kinematic correlations plots
  char classdir_n[NCLASSES][128];
  char plotdir_n[NCLASSES][NPLOTS][128];
  TDirectory * classdir[NCLASSES];
  TDirectory * plotdir[NCLASSES][NPLOTS];
  TDirectory * rtypedir;
  char rtypedir_name[256];
  char new_array_name[NCLASSES][NPLOTS][NTRIGS][256];
  char pdfname[1024];
  char pdfnamel[1024];
  char pdfnamer[1024];
  char mkpdfdir[1024];
  TCanvas * canv = new TCanvas("canv","canv",200,200);
  gStyle->SetOptStat(0);
  for(int c=0; c<NCLASSES; c++) {
    if(c!=whichClass) continue;
    sprintf(classdir_n[c],"%s",classname[c][0][0]);
    classdir[c] = outfile->mkdir(classdir_n[c]);
    classdir[c]->cd();

    // first write rdists
    for(int p=0; p<NRTYPES; p++) {
      sprintf(rtypedir_name,"%s_rdist_%s",classname[c][0][0],rdist_types[p].Data());
      rtypedir = classdir[c]->mkdir(rtypedir_name);
      rtypedir->cd();
      for(int t=0; t<NTRIGS; t++) {
        new_rarr[c][p][t]->Write(new_rarr_name[c][p][t],TObject::kSingleKey);
      };
      classdir[c]->cd();
    };



    for(int p=0; p<NPLOTS; p++ ) {
      sprintf(plotdir_n[c][p],"%s_%s",classdir_n[c],plotname[c][p][0]);
      plotdir[c][p] = classdir[c]->mkdir(plotdir_n[c][p]);
      plotdir[c][p]->cd();

      for(int t=0; t<NTRIGS; t++) {
        sprintf(new_array_name[c][p][t],"%s_%s_%s",
          classname[c][p][t],
          plotname[c][p][t],
          trigname[c][p][t]);
        new_array[c][p][t]->Write(new_array_name[c][p][t],TObject::kSingleKey);

      };
      classdir[c]->cd();
    };
    outfile->cd();
  };

};
