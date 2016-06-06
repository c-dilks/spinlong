// combines diagsetfiles
//

void add_diag() {
  //Bool_t printPDFs=false; // DEPRECATED; moved to print_diag.C
  gROOT->Reset();

  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * LT = new LevelTwo(RD->env);
  EventClass * ev = new EventClass(RD->env);



  // build array of diagset/*.root TFile pointers
  // and get set-names
  const Int_t MAX_NUM_FILES=200;
  TFile * diag_file[MAX_NUM_FILES]; 
  Int_t diag_file_cnt=0;
  gROOT->ProcessLine(".! ls diagset/diag*.root | sort > toa_files.txt");
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

      // fgets reads in "returns"; this hack gets rid of them
      sscanf(filename[diag_file_cnt],"%s",filename[diag_file_cnt]);
      sscanf(filename[diag_file_cnt],"diagset/diagset%[^.].root",setname[diag_file_cnt]);
      

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
  const Int_t NPLOTS = NPLOTS_tmp;
  Int_t NTRIGS_tmp = LT->N;
  const Int_t NTRIGS = NTRIGS_tmp;
  printf("NCLASSES=%d  NPLOTS=%d  NTRIGS=%d  NSETS=%d  NMASSES=%d\n",
    NCLASSES,NPLOTS,NTRIGS,NSETS,NMASSES);

  /////////////////////////////////////
  


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


  
  ////////////////////////////////
  // DIAGSET FILE LOOP
  ////////////////////////////////
  
  for(int s=0; s<5; s++) {
  //for(int s=0; s<NSETS; s++) {
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


  TFile * outfile = new TFile("diagset/all.root","RECREATE");

  // write trig_dist
  trig_dist_arr->Write("trig_dist_arr",TObject::kSingleKey);

  // write overlap_matrices
  TDirectory * om_dir = outfile->mkdir("overlap_matrices");
  TDirectory * om_classdir[NCLASSES];
  TString om_classname;

  for(int k=0; k<NCLASSES; k++ ) {
    om_dir->cd();
    om_classname = Form("%s_overlap",ev->Name(k));
    om_classdir[k] = om_dir->mkdir(om_classname.Data());
    om_classdir[k]->cd();
    mix_arr[k]->Write(mix_arr_name[k].Data(),TObject::kSingleKey);
  };
  outfile->cd();


  // build directory tree and store Tobjarrays of kinematic correlations plots
  char classdir_n[NCLASSES][128];
  char plotdir_n[NCLASSES][NPLOTS][128];
  TDirectory * classdir[NCLASSES];
  TDirectory * plotdir[NCLASSES][NPLOTS];
  char new_array_name[NCLASSES][NPLOTS][NTRIGS][256];
  char pdfname[1024];
  char pdfnamel[1024];
  char pdfnamer[1024];
  char mkpdfdir[1024];
  TCanvas * canv = new TCanvas("canv","canv",200,200);
  gStyle->SetOptStat(0);
  for(int c=0; c<NCLASSES; c++) {
    sprintf(classdir_n[c],"%s",classname[c][0][0]);
    classdir[c] = outfile->mkdir(classdir_n[c]);
    classdir[c]->cd();

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

  // write en- and pt-dependent mass dists
  TDirectory * massdir = outfile->mkdir("mass");
  TDirectory * masskindir[2*NMASSES];
  TString kinname;
  TString masskindir_name;
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
