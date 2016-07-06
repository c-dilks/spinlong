// TObjArray Adder
// 
// combines and organises phi distributions from TObjArrays
// of phiset files
//
// -- reads rtree.root


void toa_add(Bool_t printPDFs=false, Int_t FILTER=0)
{
  gROOT->Reset();

  
  // get bins from environment
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  EventClass * ev = new EventClass(RD->env);
  Int_t phi_bins0 = RD->env->PhiBins; const Int_t phi_bins = phi_bins0;
  Int_t eta_bins0 = RD->env->EtaBins; const Int_t eta_bins = eta_bins0;
  Int_t en_bins0 = RD->env->EnBins; const Int_t en_bins = en_bins0;
  Int_t pt_bins0 = RD->env->PtBins; const Int_t pt_bins = pt_bins0;

  // event classes
  Int_t N_CLASS_tmp = ev->N;
  const Int_t N_CLASS = N_CLASS_tmp;


  // build array of $PHISET_DIR/*.root TFile pointers
  const Int_t MAX_NUM_FILES=200;
  TFile * phi_file[MAX_NUM_FILES]; 
  Int_t phi_file_cnt=0;
  char lscommand[1024];
  sprintf(lscommand,".! ls %s/phi*.root > toa_files.txt",RD->env->phiset_dir);
  gROOT->ProcessLine(lscommand);
  const Int_t filename_buffer=64;
  char filename[MAX_NUM_FILES][filename_buffer];
  char temp[filename_buffer];
  FILE * toa_files;
  toa_files = fopen("toa_files.txt","r");
  if(toa_files==NULL)
  {
    fprintf(stderr,"Error opening toa_files.txt\n");
    return;
  }
  else
  {
    while(!feof(toa_files))
    {
      fgets(filename[phi_file_cnt],filename_buffer,toa_files);

      // fgets reads in "returns"; this hack gets rid of them
      sscanf(filename[phi_file_cnt],"%s",filename[phi_file_cnt]);

      if(strcmp(filename[phi_file_cnt],""))
      {
        printf("%d: %s\n",phi_file_cnt,filename[phi_file_cnt]);
        phi_file[phi_file_cnt] = new TFile(filename[phi_file_cnt],"READ");
        phi_file_cnt++;
      };
    };
  };
  const Int_t NFILES=phi_file_cnt;
  gROOT->ProcessLine(".! rm toa_files.txt");



  // build run number hash table
  TFile * rtree_file = new TFile(RD->env->rtree_file,"READ");
  TTree * rtree = (TTree*) rtree_file->Get("rellum");
  Int_t index,runnum;
  Int_t rtree_ent_tmp = rtree->GetEntries();
  const Int_t rtree_ent = rtree_ent_tmp;
  Int_t runnum_arr[rtree_ent];
  rtree->SetBranchAddress("i",&index);
  rtree->SetBranchAddress("runnum",&runnum);
  for(Int_t i=0; i<rtree->GetEntries(); i++)
  {
    rtree->GetEntry(i);
    if(i+1 == index) runnum_arr[i] = runnum;
    else 
    {
      fprintf(stderr,"ERROR: rtree file problem\n");
      return;
    };
    // FILTER OUT RUNS (FOR PLOTTING WDIST BEFORE & AFTER MAJOR CHANGES)
    if(FILTER==1 && runnum > 16078015) runnum_arr[i]=0; // choose period before DSM update
    if(FILTER==2 && runnum < 16078015) runnum_arr[i]=0; // choose after before DSM update
  };


  // phi_dist [10*event_class + spin] [eta] [pt] [energy] [run index - 1]
  TH1D * phi_dist[10*N_CLASS+4][eta_bins][pt_bins][en_bins][rtree_ent];
  TH1D * pt_wdist[N_CLASS][eta_bins][en_bins][rtree_ent];
  TH1D * en_wdist[N_CLASS][eta_bins][pt_bins][rtree_ent];
  TH1D * mm_wdist[N_CLASS][eta_bins][pt_bins][en_bins][rtree_ent];
  
  // infile_phi_arr [10*event_class + spin] [eta] [pt] [energy] [phi file]
  TObjArray * infile_phi_arr[10*N_CLASS+4][eta_bins][pt_bins][en_bins][NFILES];
  TObjArray * infile_pt_wdist_arr[N_CLASS][eta_bins][en_bins][NFILES];
  TObjArray * infile_en_wdist_arr[N_CLASS][eta_bins][pt_bins][NFILES];
  TObjArray * infile_mm_wdist_arr[N_CLASS][eta_bins][pt_bins][en_bins][NFILES];
  char infile_phi_arr_n[10*N_CLASS+4][eta_bins][pt_bins][en_bins][200];
  char infile_pt_wdist_arr_n[N_CLASS][eta_bins][en_bins][200];
  char infile_en_wdist_arr_n[N_CLASS][eta_bins][pt_bins][200];
  char infile_mm_wdist_arr_n[N_CLASS][eta_bins][pt_bins][en_bins][200];
  TString tmpstr;

  Int_t inrun;
  Bool_t filter_phi[N_CLASS][rtree_ent];
  for(Int_t c=0; c<N_CLASS; c++)
  {
    for(Int_t rr=0; rr<rtree_ent; rr++) 
    {
      filter_phi[c][rr]=false;
    };
  };

  for(Int_t f=0; f<NFILES; f++)
  {
    phi_file[f]->cd(); // focus on next TFile
    // organise phi dists
    for(Int_t c=0; c<N_CLASS; c++)
    {
      for(Int_t s=0; s<4; s++)
      {
        for(Int_t g=0; g<eta_bins; g++)
        {
          for(Int_t p=0; p<pt_bins; p++)
          {
            for(Int_t e=0; e<en_bins; e++)
            {
              // set up TObjArray names to be read (only needs one execution)
              if(f==0)
              {
                sprintf(infile_phi_arr_n[10*c+s][g][p][e],"/%s/phi_dist_%s_s%d_g%d_p%d_e%d",
                  ev->Name(c),ev->Name(c),s,g,p,e);
              };
              // read TObjArrays
              infile_phi_arr[10*c+s][g][p][e][f] = (TObjArray*)phi_file[f]->Get(infile_phi_arr_n[10*c+s][g][p][e]);

              // loop through  TObjArrays
              for(Int_t o=0; o<infile_phi_arr[10*c+s][g][p][e][f]->GetEntries(); o++)
              {
                // get run number "inrun" (last 8 characters)
                tmpstr = TString(infile_phi_arr[10*c+s][g][p][e][f]->At(o)->GetName());
                sscanf(tmpstr(tmpstr.Length()-8,tmpstr.Length()).Data(),"%d",&inrun);

                // linear hash --> typcast phi_dist's
                for(Int_t h=0; h<rtree_ent; h++)
                {
                  if(inrun == runnum_arr[h])
                  {
                    phi_dist[10*c+s][g][p][e][h] = (TH1D*) infile_phi_arr[10*c+s][g][p][e][f]->At(o);
                    filter_phi[c][h]=true;
                  };
                };
              };
            };
          };
        };
      };
    };

    // organise wdists
    for(Int_t c=0; c<N_CLASS; c++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          if(f==0) sprintf(infile_pt_wdist_arr_n[c][g][e],"pt_wdist_%s_g%d_e%d",ev->Name(c),g,e);
          infile_pt_wdist_arr[c][g][e][f] = (TObjArray*)phi_file[f]->Get(infile_pt_wdist_arr_n[c][g][e]);
          for(Int_t o=0; o<infile_pt_wdist_arr[c][g][e][f]->GetEntries(); o++)
          {
            tmpstr = TString(infile_pt_wdist_arr[c][g][e][f]->At(o)->GetName());
            sscanf(tmpstr(tmpstr.Length()-8,tmpstr.Length()).Data(),"%d",&inrun);
            for(Int_t h=0; h<rtree_ent; h++)
            {
              if(inrun == runnum_arr[h])
              {
                pt_wdist[c][g][e][h] = (TH1D*) infile_pt_wdist_arr[c][g][e][f]->At(o);
              };
            };
          };
        };
        for(Int_t p=0; p<pt_bins; p++)
        {
          if(f==0) sprintf(infile_en_wdist_arr_n[c][g][p],"en_wdist_%s_g%d_p%d",ev->Name(c),g,p);
          infile_en_wdist_arr[c][g][p][f] = (TObjArray*)phi_file[f]->Get(infile_en_wdist_arr_n[c][g][p]);
          for(Int_t o=0; o<infile_en_wdist_arr[c][g][p][f]->GetEntries(); o++)
          {
            tmpstr = TString(infile_en_wdist_arr[c][g][p][f]->At(o)->GetName());
            sscanf(tmpstr(tmpstr.Length()-8,tmpstr.Length()).Data(),"%d",&inrun);
            for(Int_t h=0; h<rtree_ent; h++)
            {
              if(inrun == runnum_arr[h])
              {
                en_wdist[c][g][p][h] = (TH1D*) infile_en_wdist_arr[c][g][p][f]->At(o);
              };
            };
          };
        };
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            if(f==0) sprintf(infile_mm_wdist_arr_n[c][g][p][e],"mm_wdist_%s_g%d_p%d_e%d",ev->Name(c),g,p,e);
            infile_mm_wdist_arr[c][g][p][e][f] = (TObjArray*)phi_file[f]->Get(infile_mm_wdist_arr_n[c][g][p][e]);
            for(Int_t o=0; o<infile_mm_wdist_arr[c][g][p][e][f]->GetEntries(); o++)
            {
              tmpstr = TString(infile_mm_wdist_arr[c][g][p][e][f]->At(o)->GetName());
              sscanf(tmpstr(tmpstr.Length()-8,tmpstr.Length()).Data(),"%d",&inrun);
              for(Int_t h=0; h<rtree_ent; h++)
              {
                if(inrun == runnum_arr[h])
                {
                  mm_wdist[c][g][p][e][h] = (TH1D*) infile_mm_wdist_arr[c][g][p][e][f]->At(o);
                };
              };
            };
          };
        };
      };
    };
  };


  // build final TObjArrays, one for each kinematic/geometric bin
  char outfile_name[1024];
  sprintf(outfile_name,"%s/all.root",RD->env->phiset_dir);
  TFile * outfile = new TFile(outfile_name,"RECREATE");
  outfile->cd();
  for(Int_t c=0; c<N_CLASS; c++) outfile->mkdir(ev->Name(c));
  TObjArray * combined_phi_array[10*N_CLASS+4][eta_bins][pt_bins][en_bins];
  TObjArray * combined_pt_wdist_array[N_CLASS][eta_bins][en_bins];
  TObjArray * combined_en_wdist_array[N_CLASS][eta_bins][pt_bins];
  TObjArray * combined_mm_wdist_array[N_CLASS][eta_bins][pt_bins][en_bins];
  char combined_phi_array_n[10*N_CLASS+4][eta_bins][pt_bins][en_bins][200];
  char combined_pt_wdist_array_n[N_CLASS][eta_bins][en_bins][200];
  char combined_en_wdist_array_n[N_CLASS][eta_bins][pt_bins][200];
  char combined_mm_wdist_array_n[N_CLASS][eta_bins][pt_bins][en_bins][200];

  printf("--------------------------------------------------------\n");

  for(Int_t c=0; c<N_CLASS; c++)
  {
    for(Int_t s=0; s<4; s++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            combined_phi_array[10*c+s][g][p][e] = new TObjArray();

            sprintf(combined_phi_array_n[10*c+s][g][p][e],"phi_dist_%s_s%d_g%d_p%d_e%d",
              ev->Name(c),s,g,p,e);

            for(Int_t r=0; r<rtree_ent; r++)
            {
              if(phi_dist[10*c+s][g][p][e][r]!=NULL && filter_phi[c][r]==1)
                combined_phi_array[10*c+s][g][p][e]->AddLast(phi_dist[10*c+s][g][p][e][r]);
            };
            outfile->cd();
            outfile->cd(ev->Name(c)); 
            combined_phi_array[10*c+s][g][p][e]->Write(combined_phi_array_n[10*c+s][g][p][e],
              TObject::kSingleKey);
          };
        };
      };
    };
  };
  for(Int_t c=0; c<N_CLASS; c++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        combined_pt_wdist_array[c][g][e] = new TObjArray();
        sprintf(combined_pt_wdist_array_n[c][g][e],"pt_wdist_%s_g%d_e%d",ev->Name(c),g,e);
        for(Int_t r=0; r<rtree_ent; r++)
        {
          if(pt_wdist[c][g][e][r]!=NULL && filter_phi[c][r]==1)
            combined_pt_wdist_array[c][g][e]->AddLast(pt_wdist[c][g][e][r]);
        };
        outfile->cd(); 
        combined_pt_wdist_array[c][g][e]->Write(combined_pt_wdist_array_n[c][g][e],TObject::kSingleKey);
      };
      for(Int_t p=0; p<pt_bins; p++)
      {
        combined_en_wdist_array[c][g][p] = new TObjArray();
        sprintf(combined_en_wdist_array_n[c][g][p],"en_wdist_%s_g%d_p%d",ev->Name(c),g,p);
        for(Int_t r=0; r<rtree_ent; r++)
        {
          if(en_wdist[c][g][p][r]!=NULL && filter_phi[c][r]==1)
            combined_en_wdist_array[c][g][p]->AddLast(en_wdist[c][g][p][r]);
        };
        outfile->cd(); 
        combined_en_wdist_array[c][g][p]->Write(combined_en_wdist_array_n[c][g][p],TObject::kSingleKey);
      };
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          combined_mm_wdist_array[c][g][p][e] = new TObjArray();
          sprintf(combined_mm_wdist_array_n[c][g][p][e],"mm_wdist_%s_g%d_p%d_e%d",ev->Name(c),g,p,e);
          for(Int_t r=0; r<rtree_ent; r++)
          {
            if(mm_wdist[c][g][p][e][r]!=NULL && filter_phi[c][r]==1)
              combined_mm_wdist_array[c][g][p][e]->AddLast(mm_wdist[c][g][p][e][r]);
          };
          outfile->cd();
          combined_mm_wdist_array[c][g][p][e]->Write(combined_mm_wdist_array_n[c][g][p][e],TObject::kSingleKey);
        };
      };
    };
  };




  // initialise tot wdists
  outfile->cd();
  TH1D * pt_wdist_tot[N_CLASS][eta_bins][en_bins];
  TH1D * en_wdist_tot[N_CLASS][eta_bins][pt_bins];
  TH1D * mm_wdist_tot[N_CLASS][eta_bins][pt_bins][en_bins];
  char pt_wdist_n[N_CLASS][eta_bins][en_bins][64];
  char en_wdist_n[N_CLASS][eta_bins][pt_bins][64];
  char mm_wdist_n[N_CLASS][eta_bins][pt_bins][en_bins][64];
  Int_t NWBINS = ((TH1D*)(combined_pt_wdist_array[0][0][0]->At(0)))->GetNbinsX();
  for(Int_t g=0; g<eta_bins; g++)
  {
    for(Int_t e=0; e<en_bins; e++)
    {
      for(Int_t c=0; c<N_CLASS; c++)
      {
        sprintf(pt_wdist_n[c][g][e],"pt_wdist_tot_%s_g%d_e%d",ev->Name(c),g,e);
        pt_wdist_tot[c][g][e] = new TH1D(pt_wdist_n[c][g][e],pt_wdist_n[c][g][e],NWBINS,RD->env->PtLow,RD->env->PtHigh);
      };
    };
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t c=0; c<N_CLASS; c++)
      {
        sprintf(en_wdist_n[c][g][p],"en_wdist_tot_%s_g%d_p%d",ev->Name(c),g,p);
        en_wdist_tot[c][g][p] = new TH1D(en_wdist_n[c][g][p],en_wdist_n[c][g][p],NWBINS,RD->env->EnLow,RD->env->EnHigh);
      };
    };
    for(Int_t p=0; p<pt_bins; p++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        for(Int_t c=0; c<N_CLASS; c++)
        {
          sprintf(mm_wdist_n[c][g][p][e],"mm_wdist_tot_%s_g%d_p%d_e%d",ev->Name(c),g,p,e);
          mm_wdist_tot[c][g][p][e] = new TH1D(mm_wdist_n[c][g][p][e],mm_wdist_n[c][g][p][e],NWBINS,0,1);
        };
      };
    };
  };


  // fill tot wdists
  Double_t bc,bc_old;
  for(Int_t c=0; c<N_CLASS; c++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        for(Int_t b=1; b<=NWBINS; b++)
        {
          for(Int_t o=0; o<combined_pt_wdist_array[c][g][e]->GetEntries(); o++)
          {
            bc = ((TH1D*)(combined_pt_wdist_array[c][g][e]->At(o)))->GetBinContent(b);
            bc_old = pt_wdist_tot[c][g][e]->GetBinContent(b);
            pt_wdist_tot[c][g][e]->SetBinContent(b,bc+bc_old);
          };
        };
      };
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t b=1; b<=NWBINS; b++)
        {
          for(Int_t o=0; o<combined_en_wdist_array[c][g][p]->GetEntries(); o++)
          {
            bc = ((TH1D*)(combined_en_wdist_array[c][g][p]->At(o)))->GetBinContent(b);
            bc_old = en_wdist_tot[c][g][p]->GetBinContent(b);
            en_wdist_tot[c][g][p]->SetBinContent(b,bc+bc_old);
          };
        };
      };
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          for(Int_t b=1; b<=NWBINS; b++)
          {
            for(Int_t o=0; o<combined_mm_wdist_array[c][g][p][e]->GetEntries(); o++)
            {
              bc = ((TH1D*)(combined_mm_wdist_array[c][g][p][e]->At(o)))->GetBinContent(b);
              bc_old = mm_wdist_tot[c][g][p][e]->GetBinContent(b);
              mm_wdist_tot[c][g][p][e]->SetBinContent(b,bc+bc_old);
            };
          };
        };
      };
    };
  };

  // write tot wdists
  for(Int_t c=0; c<N_CLASS; c++)
  {
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t e=0; e<en_bins; e++) pt_wdist_tot[c][g][e]->Write();
      for(Int_t p=0; p<pt_bins; p++) en_wdist_tot[c][g][p]->Write();
      for(Int_t p=0; p<pt_bins; p++) for(Int_t e=0; e<en_bins; e++) mm_wdist_tot[c][g][p][e]->Write();
    };
  };


  printf("%s written\n",outfile_name);



  // print wdists for each phiset file
  if(printPDFs)
  {
    char pt_wdist_pdf[N_CLASS][eta_bins][en_bins][64];
    char pt_wdist_pdfl[N_CLASS][eta_bins][en_bins][64];
    char pt_wdist_pdfr[N_CLASS][eta_bins][en_bins][64];
    char en_wdist_pdf[N_CLASS][eta_bins][pt_bins][64];
    char en_wdist_pdfl[N_CLASS][eta_bins][pt_bins][64];
    char en_wdist_pdfr[N_CLASS][eta_bins][pt_bins][64];
    char mm_wdist_pdf[N_CLASS][eta_bins][pt_bins][en_bins][64];
    char mm_wdist_pdfl[N_CLASS][eta_bins][pt_bins][en_bins][64];
    char mm_wdist_pdfr[N_CLASS][eta_bins][pt_bins][en_bins][64];
    char wdist_dir[16]; strcpy(wdist_dir,"wdist_pdfs");
    char clear_pdfs[32]; sprintf(clear_pdfs,"rm -v %s/*.pdf",wdist_dir);
    system(clear_pdfs);
    for(Int_t g=0; g<eta_bins; g++)
    {
      for(Int_t e=0; e<en_bins; e++)
      {
        for(Int_t c=0; c<N_CLASS; c++)
        {
          sprintf(pt_wdist_pdf[c][g][e],"%s/pt_wdist_%s_g%d_e%d.pdf",wdist_dir,ev->Name(c),g,e);
          sprintf(pt_wdist_pdfl[c][g][e],"%s(",pt_wdist_pdf[c][g][e]);
          sprintf(pt_wdist_pdfr[c][g][e],"%s)",pt_wdist_pdf[c][g][e]);
        };
      };
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t c=0; c<N_CLASS; c++)
        {
          sprintf(en_wdist_pdf[c][g][p],"%s/en_wdist_%s_g%d_p%d.pdf",wdist_dir,ev->Name(c),g,p);
          sprintf(en_wdist_pdfl[c][g][p],"%s(",en_wdist_pdf[c][g][p]);
          sprintf(en_wdist_pdfr[c][g][p],"%s)",en_wdist_pdf[c][g][p]);
        };
      };
      for(Int_t p=0; p<pt_bins; p++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          for(Int_t c=0; c<N_CLASS; c++)
          {
            sprintf(mm_wdist_pdf[c][g][p][e],"%s/mm_wdist_%s_g%d_p%d_e%d.pdf",wdist_dir,ev->Name(c),g,p,e);
            sprintf(mm_wdist_pdfl[c][g][p][e],"%s(",mm_wdist_pdf[c][g][p][e]);
            sprintf(mm_wdist_pdfr[c][g][p][e],"%s)",mm_wdist_pdf[c][g][p][e]);
          };
        };
      };
    };
    TCanvas * cc = new TCanvas("cc","cc",700,500);
    for(Int_t c=0; c<N_CLASS; c++)
    {
      for(Int_t g=0; g<eta_bins; g++)
      {
        for(Int_t e=0; e<en_bins; e++)
        {
          for(Int_t o=0; o<combined_pt_wdist_array[c][g][e]->GetEntries(); o++)
          {
            //cc->SetLogy();
            ((TH1D*)(combined_pt_wdist_array[c][g][e]->At(o)))->Draw();
            if(o==0) cc->Print(pt_wdist_pdfl[c][g][e],"pdf");
            else if(o+1==combined_pt_wdist_array[c][g][e]->GetEntries()) cc->Print(pt_wdist_pdfr[c][g][e],"pdf");
            else cc->Print(pt_wdist_pdf[c][g][e],"pdf");
            cc->Clear();
          };
        };
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t o=0; o<combined_en_wdist_array[c][g][p]->GetEntries(); o++)
          {
            //cc->SetLogy();
            ((TH1D*)(combined_en_wdist_array[c][g][p]->At(o)))->Draw();
            if(o==0) cc->Print(en_wdist_pdfl[c][g][p],"pdf");
            else if(o+1==combined_en_wdist_array[c][g][p]->GetEntries()) cc->Print(en_wdist_pdfr[c][g][p],"pdf");
            else cc->Print(en_wdist_pdf[c][g][p],"pdf");
            cc->Clear();
          };
        };
        for(Int_t p=0; p<pt_bins; p++)
        {
          for(Int_t e=0; e<en_bins; e++)
          {
            for(Int_t o=0; o<combined_mm_wdist_array[c][g][p][e]->GetEntries(); o++)
            {
              //cc->SetLogy();
              ((TH1D*)(combined_mm_wdist_array[c][g][p][e]->At(o)))->Draw();
              if(o==0) cc->Print(mm_wdist_pdfl[c][g][p][e],"pdf");
              else if(o+1==combined_mm_wdist_array[c][g][p][e]->GetEntries()) cc->Print(mm_wdist_pdfr[c][g][p][e],"pdf");
              else cc->Print(mm_wdist_pdf[c][g][p][e],"pdf");
            };
          };
        };
      };
    };
  };
}
