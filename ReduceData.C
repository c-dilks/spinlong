// reads TwoTr from Output file and extracts
// only the essential subset for spin analysis,
// adding on the relative luminosities for each run
// -- see code under comment "reduction cut"
//    -- (e.g., N12==2 events with  specific E and Pt range)
// -- writes output to redset/
//
// 145ka

void ReduceData(const char * filename="Outputset080Ba.root")
{
  // load polarization and rellum data
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  LevelTwo * T = new LevelTwo(RD->env);


  char root_file[256];
  sprintf(root_file,"%s/%s",RD->env->output_dir,filename);
  printf("opening %s\n",root_file);
  TFile * infile = new TFile(root_file,"READ");
  TTree * twotr = (TTree*) infile->Get("TwoTr");

  char outname[256];
  sprintf(outname,"%s/Red%s",RD->env->redset_dir,filename);
  printf("reducing to data set %s\n",outname);
  TFile * outfile = new TFile(outname,"RECREATE");
  TTree * str = new TTree("str","str");

  Int_t ent = twotr->GetEntries();
  printf(" TwoTr->GetEntries() = %d\n",ent);

  Int_t ClIndex,Nclust,Evnum;
  Float_t M12,N12,E12,Phi,Eta,Pt,Z;
  Int_t spin,TrigBits,runnum,fill,Bunchid7bit;
  Bool_t kicked_str;
  Int_t blue_str,yell_str,spin_str;
  Float_t b_pol_str,y_pol_str;
  Float_t R_bbc[10];
  Float_t R_zdc[10];
  Float_t R_vpd[10];
  Float_t R_bbc_err[10];
  Float_t R_zdc_err[10];
  Float_t R_vpd_err[10];
  Bool_t isConsistent;
  Float_t b_pol;
  Float_t y_pol;
  Float_t b_pol_err;
  Float_t y_pol_err;
  Int_t pattern;
  Int_t L2sum[2];
  UInt_t lastdsm[8];

  twotr->SetBranchAddress("ClIndex",&ClIndex);
  twotr->SetBranchAddress("Nclust",&Nclust);
  twotr->SetBranchAddress("Evnum",&Evnum);
  twotr->SetBranchAddress("M12",&M12);
  twotr->SetBranchAddress("N12",&N12);
  twotr->SetBranchAddress("E12",&E12);
  twotr->SetBranchAddress("Z",&Z);
  twotr->SetBranchAddress("Phi",&Phi);
  twotr->SetBranchAddress("Eta",&Eta);
  twotr->SetBranchAddress("Pt",&Pt);
  twotr->SetBranchAddress("spin",&spin);
  //twotr->SetBranchAddress("TrigBits",&TrigBits);
  twotr->SetBranchAddress("Rnum",&runnum);
  twotr->SetBranchAddress("Bunchid7bit",&Bunchid7bit);
  twotr->SetBranchAddress("L2sum",L2sum);
  twotr->SetBranchAddress("lastdsm",lastdsm);

  str->Branch("ClIndex",&ClIndex,"ClIndex/I");
  str->Branch("Nclust",&Nclust,"Nclust/I");
  str->Branch("Evnum",&Evnum,"Evnum/I");
  str->Branch("M12",&M12,"M12/F");
  str->Branch("N12",&N12,"N12/F");
  str->Branch("E12",&E12,"E12/F");
  str->Branch("Z",&Z,"Z/F");
  str->Branch("Phi",&Phi,"Phi/F");
  str->Branch("Eta",&Eta,"Eta/F");
  str->Branch("Pt",&Pt,"Pt/F");
  str->Branch("spin",&spin_str,"spin/I"); // spinbit; set to 40 if kicked or abort/empty
  str->Branch("blue",&blue_str,"blue/I");
  str->Branch("yell",&yell_str,"yell/I");
  str->Branch("kicked",&kicked_str,"kicked/O");
  str->Branch("Bunchid7bit",&Bunchid7bit,"Bunchid7bit/I");
  //str->Branch("TrigBits",&TrigBits,"TrigBits/I");
  str->Branch("runnum",&runnum,"runnum/I");
  str->Branch("fill",&fill,"fill/I");
  str->Branch("L2sum",L2sum,"L2sum[2]/I");
  str->Branch("lastdsm",lastdsm,"lastdsm[8]/i");

  /*
  char R_bbc_n[10][32];
  char R_zdc_n[10][32];
  char R_vpd_n[10][32];
  char R_bbc_err_n[10][32];
  char R_zdc_err_n[10][32];
  char R_vpd_err_n[10][32];
  char R_bbc_d[10][32];
  char R_zdc_d[10][32];
  char R_vpd_d[10][32];
  char R_bbc_err_d[10][32];
  char R_zdc_err_d[10][32];
  char R_vpd_err_d[10][32];
  for(Int_t r=1; r<10; r++)
  {
    sprintf(R_bbc_n[r],"R%d_bbc",r);
    sprintf(R_zdc_n[r],"R%d_zdc",r);
    sprintf(R_vpd_n[r],"R%d_vpd",r);
    sprintf(R_bbc_err_n[r],"R%d_bbc_err",r);
    sprintf(R_zdc_err_n[r],"R%d_zdc_err",r);
    sprintf(R_vpd_err_n[r],"R%d_vpd_err",r);

    sprintf(R_bbc_d[r],"R%d_bbc/F",r);
    sprintf(R_zdc_d[r],"R%d_zdc/F",r);
    sprintf(R_vpd_d[r],"R%d_vpd/F",r);
    sprintf(R_bbc_err_d[r],"R%d_bbc_err/F",r);
    sprintf(R_zdc_err_d[r],"R%d_zdc_err/F",r);
    sprintf(R_vpd_err_d[r],"R%d_vpd_err/F",r);

    str->Branch(R_bbc_n[r],&(R_bbc[r]),R_bbc_d[r]);
    str->Branch(R_bbc_err_n[r],&(R_bbc_err[r]),R_bbc_err_d[r]);
    str->Branch(R_zdc_n[r],&(R_zdc[r]),R_zdc_d[r]);
    str->Branch(R_zdc_err_n[r],&(R_zdc_err[r]),R_zdc_err_d[r]);
    str->Branch(R_vpd_n[r],&(R_vpd[r]),R_vpd_d[r]);
    str->Branch(R_vpd_err_n[r],&(R_vpd_err[r]),R_vpd_err_d[r]);
  };
  */

  str->Branch("isConsistent",&isConsistent,"isConsistent/O");
  str->Branch("b_pol",&b_pol,"b_pol/F");
  str->Branch("y_pol",&y_pol,"y_pol/F");
  str->Branch("b_pol_err",&b_pol_err,"b_pol_err/F");
  str->Branch("y_pol_err",&y_pol_err,"y_pol_err/F");
  str->Branch("pattern",&pattern,"pattern/I");

  
  // reduction loop
  Int_t runnum_tmp=0;
  for(Int_t q=0; q<twotr->GetEntries(); q++)
  {
    if(q%1000 == 0) printf("%d/%d\n",q+1,twotr->GetEntries());

    twotr->GetEntry(q);

    T->runnum = runnum;
    for(int x=0; x<2; x++) T->L2sum[x] = L2sum[x];
    for(int x=0; x<8; x++) T->tcu->lastdsm[x] = lastdsm[x];


    // reduction cut
    if(/*M12>=0 &&*/ T->Fired("All"))
    {
      if(runnum!=runnum_tmp)
      {
        /*
        for(Int_t r=1; r<10; r++)
        {
          R_bbc[r] = RD->Rellum(runnum,r,"bbc");
          R_zdc[r] = RD->Rellum(runnum,r,"zdc");
          R_vpd[r] = RD->Rellum(runnum,r,"vpd");

          R_bbc_err[r] = RD->RellumErr(runnum,r,"bbc");
          R_zdc_err[r] = RD->RellumErr(runnum,r,"zdc");
          R_vpd_err[r] = RD->RellumErr(runnum,r,"vpd");
        };
        */

        isConsistent = RD->RellumConsistent(runnum);
        b_pol = RD->BluePol(runnum);
        y_pol = RD->YellPol(runnum);
        b_pol_err = RD->BluePolErr(runnum);
        y_pol_err = RD->YellPolErr(runnum);
        pattern = RD->Pattern(runnum);
        fill = RD->GetFill(runnum);

        runnum_tmp = runnum;
      };

      // BXING SHIFT CORRECTION (for run 12 only)
      if(runnum/1000000==13) Bunchid7bit = (Bunchid7bit+1)%120;

      blue_str = RD->BlueSpin(runnum,Bunchid7bit);
      yell_str = RD->YellSpin(runnum,Bunchid7bit);
      kicked_str = RD->Kicked(runnum,Bunchid7bit);


      // add entry to tree iff there exists rellum and polarization data
      // AND if relative luminosity passed consistency check
      if(fill!=0 && b_pol>0 && y_pol>0 && isConsistent)
      {
        if(blue_str==-1 && yell_str==-1) spin_str = 0;
        else if(blue_str==-1 && yell_str==1) spin_str = 1;
        else if(blue_str==1 && yell_str==-1) spin_str = 2;
        else if(blue_str==1 && yell_str==1) spin_str = 3;
        else spin_str=40;
        if(kicked_str)
        {
          spin=40;
          blue_str=0;
          yell_str=0;
        };

        str->Fill();
      };
    };
  };

  str->Write("str");
};
