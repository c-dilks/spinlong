// takes pi0 sample in redset and builds bXing distribution (number of pi0s)
// for each spin bit, and a sum
//
// this was written to see if there is any structure in the bXing distributions

void BxingDistPi0(Bool_t useRedSet=false)
{
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;
  TString filenames;

  if(useRedSet) filenames = TString(RD->env->redset_dir)+"/Red*.root";
  else filenames = TString(RD->env->output_dir)+"/Outputset*.root";

  Int_t runnum;
  TChain * tc;
  
  if(useRedSet) {
    tc = new TChain("str");
    tc->Add(filenames.Data());
    tc->SetBranchAddress("runnum",&runnum);
  }
  else {
    tc = new TChain("TwoTr");
    tc->Add(filenames.Data());
    tc->SetBranchAddress("Rnum",&runnum);
  };
  
  const Int_t NRUNS = 2000; // assume 600 max runs
  TH1D * bxing_dist[NRUNS][5];
  char bxing_dist_n[NRUNS][5][64];
  char bxing_dist_t[NRUNS][5][128];

  // set branch addresses
  Bool_t kicked;
  Float_t M12,E12,Z;
  Int_t Bunchid7bit;
  //Int_t TrigBits;
  Int_t blue;
  Int_t yell;
  Float_t N12;
  Float_t Pt;
  //tc->SetBranchAddress("kicked",&kicked);
  tc->SetBranchAddress("M12",&M12);
  tc->SetBranchAddress("E12",&E12);
  tc->SetBranchAddress("N12",&N12);
  tc->SetBranchAddress("Z",&Z);
  //tc->SetBranchAddress("TrigBits",&TrigBits);
  tc->SetBranchAddress("Bunchid7bit",&Bunchid7bit);
  //tc->SetBranchAddress("blue",&blue);
  //tc->SetBranchAddress("yell",&yell);
  tc->SetBranchAddress("Pt",&Pt);



  Int_t iter=-1;
  Int_t runnum_tmp=0;
  Int_t spin_set;
  system("touch RUNLIST.dat; rm RUNLIST.dat; touch RUNLIST.dat");
  //for(Int_t i=0; i<900000; i++) {
  for(Int_t i=0; i<tc->GetEntries(); i++) {
    tc->GetEntry(i);
    //printf("M12=%f Z=%f Pt=%f TrigBits=%d N12=%d\n",M12,Z,Pt,TrigBits,N12);
    if(i%1000000==0) printf("%d/%d\n",i,tc->GetEntries());
    if(runnum!=runnum_tmp)
    {
      iter++;
      runnum_tmp = runnum;
      for(Int_t s=0; s<5; s++)
      {
        sprintf(bxing_dist_n[iter][s],"bxing_dist_%d_%d",iter,s);
        sprintf(bxing_dist_t[iter][s],
          "pi0 bXing dist (Grn:-- Orn:-+ Red:+- Blue:++ Blk:all) -- run %d",runnum);
        bxing_dist[iter][s] = new TH1D(bxing_dist_n[iter][s],bxing_dist_t[iter][s],120,0,120);
        if(s==0)
        {
          gSystem->RedirectOutput("RUNLIST.dat");
          printf("%d\n",runnum);
          gSystem->RedirectOutput(0);
        };
      };
    };
    /*
    spin_set = -1;
    if(blue==-1 && yell==-1) spin_set = 0;
    else if(blue==-1 && yell==1) spin_set = 1;
    else if(blue==1 && yell==-1) spin_set = 2;
    else if(blue==1 && yell==1) spin_set = 3;
    if(spin_set>=0)
    {
    */
      if(fabs(M12-0.135)<0.1 && Z<0.8 && Pt<15 && fabs(N12-2)<0.001)
      {
        //bxing_dist[iter][spin_set]->Fill(Bunchid7bit);
        bxing_dist[iter][4]->Fill(Bunchid7bit);
      };
    //};
  };


  // set formatting and draw
  TCanvas * canv = new TCanvas("canv","canv",1100,500);
  TString pdfname,pdfnamel,pdfnamer;
  if(useRedSet) pdfname="pi0_bx_dist_via_redset.pdf";
  else pdfname="pi0_bx_dist_via_twotr.pdf";
  pdfnamel = pdfname+"(";
  pdfnamer = pdfname+")";
  gStyle->SetOptStat(0);
  for(Int_t i=0; i<=iter; i++)
  {
    bxing_dist[i][0]->SetLineColor(kGreen+2);
    bxing_dist[i][1]->SetLineColor(kOrange+7);
    bxing_dist[i][2]->SetLineColor(kRed);
    bxing_dist[i][3]->SetLineColor(kBlue);
    bxing_dist[i][4]->SetLineColor(kBlack);
    canv->SetGrid(1,1);
    bxing_dist[i][4]->Draw();
    //for(Int_t s=0; s<4; s++) bxing_dist[i][s]->Draw("same");
    if(i==0) canv->Print(pdfnamel.Data(),"pdf");
    else if(i==iter) canv->Print(pdfnamer.Data(),"pdf");
    else canv->Print(pdfname.Data(),"pdf");
    canv->Clear();
  }
}



