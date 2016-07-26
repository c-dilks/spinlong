void TestRunInfo(Int_t runnum = 13080001) {
  gSystem->Load("src/RunInfo.so");
  RunInfo * RD = new RunInfo();
  if(!(RD->env->success)) return;

  printf("Index=%d\n",RD->Index(runnum));
  printf("IndexCounts=%d\n",RD->IndexCounts(runnum));
  printf("IndexPol=%d\n",RD->IndexPol(runnum));
  printf("Fill=%d\n",RD->GetFill(runnum));

  printf("\nVPD RELATIVE LUMINOSITY\n");
  for(int r=1; r<10; r++) {
    printf(" R%d = %f +/- %f  (%s)\n",
      r,
      RD->Rellum(runnum,r,"vpd"),
      RD->RellumErr(runnum,r,"vpd"),
      RD->RellumConsistent(runnum) ? "is consistent":"is not consistent"
    );
  };

  printf("\nZDC RELATIVE LUMINOSITY\n");
  for(int r=1; r<10; r++) {
    printf(" R%d = %f +/- %f  (%s)\n",
      r,
      RD->Rellum(runnum,r,"zdc"),
      RD->RellumErr(runnum,r,"zdc"),
      RD->RellumConsistent(runnum) ? "is consistent":"is not consistent"
    );
  };

  printf("\nPOLARIZATION\n");
  printf(" b_pol = %f +/- %f\n",RD->BluePol(runnum),RD->BluePolErr(runnum));
  printf(" y_pol = %f +/- %f\n",RD->YellPol(runnum),RD->YellPolErr(runnum));

  printf("\nSPIN PATTERN\n");
  printf(" pattern = %d\n",RD->Pattern(runnum));
  for(int b=0; b<120; b++) {
    printf("  [bx%3d] [b%2d] [y%2d] (%s)\n",
      b,
      RD->BlueSpin(runnum,b),
      RD->YellSpin(runnum,b),
      RD->Kicked(runnum,b) ? "kicked":"not kicked"
    );
  };
};
