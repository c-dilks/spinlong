// generates naive mass_cuts.dat, with just 
// simple mass cuts of pi0_mass +/- width

void InitMassCuts(Float_t width=0.1) {
  const Float_t M = 0.135;
  gSystem->RedirectOutput("mass_cuts.dat","w");
  for(Int_t E=10; E<100; E+=10) {
    printf("%f %f %f %f %f\n",(Float_t)E,
                              (Float_t)(E+10),
                              M-width,
                              M,
                              M+width
    );
  };
  gSystem->RedirectOutput(0);

  printf("\nmass_cuts.dat:\n------------------------\n");
  system("cat mass_cuts.dat");
};
