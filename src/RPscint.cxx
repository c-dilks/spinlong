#include "RPscint.h"

ClassImp(RPscint)

namespace
{
  enum ew_enum {kE,kW};
  enum io_enum {kI,kO};
  enum ud_enum {kU,kD};
  enum ns_enum {kN,kS};
  
  const TString sEW[2] = {"E","W"};
  const TString sIO[2] = {"I","O"};
  const TString sUD[2] = {"U","D"};
  const TString sNS[2] = {"N","S"};

  // ADC threshold and TAC bounds which were used for the RP
  // triggers and TCU bits
  const Double_t ADC_TRIGGER_THRESH = 100;
  const Double_t TAC_TRIGGER_MIN = 100;
  const Double_t TAC_TRIGGER_MAX = 2500;
};


RPscint::RPscint()
{
  //for(ch=0; ch<16; ch++) printf("%s %d\n",RPname(ch).Data(),ch);

  ResetBranches();
  ResetBits();

  // MIP region ADC thresholds (naive, determined visually from code in ../RP)
  MIPthresh[0]  = 230;
  MIPthresh[1]  = 230;
  MIPthresh[2]  = 230;
  MIPthresh[3]  = 240;
  MIPthresh[4]  = 230;
  MIPthresh[5]  = 160;
  MIPthresh[6]  = 250;
  MIPthresh[7]  = 250;
  MIPthresh[8]  = 250;
  MIPthresh[9]  = 230;
  MIPthresh[10] = 220;
  MIPthresh[11] = 240;
  MIPthresh[12] = 300;
  MIPthresh[13] = 125;
  MIPthresh[14] = 240;
  MIPthresh[15] = 260;

  // TAC shifts
  // these were the shifts of the TAC distributions applied in the PullRP class, which
  // read the RP QTs; in this RPscint class, they are used to shift the TAC dists back to where
  // the were so that the actual TAC value read in the QT can be compared to the 
  // TAC_TRIGGER_MIN and TAC_TRIGGER_MAX values, which are what were applied to the TCU
  // bits
  // TAC in OFile = TAC from QT - TACshift
  // --> TAC from QT = TAC in OFile + TACshift
  TACshift[0]  = 553;
  TACshift[1]  = 537;
  TACshift[2]  = 528;
  TACshift[3]  = 600;
  TACshift[4]  = 628;
  TACshift[5]  = 446;
  TACshift[6]  = 491;
  TACshift[7]  = 521;
  TACshift[8]  = 338;
  TACshift[9]  = 290;
  TACshift[10] = 440;
  TACshift[11] = 163;
  TACshift[12] = 704;
  TACshift[13] = 768;
  TACshift[14] = 642;
  TACshift[15] = 757;
};

void RPscint::Process()
{
  ResetBits();

  if(!(N[kE]==0 && N[kW]==0))
  {
    // loop through hit channels
    for(ew=0; ew<2; ew++)
    {
      for(q=0; q<N[ew]; q++)
      {
        ADCtmp[Idx[ew][q]] = ADC[ew][q];

        // re-shift TAC back to value recorded in QT (for comparing
        // to thresholds, see comments in constructor)
        TACtmp[Idx[ew][q]] = TAC[ew][q] + TACshift[Idx[ew][q]];

        if( ADCtmp[Idx[ew][q]] > ADC_TRIGGER_THRESH &&
            TACtmp[Idx[ew][q]] > TAC_TRIGGER_MIN &&
            TACtmp[Idx[ew][q]] < TAC_TRIGGER_MAX)
        {
          fired[Idx[ew][q]] = true;

          /* stg 0 track_trg same for all mipn */
          for(mipn=0; mipn<3; mipn++) 
          {
            track_trg[ew][0][mipn]=1;
            ud_track_trg[ew][iUD(Idx[ew][q])][0][mipn]=1;
          };
        };
      };
    };


    // track trigger bits
    for(ew=0; ew<2; ew++)
    {
      // loop through inner seqs
      for(udi=0; udi<2; udi++)
      {
        for(nsi=0; nsi<2; nsi++)
        {
          // loop through outer seqs if inner seq fired
          ii = EiunToIdx(ew,kI,udi,nsi);
          if(fired[ii])
          {
            mipn = (ADCtmp[ii] < MIPthresh[ii]) ? 1:2;
            for(udo=0; udo<2; udo++)
            {
              for(nso=0; nso<2; nso++)
              {
                // if outer seq fired too, check it's alignment
                // w.r.t. inner seq
                oo = EiunToIdx(ew,kO,udo,nso);
                if(fired[oo]) 
                {
                  track_trg[ew][1][0]=1;
                  track_trg[ew][1][mipn]=1;

                  // u/d I hits u/d O
                  if(udo==udi)
                  {
                    track_trg[ew][2][0]=1;
                    track_trg[ew][2][mipn]=1;
                    ud_track_trg[ew][udo][1][0]=1;
                    ud_track_trg[ew][udo][1][mipn]=1;

                    // n/s I hits n/s O
                    if(nso==nsi)
                    {
                      track_trg[ew][3][0]=1;
                      track_trg[ew][3][mipn]=1;
                      ud_track_trg[ew][udo][2][0]=1;
                      ud_track_trg[ew][udo][2][mipn]=1;
                    };
                  };
                }; // eo if outer seq fired
              };
            };
          }; // eo if innner seq fired
        };
      }; // eo inner seq loop
    }; // eo e/w loop


    // elatic and inelastic triggers
    for(stg=0; stg<3; stg++)
    {
      for(mipn=0; mipn<3; mipn++)
      {
        elastic_trg[stg][mipn] = (ud_track_trg[kE][kU][stg][mipn] && ud_track_trg[kW][kD][stg][mipn]) ||
                                 (ud_track_trg[kE][kD][stg][mipn] && ud_track_trg[kW][kU][stg][mipn]);
        inelastic_trg[stg][mipn] = (ud_track_trg[kE][kU][stg][mipn] && ud_track_trg[kW][kU][stg][mipn]) ||
                                   (ud_track_trg[kE][kD][stg][mipn] && ud_track_trg[kW][kD][stg][mipn]);
      };
    };
  }; // eo multiplity>0 cut

  return;
};

void RPscint::IdxToEiun(Int_t idx0, Int_t &ew0, Int_t &io0, Int_t &ud0, Int_t &ns0)
{
  ew0 = iEW(idx0);
  io0 = iIO(idx0);
  ud0 = iUD(idx0);
  ns0 = iNS(idx0);
};
  
Short_t RPscint::EiunToIdx(Int_t ew0, Int_t io0, Int_t ud0, Int_t ns0)
{
  return 8*ew0 + 4*io0 + 2*ud0 + ns0;
};

TString RPscint::RPname(Int_t idx0)
{
  return sEW[iEW(idx0)] + sIO[iIO(idx0)] + sUD[iUD(idx0)] + sNS[iNS(idx0)];
};


void RPscint::ResetBranches()
{
  // reset ADC, TAC, and multiplicity
  for(ew=0; ew<2; ew++)
  {
    N[ew]=0;
    for(ch=0; ch<8; ch++)
    {
      Idx[ew][ch]=0;
      ADC[ew][ch]=0;
      TAC[ew][ch]=0;
    };
  };
  vertex=0;
};


void RPscint::ResetBits()
{
  // reset bits
  for(int i=0; i<16; i++) fired[i]=false;
  for(ew=0; ew<2; ew++)
  {
    for(stg=0; stg<4; stg++)
    {
      for(mipn=0; mipn<3; mipn++)
      {
        track_trg[ew][stg][mipn]=false;
        if(stg<3)
        {
          for(ud=0; ud<2; ud++)
          {
            ud_track_trg[ew][ud][stg][mipn]=false;
          };
        };
      };
    };
  };
  for(stg=0; stg<3; stg++)
  {
    for(mipn=0; mipn<3; mipn++)
    {
      elastic_trg[stg][mipn]=false;
      inelastic_trg[stg][mipn]=false;
    };
  };

  // debug reset-------------------------
  /*
  for(int i=0; i<16; i++) fired[i]=false;
  for(ew=0; ew<2; ew++)
  {
    for(stg=0; stg<4; stg++)
    {
      for(mipn=0; mipn<3; mipn++)
      {
        printf("trac_trg[%d][%d][%d]=%d\n",ew,stg,mipn,track_trg[ew][stg][mipn]=false);
        if(stg<3)
        {
          for(ud=0; ud<2; ud++)
          {
            printf("ud_track[%d][%d][%d][%d]=%d\n",ew,ud,stg,mipn,ud_track_trg[ew][ud][stg][mipn]=false);
          };
        };
      };
    };
  };
  for(stg=0; stg<3; stg++)
  {
    for(mipn=0; mipn<3; mipn++)
    {
      printf("elastic_trg[%d][%d]=%d\n",stg,mipn,elastic_trg[stg][mipn]=false);
      printf("inelastic_trg[%d][%d]=%d\n",stg,mipn,inelastic_trg[stg][mipn]=false);
    };
  };
  */
  // debug reset-------------------------
};
