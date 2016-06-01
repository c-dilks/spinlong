// Implementation of L1 trigger bits (from TCU input)
// -- see ../TCUbits directory for further information

#include "TCUbits.h"

ClassImp(TCUbits)

using namespace std;

TCUbits::TCUbits()
{
  // open tcu bit tables
  env = new Environ();
  char tcu_tr_file[512];
  char tcuchan_tr_file[512];
  sprintf(tcu_tr_file,"%s/TCUbits/tcu.dat",env->SpinDir);
  sprintf(tcuchan_tr_file,"%s/TCUbits/tcuchan.dat",env->SpinDir);
  tcu_tr = new TTree();
  tcuchan_tr = new TTree();
  tcu_tr->ReadFile(tcu_tr_file,"dsm/C:bit/i:trigger/C");
  tcuchan_tr->ReadFile(tcuchan_tr_file,"dsm/C:tcuchan/i");
  UInt_t bit,tcuchan;
  char dsm[16];
  char trig[32];
  tcu_tr->SetBranchAddress("dsm",dsm);
  tcu_tr->SetBranchAddress("bit",&bit);
  tcu_tr->SetBranchAddress("trigger",trig);
  tcuchan_tr->SetBranchAddress("dsm",dsm);
  tcuchan_tr->SetBranchAddress("tcuchan",&tcuchan);


  // fill tcu_bit map :: trigger --> (DSM,bit)
  for(Int_t i=0; i<tcu_tr->GetEntries(); i++)
  {
    tcu_tr->GetEntry(i);
    tcu_bit.insert(pair<string,pair<string,UInt_t> >(string(trig),pair<string,UInt_t>(string(dsm),bit)));
  };

  // fill tcu_chan map :: DSM --> TCU channel
  for(Int_t i=0; i<tcuchan_tr->GetEntries(); i++)
  {
    tcuchan_tr->GetEntry(i);
    tcu_chan.insert(pair<string,UInt_t>(string(dsm),tcuchan));
  };
  debug=false;

  // initialize RP trigger definitions
  InitRPdefs();
};


// RP trigger definitions initialization
void TCUbits::InitRPdefs()
{
  Int_t ii=0;
  rp_idx.insert(pair<string,Int_t>(string("EOR"),ii++));
  rp_idx.insert(pair<string,Int_t>(string("WOR"),ii++));
  rp_idx.insert(pair<string,Int_t>(string("EXOR"),ii++));
  rp_idx.insert(pair<string,Int_t>(string("WXOR"),ii++));
  rp_idx.insert(pair<string,Int_t>(string("SDE"),ii++));
  rp_idx.insert(pair<string,Int_t>(string("SDW"),ii++));
  rp_idx.insert(pair<string,Int_t>(string("ET"),ii++));
  rp_idx.insert(pair<string,Int_t>(string("IT"),ii++));
  rp_idx.insert(pair<string,Int_t>(string("DD"),ii++));

  NRP = ii;

  map<string,Int_t>::iterator iter;
  for(iter=rp_idx.begin(); iter!=rp_idx.end(); ++iter)
  {
    rp_name.insert(pair<Int_t,string>(iter->second,iter->first));
  };
};


const char * TCUbits::RPname(Int_t idx0)
{
  string return_str;
  try { return_str = rp_name.at(idx0); }
  catch(const out_of_range& e)
  {
    fprintf(stderr,"ERROR: rp_name out of range\n");
    return "";
  };
  return return_str.data();
};

Int_t TCUbits::RPidx(char * name0)
{
  Int_t return_idx;
  try { return_idx = rp_idx.at(string(name0)); }
  catch(const out_of_range& e)
  {
    fprintf(stderr,"ERROR: rp_idx out of range\n");
    return -1;
  };
  return return_idx;
};



Bool_t TCUbits::Fired(char * trg)
{
  UInt_t tcuchan0,bit0;
  string dsm0;
  try
  {
    dsm0 = tcu_bit.at(string(trg)).first;
    bit0 = tcu_bit.at(string(trg)).second;
    tcuchan0 = tcu_chan.at(dsm0);
  }
  catch(const out_of_range& e)
  {
    fprintf(stderr,"ERROR: invalid trigger requested in TCUbits::Fired\n");
    return 0;
  };
  if(debug) printf("trigger=%s DSM=%s bit=%d TCUchan=%d\n",trg,dsm0.data(),bit0,tcuchan0);
  return (lastdsm[tcuchan0] >> bit0) & 1;
};


// TOF trigger; returns true if anything seen in TOF
Bool_t TCUbits::FiredTOF()
{
  return (Fired("TOF_UPC") ||
          Fired("TOFmult0") ||
          Fired("TOFmult1") ||
          Fired("TOFmult2") ||
          Fired("TOFmult3") ||
          Fired("TOFsector0_3") ||
          Fired("TOFsector1_4") ||
          Fired("TOFsector2_5"));
};

// BBC trigger; returns true if anything seen in BBC
Bool_t TCUbits::FiredBBC()
{
  return (Fired("BBC-E") || Fired("BBC-W"));
};



/*
Bool_t TCUbits::FiredRP(Int_t idx0)
{
  string name0;
  try { name0 = rp_name.at(idx0); }
  catch(const out_of_range& e)
  {
    fprintf(stderr,"ERROR: RP idx out of range\n");
    return 0;
  };
  return FiredRP((char*)(name0.data()));
};
*/


// RP TRIGGER LOGIC DEFINITIONS ////////////////////////////////////
/*

EOR = E1U | E1D | E2U | E2D 
WOR = W1U | W1D | W2U | W2D 
 
Elastic Trigger:    ET = [ (E1U | E2U) & (W1D | W2D) ] | [ (E1D | E2D) & (W1U | W2U) ]
Inelastic Trigger:  IT = [ (E1U | E2U) & (W1U | W2U) ] | [ (E1D | E2D) & (W1D | W2D) ] 

Single Diffraction to W: SDE = EOR & !ZDCE & !BBCE & (ZDCW | BBCW) 
Single Diffraction to E: SDW = WOR & !ZDCW & !BBCW & (ZDCE | BBCE)

Double Diffractive: DD = EOR & WOR & TOF & !BBC

modified SD, suggested by Akio
  SDE1 = EOR & !ZDCE
  SDW1 = WOR & !ZDCW & (ZDCE | BBCE)

[ In TCU inputs: EOR,WOR,ET,IT ]

*/
/*
Bool_t TCUbits::FiredRP(char * name0)
{
  if(!strcmp(name0,"N")) return true; // no RP bias
  else if(!strcmp(name0,"EOR")) return Fired("RP_EOR");
  else if(!strcmp(name0,"WOR")) return Fired("RP_WOR");

  if(!strcmp(name0,"EXOR"))
    return Fired("RP_EOR") && !(Fired("RP_WOR"));
  if(!strcmp(name0,"WXOR"))
    return Fired("RP_WOR") && !(Fired("RP_EOR"));

  else if(!strcmp(name0,"SDE"))
    return (Fired("RP_EOR") &&
           !(Fired("ZDC-E")) && !(Fired("BBC-E")) &&
            (Fired("ZDC-W") || Fired("BBC-W")));

  else if(!strcmp(name0,"SDW"))
    return (Fired("RP_WOR") &&
           !(Fired("ZDC-W")) && !(Fired("BBC-W")) &&
            (Fired("ZDC-E") || Fired("BBC-E")));

  else if(!strcmp(name0,"SDE1"))
    return (Fired("RP_EOR") &&
           !(Fired("ZDC-E")));

  else if(!strcmp(name0,"SDW1"))
    return (Fired("RP_WOR") &&
           !(Fired("ZDC-W")) &&
            (Fired("ZDC-E") || Fired("BBC-E")));

  else if(!strcmp(name0,"IT")) return Fired("RP_IT");
  else if(!strcmp(name0,"ET")) return Fired("RP_ET");

  else if(!strcmp(name0,"DD"))
    return (Fired("RP_EOR") &&
            Fired("RP_WOR") &&
            FiredTOF() &&
            !(FiredBBC()));

  else if(!strcmp(name0,"BBCE")) return Fired("BBC-E");
  else if(!strcmp(name0,"BBCW")) return Fired("BBC-W");
  else if(!strcmp(name0,"NOT_BBCE")) return !(Fired("BBC-E"));
  else if(!strcmp(name0,"NOT_BBCW")) return !(Fired("BBC-W"));
};
*/
