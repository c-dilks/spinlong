// Implementation of L1 trigger bits (from TCU input)
// -- see ../TCUbits directory for further information

#include "TCUbits.h"

ClassImp(TCUbits)

using namespace std;

TCUbits::TCUbits(Environ * env0)
{
  // open tcu bit tables
  debug=false;
  env = env0;
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
    if(debug) printf("TCUbits::TCUbits DEBUG -- tcu_bit: %s %s %d\n",dsm,trig,bit);
  };

  // fill tcu_chan map :: DSM --> TCU channel
  for(Int_t i=0; i<tcuchan_tr->GetEntries(); i++)
  {
    tcuchan_tr->GetEntry(i);
    tcu_chan.insert(pair<string,UInt_t>(string(dsm),tcuchan));
    if(debug) printf("TCUbits::TCUbits DEBUG -- tcu_chan: %s %d\n",dsm,tcuchan);
  };

};


Bool_t TCUbits::Fired(const char * trg)
{
  std::string dsm0 = WhichDSM(trg);
  UInt_t bit0 = WhichBit(trg);
  UInt_t tcuchan0 = WhichTCUchan(trg,dsm0);

  if(debug) printf("trigger=%s DSM=%s bit=%d TCUchan=%d\n",trg,dsm0.data(),bit0,tcuchan0);
  return (lastdsm[tcuchan0] >> bit0) & 1;
};


std::string TCUbits::WhichDSM(const char * trg) {
  std::string retstr;
  try { retstr = tcu_bit.at(string(trg)).first; }
  catch(const out_of_range& e) { 
    //fprintf(stderr,"ERROR: invalid trigger requested in TCUbits::WhichDSM\n");
    retstr="n/a";
  };
  if(debug) printf("TCUbits::WhichDSM DEBUG -- trg=%s dsm=%s\n",trg,retstr.data());
  return retstr;
};

Int_t TCUbits::WhichTCUchan(const char * trg, std::string dsm_0) {
  Int_t retval;
  try { retval = tcu_chan.at(dsm_0); }
  catch(const out_of_range& e) { 
    //fprintf(stderr,"ERROR: invalid trigger requested in TCUbits::WhichTCUchan\n");
    retval=-1;
  };
  return retval;
};

Int_t TCUbits::WhichBit(const char * trg) {
  Int_t retval;
  try { retval = tcu_bit.at(string(trg)).second; }
  catch(const out_of_range& e) { 
    //fprintf(stderr,"ERROR: invalid trigger requested in TCUbits::WhichBit\n");
    retval=-1;
  };
  return retval;
};
