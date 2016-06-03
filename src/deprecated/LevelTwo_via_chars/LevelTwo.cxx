// Implementation of L2 trigger daq ids, used for masking 
// L2sum in the trigger files

#include "LevelTwo.h"

ClassImp(LevelTwo)

namespace
{
  const Int_t N_TRIG = 12;
  enum trigger_enum 
  {
    kAll,
    kHT0,
    kSmBS0,
    kSmBS1,
    kLgBS0,
    kLgBS1,
    kJP1,
    kJP2,
    kDijet,
    kLED,
    kBase,
    kJP1withBarrelJP0
  };
  enum which_branch_enum
  {
    kNeither,
    kL2sum,
    kLastDSM
  };
  const Int_t int_max = 0xFFFFFFFF;
};

LevelTwo::LevelTwo(Environ * env0)
{
  N = N_TRIG;
  debug = false;
  
  env = env0;
  tcu = new TCUbits(env);


  // fms trigger names
  trigger_name.insert(std::pair<Int_t,char*>(kAll,"All")); // OR of all triggers except LED
  trigger_name.insert(std::pair<Int_t,char*>(kHT0,"HT0"));
  trigger_name.insert(std::pair<Int_t,char*>(kSmBS0,"SmBS0"));
  trigger_name.insert(std::pair<Int_t,char*>(kSmBS1,"SmBS1"));
  trigger_name.insert(std::pair<Int_t,char*>(kLgBS0,"LgBS0"));
  trigger_name.insert(std::pair<Int_t,char*>(kLgBS1,"LgBS1"));
  trigger_name.insert(std::pair<Int_t,char*>(kJP1,"JP1"));
  trigger_name.insert(std::pair<Int_t,char*>(kJP2,"JP2"));
  trigger_name.insert(std::pair<Int_t,char*>(kDijet,"Dijet"));
  trigger_name.insert(std::pair<Int_t,char*>(kLED,"LED"));
  trigger_name.insert(std::pair<Int_t,char*>(kBase,"Base"));
  trigger_name.insert(std::pair<Int_t,char*>(kJP1withBarrelJP0,"JP1BEMC"));

  takeall.insert(std::pair<std::string,Bool_t>("All",false));
  takeall.insert(std::pair<std::string,Bool_t>("HT0",false));
  takeall.insert(std::pair<std::string,Bool_t>("SmBS0",false));
  takeall.insert(std::pair<std::string,Bool_t>("SmBS1",true));
  takeall.insert(std::pair<std::string,Bool_t>("LgBS0",false));
  takeall.insert(std::pair<std::string,Bool_t>("LgBS1",true));
  takeall.insert(std::pair<std::string,Bool_t>("JP1",false));
  takeall.insert(std::pair<std::string,Bool_t>("JP2",true));
  takeall.insert(std::pair<std::string,Bool_t>("Dijet",true));
  takeall.insert(std::pair<std::string,Bool_t>("LED",true));
  takeall.insert(std::pair<std::string,Bool_t>("Base",false));
  takeall.insert(std::pair<std::string,Bool_t>("JP1BEMC",false));

  trigger_dbname.insert(std::pair<Int_t,char*>(kAll,"All")); // unused
  trigger_dbname.insert(std::pair<Int_t,char*>(kHT0,"FMSHT0"));
  trigger_dbname.insert(std::pair<Int_t,char*>(kSmBS0,"FMSSmBS0"));
  trigger_dbname.insert(std::pair<Int_t,char*>(kSmBS1,"FMSSmBS1"));
  trigger_dbname.insert(std::pair<Int_t,char*>(kLgBS0,"FMSLgBS0"));
  trigger_dbname.insert(std::pair<Int_t,char*>(kLgBS1,"FMSLgBS1"));
  trigger_dbname.insert(std::pair<Int_t,char*>(kJP1,"FMSJP1"));
  trigger_dbname.insert(std::pair<Int_t,char*>(kJP2,"FMSJP2"));
  trigger_dbname.insert(std::pair<Int_t,char*>(kDijet,"FMSDijet"));
  trigger_dbname.insert(std::pair<Int_t,char*>(kLED,"FMSLed"));
  trigger_dbname.insert(std::pair<Int_t,char*>(kBase,"FMS-base"));
  trigger_dbname.insert(std::pair<Int_t,char*>(kJP1withBarrelJP0,"FMSJP1*JP0"));


  std::map<char*,Int_t> trigger_name_idx; // trigger name --> trigger idx
  std::map<char*,Int_t> trigger_dbname_idx; // trigger dbname --> trigger idx
  for(Int_t t=0; t<N_TRIG; t++)
  {
    trigger_idx.insert(std::pair<std::string,Int_t>(std::string(trigger_name[t]),t));
    trigger_dbidx.insert(std::pair<std::string,Int_t>(std::string(trigger_dbname[t]),t));
  };



  // read trigid.dat
  id_tr = new TTree();
  char trigidfile[512];
  sprintf(trigidfile,"%s/trigid.dat",env->SpinDir);
  id_tr->ReadFile(trigidfile,"rn/I:dbidx/I:name/C");
  Int_t rn,dbidx;
  char name[32];
  id_tr->SetBranchAddress("rn",&rn);
  id_tr->SetBranchAddress("dbidx",&dbidx);
  id_tr->SetBranchAddress("name",name);
  std::map<Int_t,Long_t> trig_idx; // trigger idx --> trigger dbidx
  id_tr->GetEntry(0); 
  Int_t rn_tmp=rn;
  Int_t rprn_tmp=rn;
  Int_t idx;
  Long_t mask;
  Long_t allmask=0;
  TString name_str;
  
  // build mask map
  for(Int_t i=0; i<id_tr->GetEntries(); i++)
  {
    id_tr->GetEntry(i);
    name_str = TString(name);

    // check for FMS triggers
    if(name_str(0,3)=="FMS")
    {
      // we use try/catch since not all triggers in trigid.dat
      // are being read here...
      try { idx = trigger_dbidx.at(std::string(name)); }
      catch(const std::out_of_range& e) { continue; };


      mask = (long)1 << dbidx;
      if(idx!=kLED) allmask = allmask | mask;

      trig_idx.insert(std::pair<Int_t,Long_t>(idx,mask));
      if(rn != rn_tmp || i+1==id_tr->GetEntries())
      {
        idx = trigger_idx.at(std::string("All"));
        trig_idx.insert(std::pair<Int_t,Long_t>(idx,allmask));
        mask_map.insert(std::pair<Int_t,std::map<Int_t,Long_t> >(rn_tmp,trig_idx));
        rn_tmp=rn;
        trig_idx.clear();
        allmask=0;
      };
    }
    else if(name_str(0,3)=="FPD") continue;
    else 
      fprintf(stderr,"LevelTwo::LevelTwo -- ERROR: %s is not FMS trigger\n",
        name_str.Data());
  };
};


// fms accessors
//////////////////
Int_t LevelTwo::Index(const char * trigger0)
{
  Int_t retval;
  try { retval = trigger_idx.at(std::string(trigger0)); }
  catch(const std::out_of_range& e)
  {
    fprintf(stderr,"[+] LevelTwo::Index(\"%s\") failed\n",trigger0);
    retval=-1;
  };
  return retval;
};

Int_t LevelTwo::Mask(const char * trigger0, Int_t dword)
{
  Int_t retval;
  try { retval = (mask_map.at(runnum).at(Index(trigger0)) >> (32*dword)) & int_max; }
  catch(const std::out_of_range& e) 
  {
    fprintf(stderr,"[+] LevelTwo::Mask -- no FMS triggers for run %d\n",runnum);
    retval=0;
  };
  return retval;
};

Int_t LevelTwo::Mask(Int_t num0, Int_t dword)
{
  Int_t retval;
  try { retval = (mask_map.at(runnum).at(num0) >> (32*dword)) & int_max; }
  catch(const std::out_of_range& e) 
  {
    fprintf(stderr,"[+] LevelTwo::Mask -- no FMS triggers for run %d\n",runnum);
    retval=0;
  };
  return retval;
};

const char * LevelTwo::Name(Int_t num0)
{
  TString retchar;
  try { retchar=TString(trigger_name.at(num0)); }
  catch(const std::out_of_range& e) 
  {
    fprintf(stderr,"[+] LevelTwo::Name(%d) failed\n",num0);
    retchar="";
  };
  printf("Name %s\n",retchar.Data());
  return retchar.Data();
};


// print
void LevelTwo::PrintTrigIds() {
  printf("LevelTwo Printout for run %d\n",runnum);
  printf("------------------------------------------------------------------------\n");
  printf("idx       name  L2sum[1] L2sum[0]    DSM  chan bit  which_branch takeall\n");
  TString which_branch_name[3];
  which_branch_name[kNeither] = "neither";
  which_branch_name[kLastDSM] = "lastdsm";
  which_branch_name[kL2sum] = "L2sum";
  for(int nn=0; nn<N; nn++) {
    sprintf(trig,"%s",Name(nn));
    printf("%3d %10s  %8x %8x  %5s  %4d %3d  %12s     %3s\n",
      nn,
      trig,
      Mask(nn,1),
      Mask(nn,0),
      (tcu->WhichDSM(trig)).data(),
      tcu->WhichTCUchan(trig,tcu->WhichDSM(trig)),
      tcu->WhichBit(trig),
      which_branch_name[WhichBranch(trig)].Data(),
      IsTakeall(trig)?"yes":"no"
    );
  };
};


Int_t LevelTwo::WhichBranch(const char * name0) {
  // use L2sum if !(L2sum[1]>0)
  // use lastdsm if L2sum[1]>0 && if it's not prescaled
  // else use neither
  return Mask(Index(name0),1)>0 ?
    ((tcu->WhichBit(name0)>=0 && IsTakeall(name0)) ? kLastDSM:kNeither):kL2sum;
};


Bool_t LevelTwo::IsTakeall(const char * name0) {
  Bool_t retbool;
  try { retbool = takeall.at(std::string(name0)); }
  catch(const std::out_of_range& e) { retbool=false; };
  return retbool;
};


Bool_t LevelTwo::Fired(const char * trg) {
  Bool_t retbool;
  switch(WhichBranch(trg)) {
    case kL2sum:
      retbool = L2sum[0] & Mask(trg,0);
      break;
    case kLastDSM:
      retbool = tcu->Fired(trg);
      break;
    case kNeither:
      retbool = false;
      break;
  };
  if(debug) printf(" -- trg=%s  fired=%s\n",trg,retbool?"yes":"no");
  return retbool;
};


Bool_t LevelTwo::Fired(Int_t num1) {
  printf("Fired %s\n",Name(num1));
  return Fired(Name(num1));
};
  

void LevelTwo::PrintVars() {
  printf("runnum=%d L2sum[0]=0x%x lastdsm[5]=%s.%s.%s.%s\n",
      runnum,
      L2sum[0],
      std::bitset<4>(tcu->lastdsm[5]>>12).to_string().c_str(),
      std::bitset<4>(tcu->lastdsm[5]>>8).to_string().c_str(),
      std::bitset<4>(tcu->lastdsm[5]>>4).to_string().c_str(),
      std::bitset<4>(tcu->lastdsm[5]).to_string().c_str()
  );
};
