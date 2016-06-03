// Implementation of L2 trigger daq ids, used for masking 
// L2sum in the trigger files

#ifndef LevelTwo_
#define LevelTwo_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <map>
#include <stdexcept>
#include <bitset>

#include "TSystem.h"
#include "TObject.h"
#include "TTree.h"
#include "TString.h"

#include "Environ.h"
#include "TCUbits.h"

class LevelTwo : public TObject
{
  public:
    LevelTwo(Environ * env0);

    Int_t Mask(const char * trigger0, Int_t dword);
    Int_t Mask(Int_t num0, Int_t dword);
    Int_t Index(const char * trigger0);
    const char * Name(Int_t num0);
    void PrintTrigIds();
    Int_t WhichBranch(const char * name0);
    Bool_t IsTakeall(const char * name0);
    Bool_t Fired(const char * trg);
    Bool_t Fired(Int_t num1);
    void PrintVars();

    Int_t N; // number of FMS triggers

    //------------------------
    // EVENT VARIABLES
    Int_t L2sum[2]; // it's signed in Outputfiles for some reason
    Int_t runnum;
    //------------------------
    TCUbits * tcu;
    
    Bool_t debug;

  protected:
    Environ * env;
    TTree * id_tr;
    char trig[64];

    // MAPS:
    // here trigger idx refers to the enumerator in the namespace
    //
    std::map<Int_t, char*> trigger_name; // trigger idx --> trigger name
    std::map<Int_t, char*> trigger_dbname; // trigger idx --> trigger daq name
    std::map<std::string, Int_t> trigger_idx; // trigger name --> trigger idx
    std::map<std::string, Int_t> trigger_dbidx; // trigger name --> trigger daq index
    // mask_map[2]: (it's a nested map) 
    // key = run number, value = trig_idx map
    //                           -->trig_idx map:
    //                              key = trigger index, value = daq trigger id
    std::map<Int_t, std::map<Int_t,Long_t> > mask_map;

    std::map<std::string, Bool_t> takeall; // trigger name --> is takeall

    ClassDef(LevelTwo,1);
};

#endif
