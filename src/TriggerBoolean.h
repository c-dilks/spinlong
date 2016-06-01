// main class for defining trigger booleans using TCUbits and RPscint members

#ifndef TriggerBoolean_
#define TriggerBoolean_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <map>
#include <stdexcept>

#include "TSystem.h"
#include "TObject.h"
#include "TTree.h"
#include "TString.h"
#include "Environ.h"

#include "TCUbits.h"
#include "RPscint.h"

class TriggerBoolean : public TObject
{
  public:
    TriggerBoolean(Int_t stg1_in, Int_t stg2_in, Int_t mipn_in, Int_t use_tcu);
    const char * Name(Int_t idx0); // boolean index --> name
    Int_t Idx(char * name0); // boolean name --> index
    Bool_t Fired(Int_t idx0);
    Bool_t Fired(char * name0);
    Bool_t FiredAlternate(Int_t idx0, Int_t stg1_in, Int_t stg2_in, Int_t mipn_in, Int_t use_tcu);
    Bool_t FiredAlternate(char * name0, Int_t stg1_in, Int_t stg2_in, Int_t mipn_in, Int_t use_tcu);
    void PrintParameters() { printf("STG1=%d STG2=%d MIPN=%d USE_TCU_BITS=%d\n",STG1,STG2,MIPN,USE_TCU_BITS); };
    void Diagnostic(Int_t runnum0, Int_t event0);

    // pointers to TCUbits and RPscint instances, set in the constructor
    TCUbits * TCU;
    RPscint * RPSCI;

    Int_t NBOOL;

    Float_t BBCvertex;


  private:
    std::map<std::string, Int_t> trg_idx; // trigger boolean name to an index
    std::map<Int_t, std::string> trg_name; // trigger boolean index to name
    Bool_t EOR,WOR,IT,ET;

    Int_t STG1, STG2, MIPN, USE_TCU_BITS;

    Int_t STG1_tmp, STG2_tmp, MIPN_tmp, USE_TCU_BITS_tmp;

    ClassDef(TriggerBoolean,1);
};

#endif
